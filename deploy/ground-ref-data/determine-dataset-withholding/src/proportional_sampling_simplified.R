# ==============================================================================
# PROPORTIONAL STRATIFIED SAMPLING: SIMPLIFIED VERSION (SINGLE DATA FRAME)
# ==============================================================================
# Purpose: Select ~20% of drone footprint groups to withhold for validation
#          such that withheld plots match the catalog's distribution
#
# Method: Greedy removal, random sampling, or hybrid (random + greedy) to find 
#         groups that best match the catalog distribution while meeting target 
#         percentage constraints
#
# Hybrid Method: 
#   - Phase 1: Random search to select ~30% of plots/trees
#   - Phase 2: Greedy removal from 30% down to ~20%
#   - Combines exploration (random) with refinement (greedy)
#
# Input: Single data frame with plots and their group_id assignments
#
# Author: [Your Name]
# Date: 2025-10-24
# Updated: 2025-10-28 (added hybrid method)
# ==============================================================================

library(tidyverse)
library(furrr)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Set random seed for reproducibility
set.seed(42)

# Define factorial combinations to track (optional, for diagnostics)
FACTORIAL_COMBINATIONS <- list(
  #list(var1 = "ecoregion", var2 = "sp_comp_group"),
  list(var1 = "mean_ba_live", var2 = "sp_comp_group"),
  list(var1 = "trees_per_ha", var2 = "sp_comp_group"),
  list(var1 = "ppt", var2 = "sp_comp_group"),
  list(var1 = "mean_ba_live", var2 = "ecoregion"),
  list(var1 = "trees_per_ha", var2 = "ecoregion"),
  list(var1 = "ppt", var2 = "ecoregion"),
  list(var1 = "trees_per_ha", var2 = "mean_ba_live"),
  
  list(var1 = "mean_ba_live", var2 = "pairing_tier"),
  list(var1 = "trees_per_ha", var2 = "pairing_tier"),
  list(var1 = "ppt", var2 = "pairing_tier"),
  list(var1 = "area_ha", var2 = "pairing_tier"),
  list(var1 = "pairing_tier", var2 = "sp_comp_group"),
    list(var1 = "pairing_tier", var2 = "ecoregion"),
    list(var1 = "pairing_tier", var2 = "project_name")

)

# Continuous variables to stratify
CONTINUOUS_VARS <- c("ppt", "trees_per_ha", "mean_ba_live", "area_ha")

# Categorical variables to stratify
CATEGORICAL_VARS <- c("ecoregion", "sp_comp_group", "project_name", "pairing_tier")

# Binning parameters
TARGET_PLOTS_PER_BIN <- 10   # Target average plots per bin in full catalog
MAX_BINS <- 5                # Maximum number of quantile bins to use
N_BINS_FACTORIAL <- 3        # Number of bins for continuous vars in factorial breakdowns

# Selection parameters
TARGET_PCT <- 17             # Target percentage to withhold
MIN_PCT <- 12                # Minimum acceptable percentage (plots/trees)
MAX_PCT <- 22                # Maximum acceptable percentage (plots/trees)
MIN_GROUPS_PCT <- 15         # Minimum acceptable percentage of groups
MAX_GROUPS_PCT <- 25         # Maximum acceptable percentage of groups

# Hybrid algorithm parameters (Phase 1: random search to ~30%)
PHASE1_TARGET_PCT <- 30      # Target percentage for Phase 1 random search
PHASE1_MIN_PCT <- 28         # Minimum acceptable percentage for Phase 1
PHASE1_MAX_PCT <- 35         # Maximum acceptable percentage for Phase 1

# Algorithm selection
ALGORITHM <- "greedy"        # Options: "greedy", "random", "hybrid"

# Greedy algorithm parameters (used if ALGORITHM = "greedy")
N_RUNS <- 60                 # Number of independent runs
TOP_K_CANDIDATES <- 5        # Consider top K groups when selecting which to remove
STOCHASTIC_TEMP <- 3         # Temperature for probability weighting
GROUP_SIZE_PENALTY_EXP <- 0.4  # Exponent for group size penalty in per-plot distortion calculation
                             # Lower values (0.3-0.5) = stronger penalty, favor larger groups
                             # Higher values (0.7-1.0) = less penalty, more distribution matching
                             # Previous (1.0) was over-penalizing, causing too many small groups to be selected
                             # Use 0.5 (sqrt) as balanced starting point

# Random sampling parameters (used if ALGORITHM = "random" or "hybrid")
N_RANDOM_SAMPLES <- 50000 * 100    # Number of random combinations to try for pure "random"
N_RANDOM_SAMPLES_PHASE1 <- 10000 # * 8  # Number for Phase 1 of hybrid (warm start)

# Shared parameters
PARALLELIZE <- TRUE          # Use parallel processing (uses all available cores)
FACTORIAL_WEIGHT <- 0.5      # Weight for factorial distribution in objective (0=ignore, 1=only factorial)

# ==============================================================================
# HELPER FUNCTIONS: QUANTILE CALCULATION
# ==============================================================================

#' Calculate optimal number of quantiles for a dataset
#'
#' @param n_plots_total Total number of plots in catalog
#' @return Optimal number of quantiles (between 2 and MAX_BINS)
calculate_optimal_n_quantiles <- function(n_plots_total) {
  
  # Use full catalog size to determine bins (not the selected subset)
  max_quantiles_allowed <- floor(n_plots_total / TARGET_PLOTS_PER_BIN)
  n_quantiles <- min(max_quantiles_allowed, MAX_BINS)
  n_quantiles <- max(2, n_quantiles)
  
  return(n_quantiles)
}

#' Create quantile bins for all continuous variables
#'
#' @param plots_df Data frame of plots
#' @param n_quantiles Number of quantiles to use
#' @return List with bin breaks and labels for each variable
prepare_quantile_bins <- function(plots_df, n_quantiles) {
  
  bin_info <- list()
  
  for (var in CONTINUOUS_VARS) {
    probs <- seq(0, 1, length.out = n_quantiles + 1)
    breaks <- quantile(plots_df[[var]], probs = probs, na.rm = TRUE)
    breaks <- unique(breaks)
    
    if (length(breaks) > 2) {
      labels <- paste0("Q", 1:(length(breaks) - 1))
    } else {
      labels <- c("Low", "High")
    }
    
    bin_info[[var]] <- list(
      breaks = breaks,
      labels = labels,
      n_bins = length(labels)
    )
  }
  
  return(bin_info)
}

#' Assign plots to quantile bins
#'
#' @param plots_df Data frame of plots
#' @param bin_info Bin information from prepare_quantile_bins()
#' @return Data frame with added quantile columns
assign_quantile_bins <- function(plots_df, bin_info) {
  
  result <- plots_df
  
  for (var in CONTINUOUS_VARS) {
    breaks <- bin_info[[var]]$breaks
    labels <- bin_info[[var]]$labels
    bin_col <- paste0(var, "_bin")
    
    result[[bin_col]] <- cut(result[[var]], 
                             breaks = breaks,
                             labels = labels,
                             include.lowest = TRUE)
  }
  
  return(result)
}

# ==============================================================================
# HELPER FUNCTIONS: DISTRIBUTION CALCULATION
# ==============================================================================

#' Calculate distribution of plots and trees across bins for all variables
#'
#' @param plots_df Data frame of plots (with bin columns and n_trees)
#' @return List of distributions for each variable
calculate_distribution <- function(plots_df) {
  
  distribution <- list()
  
  # Continuous variables (using bins)
  for (var in CONTINUOUS_VARS) {
    bin_col <- paste0(var, "_bin")
    
    dist <- plots_df |>
      filter(!is.na(.data[[bin_col]])) |>
      group_by(.data[[bin_col]]) |>
      summarise(
        n_plots = n(),
        n_trees = sum(n_trees, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(
        pct_plots = n_plots / sum(n_plots) * 100,
        pct_trees = n_trees / sum(n_trees) * 100
      ) |>
      rename(bin = !!bin_col)
    
    distribution[[var]] <- dist
  }
  
  # Categorical variables
  for (var in CATEGORICAL_VARS) {
    dist <- plots_df |>
      filter(!is.na(.data[[var]])) |>
      group_by(.data[[var]]) |>
      summarise(
        n_plots = n(),
        n_trees = sum(n_trees, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(
        pct_plots = n_plots / sum(n_plots) * 100,
        pct_trees = n_trees / sum(n_trees) * 100
      ) |>
      rename(category = !!var)
    
    distribution[[var]] <- dist
  }
  
  return(distribution)
}

#' Calculate distance between two distributions
#'
#' @param dist1 Distribution 1
#' @param dist2 Distribution 2
#' @param var Variable name
#' @return Distance value (0 = perfect match)
distribution_distance <- function(dist1, dist2, var) {
  
  key_col <- if (var %in% CONTINUOUS_VARS) "bin" else "category"
  
  comparison <- dist1 |>
    full_join(dist2, by = key_col, suffix = c("_1", "_2")) |>
    mutate(
      pct_plots_1 = replace_na(pct_plots_1, 0),
      pct_plots_2 = replace_na(pct_plots_2, 0),
      abs_diff = abs(pct_plots_1 - pct_plots_2)
    )
  
  total_distance <- sum(comparison$abs_diff)
  
  return(total_distance)
}

#' Calculate overall distribution distance across all variables
#'
#' @param dist_selected Distribution of selected plots
#' @param dist_catalog Distribution of catalog plots
#' @return Average distance across all variables
overall_distribution_distance <- function(dist_selected, dist_catalog) {
  
  all_vars <- c(CONTINUOUS_VARS, CATEGORICAL_VARS)
  
  distances <- map_dbl(all_vars, function(var) {
    distribution_distance(dist_selected[[var]], dist_catalog[[var]], var)
  })
  
  return(mean(distances))
}

#' Calculate factorial distribution distance (2-way combinations)
#'
#' @param factorial_selected Factorial distribution of selected plots
#' @param factorial_catalog Factorial distribution of catalog plots
#' @return Average distance across all factorial combinations
factorial_distribution_distance <- function(factorial_selected, factorial_catalog) {
  
  if (length(FACTORIAL_COMBINATIONS) == 0) {
    return(0)  # No factorial combinations to evaluate
  }
  
  distances <- map_dbl(names(factorial_catalog), function(combo_name) {
    
    # Get distributions for this combination
    catalog_dist <- factorial_catalog[[combo_name]]$distribution
    selected_dist <- factorial_selected[[combo_name]]$distribution
    
    # Combine and calculate percentage differences
    comparison <- catalog_dist |>
      rename(
        catalog_n_plots = n_plots,
        catalog_n_trees = n_trees
      ) |>
      full_join(
        selected_dist |>
          rename(
            selected_n_plots = n_plots,
            selected_n_trees = n_trees
          ),
        by = c("var1_value", "var2_value")
      ) |>
      mutate(
        catalog_n_plots = replace_na(catalog_n_plots, 0),
        selected_n_plots = replace_na(selected_n_plots, 0),
        catalog_pct_plots = catalog_n_plots / sum(catalog_n_plots) * 100,
        selected_pct_plots = selected_n_plots / sum(selected_n_plots) * 100,
        abs_diff = abs(catalog_pct_plots - selected_pct_plots)
      )
    
    total_distance <- sum(comparison$abs_diff)
    return(total_distance)
  })
  
  return(mean(distances))
}

#' Calculate combined distribution distance (single vars + factorial)
#'
#' @param dist_selected Single-variable distribution of selected plots
#' @param dist_catalog Single-variable distribution of catalog plots
#' @param factorial_selected Factorial distribution of selected plots (optional)
#' @param factorial_catalog Factorial distribution of catalog plots (optional)
#' @param factorial_weight Weight for factorial distance (default 0.5)
#' @return Combined distance score
combined_distribution_distance <- function(dist_selected, dist_catalog,
                                          factorial_selected = NULL, 
                                          factorial_catalog = NULL,
                                          factorial_weight = 0.5) {
  
  # Single-variable distance
  single_var_distance <- overall_distribution_distance(dist_selected, dist_catalog)
  
  # If no factorial distributions provided, return single-var distance
  if (is.null(factorial_selected) || is.null(factorial_catalog)) {
    return(single_var_distance)
  }
  
  # Factorial distance
  factorial_distance <- factorial_distribution_distance(factorial_selected, factorial_catalog)
  
  # Weighted combination
  combined_distance <- (1 - factorial_weight) * single_var_distance + 
                       factorial_weight * factorial_distance
  
  return(combined_distance)
}

# ==============================================================================
# HELPER FUNCTIONS: METRICS AND CONSTRAINTS
# ==============================================================================

#' Calculate current sampling metrics
#'
#' @param selected_groups Vector of selected group IDs
#' @param plots_df Data frame of plots with group_id column
#' @param total_groups Total number of groups
#' @param total_plots Total number of plots
#' @param total_trees Total number of trees
#' @return List with counts and percentages
calculate_metrics <- function(selected_groups, plots_df,
                             total_groups, total_plots, total_trees) {
  
  if (length(selected_groups) == 0) {
    return(list(
      n_groups = 0, pct_groups = 0,
      n_plots = 0, pct_plots = 0,
      n_trees = 0, pct_trees = 0
    ))
  }
  
  # Get plots in selected groups (simple!)
  selected_plots <- plots_df |>
    filter(group_id %in% selected_groups)
  
  # Calculate metrics
  list(
    n_groups = length(selected_groups),
    pct_groups = (length(selected_groups) / total_groups) * 100,
    n_plots = nrow(selected_plots),
    pct_plots = (nrow(selected_plots) / total_plots) * 100,
    n_trees = sum(selected_plots$n_trees, na.rm = TRUE),
    pct_trees = (sum(selected_plots$n_trees, na.rm = TRUE) / total_trees) * 100
  )
}

#' Check if metrics satisfy constraints
#'
#' @param metrics List from calculate_metrics()
#' @return TRUE if all constraints satisfied
check_constraints <- function(metrics) {
  
  groups_ok <- metrics$pct_groups >= MIN_GROUPS_PCT & metrics$pct_groups <= MAX_GROUPS_PCT
  plot_ok <- metrics$pct_plots >= MIN_PCT & metrics$pct_plots <= MAX_PCT
  tree_ok <- metrics$pct_trees >= MIN_PCT & metrics$pct_trees <= MAX_PCT
  
  return(groups_ok & plot_ok & tree_ok)
}

#' Calculate distance from target (20%)
#'
#' @param metrics List from calculate_metrics()
#' @return Distance score (lower is better)
calculate_target_distance <- function(metrics) {
  
  distance <- abs(metrics$pct_plots - TARGET_PCT) +
              abs(metrics$pct_trees - TARGET_PCT)
  
  return(distance)
}

# ==============================================================================
# ALGORITHM: RANDOM SAMPLING
# ==============================================================================

#' Select groups by random sampling
#'
#' @param plots_df Data frame of plots with group_id
#' @param catalog_distribution Catalog distribution
#' @param catalog_factorial Catalog factorial distribution
#' @param total_groups Total number of groups
#' @param total_plots Total number of plots
#' @param total_trees Total number of trees
#' @param required_groups Vector of group_ids that must be included
#' @param n_samples Number of random samples to try
#' @param min_pct Minimum acceptable percentage (default: MIN_PCT)
#' @param max_pct Maximum acceptable percentage (default: MAX_PCT)
#' @param target_pct Target percentage (default: TARGET_PCT)
#' @param previous_selected_plots Data frame of plots already selected in previous tiers (optional)
#' @return Data frame of valid solutions
random_sampling <- function(plots_df, catalog_distribution, catalog_factorial,
                           total_groups, total_plots, total_trees,
                           required_groups = c(), n_samples = N_RANDOM_SAMPLES,
                           min_pct = MIN_PCT, max_pct = MAX_PCT, 
                           target_pct = TARGET_PCT,
                           previous_selected_plots = NULL) {
  
  cat(sprintf("Random sampling: Generating %d random combinations...\n", n_samples))
  
  # Check if we're filling gaps from previous selections
  filling_gaps <- !is.null(previous_selected_plots)
  if (filling_gaps) {
    cat(sprintf("  Gap-filling mode: Including %d previously selected plots in distribution calculation\n",
                nrow(previous_selected_plots)))
  }
  
  all_groups <- unique(plots_df$group_id)
  removable_groups <- setdiff(all_groups, required_groups)
  n_required <- length(required_groups)
  
  # Calculate target number of groups based on min_pct and max_pct of plots
  # Average plots per group to inform the lower bound calculation
  avg_plots_per_group <- total_plots / total_groups
  target_n_groups_upper <- round(total_plots * max_pct / 100)
  target_n_groups_lower <- round((total_plots * min_pct / 100) / (avg_plots_per_group * 2))
  
  # Adjust bounds based on what's available (required + removable)
  max_n_groups <- min(target_n_groups_upper, n_required + length(removable_groups))
  min_n_groups <- max(n_required, round(target_n_groups_lower))
  
  # Range of removable groups to add
  min_removable <- max(0, min_n_groups - n_required)
  max_removable <- min(max_n_groups - n_required, length(removable_groups))
  
  # Pre-calculate group plot and tree counts for faster metric calculation
  group_stats <- plots_df |>
    group_by(group_id) |>
    summarise(
      n_plots = n(),
      n_trees = sum(n_trees, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat(sprintf("  Target groups: %d-%d (based on %.1f%%-%.1f%% of %d plots)\n",
              target_n_groups_lower, target_n_groups_upper,
              min_pct, max_pct, total_plots))
  cat(sprintf("  Sampling %d-%d groups (%d required + %d-%d removable from %d available)\n",
              min_n_groups, max_n_groups, n_required, 
              min_removable, max_removable, length(removable_groups)))
  
  # ============================================================================
  # Step 1: Generate all candidate group selections (not parallelized)
  # ============================================================================
  
  cat("  Step 1: Generating candidate group selections...\n")
  
  set.seed(42)
  candidates <- vector("list", n_samples)
  
  for (i in 1:n_samples) {
    # Set sample-specific seed for reproducibility
    set.seed(42 + i)
    
    # Randomly vary the number of removable groups to add
    n_removable_to_select <- sample(min_removable:max_removable, size = 1)
    
    # Randomly sample groups
    if (n_removable_to_select > 0) {
      sampled_removable <- sample(removable_groups, size = n_removable_to_select, replace = FALSE)
      selected_groups <- c(required_groups, sampled_removable)
    } else {
      selected_groups <- required_groups
    }
    
    candidates[[i]] <- list(
      sample_id = i,
      selected_groups = selected_groups
    )
  }
  
  cat(sprintf("  Generated %d candidate selections\n", length(candidates)))
  
  # ============================================================================
  # Step 2: Filter to valid candidates (parallelized by chunks)
  # ============================================================================
  
  cat("  Step 2: Filtering to valid candidates (checking plot/tree percentages)...\n")
  
  # Set up parallel processing for filtering
  if (PARALLELIZE) {
    n_cores <- availableCores()
    cat(sprintf("  Using %d cores for parallel filtering\n", n_cores))
    plan(multisession, workers = n_cores)
    
    # Split candidates into chunks (one per core)
    chunk_size <- ceiling(length(candidates) / n_cores)
    candidate_chunks <- split(candidates, ceiling(seq_along(candidates) / chunk_size))
    cat(sprintf("  Processing %d candidates in %d chunks of ~%d each\n", 
                length(candidates), length(candidate_chunks), chunk_size))
  } else {
    cat("  Running sequentially (PARALLELIZE = FALSE)\n")
    plan(sequential)
    candidate_chunks <- list(candidates)  # Single chunk for sequential
  }
  
  # Process each chunk in parallel
  valid_candidates_by_chunk <- future_map(candidate_chunks, function(chunk) {
    # Process all candidates in this chunk
    valid_in_chunk <- list()
    
    for (candidate in chunk) {
      # Calculate metrics using pre-computed group stats (fast - just summing)
      selected_stats <- group_stats |> filter(group_id %in% candidate$selected_groups)
      n_groups <- length(candidate$selected_groups)
      n_plots <- sum(selected_stats$n_plots)
      n_trees <- sum(selected_stats$n_trees)
      pct_groups <- (n_groups / total_groups) * 100
      pct_plots <- (n_plots / total_plots) * 100
      pct_trees <- (n_trees / total_trees) * 100
      
      # Check if in valid range (including group percentage)
      if (pct_groups >= MIN_GROUPS_PCT && pct_groups <= MAX_GROUPS_PCT &&
        pct_plots >= min_pct && pct_plots <= max_pct &&
        pct_trees >= min_pct && pct_trees <= max_pct) {
        valid_in_chunk[[length(valid_in_chunk) + 1]] <- list(
          sample_id = candidate$sample_id,
          selected_groups = candidate$selected_groups,
          n_groups = n_groups,
          pct_groups = pct_groups,
          n_plots = n_plots,
          pct_plots = pct_plots,
          n_trees = n_trees,
          pct_trees = pct_trees
        )
      }
      
      gc()

    }
    
    return(valid_in_chunk)
  }, .options = furrr_options(seed = TRUE))
  
  # Combine results from all chunks
  valid_candidates <- unlist(valid_candidates_by_chunk, recursive = FALSE)
  
  # Retain only the unique valid candidates
  valid_candidates <- unique(valid_candidates)
  
  # Reset to sequential
  plan(sequential)
  
  n_valid <- length(valid_candidates)
  cat(sprintf("  Found %d valid candidates (%.1f%% of total)\n",
              n_valid, 100 * n_valid / n_samples))
  
  if (n_valid == 0) {
    warning("No valid candidates found! Returning empty result.")
    return(tibble())
  }
  
  # ============================================================================
  # Step 3: Calculate distribution distances for valid candidates (parallelized)
  # ============================================================================
  
  cat(sprintf("  Step 3: Computing distribution distances for %d valid candidates...\n", n_valid))
  
  # Set up parallel processing
  if (PARALLELIZE) {
    n_cores <- availableCores()
    cat(sprintf("  Using %d cores for parallel distance calculation\n", n_cores))
    plan(multisession, workers = n_cores)
  } else {
    cat("  Running sequentially (PARALLELIZE = FALSE)\n")
    plan(sequential)
  }
  
  # Calculate distribution distances in parallel
  solutions <- future_map_dfr(valid_candidates, function(candidate) {
    
    # Get plots from current tier selection
    selected_plots <- plots_df |> filter(group_id %in% candidate$selected_groups)
    
    # If filling gaps, combine with previous selections for distribution calculation
    if (filling_gaps) {
      combined_plots <- bind_rows(previous_selected_plots, selected_plots)
      selected_dist <- calculate_distribution(combined_plots)
      selected_factorial <- calculate_factorial_distribution(combined_plots, reference_df = plots_df)
    } else {
      selected_dist <- calculate_distribution(selected_plots)
      selected_factorial <- calculate_factorial_distribution(selected_plots, reference_df = plots_df)
    }
    
    # Calculate distribution distance
    dist_distance <- combined_distribution_distance(selected_dist, catalog_distribution,
                                                    selected_factorial, catalog_factorial,
                                                    factorial_weight = FACTORIAL_WEIGHT)
    
    # Calculate target distance
    target_distance <- abs(candidate$pct_plots - target_pct) + 
                      abs(candidate$pct_trees - target_pct)
    
    # Return solution
    tibble(
      sample_id = candidate$sample_id,
      n_groups = candidate$n_groups,
      n_plots = candidate$n_plots,
      pct_plots = candidate$pct_plots,
      n_trees = candidate$n_trees,
      pct_trees = candidate$pct_trees,
      dist_distance = dist_distance,
      target_distance = target_distance,
      selected_groups = list(candidate$selected_groups)
    )
  }, .options = furrr_options(seed = TRUE))
  
  # Reset to sequential
  plan(sequential)
  
  cat(sprintf("  Found %d valid solutions (%.1f%% of samples)\n",
              nrow(solutions), 100 * nrow(solutions) / n_samples))
  
  return(solutions)
}

# ==============================================================================
# ALGORITHM: HYBRID (RANDOM + GREEDY)
# ==============================================================================

#' Select groups by hybrid random + greedy approach
#'
#' @param plots_df Data frame of plots with group_id
#' @param catalog_distribution Catalog distribution
#' @param catalog_factorial Catalog factorial distribution
#' @param bin_info Bin information
#' @param total_groups Total number of groups
#' @param total_plots Total number of plots
#' @param total_trees Total number of trees
#' @param required_groups Vector of group_ids that must be included
#' @param n_runs Number of independent hybrid runs
#' @param previous_selected_plots Data frame of plots already selected in previous tiers (optional)
#' @return Data frame of valid solutions
hybrid_sampling <- function(plots_df, catalog_distribution, catalog_factorial, bin_info,
                           total_groups, total_plots, total_trees,
                           required_groups = c(), n_runs = N_RUNS,
                           previous_selected_plots = NULL) {
  
  cat(sprintf("Hybrid sampling: Running %d independent hybrid runs...\n", n_runs))
  cat("  Phase 1: Random search to ~30%\n")
  cat("  Phase 2: Greedy removal to ~20%\n\n")
  
  # ============================================================================
  # PHASE 1: Random search (one search for all runs, parallelized internally)
  # ============================================================================
  
  cat("=== PHASE 1: Random search to ~30% ===\n")
  cat(sprintf("Generating %d random combinations...\n", N_RANDOM_SAMPLES_PHASE1))
  
  phase1_solutions <- random_sampling(
    plots_df = plots_df,
    catalog_distribution = catalog_distribution,
    catalog_factorial = catalog_factorial,
    total_groups = total_groups,
    total_plots = total_plots,
    total_trees = total_trees,
    required_groups = required_groups,
    n_samples = N_RANDOM_SAMPLES_PHASE1,
    min_pct = PHASE1_MIN_PCT,
    max_pct = PHASE1_MAX_PCT,
    target_pct = PHASE1_TARGET_PCT,
    previous_selected_plots = previous_selected_plots
  )
  
  # Note: random_sampling already filters internally using the Phase 1 constraints
  # So phase1_solutions already contains only valid Phase 1 solutions
  
  if (nrow(phase1_solutions) == 0) {
    stop("No valid Phase 1 solutions found! Try adjusting Phase 1 constraints.")
  }
  
  cat(sprintf("\nPhase 1 complete: Found %d valid solutions in range [%d%%, %d%%]\n",
              nrow(phase1_solutions), PHASE1_MIN_PCT, PHASE1_MAX_PCT))
  
  # Select THE BEST solution based on distribution distance
  phase1_best <- phase1_solutions |>
    arrange(dist_distance) |>
    slice(1)
  
  cat(sprintf("Best Phase 1 solution selected:\n"))
  cat(sprintf("  Groups: %d (%.1f%% plots, %.1f%% trees)\n",
              phase1_best$n_groups,
              phase1_best$pct_plots,
              phase1_best$pct_trees))
  cat(sprintf("  Distribution distance: %.2f\n\n",
              phase1_best$dist_distance))
  
  # Extract the initial groups for Phase 2
  initial_groups <- phase1_best$selected_groups[[1]]
  
  # ============================================================================
  # PHASE 2: Greedy removal (N_RUNS times from same Phase 1 solution)
  # ============================================================================
  
  cat(sprintf("=== PHASE 2: Greedy removal to ~20%% (%d parallel runs) ===\n", n_runs))
  cat(sprintf("All runs will start from the same Phase 1 solution\n"))
  cat(sprintf("Starting point: %d groups (%.1f%% plots, %.1f%% trees)\n\n",
              phase1_best$n_groups,
              phase1_best$pct_plots,
              phase1_best$pct_trees))
  
  # Set up parallel processing for Phase 2
  if (PARALLELIZE) {
    n_cores <- availableCores()
    cat(sprintf("  Using parallel processing with %d cores\n\n", n_cores))
    plan(multisession, workers = n_cores)
  } else {
    cat("  Running sequentially (PARALLELIZE = FALSE)\n\n")
    plan(sequential)
  }
  
  # Run N_RUNS greedy removals in parallel, all starting from same Phase 1 solution
  all_states <- future_map(1:n_runs, function(run) {
    
    # Run greedy removal starting from the best Phase 1 solution
    # Each run will be stochastic due to the random selection in greedy_removal
    states <- greedy_removal(
      plots_df = plots_df,
      catalog_distribution = catalog_distribution,
      catalog_factorial = catalog_factorial,
      bin_info = bin_info,
      total_groups = total_groups,
      total_plots = total_plots,
      total_trees = total_trees,
      run_id = run,
      required_groups = required_groups,
      initial_groups = initial_groups,
      previous_selected_plots = previous_selected_plots
    )
    
    # Add Phase 1 metadata
    states$run <- run
    states$phase1_n_groups <- phase1_best$n_groups
    states$phase1_pct_plots <- phase1_best$pct_plots
    states$phase1_pct_trees <- phase1_best$pct_trees
    states$phase1_dist_distance <- phase1_best$dist_distance
    
    return(states)
    
  }, .options = furrr_options(seed = TRUE))
  
  # Reset to sequential after parallel work
  plan(sequential)
  
  cat(sprintf("\n=== Phase 2 complete: %d greedy runs finished ===\n", length(all_states)))
  
  # Combine all states
  all_states_df <- bind_rows(all_states)
  
  # Report summary statistics across all runs
  final_states <- all_states_df |>
    group_by(run) |>
    slice(n()) |>
    ungroup()
  
  cat(sprintf("Final results across %d runs:\n", nrow(final_states)))
  cat(sprintf("  Plots: %.1f%% - %.1f%% (mean: %.1f%%)\n",
              min(final_states$pct_plots),
              max(final_states$pct_plots),
              mean(final_states$pct_plots)))
  cat(sprintf("  Trees: %.1f%% - %.1f%% (mean: %.1f%%)\n",
              min(final_states$pct_trees),
              max(final_states$pct_trees),
              mean(final_states$pct_trees)))
  cat(sprintf("  Distribution distance: %.2f - %.2f (mean: %.2f)\n\n",
              min(final_states$dist_distance),
              max(final_states$dist_distance),
              mean(final_states$dist_distance)))
  
  return(all_states_df)
}

# ==============================================================================
# ALGORITHM: GREEDY REMOVAL
# ==============================================================================

#' Run one iteration of greedy removal
#'
#' @param plots_df Data frame of plots (with bin columns and group_id)
#' @param catalog_distribution Distribution of full catalog
#' @param bin_info Bin information
#' @param total_groups Total number of groups
#' @param total_plots Total number of plots
#' @param total_trees Total number of trees
#' @param run_id Run identifier for random seed
#' @param required_groups Vector of group_ids that cannot be removed
#' @param initial_groups Optional vector of group_ids to start with (default NULL = all groups)
#' @param previous_selected_plots Data frame of plots already selected in previous tiers (optional)
#' @return Data frame of states at each iteration
greedy_removal <- function(plots_df, catalog_distribution, catalog_factorial, bin_info,
                          total_groups, total_plots, total_trees, run_id,
                          required_groups = c(), initial_groups = NULL,
                          previous_selected_plots = NULL) {
  
  # Set run-specific seed
  set.seed(42 + run_id)
  
  # Check if we're filling gaps from previous selections
  filling_gaps <- !is.null(previous_selected_plots)
  
  # Initialize: Use provided initial groups or ALL groups
  all_groups <- unique(plots_df$group_id)
  if (!is.null(initial_groups)) {
    selected_groups <- initial_groups
  } else {
    selected_groups <- all_groups
  }
  
  # Identify removable groups (all except required)
  removable_groups <- setdiff(selected_groups, required_groups)
  
  # Pre-calculate group statistics (plots and trees per group) - OPTIMIZATION
  group_stats <- plots_df |>
    group_by(group_id) |>
    summarise(
      n_plots = n(),
      n_trees = sum(n_trees, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create a named vector for fast lookup - OPTIMIZATION
  group_n_plots <- setNames(group_stats$n_plots, group_stats$group_id)
  group_n_trees <- setNames(group_stats$n_trees, group_stats$group_id)
  
  # Storage for state history
  state_history <- list()
  iteration <- 1
  
  # Calculate and save initial state (100% selected)
  metrics <- calculate_metrics(selected_groups, plots_df,
                               total_groups, total_plots, total_trees)
  selected_plots <- plots_df |> filter(group_id %in% selected_groups)
  
  # If filling gaps, combine with previous selections for distribution calculation
  if (filling_gaps) {
    combined_plots <- bind_rows(previous_selected_plots, selected_plots)
    selected_dist <- calculate_distribution(combined_plots)
    selected_factorial <- calculate_factorial_distribution(combined_plots, reference_df = plots_df)
  } else {
    selected_dist <- calculate_distribution(selected_plots)
    selected_factorial <- calculate_factorial_distribution(selected_plots, reference_df = plots_df)
  }
  
  dist_distance <- combined_distribution_distance(selected_dist, catalog_distribution,
                                                  selected_factorial, catalog_factorial,
                                                  factorial_weight = FACTORIAL_WEIGHT)
  
  state_history[[iteration]] <- list(
    iteration = iteration,
    selected_groups = selected_groups,
    metrics = metrics,
    dist_distance = dist_distance
  )
  
  # Greedy removal loop
  while (metrics$pct_groups > MAX_GROUPS_PCT || 
         metrics$pct_plots > MAX_PCT || 
         metrics$pct_trees > MAX_PCT) {
    
    # Pre-allocate candidate_scores for efficiency - OPTIMIZATION
    n_removable <- length(removable_groups)
    candidate_scores <- tibble(
      group_id = rep(NA_real_, n_removable),
      dist_distance = rep(NA_real_, n_removable),
      per_plot_distortion = rep(NA_real_, n_removable),
      n_plots = rep(NA_integer_, n_removable),
      test_pct_plots = rep(NA_real_, n_removable),
      test_pct_trees = rep(NA_real_, n_removable)
    )
    
    valid_idx <- 0  # Track number of valid candidates
    
    for (candidate_group in removable_groups) {
      
      # Fast lookup of group stats - OPTIMIZATION
      candidate_n_plots <- group_n_plots[as.character(candidate_group)]
      candidate_n_trees <- group_n_trees[as.character(candidate_group)]
      
      # Quick check: Can we remove this group without violating constraints? - OPTIMIZATION
      # Calculate what metrics would be if we remove this group (incremental update)
      test_n_groups <- length(selected_groups) - 1
      test_n_plots <- metrics$n_plots - candidate_n_plots
      test_n_trees <- metrics$n_trees - candidate_n_trees
      test_pct_groups <- (test_n_groups / total_groups) * 100
      test_pct_plots <- (test_n_plots / total_plots) * 100
      test_pct_trees <- (test_n_trees / total_trees) * 100
      
      # Check if removal would violate minimum constraints - EARLY EXIT
      if (test_pct_groups < MIN_GROUPS_PCT ||
          test_pct_plots < MIN_PCT || 
          test_pct_trees < MIN_PCT) {
        next  # Skip - can't remove without violating constraints
      }
      
      # Only calculate expensive distribution if constraints are satisfied - OPTIMIZATION
      test_groups <- setdiff(selected_groups, candidate_group)
      test_plots <- plots_df |> filter(group_id %in% test_groups)
      
      # If filling gaps, combine with previous selections for distribution calculation
      if (filling_gaps) {
        combined_test_plots <- bind_rows(previous_selected_plots, test_plots)
        test_dist <- calculate_distribution(combined_test_plots)
        test_factorial <- calculate_factorial_distribution(combined_test_plots, reference_df = plots_df)
      } else {
        test_dist <- calculate_distribution(test_plots)
        test_factorial <- calculate_factorial_distribution(test_plots, reference_df = plots_df)
      }
      
      test_dist_distance <- combined_distribution_distance(test_dist, catalog_distribution,
                                                           test_factorial, catalog_factorial,
                                                           factorial_weight = FACTORIAL_WEIGHT)
      
      # Normalize by group size with tunable penalty exponent
      per_plot_distortion <- test_dist_distance / (candidate_n_plots ^ GROUP_SIZE_PENALTY_EXP)
      
      # Store candidate in pre-allocated tibble - OPTIMIZATION
      valid_idx <- valid_idx + 1
      candidate_scores$group_id[valid_idx] <- candidate_group
      candidate_scores$dist_distance[valid_idx] <- test_dist_distance
      candidate_scores$per_plot_distortion[valid_idx] <- per_plot_distortion
      candidate_scores$n_plots[valid_idx] <- candidate_n_plots
      candidate_scores$test_pct_plots[valid_idx] <- test_pct_plots
      candidate_scores$test_pct_trees[valid_idx] <- test_pct_trees
    }
    
    # Remove unused rows from pre-allocated tibble - OPTIMIZATION
    if (valid_idx == 0) {
      break  # No valid removals possible
    }
    candidate_scores <- candidate_scores |> slice(1:valid_idx)
    
    # Sort by per-plot distortion (ascending)
    candidate_scores <- candidate_scores |>
      arrange(per_plot_distortion)
    
    # Stochastic selection from top K candidates
    n_candidates <- min(TOP_K_CANDIDATES, nrow(candidate_scores))
    top_candidates <- candidate_scores |>
      slice(1:n_candidates)
    
    # Calculate selection probabilities
    weights <- exp(-((1:n_candidates) - 1) / STOCHASTIC_TEMP)
    weights <- weights / sum(weights)
    
    # Select group to remove
    selected_idx <- sample(1:n_candidates, size = 1, prob = weights)
    group_to_remove <- top_candidates$group_id[selected_idx]
    
    # Remove the selected group
    selected_groups <- setdiff(selected_groups, group_to_remove)
    removable_groups <- setdiff(removable_groups, group_to_remove)
    
    # Update metrics and save state
    iteration <- iteration + 1
    metrics <- calculate_metrics(selected_groups, plots_df,
                                total_groups, total_plots, total_trees)
    selected_plots <- plots_df |> filter(group_id %in% selected_groups)
    
    # If filling gaps, combine with previous selections for distribution calculation
    if (filling_gaps) {
      combined_plots <- bind_rows(previous_selected_plots, selected_plots)
      selected_dist <- calculate_distribution(combined_plots)
      selected_factorial <- calculate_factorial_distribution(combined_plots, reference_df = plots_df)
    } else {
      selected_dist <- calculate_distribution(selected_plots)
      selected_factorial <- calculate_factorial_distribution(selected_plots, reference_df = plots_df)
    }
    
    dist_distance <- combined_distribution_distance(selected_dist, catalog_distribution,
                                                    selected_factorial, catalog_factorial,
                                                    factorial_weight = FACTORIAL_WEIGHT)
    
    state_history[[iteration]] <- list(
      iteration = iteration,
      selected_groups = selected_groups,
      metrics = metrics,
      dist_distance = dist_distance
    )
  }
  
  # Convert state history to data frame
  states_df <- map_dfr(state_history, function(state) {
    tibble(
      iteration = state$iteration,
      n_groups = length(state$selected_groups),
      n_plots = state$metrics$n_plots,
      pct_plots = state$metrics$pct_plots,
      n_trees = state$metrics$n_trees,
      pct_trees = state$metrics$pct_trees,
      dist_distance = state$dist_distance,
      selected_groups = list(state$selected_groups)
    )
  })
  
  return(states_df)
}

# ==============================================================================
# DIAGNOSTIC FUNCTIONS
# ==============================================================================

#' Create detailed distribution comparison
#'
#' @param catalog_dist Distribution of full catalog
#' @param selected_dist Distribution of selected plots
#' @return List of comparison tables
create_distribution_comparison <- function(catalog_dist, selected_dist) {
  
  comparisons <- list()
  
  # Continuous variables
  for (var in CONTINUOUS_VARS) {
    comparison <- catalog_dist[[var]] |>
      rename(
        catalog_n_plots = n_plots,
        catalog_n_trees = n_trees,
        catalog_pct_plots = pct_plots,
        catalog_pct_trees = pct_trees
      ) |>
      full_join(
        selected_dist[[var]] |>
          rename(
            selected_n_plots = n_plots,
            selected_n_trees = n_trees,
            selected_pct_plots = pct_plots,
            selected_pct_trees = pct_trees
          ),
        by = "bin"
      ) |>
      mutate(
        catalog_n_plots = replace_na(catalog_n_plots, 0),
        catalog_n_trees = replace_na(catalog_n_trees, 0),
        catalog_pct_plots = replace_na(catalog_pct_plots, 0),
        catalog_pct_trees = replace_na(catalog_pct_trees, 0),
        selected_n_plots = replace_na(selected_n_plots, 0),
        selected_n_trees = replace_na(selected_n_trees, 0),
        selected_pct_plots = replace_na(selected_pct_plots, 0),
        selected_pct_trees = replace_na(selected_pct_trees, 0),
        diff_pct_plots = selected_pct_plots - catalog_pct_plots,
        diff_pct_trees = selected_pct_trees - catalog_pct_trees,
        abs_diff_plots = abs(diff_pct_plots),
        abs_diff_trees = abs(diff_pct_trees)
      ) |>
      arrange(bin)
    
    comparisons[[var]] <- comparison
  }
  
  # Categorical variables
  for (var in CATEGORICAL_VARS) {
    comparison <- catalog_dist[[var]] |>
      rename(
        catalog_n_plots = n_plots,
        catalog_n_trees = n_trees,
        catalog_pct_plots = pct_plots,
        catalog_pct_trees = pct_trees
      ) |>
      full_join(
        selected_dist[[var]] |>
          rename(
            selected_n_plots = n_plots,
            selected_n_trees = n_trees,
            selected_pct_plots = pct_plots,
            selected_pct_trees = pct_trees
          ),
        by = "category"
      ) |>
      mutate(
        catalog_n_plots = replace_na(catalog_n_plots, 0),
        catalog_n_trees = replace_na(catalog_n_trees, 0),
        catalog_pct_plots = replace_na(catalog_pct_plots, 0),
        catalog_pct_trees = replace_na(catalog_pct_trees, 0),
        selected_n_plots = replace_na(selected_n_plots, 0),
        selected_n_trees = replace_na(selected_n_trees, 0),
        selected_pct_plots = replace_na(selected_pct_plots, 0),
        selected_pct_trees = replace_na(selected_pct_trees, 0),
        diff_pct_plots = selected_pct_plots - catalog_pct_plots,
        diff_pct_trees = selected_pct_trees - catalog_pct_trees,
        abs_diff_plots = abs(diff_pct_plots),
        abs_diff_trees = abs(diff_pct_trees)
      ) |>
      arrange(desc(catalog_pct_plots))
    
    comparisons[[var]] <- comparison
  }
  
  return(comparisons)
}

#' Create summary statistics comparison
#'
#' @param catalog_plots All plots
#' @param selected_plots Selected plots
#' @return Data frame comparing summary statistics
create_summary_stats_comparison <- function(catalog_plots, selected_plots) {
  
  stats <- map_dfr(CONTINUOUS_VARS, function(var) {
    tibble(
      variable = var,
      catalog_mean = mean(catalog_plots[[var]], na.rm = TRUE),
      selected_mean = mean(selected_plots[[var]], na.rm = TRUE),
      catalog_sd = sd(catalog_plots[[var]], na.rm = TRUE),
      selected_sd = sd(selected_plots[[var]], na.rm = TRUE),
      catalog_min = min(catalog_plots[[var]], na.rm = TRUE),
      selected_min = min(selected_plots[[var]], na.rm = TRUE),
      catalog_max = max(catalog_plots[[var]], na.rm = TRUE),
      selected_max = max(selected_plots[[var]], na.rm = TRUE)
    ) |>
      mutate(
        mean_diff = selected_mean - catalog_mean,
        sd_diff = selected_sd - catalog_sd,
        mean_pct_diff = (mean_diff / catalog_mean) * 100,
        sd_pct_diff = (sd_diff / catalog_sd) * 100
      )
  })
  
  return(stats)
}

#' Calculate factorial distribution (2-way combinations)
#'
#' @param plots_df Data frame of plots (with n_trees and stratification variables)
#' @param reference_df Optional reference data frame to use for calculating bin breaks (for continuous vars)
#' @return List of factorial distributions for each combination
calculate_factorial_distribution <- function(plots_df, reference_df = NULL) {
  
  # If no reference provided, use the data itself
  if (is.null(reference_df)) {
    reference_df <- plots_df
  }
  
  factorial_dists <- list()
  
  for (combo in FACTORIAL_COMBINATIONS) {
    var1 <- combo$var1
    var2 <- combo$var2
    
    # Create a working copy for this combination
    plots_working <- plots_df
    
    # For continuous variables, create factorial-specific bins
    if (var1 %in% CONTINUOUS_VARS) {
      probs <- seq(0, 1, length.out = N_BINS_FACTORIAL + 1)
      breaks <- quantile(reference_df[[var1]], probs = probs, na.rm = TRUE)
      breaks <- unique(breaks)
      
      if (length(breaks) > 1) {
        plots_working[[paste0(var1, "_factorial_bin")]] <- cut(
          plots_working[[var1]],
          breaks = breaks,
          include.lowest = TRUE,
          dig.lab = 3
        )
      } else {
        plots_working[[paste0(var1, "_factorial_bin")]] <- as.factor(paste0(var1, "_all"))
      }
      col1 <- paste0(var1, "_factorial_bin")
    } else {
      col1 <- var1
    }
    
    if (var2 %in% CONTINUOUS_VARS) {
      probs <- seq(0, 1, length.out = N_BINS_FACTORIAL + 1)
      breaks <- quantile(reference_df[[var2]], probs = probs, na.rm = TRUE)
      breaks <- unique(breaks)
      
      if (length(breaks) > 1) {
        plots_working[[paste0(var2, "_factorial_bin")]] <- cut(
          plots_working[[var2]],
          breaks = breaks,
          include.lowest = TRUE,
          dig.lab = 3
        )
      } else {
        plots_working[[paste0(var2, "_factorial_bin")]] <- as.factor(paste0(var2, "_all"))
      }
      col2 <- paste0(var2, "_factorial_bin")
    } else {
      col2 <- var2
    }
    
    # Calculate 2-way distribution
    dist <- plots_working |>
      filter(!is.na(.data[[col1]]) & !is.na(.data[[col2]])) |>
      group_by(.data[[col1]], .data[[col2]]) |>
      summarise(
        n_plots = n(),
        n_trees = sum(n_trees, na.rm = TRUE),
        .groups = "drop"
      ) |>
      rename(
        var1_value = !!col1,
        var2_value = !!col2
      )
    
    factorial_dists[[paste(var1, var2, sep = "_x_")]] <- list(
      var1 = var1,
      var2 = var2,
      distribution = dist
    )
  }
  
  return(factorial_dists)
}

#' Create factorial distribution visualizations
#'
#' @param catalog_factorial Factorial distribution of catalog
#' @param selected_factorial Factorial distribution of selected plots
#' @return List of ggplot objects (heatmaps)
create_factorial_plots <- function(catalog_factorial, selected_factorial) {
  
  plots <- list()
  
  for (combo_name in names(catalog_factorial)) {
    combo_info <- catalog_factorial[[combo_name]]
    var1 <- combo_info$var1
    var2 <- combo_info$var2
    
    # Combine catalog and selected distributions
    combined <- catalog_factorial[[combo_name]]$distribution |>
      rename(
        catalog_n_plots = n_plots,
        catalog_n_trees = n_trees
      ) |>
      full_join(
        selected_factorial[[combo_name]]$distribution |>
          rename(
            selected_n_plots = n_plots,
            selected_n_trees = n_trees
          ),
        by = c("var1_value", "var2_value")
      ) |>
      mutate(
        catalog_n_plots = replace_na(catalog_n_plots, 0),
        catalog_n_trees = replace_na(catalog_n_trees, 0),
        selected_n_plots = replace_na(selected_n_plots, 0),
        selected_n_trees = replace_na(selected_n_trees, 0),
        # Calculate selection percentages
        pct_plots_selected = ifelse(catalog_n_plots > 0, 
                                     100 * selected_n_plots / catalog_n_plots, 
                                     0),
        pct_trees_selected = ifelse(catalog_n_trees > 0,
                                     100 * selected_n_trees / catalog_n_trees,
                                     0)
      )
    
    # Create heatmap for plots
    p_plots <- ggplot(combined, aes(x = var1_value, y = var2_value)) +
      geom_tile(aes(fill = pct_plots_selected), color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("C:%d\nS:%d\n%.0f%%", 
                                     catalog_n_plots, 
                                     selected_n_plots,
                                     pct_plots_selected)),
                size = 3, lineheight = 0.8) +
      scale_fill_gradient2(low = "#d73027", mid = "#fee08b", high = "#1a9850",
                          midpoint = 20,
                          limits = c(0, 100),
                          name = "% Selected") +
      labs(
        title = paste("Factorial Distribution:", var1, "×", var2, "(Plots)"),
        subtitle = "Colors show % of catalog selected in each bin",
        x = var1,
        y = var2
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold")
      )
    
    # Create heatmap for trees
    p_trees <- ggplot(combined, aes(x = var1_value, y = var2_value)) +
      geom_tile(aes(fill = pct_trees_selected), color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("C:%d\nS:%d\n%.0f%%", 
                                     catalog_n_trees, 
                                     selected_n_trees,
                                     pct_trees_selected)),
                size = 3, lineheight = 0.8) +
      scale_fill_gradient2(low = "#d73027", mid = "#fee08b", high = "#1a9850",
                          midpoint = 20,
                          limits = c(0, 100),
                          name = "% Selected") +
      labs(
        title = paste("Factorial Distribution:", var1, "×", var2, "(Trees)"),
        subtitle = "Colors show % of catalog selected in each bin",
        x = var1,
        y = var2
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold")
      )
    
    plots[[paste0(combo_name, "_plots")]] <- p_plots
    plots[[paste0(combo_name, "_trees")]] <- p_trees
  }
  
  return(plots)
}

#' Create dual factorial distribution visualizations (reference + tier)
#'
#' @param reference_factorial Factorial distribution of reference catalog
#' @param tier_factorial Factorial distribution of tier catalog
#' @param selected_factorial Factorial distribution of selected plots
#' @return List of ggplot objects (dual heatmaps side by side)
create_dual_factorial_plots <- function(reference_factorial, tier_factorial, selected_factorial) {
  
  library(patchwork)
  
  plots <- list()
  
  for (combo_name in names(reference_factorial)) {
    combo_info <- reference_factorial[[combo_name]]
    var1 <- combo_info$var1
    var2 <- combo_info$var2
    
    # === REFERENCE CATALOG ===
    combined_ref <- reference_factorial[[combo_name]]$distribution |>
      rename(
        catalog_n_plots = n_plots,
        catalog_n_trees = n_trees
      ) |>
      full_join(
        selected_factorial[[combo_name]]$distribution |>
          rename(
            selected_n_plots = n_plots,
            selected_n_trees = n_trees
          ),
        by = c("var1_value", "var2_value")
      ) |>
      mutate(
        catalog_n_plots = replace_na(catalog_n_plots, 0),
        catalog_n_trees = replace_na(catalog_n_trees, 0),
        selected_n_plots = replace_na(selected_n_plots, 0),
        selected_n_trees = replace_na(selected_n_trees, 0),
        pct_plots_selected_raw = ifelse(catalog_n_plots > 0, 
                                     100 * selected_n_plots / catalog_n_plots, 
                                     0),
        pct_trees_selected_raw = ifelse(catalog_n_trees > 0,
                                     100 * selected_n_trees / catalog_n_trees,
                                     0),
        # Cap at 40 for color scale
        pct_plots_selected = pmin(pct_plots_selected_raw, 40),
        pct_trees_selected = pmin(pct_trees_selected_raw, 40),
        # Labels show "40+" for capped values
        pct_plots_label = ifelse(pct_plots_selected_raw > 40, "40+", sprintf("%.0f%%", pct_plots_selected_raw)),
        pct_trees_label = ifelse(pct_trees_selected_raw > 40, "40+", sprintf("%.0f%%", pct_trees_selected_raw))
      )
    
    # === TIER CATALOG ===
    combined_tier <- tier_factorial[[combo_name]]$distribution |>
      rename(
        catalog_n_plots = n_plots,
        catalog_n_trees = n_trees
      ) |>
      full_join(
        selected_factorial[[combo_name]]$distribution |>
          rename(
            selected_n_plots = n_plots,
            selected_n_trees = n_trees
          ),
        by = c("var1_value", "var2_value")
      ) |>
      mutate(
        catalog_n_plots = replace_na(catalog_n_plots, 0),
        catalog_n_trees = replace_na(catalog_n_trees, 0),
        selected_n_plots = replace_na(selected_n_plots, 0),
        selected_n_trees = replace_na(selected_n_trees, 0),
        pct_plots_selected_raw = ifelse(catalog_n_plots > 0, 
                                     100 * selected_n_plots / catalog_n_plots, 
                                     0),
        pct_trees_selected_raw = ifelse(catalog_n_trees > 0,
                                     100 * selected_n_trees / catalog_n_trees,
                                     0),
        # Cap at 40 for color scale
        pct_plots_selected = pmin(pct_plots_selected_raw, 40),
        pct_trees_selected = pmin(pct_trees_selected_raw, 40),
        # Labels show "40+" for capped values
        pct_plots_label = ifelse(pct_plots_selected_raw > 40, "40+", sprintf("%.0f%%", pct_plots_selected_raw)),
        pct_trees_label = ifelse(pct_trees_selected_raw > 40, "40+", sprintf("%.0f%%", pct_trees_selected_raw))
      )
    
    # Create reference heatmap for plots
    p_ref_plots <- ggplot(combined_ref, aes(x = var1_value, y = var2_value)) +
      geom_tile(aes(fill = pct_plots_selected), color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("C:%d\nS:%d\n%s", 
                                     catalog_n_plots, 
                                     selected_n_plots,
                                     pct_plots_label)),
                size = 2.5, lineheight = 0.8) +
      scale_fill_gradient2(low = "#d73027", mid = "#fee08b", high = "#1a9850",
                          midpoint = 20,
                          limits = c(0, 40),
                          name = "% Selected") +
      labs(
        title = "Reference (All Tiers)",
        subtitle = "Optimized for this distribution",
        x = var1,
        y = var2
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 10)
      )
    
    # Create tier heatmap for plots
    p_tier_plots <- ggplot(combined_tier, aes(x = var1_value, y = var2_value)) +
      geom_tile(aes(fill = pct_plots_selected), color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("C:%d\nS:%d\n%s", 
                                     catalog_n_plots, 
                                     selected_n_plots,
                                     pct_plots_label)),
                size = 2.5, lineheight = 0.8) +
      scale_fill_gradient2(low = "#d73027", mid = "#fee08b", high = "#1a9850",
                          midpoint = 20,
                          limits = c(0, 40),
                          name = "% Selected") +
      labs(
        title = "Current Tier Only",
        subtitle = "Actual within-tier distribution",
        x = var1,
        y = var2
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 10)
      )
    
    # Combine plots side by side
    p_plots_combined <- p_ref_plots + p_tier_plots +
      plot_annotation(
        title = paste("Factorial Distribution:", var1, "×", var2, "(Plots)"),
        theme = theme(plot.title = element_text(face = "bold", size = 12))
      )
    
    # Create reference heatmap for trees
    p_ref_trees <- ggplot(combined_ref, aes(x = var1_value, y = var2_value)) +
      geom_tile(aes(fill = pct_trees_selected), color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("C:%d\nS:%d\n%s", 
                                     catalog_n_trees, 
                                     selected_n_trees,
                                     pct_trees_label)),
                size = 2.5, lineheight = 0.8) +
      scale_fill_gradient2(low = "#d73027", mid = "#fee08b", high = "#1a9850",
                          midpoint = 20,
                          limits = c(0, 40),
                          name = "% Selected") +
      labs(
        title = "Reference (All Tiers)",
        subtitle = "Optimized for this distribution",
        x = var1,
        y = var2
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 10)
      )
    
    # Create tier heatmap for trees
    p_tier_trees <- ggplot(combined_tier, aes(x = var1_value, y = var2_value)) +
      geom_tile(aes(fill = pct_trees_selected), color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("C:%d\nS:%d\n%s", 
                                     catalog_n_trees, 
                                     selected_n_trees,
                                     pct_trees_label)),
                size = 2.5, lineheight = 0.8) +
      scale_fill_gradient2(low = "#d73027", mid = "#fee08b", high = "#1a9850",
                          midpoint = 20,
                          limits = c(0, 40),
                          name = "% Selected") +
      labs(
        title = "Current Tier Only",
        subtitle = "Actual within-tier distribution",
        x = var1,
        y = var2
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 10)
      )
    
    # Combine plots side by side
    p_trees_combined <- p_ref_trees + p_tier_trees +
      plot_annotation(
        title = paste("Factorial Distribution:", var1, "×", var2, "(Trees)"),
        theme = theme(plot.title = element_text(face = "bold", size = 12))
      )
    
    plots[[paste0(combo_name, "_plots")]] <- p_plots_combined
    plots[[paste0(combo_name, "_trees")]] <- p_trees_combined
  }
  
  return(plots)
}

# ==============================================================================
# MAIN EXECUTION FUNCTION
# ==============================================================================

#' Main function to select withheld groups
#'
#' @param plots_df Data frame with columns: plot_id, group_id, attributes, n_trees
#' @param required_groups Vector of group_ids that must be in the withheld set (for hierarchical stratification)
#' @param reference_distribution Optional: Pre-calculated distribution to match (instead of plots_df's distribution)
#' @param reference_factorial Optional: Pre-calculated factorial distribution to match
#' @param reference_plots_df Optional: Full reference dataset for binning consistency (used when reference_distribution is provided)
#' @param previous_selected_plots Optional: Plots already selected in previous tiers (for gap-filling mode)
#' @return List containing withheld and training plots with diagnostics
select_withheld_groups <- function(plots_df, required_groups = c(),
                                  reference_distribution = NULL,
                                  reference_factorial = NULL,
                                  reference_plots_df = NULL,
                                  previous_selected_plots = NULL) {
  
  cat("Starting proportional stratified group selection (SIMPLIFIED)...\n")
  cat(sprintf("Method: %s\n", 
              if (ALGORITHM == "greedy") {
                "Greedy removal with stochastic selection"
              } else if (ALGORITHM == "random") {
                "Random sampling"
              } else if (ALGORITHM == "hybrid") {
                "Hybrid (Random search to 30% + Greedy removal to 20%)"
              } else {
                ALGORITHM
              }))
  
  # Check if using reference distribution
  use_reference_distribution <- !is.null(reference_distribution)
  if (use_reference_distribution) {
    cat("Using provided reference distribution (cross-tier matching)\n")
  } else {
    cat("Using within-tier distribution (tier-specific matching)\n")
  }
  
  # Check if using gap-filling mode
  filling_gaps <- !is.null(previous_selected_plots)
  if (filling_gaps) {
    cat(sprintf("Gap-filling mode: Considering %d previously selected plots in optimization\n",
                nrow(previous_selected_plots)))
  }
  cat("\n")
  
  # ============================================================================
  # STEP 1: Data preparation
  # ============================================================================
  
  cat("Step 1: Preparing data and calculating bins...\n")
  
  # Calculate totals
  total_groups <- n_distinct(plots_df$group_id)
  total_plots <- nrow(plots_df)
  total_trees <- sum(plots_df$n_trees, na.rm = TRUE)
  
  cat(sprintf("  Total groups: %d\n", total_groups))
  cat(sprintf("  Total plots: %d\n", total_plots))
  cat(sprintf("  Total trees: %d\n", total_trees))
  
  # Determine which dataset to use for binning
  if (use_reference_distribution) {
    if (is.null(reference_plots_df)) {
      stop("reference_plots_df must be provided when using reference_distribution")
    }
    binning_df <- reference_plots_df
    cat("  Using reference dataset for bin calculation\n")
  } else {
    binning_df <- plots_df
    cat("  Using current tier for bin calculation\n")
  }
  
  # Calculate optimal number of quantiles
  n_quantiles <- calculate_optimal_n_quantiles(nrow(binning_df))
  cat(sprintf("  Using %d quantiles for continuous variables\n", n_quantiles))
  cat(sprintf("  Target: ≥%d plots per bin in catalog\n", TARGET_PLOTS_PER_BIN))
  
  # Prepare quantile bins (using reference data if provided)
  bin_info <- prepare_quantile_bins(binning_df, n_quantiles)
  plots_df_binned <- assign_quantile_bins(plots_df, bin_info)
  
  # Calculate catalog distribution (baseline to match)
  if (use_reference_distribution) {
    catalog_distribution <- reference_distribution
    catalog_factorial_ref <- reference_factorial
    cat("  Using provided reference distribution\n")
  } else {
    catalog_distribution <- calculate_distribution(plots_df_binned)
    catalog_factorial_ref <- calculate_factorial_distribution(plots_df_binned)
    cat("  Calculated catalog distribution from current tier\n")
  }
  
  # For diagnostics, also calculate tier-specific distribution if using reference
  if (use_reference_distribution) {
    tier_distribution <- calculate_distribution(plots_df_binned)
    # Use reference dataset for bin calculation to match reference bins in dual plots
    tier_factorial <- calculate_factorial_distribution(plots_df_binned, reference_df = reference_plots_df)
  } else {
    tier_distribution <- NULL
    tier_factorial <- NULL
  }
  
  cat("  Catalog distribution prepared\n\n")
  
  # ============================================================================
  # STEP 1.5: Handle required groups (for hierarchical stratification)
  # ============================================================================
  
  # Filter required groups to those actually present in this dataset
  required_groups <- intersect(required_groups, unique(plots_df$group_id))
  
  if (length(required_groups) > 0) {
    cat(sprintf("Step 1.5: Checking required groups (%d groups)...\n", 
                length(required_groups)))
    
    # Calculate metrics for required groups
    required_metrics <- calculate_metrics(required_groups, plots_df_binned,
                                          total_groups, total_plots, total_trees)
    
    cat(sprintf("  Required groups contain: %d plots (%.1f%%), %d trees (%.1f%%)\n",
                required_metrics$n_plots, required_metrics$pct_plots,
                required_metrics$n_trees, required_metrics$pct_trees))
    
    # Check if required groups already meet or exceed targets
    if (required_metrics$pct_plots >= MIN_PCT && 
        required_metrics$pct_plots <= MAX_PCT &&
        required_metrics$pct_trees >= MIN_PCT && 
        required_metrics$pct_trees <= MAX_PCT) {
      
      # Already in acceptable range! Just use required groups
      cat(sprintf("  ✓ Required groups already in acceptable range [%d%%, %d%%]\n",
                  MIN_PCT, MAX_PCT))
      cat("  Using required groups as final selection (no additional stratification needed)\n\n")
      
      # Calculate distribution for reporting
      required_plots <- plots_df_binned |> filter(group_id %in% required_groups)
      required_dist <- calculate_distribution(required_plots)
      dist_distance <- overall_distribution_distance(required_dist, catalog_distribution)
      
      # Create simple results object
      training_groups <- setdiff(unique(plots_df$group_id), required_groups)
      training_plots <- plots_df_binned |> filter(group_id %in% training_groups)
      
      # Calculate factorial distributions
      catalog_factorial <- calculate_factorial_distribution(plots_df_binned)
      required_factorial <- calculate_factorial_distribution(required_plots, reference_df = plots_df_binned)
      
      results <- list(
        withheld_group_ids = required_groups,
        training_group_ids = training_groups,
        withheld_plots = required_plots,
        training_plots = training_plots,
        diagnostics = list(
          summary = tibble(
            metric = c("Groups", "Plots", "Trees"),
            n_withheld = c(length(required_groups), 
                          nrow(required_plots),
                          sum(required_plots$n_trees, na.rm = TRUE)),
            n_total = c(total_groups, total_plots, total_trees),
            pct_withheld = c((length(required_groups) / total_groups) * 100,
                            required_metrics$pct_plots,
                            required_metrics$pct_trees),
            distance_from_target = abs(pct_withheld - TARGET_PCT)
          ),
          distribution_comparison = create_distribution_comparison(catalog_distribution, 
                                                                  required_dist),
          summary_stats = create_summary_stats_comparison(plots_df, required_plots),
          factorial_distributions = list(
            catalog = catalog_factorial,
            withheld = required_factorial
          ),
          factorial_plots = create_factorial_plots(catalog_factorial, required_factorial),
          target_distance = abs(required_metrics$pct_plots - TARGET_PCT) + 
                           abs(required_metrics$pct_trees - TARGET_PCT),
          distribution_distance = dist_distance
        ),
        config = list(
          n_quantiles = n_quantiles,
          target_pct = TARGET_PCT,
          min_pct = MIN_PCT,
          max_pct = MAX_PCT,
          n_runs = 0,
          note = "Used required groups only (already in range)"
        )
      )
      
      cat("Selection complete (using required groups)!\n\n")
      return(results)
      
    } else if (required_metrics$pct_plots > MAX_PCT || 
               required_metrics$pct_trees > MAX_PCT) {
      
      # Required groups exceed maximum constraint
      cat(sprintf("  ⚠ Warning: Required groups exceed MAX_PCT (%.1f%% plots, %.1f%% trees)\n",
                  required_metrics$pct_plots, required_metrics$pct_trees))
      cat(sprintf("  Accepting as-is since these are required from higher hierarchy levels\n\n"))
      
      # Calculate distribution for reporting
      required_plots <- plots_df_binned |> filter(group_id %in% required_groups)
      required_dist <- calculate_distribution(required_plots)
      dist_distance <- overall_distribution_distance(required_dist, catalog_distribution)
      
      # Create results object
      training_groups <- setdiff(unique(plots_df$group_id), required_groups)
      training_plots <- plots_df_binned |> filter(group_id %in% training_groups)
      
      # Calculate factorial distributions
      catalog_factorial <- calculate_factorial_distribution(plots_df_binned)
      required_factorial <- calculate_factorial_distribution(required_plots, reference_df = plots_df_binned)
      
      results <- list(
        withheld_group_ids = required_groups,
        training_group_ids = training_groups,
        withheld_plots = required_plots,
        training_plots = training_plots,
        diagnostics = list(
          summary = tibble(
            metric = c("Groups", "Plots", "Trees"),
            n_withheld = c(length(required_groups), 
                          nrow(required_plots),
                          sum(required_plots$n_trees, na.rm = TRUE)),
            n_total = c(total_groups, total_plots, total_trees),
            pct_withheld = c((length(required_groups) / total_groups) * 100,
                            required_metrics$pct_plots,
                            required_metrics$pct_trees),
            distance_from_target = abs(pct_withheld - TARGET_PCT)
          ),
          distribution_comparison = create_distribution_comparison(catalog_distribution, 
                                                                  required_dist),
          summary_stats = create_summary_stats_comparison(plots_df, required_plots),
          factorial_distributions = list(
            catalog = catalog_factorial,
            withheld = required_factorial
          ),
          factorial_plots = create_factorial_plots(catalog_factorial, required_factorial),
          target_distance = abs(required_metrics$pct_plots - TARGET_PCT) + 
                           abs(required_metrics$pct_trees - TARGET_PCT),
          distribution_distance = dist_distance
        ),
        config = list(
          n_quantiles = n_quantiles,
          target_pct = TARGET_PCT,
          min_pct = MIN_PCT,
          max_pct = MAX_PCT,
          n_runs = 0,
          note = "Over-constrained by required groups"
        )
      )
      
      cat("Selection complete (over-constrained by required groups)!\n\n")
      return(results)
      
    } else {
      # Required groups below MIN_PCT - will proceed with algorithm
      cat(sprintf("  Required groups below target. Will proceed with %s algorithm\n", ALGORITHM))
      cat(sprintf("  (Required groups locked in, cannot be removed)\n\n"))
    }
  }
  
  # ============================================================================
  # STEP 2: Run algorithm (greedy, random, or hybrid)
  # ============================================================================
  
  if (ALGORITHM == "greedy") {
    
    cat(sprintf("Step 2: Running %d independent greedy removal runs...\n", N_RUNS))
    
    # Set up parallel processing
    if (PARALLELIZE) {
      n_cores <- availableCores()
      cat(sprintf("  Using parallel processing with %d cores\n", n_cores))
      plan(multisession, workers = n_cores)
    } else {
      cat("  Running sequentially (PARALLELIZE = FALSE)\n")
      plan(sequential)
    }
    
    # Run greedy removal in parallel
    all_states <- future_map(1:N_RUNS, function(run) {
      
      states <- greedy_removal(
        plots_df = plots_df_binned,
        catalog_distribution = catalog_distribution,
        catalog_factorial = catalog_factorial_ref,
        bin_info = bin_info,
        total_groups = total_groups,
        total_plots = total_plots,
        total_trees = total_trees,
        run_id = run,
        required_groups = required_groups,
        previous_selected_plots = previous_selected_plots
      )
      
      states$run <- run
      return(states)
    }, .options = furrr_options(seed = TRUE))
    
    # Reset to sequential after parallel work
    plan(sequential)
    
    cat(sprintf("  Completed all %d runs\n\n", N_RUNS))
    
    # Combine all states
    all_states_df <- bind_rows(all_states)
    
    # Filter to valid states
    valid_states <- all_states_df |>
      filter(
        pct_plots >= MIN_PCT & pct_plots <= MAX_PCT,
        pct_trees >= MIN_PCT & pct_trees <= MAX_PCT
      )
    
  } else if (ALGORITHM == "random") {
    
    cat("Step 2: Running random sampling...\n")
    
    # Run random sampling
    valid_states <- random_sampling(
      plots_df = plots_df_binned,
      catalog_distribution = catalog_distribution,
      catalog_factorial = catalog_factorial_ref,
      total_groups = total_groups,
      total_plots = total_plots,
      total_trees = total_trees,
      required_groups = required_groups,
      n_samples = N_RANDOM_SAMPLES,
      previous_selected_plots = previous_selected_plots
    )
    
    cat("\n")
    
  } else if (ALGORITHM == "hybrid") {
    
    cat(sprintf("Step 2: Running %d independent hybrid runs...\n", N_RUNS))
    
    # Run hybrid sampling
    all_states_df <- hybrid_sampling(
      plots_df = plots_df_binned,
      catalog_distribution = catalog_distribution,
      catalog_factorial = catalog_factorial_ref,
      bin_info = bin_info,
      total_groups = total_groups,
      total_plots = total_plots,
      total_trees = total_trees,
      required_groups = required_groups,
      n_runs = N_RUNS,
      previous_selected_plots = previous_selected_plots
    )
    
    # Filter to valid states (final 20% target)
    valid_states <- all_states_df |>
      filter(
        pct_plots >= MIN_PCT & pct_plots <= MAX_PCT,
        pct_trees >= MIN_PCT & pct_trees <= MAX_PCT
      )
    
  } else {
    stop(sprintf("Unknown algorithm: %s. Must be 'greedy', 'random', or 'hybrid'", ALGORITHM))
  }
  
  # ============================================================================
  # STEP 3: Find best solution
  # ============================================================================
  
  cat("Step 3: Selecting best solution...\n")
  
  cat(sprintf("  Found %d valid solutions\n", nrow(valid_states)))
  
  if (nrow(valid_states) == 0) {
    stop("No valid solutions found! Try adjusting MIN_PCT/MAX_PCT constraints.")
  }
  
  # Calculate target distance (if not already present from random sampling)
  if (!"target_distance" %in% names(valid_states)) {
    valid_states <- valid_states |>
      mutate(
        target_distance = abs(pct_plots - TARGET_PCT) +
                          abs(pct_trees - TARGET_PCT)
      )
  }
  
  # Select best solution
  best_solution <- valid_states |>
    arrange(dist_distance, target_distance) |>
    slice(1)
  
  withheld_groups <- best_solution$selected_groups[[1]]
  
  # Print solution info based on algorithm
  if (ALGORITHM == "greedy") {
    cat(sprintf("  Best solution from run %d, iteration %d\n", 
                best_solution$run, best_solution$iteration))
  } else if (ALGORITHM == "random") {
    cat(sprintf("  Best solution: sample %d\n", best_solution$sample_id))
  } else if (ALGORITHM == "hybrid") {
    cat(sprintf("  Best solution from run %d, iteration %d\n", 
                best_solution$run, best_solution$iteration))
    cat(sprintf("  Phase 1 start: %.1f%% plots, %.1f%% trees (dist=%.2f)\n",
                best_solution$phase1_pct_plots, 
                best_solution$phase1_pct_trees,
                best_solution$phase1_dist_distance))
  }
  cat(sprintf("  Target distance: %.2f\n", best_solution$target_distance))
  cat(sprintf("  Distribution distance: %.2f\n\n", best_solution$dist_distance))
  
  # ============================================================================
  # STEP 4: Extract withheld and training sets
  # ============================================================================
  
  cat("Step 4: Extracting final selection...\n")
  
  # Withheld set
  withheld_plots <- plots_df_binned |>
    filter(group_id %in% withheld_groups)
  
  # Training set
  training_groups <- setdiff(unique(plots_df$group_id), withheld_groups)
  training_plots <- plots_df_binned |>
    filter(group_id %in% training_groups)
  
  cat(sprintf("  Withheld: %d groups, %d plots\n",
              length(withheld_groups),
              nrow(withheld_plots)))
  cat(sprintf("  Training: %d groups, %d plots\n\n",
              length(training_groups),
              nrow(training_plots)))
  
  # ============================================================================
  # STEP 5: Generate diagnostics
  # ============================================================================
  
  cat("Step 5: Generating diagnostics...\n")
  
  diagnostics <- list()
  
  # Summary metrics
  diagnostics$summary <- tibble(
    metric = c("Groups", "Plots", "Trees"),
    n_withheld = c(
      length(withheld_groups),
      nrow(withheld_plots),
      sum(withheld_plots$n_trees, na.rm = TRUE)
    ),
    n_total = c(total_groups, total_plots, total_trees),
    pct_withheld = c(
      (length(withheld_groups) / total_groups) * 100,
      best_solution$pct_plots,
      best_solution$pct_trees
    ),
    distance_from_target = abs(pct_withheld - TARGET_PCT)
  )
  
  # Calculate withheld distribution
  withheld_distribution <- calculate_distribution(withheld_plots)
  
  # Calculate factorial distribution using appropriate reference
  if (use_reference_distribution) {
    # Use reference dataset for bin calculation to match reference catalog bins
    withheld_factorial <- calculate_factorial_distribution(withheld_plots, reference_df = reference_plots_df)
  } else {
    # Use current tier data for bin calculation
    withheld_factorial <- calculate_factorial_distribution(withheld_plots, reference_df = plots_df_binned)
  }
  
  # Generate diagnostics based on whether using reference distribution
  if (use_reference_distribution) {
    cat("  Generating dual diagnostics (reference + tier-specific)...\n")
    
    # Reference distribution comparisons (what we optimized for)
    diagnostics$distribution_comparison_reference <- create_distribution_comparison(
      catalog_distribution,
      withheld_distribution
    )
    
    # Tier-specific distribution comparisons (for comparison)
    diagnostics$distribution_comparison_tier <- create_distribution_comparison(
      tier_distribution,
      withheld_distribution
    )
    
    # Summary statistics (use reference if provided, otherwise tier)
    if (!is.null(reference_plots_df)) {
      diagnostics$summary_stats_reference <- create_summary_stats_comparison(
        reference_plots_df,
        withheld_plots
      )
    }
    diagnostics$summary_stats_tier <- create_summary_stats_comparison(
      plots_df,
      withheld_plots
    )
    
    # Factorial distributions
    diagnostics$factorial_distributions <- list(
      reference_catalog = catalog_factorial_ref,
      tier_catalog = tier_factorial,
      withheld = withheld_factorial
    )
    
    # Create dual factorial plots (side by side)
    diagnostics$factorial_plots <- create_dual_factorial_plots(
      catalog_factorial_ref,
      tier_factorial,
      withheld_factorial
    )
    
    # Overall scores (based on reference distribution)
    diagnostics$target_distance <- best_solution$target_distance
    diagnostics$distribution_distance_reference <- best_solution$dist_distance
    diagnostics$distribution_distance_tier <- overall_distribution_distance(withheld_distribution, tier_distribution)
    
  } else {
    cat("  Generating standard diagnostics (tier-specific only)...\n")
    
    # Standard distribution comparisons
    diagnostics$distribution_comparison <- create_distribution_comparison(
      catalog_distribution,
      withheld_distribution
    )
    
    # Summary statistics comparison
    diagnostics$summary_stats <- create_summary_stats_comparison(
      plots_df,
      withheld_plots
    )
    
    # Factorial distributions and plots
    catalog_factorial <- calculate_factorial_distribution(plots_df_binned)
    diagnostics$factorial_distributions <- list(
      catalog = catalog_factorial,
      withheld = withheld_factorial
    )
    diagnostics$factorial_plots <- create_factorial_plots(
      catalog_factorial,
      withheld_factorial
    )
    
    # Overall scores
    diagnostics$target_distance <- best_solution$target_distance
    diagnostics$distribution_distance <- best_solution$dist_distance
  }
  
  cat("  Diagnostics complete\n\n")
  
  # ============================================================================
  # Return results
  # ============================================================================
  
  results <- list(
    # Primary outputs
    withheld_group_ids = withheld_groups,
    training_group_ids = training_groups,
    
    # Detailed data frames
    withheld_plots = withheld_plots,
    training_plots = training_plots,
    
    # Diagnostics
    diagnostics = diagnostics,
    
    # Full state history (greedy and hybrid only)
    all_states = if (ALGORITHM %in% c("greedy", "hybrid")) all_states_df else NULL,
    best_solution = best_solution,
    
    # Configuration used
    config = list(
      algorithm = ALGORITHM,
      n_quantiles = n_quantiles,
      target_pct = TARGET_PCT,
      min_pct = MIN_PCT,
      max_pct = MAX_PCT,
      n_runs = if (ALGORITHM %in% c("greedy", "hybrid")) N_RUNS else NULL,
      n_samples = if (ALGORITHM == "random") N_RANDOM_SAMPLES else NULL,
      phase1_target_pct = if (ALGORITHM == "hybrid") PHASE1_TARGET_PCT else NULL,
      phase1_min_pct = if (ALGORITHM == "hybrid") PHASE1_MIN_PCT else NULL,
      phase1_max_pct = if (ALGORITHM == "hybrid") PHASE1_MAX_PCT else NULL,
      phase1_n_samples = if (ALGORITHM == "hybrid") N_RANDOM_SAMPLES_PHASE1 else NULL,
      factorial_weight = FACTORIAL_WEIGHT,
      use_reference_distribution = use_reference_distribution
    )
  )
  
  cat("Selection complete!\n\n")
  
  return(results)
}

# ==============================================================================
# REPORTING FUNCTION
# ==============================================================================

#' Print comprehensive selection report
#'
#' @param results Output from select_withheld_groups()
print_selection_report <- function(results) {
  
  cat("=" |> str_pad(width = 80, side = "both", pad = "="), "\n")
  cat("PROPORTIONAL STRATIFIED SAMPLING REPORT (SIMPLIFIED)\n")
  cat("=" |> str_pad(width = 80, side = "both", pad = "="), "\n\n")
  
  # Check if using reference distribution
  use_reference <- results$config$use_reference_distribution
  
  if (use_reference) {
    cat("MODE: Cross-tier stratification (matching reference distribution)\n\n")
  } else {
    cat("MODE: Within-tier stratification\n\n")
  }
  
  # Summary metrics
  cat("SUMMARY METRICS\n")
  cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
  print(results$diagnostics$summary, n = Inf)
  cat("\n")
  
  cat(sprintf("Target distance from 20%%: %.2f\n", 
              results$diagnostics$target_distance))
  
  if (use_reference) {
    cat(sprintf("Distribution distance (reference): %.2f\n", 
                results$diagnostics$distribution_distance_reference))
    cat(sprintf("Distribution distance (tier):      %.2f\n\n", 
                results$diagnostics$distribution_distance_tier))
  } else {
    cat(sprintf("Distribution distance: %.2f\n\n", 
                results$diagnostics$distribution_distance))
  }
  
  # Determine which distribution comparison to show
  dist_comp_key <- if (use_reference) "distribution_comparison_reference" else "distribution_comparison"
  summary_stats_key <- if (use_reference) "summary_stats_reference" else "summary_stats"
  
  # Distribution comparisons for continuous variables
  cat("CONTINUOUS VARIABLE DISTRIBUTIONS\n")
  if (use_reference) {
    cat("(Showing REFERENCE distribution comparison - what was optimized)\n")
  }
  cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
  for (var in CONTINUOUS_VARS) {
    cat(sprintf("\n%s:\n", toupper(var)))
    print(results$diagnostics[[dist_comp_key]][[var]] |> 
            select(bin, 
                   catalog_n_plots, catalog_pct_plots, 
                   catalog_n_trees, catalog_pct_trees,
                   selected_n_plots, selected_pct_plots,
                   selected_n_trees, selected_pct_trees,
                   diff_pct_plots, diff_pct_trees), 
          n = Inf)
    total_distance_plots <- sum(results$diagnostics[[dist_comp_key]][[var]]$abs_diff_plots)
    total_distance_trees <- sum(results$diagnostics[[dist_comp_key]][[var]]$abs_diff_trees)
    cat(sprintf("  Total absolute difference (plots): %.2f\n", total_distance_plots))
    cat(sprintf("  Total absolute difference (trees): %.2f\n", total_distance_trees))
  }
  cat("\n")
  
  # Distribution comparisons for categorical variables
  cat("CATEGORICAL VARIABLE DISTRIBUTIONS\n")
  if (use_reference) {
    cat("(Showing REFERENCE distribution comparison - what was optimized)\n")
  }
  cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
  for (var in CATEGORICAL_VARS) {
    cat(sprintf("\n%s:\n", toupper(var)))
    print(results$diagnostics[[dist_comp_key]][[var]] |>
            select(category,
                   catalog_n_plots, catalog_pct_plots,
                   catalog_n_trees, catalog_pct_trees,
                   selected_n_plots, selected_pct_plots,
                   selected_n_trees, selected_pct_trees,
                   diff_pct_plots, diff_pct_trees),
          n = Inf)
    total_distance_plots <- sum(results$diagnostics[[dist_comp_key]][[var]]$abs_diff_plots)
    total_distance_trees <- sum(results$diagnostics[[dist_comp_key]][[var]]$abs_diff_trees)
    cat(sprintf("  Total absolute difference (plots): %.2f\n", total_distance_plots))
    cat(sprintf("  Total absolute difference (trees): %.2f\n", total_distance_trees))
  }
  cat("\n")
  
  # Summary statistics
  cat("SUMMARY STATISTICS COMPARISON\n")
  if (use_reference) {
    cat("(Showing REFERENCE comparison)\n")
  }
  cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
  
  if (!is.null(results$diagnostics[[summary_stats_key]])) {
    print(results$diagnostics[[summary_stats_key]] |>
            select(variable, catalog_mean, selected_mean, mean_pct_diff, 
                   catalog_sd, selected_sd, sd_pct_diff),
          n = Inf)
  }
  cat("\n")
  
  # If using reference, also show tier-specific comparison
  if (use_reference) {
    cat("=" |> str_pad(width = 80, side = "both", pad = "="), "\n")
    cat("TIER-SPECIFIC COMPARISON (for reference only)\n")
    cat("=" |> str_pad(width = 80, side = "both", pad = "="), "\n\n")
    
    cat("CONTINUOUS VARIABLE DISTRIBUTIONS (TIER)\n")
    cat("(Showing how well selection captures within-tier distribution)\n")
    cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
    for (var in CONTINUOUS_VARS) {
      cat(sprintf("\n%s:\n", toupper(var)))
      print(results$diagnostics$distribution_comparison_tier[[var]] |> 
            select(bin, 
                   catalog_n_plots, catalog_pct_plots, 
                   catalog_n_trees, catalog_pct_trees,
                   selected_n_plots, selected_pct_plots,
                   selected_n_trees, selected_pct_trees,
                   diff_pct_plots, diff_pct_trees), 
          n = Inf)
      total_distance_plots <- sum(results$diagnostics$distribution_comparison_tier[[var]]$abs_diff_plots)
      total_distance_trees <- sum(results$diagnostics$distribution_comparison_tier[[var]]$abs_diff_trees)
      cat(sprintf("  Total absolute difference (plots): %.2f\n", total_distance_plots))
      cat(sprintf("  Total absolute difference (trees): %.2f\n", total_distance_trees))
    }
    cat("\n")
    
    cat("CATEGORICAL VARIABLE DISTRIBUTIONS (TIER)\n")
    cat("(Showing how well selection captures within-tier distribution)\n")
    cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
    for (var in CATEGORICAL_VARS) {
      cat(sprintf("\n%s:\n", toupper(var)))
      print(results$diagnostics$distribution_comparison_tier[[var]] |>
            select(category,
                   catalog_n_plots, catalog_pct_plots,
                   catalog_n_trees, catalog_pct_trees,
                   selected_n_plots, selected_pct_plots,
                   selected_n_trees, selected_pct_trees,
                   diff_pct_plots, diff_pct_trees),
          n = Inf)
      total_distance_plots <- sum(results$diagnostics$distribution_comparison_tier[[var]]$abs_diff_plots)
      total_distance_trees <- sum(results$diagnostics$distribution_comparison_tier[[var]]$abs_diff_trees)
      cat(sprintf("  Total absolute difference (plots): %.2f\n", total_distance_plots))
      cat(sprintf("  Total absolute difference (trees): %.2f\n", total_distance_trees))
    }
    cat("\n")
    
    # Tier-specific summary statistics
    cat("SUMMARY STATISTICS COMPARISON (TIER)\n")
    cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
    if (!is.null(results$diagnostics$summary_stats_tier)) {
      print(results$diagnostics$summary_stats_tier |>
              select(variable, catalog_mean, selected_mean, mean_pct_diff, 
                     catalog_sd, selected_sd, sd_pct_diff),
            n = Inf)
    }
    cat("\n")
  }
  
  cat("=" |> str_pad(width = 80, side = "both", pad = "="), "\n\n")
  
  # List withheld groups
  cat("Withheld group IDs:\n")
  cat(paste(results$withheld_group_ids, collapse = ", "))
  cat("\n\n")
  
  cat("Configuration used:\n")
  cat(sprintf("  Algorithm: %s\n", results$config$algorithm))
  cat(sprintf("  Quantiles: %d\n", results$config$n_quantiles))
  cat(sprintf("  Target: %d%%, Range: [%d%%, %d%%]\n", 
              results$config$target_pct, results$config$min_pct, 
              results$config$max_pct))
  if (!is.null(results$config$n_runs)) {
    cat(sprintf("  Runs: %d\n", results$config$n_runs))
  }
  cat(sprintf("  Reference distribution: %s\n", 
              ifelse(use_reference, "Yes (cross-tier)", "No (within-tier)")))
  cat("\n")
}

#' Save factorial distribution plots to files
#'
#' @param results Output from select_withheld_groups()
#' @param output_dir Directory to save plots (default: current directory)
#' @param width Plot width in inches (default: 8)
#' @param height Plot height in inches (default: 6)
save_factorial_plots <- function(results, output_dir = ".", width = 8, height = 6) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("Saving factorial distribution plots...\n")
  
  for (plot_name in names(results$diagnostics$factorial_plots)) {
    plot_obj <- results$diagnostics$factorial_plots[[plot_name]]
    filename <- file.path(output_dir, paste0(plot_name, ".png"))
    
    ggsave(
      filename = filename,
      plot = plot_obj,
      width = width,
      height = height,
      dpi = 300
    )
    
    cat(sprintf("  Saved: %s\n", filename))
  }
  
  cat("All plots saved!\n")
}

# ==============================================================================
# HELPER: PREPARE REFERENCE DISTRIBUTIONS
# ==============================================================================

#' Prepare reference distributions for cross-tier stratification
#'
#' This function prepares the reference distributions needed for hierarchical
#' stratification across multiple tiers while maintaining overall distribution.
#'
#' @param all_plots_df Data frame containing ALL plots across all tiers
#' @return List with reference_distribution, reference_factorial, and binned data
#' @export
prepare_reference_distributions <- function(all_plots_df) {
  
  cat("Preparing reference distributions from full dataset...\n")
  cat(sprintf("  Total plots: %d\n", nrow(all_plots_df)))
  
  # Calculate optimal number of quantiles
  n_quantiles <- calculate_optimal_n_quantiles(nrow(all_plots_df))
  cat(sprintf("  Using %d quantiles for continuous variables\n", n_quantiles))
  
  # Prepare quantile bins
  bin_info <- prepare_quantile_bins(all_plots_df, n_quantiles)
  all_plots_binned <- assign_quantile_bins(all_plots_df, bin_info)
  
  # Calculate distributions
  reference_distribution <- calculate_distribution(all_plots_binned)
  reference_factorial <- calculate_factorial_distribution(all_plots_binned)
  
  cat("  Reference distributions calculated\n\n")
  
  return(list(
    reference_distribution = reference_distribution,
    reference_factorial = reference_factorial,
    reference_plots_df = all_plots_binned,
    bin_info = bin_info,
    n_quantiles = n_quantiles
  ))
}

# ==============================================================================
# HELPER: EVALUATE COMBINED SELECTION ACROSS TIERS
# ==============================================================================

#' Evaluate combined selection across multiple tiers
#'
#' This function generates a comprehensive distribution comparison report for
#' a combined set of selected groups across all tiers, evaluated against the
#' reference distribution.
#'
#' @param all_plots_df Data frame containing ALL plots across all tiers
#' @param selected_group_ids Vector of all selected group IDs (across all tiers)
#' @param reference_distributions List from prepare_reference_distributions() (optional, will calculate if not provided)
#' @return List with diagnostics comparing combined selection to reference
#' @export
evaluate_combined_selection <- function(all_plots_df, selected_group_ids, 
                                        reference_distributions = NULL) {
  
  cat("=" |> str_pad(width = 80, side = "both", pad = "="), "\n")
  cat("COMBINED SELECTION EVALUATION (ACROSS ALL TIERS)\n")
  cat("=" |> str_pad(width = 80, side = "both", pad = "="), "\n\n")
  
  # Prepare reference distributions if not provided
  if (is.null(reference_distributions)) {
    cat("Calculating reference distributions...\n")
    reference_distributions <- prepare_reference_distributions(all_plots_df)
  }
  
  # Extract reference components
  reference_distribution <- reference_distributions$reference_distribution
  reference_factorial <- reference_distributions$reference_factorial
  reference_plots_df <- reference_distributions$reference_plots_df
  
  # Calculate totals
  total_groups <- n_distinct(all_plots_df$group_id)
  total_plots <- nrow(all_plots_df)
  total_trees <- sum(all_plots_df$n_trees, na.rm = TRUE)
  
  cat(sprintf("Total catalog: %d groups, %d plots, %d trees\n",
              total_groups, total_plots, total_trees))
  cat(sprintf("Selected: %d groups\n\n", length(selected_group_ids)))
  
  # Filter to selected plots
  selected_plots <- reference_plots_df |>
    filter(group_id %in% selected_group_ids)
  
  # Calculate metrics
  n_selected_plots <- nrow(selected_plots)
  n_selected_trees <- sum(selected_plots$n_trees, na.rm = TRUE)
  pct_plots <- 100 * n_selected_plots / total_plots
  pct_trees <- 100 * n_selected_trees / total_trees
  
  cat("SUMMARY METRICS\n")
  cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
  summary_table <- tibble(
    metric = c("Groups", "Plots", "Trees"),
    n_selected = c(length(selected_group_ids), n_selected_plots, n_selected_trees),
    n_total = c(total_groups, total_plots, total_trees),
    pct_selected = c(100 * length(selected_group_ids) / total_groups,
                     pct_plots, pct_trees),
    distance_from_target = abs(pct_selected - TARGET_PCT)
  )
  print(summary_table, n = Inf)
  cat("\n")
  
  # Calculate distributions
  selected_distribution <- calculate_distribution(selected_plots)
  selected_factorial <- calculate_factorial_distribution(selected_plots, 
                                                         reference_df = reference_plots_df)
  
  # Distribution distance
  dist_distance <- combined_distribution_distance(
    selected_distribution, 
    reference_distribution,
    selected_factorial, 
    reference_factorial,
    factorial_weight = FACTORIAL_WEIGHT
  )
  
  cat(sprintf("Distribution distance: %.2f\n", dist_distance))
  cat(sprintf("Target distance from %.0f%%: %.2f\n\n", 
              TARGET_PCT, abs(pct_plots - TARGET_PCT) + abs(pct_trees - TARGET_PCT)))
  
  # Create detailed comparisons
  cat("Generating detailed distribution comparisons...\n")
  
  distribution_comparison <- create_distribution_comparison(
    reference_distribution,
    selected_distribution
  )
  
  summary_stats <- create_summary_stats_comparison(
    all_plots_df,
    selected_plots
  )
  
  factorial_plots <- create_factorial_plots(
    reference_factorial,
    selected_factorial
  )
  
  cat("Evaluation complete!\n\n")
  
  # Return results
  results <- list(
    # Summary
    summary = summary_table,
    n_groups = length(selected_group_ids),
    n_plots = n_selected_plots,
    n_trees = n_selected_trees,
    pct_plots = pct_plots,
    pct_trees = pct_trees,
    
    # Group IDs
    selected_group_ids = selected_group_ids,
    
    # Data frames
    selected_plots = selected_plots,
    reference_plots = reference_plots_df,
    
    # Diagnostics
    diagnostics = list(
      summary = summary_table,
      distribution_comparison = distribution_comparison,
      summary_stats = summary_stats,
      factorial_distributions = list(
        catalog = reference_factorial,
        withheld = selected_factorial
      ),
      factorial_plots = factorial_plots,
      target_distance = abs(pct_plots - TARGET_PCT) + abs(pct_trees - TARGET_PCT),
      distribution_distance = dist_distance
    ),
    
    # Configuration
    config = list(
      target_pct = TARGET_PCT,
      factorial_weight = FACTORIAL_WEIGHT,
      evaluation_type = "combined_across_tiers"
    )
  )
  
  return(results)
}

#' Print evaluation report for combined selection
#'
#' @param eval_results Output from evaluate_combined_selection()
#' @export
print_combined_evaluation_report <- function(eval_results) {
  
  cat("=" |> str_pad(width = 80, side = "both", pad = "="), "\n")
  cat("COMBINED SELECTION EVALUATION REPORT\n")
  cat("=" |> str_pad(width = 80, side = "both", pad = "="), "\n\n")
  
  # Summary metrics
  cat("SUMMARY METRICS\n")
  cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
  print(eval_results$diagnostics$summary, n = Inf)
  cat("\n")
  
  cat(sprintf("Target distance from %.0f%%: %.2f\n", 
              eval_results$config$target_pct,
              eval_results$diagnostics$target_distance))
  cat(sprintf("Distribution distance: %.2f\n\n", 
              eval_results$diagnostics$distribution_distance))
  
  # Distribution comparisons for continuous variables
  cat("CONTINUOUS VARIABLE DISTRIBUTIONS\n")
  cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
  for (var in CONTINUOUS_VARS) {
    cat(sprintf("\n%s:\n", toupper(var)))
    print(eval_results$diagnostics$distribution_comparison[[var]] |> 
            select(bin, 
                   catalog_n_plots, catalog_pct_plots, 
                   catalog_n_trees, catalog_pct_trees,
                   selected_n_plots, selected_pct_plots,
                   selected_n_trees, selected_pct_trees,
                   diff_pct_plots, diff_pct_trees), 
          n = Inf)
    total_distance_plots <- sum(eval_results$diagnostics$distribution_comparison[[var]]$abs_diff_plots)
    total_distance_trees <- sum(eval_results$diagnostics$distribution_comparison[[var]]$abs_diff_trees)
    cat(sprintf("  Total absolute difference (plots): %.2f\n", total_distance_plots))
    cat(sprintf("  Total absolute difference (trees): %.2f\n", total_distance_trees))
  }
  cat("\n")
  
  # Distribution comparisons for categorical variables
  cat("CATEGORICAL VARIABLE DISTRIBUTIONS\n")
  cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
  for (var in CATEGORICAL_VARS) {
    cat(sprintf("\n%s:\n", toupper(var)))
    print(eval_results$diagnostics$distribution_comparison[[var]] |>
            select(category,
                   catalog_n_plots, catalog_pct_plots,
                   catalog_n_trees, catalog_pct_trees,
                   selected_n_plots, selected_pct_plots,
                   selected_n_trees, selected_pct_trees,
                   diff_pct_plots, diff_pct_trees),
          n = Inf)
    total_distance_plots <- sum(eval_results$diagnostics$distribution_comparison[[var]]$abs_diff_plots)
    total_distance_trees <- sum(eval_results$diagnostics$distribution_comparison[[var]]$abs_diff_trees)
    cat(sprintf("  Total absolute difference (plots): %.2f\n", total_distance_plots))
    cat(sprintf("  Total absolute difference (trees): %.2f\n", total_distance_trees))
  }
  cat("\n")
  
  # Summary statistics
  cat("SUMMARY STATISTICS COMPARISON\n")
  cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
  print(eval_results$diagnostics$summary_stats |>
          select(variable, catalog_mean, selected_mean, mean_pct_diff, 
                 catalog_sd, selected_sd, sd_pct_diff),
        n = Inf)
  cat("\n")
  
  cat("=" |> str_pad(width = 80, side = "both", pad = "="), "\n\n")
  
  # List selected groups
  cat(sprintf("Selected %d group IDs:\n", length(eval_results$selected_group_ids)))
  cat(paste(eval_results$selected_group_ids, collapse = ", "))
  cat("\n\n")
}

# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================

# Uncomment and modify to run with your data:

# # ============================================================================
# # EXAMPLE 1: Standard within-tier stratification
# # ============================================================================
# plots <- read_csv("path/to/plots.csv")
# # Required columns: plot_id, group_id, temp, ppt, ecoregion, mean_ba_live,
# #                   trees_per_ha, pct_dead, sp_comp_group, n_trees
# 
# # Run selection
# results <- select_withheld_groups(plots_df = plots)
# 
# # Print report
# print_selection_report(results)
# 
# # Save results
# write_csv(results$withheld_plots, "withheld_plots.csv")
# write_csv(results$training_plots, "training_plots.csv")
# write_csv(results$diagnostics$summary, "selection_summary.csv")
# 
# # Save factorial distribution plots
# save_factorial_plots(results, output_dir = "factorial_plots")
# 
# # Or view individual plots interactively
# print(results$diagnostics$factorial_plots$ecoregion_x_sp_comp_group_plots)
# print(results$diagnostics$factorial_plots$ecoregion_x_sp_comp_group_trees)
#
# # ============================================================================
# # EXAMPLE 2: Cross-tier stratification (hierarchical with shared reference)
# # ============================================================================
# # Load full dataset
# all_data <- read_csv("path/to/all_plots.csv")
# 
# # Prepare reference distributions from ALL data
# ref <- prepare_reference_distributions(all_data)
# 
# # Tier 1: aligned_paired_drone_footprint
# tier1_plots <- all_data |> filter(pairing_tier == "aligned_paired_drone_footprint")
# 
# res1 <- select_withheld_groups(
#   plots_df = tier1_plots,
#   required_groups = c(),
#   reference_distribution = ref$reference_distribution,
#   reference_factorial = ref$reference_factorial,
#   reference_plots_df = ref$reference_plots_df
# )
# 
# groups_withheld1 <- res1$withheld_plots$group_id
# print_selection_report(res1)
# 
# # Tier 2: paired_drone_footprint
# tier2_plots <- all_data |> filter(pairing_tier == "paired_drone_footprint")
# 
# res2 <- select_withheld_groups(
#   plots_df = tier2_plots,
#   required_groups = groups_withheld1,  # Include tier 1 selections
#   reference_distribution = ref$reference_distribution,
#   reference_factorial = ref$reference_factorial,
#   reference_plots_df = ref$reference_plots_df
# )
# 
# groups_withheld2 <- unique(c(groups_withheld1, res2$withheld_plots$group_id))
# print_selection_report(res2)
# 
# # View dual plots showing both reference and tier-specific matching
# print(res2$diagnostics$factorial_plots$mean_ba_live_x_sp_comp_group_plots)
# 
# # Continue for additional tiers...

