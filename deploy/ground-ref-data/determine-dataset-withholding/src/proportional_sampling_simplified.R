# ==============================================================================
# PROPORTIONAL STRATIFIED SAMPLING: SIMPLIFIED VERSION (SINGLE DATA FRAME)
# ==============================================================================
# Purpose: Select ~20% of drone footprint groups to withhold for validation
#          such that withheld plots match the catalog's distribution
#
# Method: Greedy removal - start with all groups, remove those that minimally
#         disturb the distribution (normalized per plot)
#
# Input: Single data frame with plots and their group_id assignments
#
# Author: [Your Name]
# Date: 2025-10-24
# ==============================================================================

library(tidyverse)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Set random seed for reproducibility
set.seed(42)

# Define factorial combinations to track (optional, for diagnostics)
FACTORIAL_COMBINATIONS <- list(
  list(var1 = "ecoregion", var2 = "sp_comp_group"),
  list(var1 = "mean_ba_live", var2 = "sp_comp_group"),
  list(var1 = "trees_per_ha", var2 = "ecoregion")
)

# Continuous variables to stratify
CONTINUOUS_VARS <- c("ppt", "trees_per_ha", "mean_ba_live", "area_ha")

# Categorical variables to stratify
CATEGORICAL_VARS <- c("ecoregion", "sp_comp_group")

# Binning parameters
TARGET_PLOTS_PER_BIN <- 5   # Target average plots per bin in full catalog
MAX_BINS <- 5                # Maximum number of quantile bins to use

# Selection parameters
TARGET_PCT <- 20             # Target percentage to withhold
MIN_PCT <- 15                # Minimum acceptable percentage
MAX_PCT <- 25                # Maximum acceptable percentage

# Stochastic selection parameters
N_RUNS <- 30                 # Number of independent runs
TOP_K_CANDIDATES <- 3        # Consider top K groups when selecting which to remove
STOCHASTIC_TEMP <- 2.0       # Temperature for probability weighting

# ==============================================================================
# HELPER FUNCTIONS: QUANTILE CALCULATION
# ==============================================================================

#' Calculate optimal number of quantiles for a dataset
#'
#' @param n_plots_total Total number of plots in catalog
#' @param selection_pct Percentage to be selected (default 20%)
#' @return Optimal number of quantiles (between 2 and MAX_BINS)
calculate_optimal_n_quantiles <- function(n_plots_total, 
                                          selection_pct = TARGET_PCT) {
  
  n_plots_selected <- n_plots_total * (selection_pct / 100)
  max_quantiles_allowed <- floor(n_plots_selected / TARGET_PLOTS_PER_BIN)
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

#' Calculate distribution of plots across bins for all variables
#'
#' @param plots_df Data frame of plots (with bin columns)
#' @return List of distributions for each variable
calculate_distribution <- function(plots_df) {
  
  distribution <- list()
  
  # Continuous variables (using bins)
  for (var in CONTINUOUS_VARS) {
    bin_col <- paste0(var, "_bin")
    
    dist <- plots_df |>
      filter(!is.na(.data[[bin_col]])) |>
      count(.data[[bin_col]], name = "n") |>
      mutate(pct = n / sum(n) * 100) |>
      rename(bin = !!bin_col)
    
    distribution[[var]] <- dist
  }
  
  # Categorical variables
  for (var in CATEGORICAL_VARS) {
    dist <- plots_df |>
      filter(!is.na(.data[[var]])) |>
      count(.data[[var]], name = "n") |>
      mutate(pct = n / sum(n) * 100) |>
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
      pct_1 = replace_na(pct_1, 0),
      pct_2 = replace_na(pct_2, 0),
      abs_diff = abs(pct_1 - pct_2)
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
      n_groups = 0,
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
  
  plot_ok <- metrics$pct_plots >= MIN_PCT & metrics$pct_plots <= MAX_PCT
  tree_ok <- metrics$pct_trees >= MIN_PCT & metrics$pct_trees <= MAX_PCT
  
  return(plot_ok & tree_ok)
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
# CORE ALGORITHM: GREEDY REMOVAL
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
#' @return Data frame of states at each iteration
greedy_removal <- function(plots_df, catalog_distribution, bin_info,
                          total_groups, total_plots, total_trees, run_id,
                          required_groups = c()) {
  
  # Set run-specific seed
  set.seed(42 + run_id)
  
  # Initialize: ALL groups selected
  all_groups <- unique(plots_df$group_id)
  selected_groups <- all_groups
  
  # Identify removable groups (all except required)
  removable_groups <- setdiff(all_groups, required_groups)
  
  # Pre-calculate number of plots in each group
  group_plot_counts <- plots_df |>
    count(group_id, name = "n_plots")
  
  # Storage for state history
  state_history <- list()
  iteration <- 1
  
  # Calculate and save initial state (100% selected)
  metrics <- calculate_metrics(selected_groups, plots_df,
                               total_groups, total_plots, total_trees)
  selected_plots <- plots_df |> filter(group_id %in% selected_groups)
  selected_dist <- calculate_distribution(selected_plots)
  dist_distance <- overall_distribution_distance(selected_dist, catalog_distribution)
  
  state_history[[iteration]] <- list(
    iteration = iteration,
    selected_groups = selected_groups,
    metrics = metrics,
    dist_distance = dist_distance
  )
  
  # Greedy removal loop
  while (metrics$pct_plots > MAX_PCT || metrics$pct_trees > MAX_PCT) {
    
    # Evaluate each currently-selected REMOVABLE group (not required groups!)
    candidate_scores <- tibble(
      group_id = numeric(),
      dist_distance = numeric(),
      per_plot_distortion = numeric(),
      n_plots = numeric(),
      test_pct_plots = numeric(),
      test_pct_trees = numeric()
    )
    
    for (candidate_group in removable_groups) {
      
      # Test removing this group
      test_groups <- setdiff(selected_groups, candidate_group)
      test_metrics <- calculate_metrics(test_groups, plots_df,
                                       total_groups, total_plots, total_trees)
      
      # Check if removal would violate minimum constraints
      if (test_metrics$pct_plots < MIN_PCT || test_metrics$pct_trees < MIN_PCT) {
        next  # Skip - can't remove without violating constraints
      }
      
      # Calculate distribution distance if we remove this group
      test_plots <- plots_df |> filter(group_id %in% test_groups)
      test_dist <- calculate_distribution(test_plots)
      test_dist_distance <- overall_distribution_distance(test_dist, 
                                                         catalog_distribution)
      
      # Normalize by number of plots in this group
      n_plots_in_group <- group_plot_counts |>
        filter(group_id == candidate_group) |>
        pull(n_plots)
      
      per_plot_distortion <- test_dist_distance / n_plots_in_group
      
      # Store candidate
      candidate_scores <- candidate_scores |>
        add_row(
          group_id = candidate_group,
          dist_distance = test_dist_distance,
          per_plot_distortion = per_plot_distortion,
          n_plots = n_plots_in_group,
          test_pct_plots = test_metrics$pct_plots,
          test_pct_trees = test_metrics$pct_trees
        )
    }
    
    # Check if we have any valid candidates
    if (nrow(candidate_scores) == 0) {
      break  # No valid removals possible
    }
    
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
    selected_dist <- calculate_distribution(selected_plots)
    dist_distance <- overall_distribution_distance(selected_dist, 
                                                  catalog_distribution)
    
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
      rename(catalog_pct = pct, catalog_n = n) |>
      full_join(
        selected_dist[[var]] |>
          rename(selected_pct = pct, selected_n = n),
        by = "bin"
      ) |>
      mutate(
        catalog_pct = replace_na(catalog_pct, 0),
        selected_pct = replace_na(selected_pct, 0),
        difference = selected_pct - catalog_pct,
        abs_difference = abs(difference)
      ) |>
      arrange(bin)
    
    comparisons[[var]] <- comparison
  }
  
  # Categorical variables
  for (var in CATEGORICAL_VARS) {
    comparison <- catalog_dist[[var]] |>
      rename(catalog_pct = pct, catalog_n = n) |>
      full_join(
        selected_dist[[var]] |>
          rename(selected_pct = pct, selected_n = n),
        by = "category"
      ) |>
      mutate(
        catalog_pct = replace_na(catalog_pct, 0),
        selected_pct = replace_na(selected_pct, 0),
        difference = selected_pct - catalog_pct,
        abs_difference = abs(difference)
      ) |>
      arrange(desc(catalog_pct))
    
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

# ==============================================================================
# MAIN EXECUTION FUNCTION
# ==============================================================================

#' Main function to select withheld groups using greedy removal
#'
#' @param plots_df Data frame with columns: plot_id, group_id, attributes, n_trees
#' @param required_groups Vector of group_ids that must be in the withheld set (for hierarchical stratification)
#' @return List containing withheld and training plots with diagnostics
select_withheld_groups <- function(plots_df, required_groups = c()) {
  
  cat("Starting proportional stratified group selection (SIMPLIFIED)...\n")
  cat("Method: Greedy removal with stochastic selection\n\n")
  
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
  
  # Calculate optimal number of quantiles
  n_quantiles <- calculate_optimal_n_quantiles(total_plots, TARGET_PCT)
  cat(sprintf("  Using %d quantiles for continuous variables\n", n_quantiles))
  cat(sprintf("  Target: ≥%d plots per bin in selected set\n", TARGET_PLOTS_PER_BIN))
  
  # Prepare quantile bins
  bin_info <- prepare_quantile_bins(plots_df, n_quantiles)
  plots_df_binned <- assign_quantile_bins(plots_df, bin_info)
  
  # Calculate catalog distribution (baseline)
  catalog_distribution <- calculate_distribution(plots_df_binned)
  
  cat("  Catalog distribution calculated\n\n")
  
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
      # Required groups below MIN_PCT - will proceed with greedy removal
      cat(sprintf("  Required groups below target. Will proceed with greedy removal\n"))
      cat(sprintf("  (Required groups locked in, cannot be removed)\n\n"))
    }
  }
  
  # ============================================================================
  # STEP 2: Run greedy removal with multiple runs
  # ============================================================================
  
  cat(sprintf("Step 2: Running %d independent greedy removal runs...\n", N_RUNS))
  
  all_states <- list()
  
  for (run in 1:N_RUNS) {
    
    cat(sprintf("  Beginning %d/%d runs...\n", run, N_RUNS))
    
    states <- greedy_removal(
      plots_df = plots_df_binned,
      catalog_distribution = catalog_distribution,
      bin_info = bin_info,
      total_groups = total_groups,
      total_plots = total_plots,
      total_trees = total_trees,
      run_id = run,
      required_groups = required_groups
    )
    
    states$run <- run
    all_states[[run]] <- states
  }
  
  cat(sprintf("  Completed all %d runs\n\n", N_RUNS))
  
  # Combine all states
  all_states_df <- bind_rows(all_states)
  
  # ============================================================================
  # STEP 3: Find best solution
  # ============================================================================
  
  cat("Step 3: Selecting best solution...\n")
  
  # Filter to valid states
  valid_states <- all_states_df |>
    filter(
      pct_plots >= MIN_PCT & pct_plots <= MAX_PCT,
      pct_trees >= MIN_PCT & pct_trees <= MAX_PCT
    )
  
  cat(sprintf("  Found %d valid solutions across all runs\n", nrow(valid_states)))
  
  if (nrow(valid_states) == 0) {
    stop("No valid solutions found! Try adjusting MIN_PCT/MAX_PCT constraints.")
  }
  
  # Calculate target distance
  valid_states <- valid_states |>
    mutate(
      target_distance = abs(pct_plots - TARGET_PCT) +
                        abs(pct_trees - TARGET_PCT)
    )
  
  # Select best solution
  best_solution <- valid_states |>
    arrange(target_distance, dist_distance) |>
    slice(1)
  
  withheld_groups <- best_solution$selected_groups[[1]]
  
  cat(sprintf("  Best solution from run %d, iteration %d\n", 
              best_solution$run, best_solution$iteration))
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
  
  cat("Step 5: Generating diagnostics...\n\n")
  
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
  
  # Distribution comparisons
  withheld_distribution <- calculate_distribution(withheld_plots)
  diagnostics$distribution_comparison <- create_distribution_comparison(
    catalog_distribution, 
    withheld_distribution
  )
  
  # Summary statistics comparison
  diagnostics$summary_stats <- create_summary_stats_comparison(
    plots_df, 
    withheld_plots
  )
  
  # Overall scores
  diagnostics$target_distance <- best_solution$target_distance
  diagnostics$distribution_distance <- best_solution$dist_distance
  
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
    
    # Full state history
    all_states = all_states_df,
    best_solution = best_solution,
    
    # Configuration used
    config = list(
      n_quantiles = n_quantiles,
      target_pct = TARGET_PCT,
      min_pct = MIN_PCT,
      max_pct = MAX_PCT,
      n_runs = N_RUNS
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
  
  # Summary metrics
  cat("SUMMARY METRICS\n")
  cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
  print(results$diagnostics$summary, n = Inf)
  cat("\n")
  
  cat(sprintf("Target distance from 20%%: %.2f\n", 
              results$diagnostics$target_distance))
  cat(sprintf("Distribution distance: %.2f\n\n", 
              results$diagnostics$distribution_distance))
  
  # Distribution comparisons for continuous variables
  cat("CONTINUOUS VARIABLE DISTRIBUTIONS\n")
  cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
  for (var in CONTINUOUS_VARS) {
    cat(sprintf("\n%s:\n", toupper(var)))
    print(results$diagnostics$distribution_comparison[[var]] |> 
            select(bin, catalog_pct, selected_pct, difference), 
          n = Inf)
    total_distance <- sum(results$diagnostics$distribution_comparison[[var]]$abs_difference)
    cat(sprintf("  Total absolute difference: %.2f\n", total_distance))
  }
  cat("\n")
  
  # Distribution comparisons for categorical variables
  cat("CATEGORICAL VARIABLE DISTRIBUTIONS\n")
  cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
  for (var in CATEGORICAL_VARS) {
    cat(sprintf("\n%s:\n", toupper(var)))
    print(results$diagnostics$distribution_comparison[[var]] |>
            select(category, catalog_pct, selected_pct, difference),
          n = Inf)
    total_distance <- sum(results$diagnostics$distribution_comparison[[var]]$abs_difference)
    cat(sprintf("  Total absolute difference: %.2f\n", total_distance))
  }
  cat("\n")
  
  # Summary statistics
  cat("SUMMARY STATISTICS COMPARISON\n")
  cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
  print(results$diagnostics$summary_stats |>
          select(variable, catalog_mean, selected_mean, mean_pct_diff, 
                catalog_sd, selected_sd, sd_pct_diff),
        n = Inf)
  cat("\n")
  
  cat("=" |> str_pad(width = 80, side = "both", pad = "="), "\n\n")
  
  # List withheld groups
  cat("Withheld group IDs:\n")
  cat(paste(results$withheld_group_ids, collapse = ", "))
  cat("\n\n")
  
  cat("Configuration used:\n")
  cat(sprintf("  Quantiles: %d\n", results$config$n_quantiles))
  cat(sprintf("  Target: %d%%, Range: [%d%%, %d%%]\n", 
              results$config$target_pct, results$config$min_pct, 
              results$config$max_pct))
  cat(sprintf("  Runs: %d\n\n", results$config$n_runs))
}

# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================

# Uncomment and modify to run with your data:

# # Load your data (single data frame!)
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
