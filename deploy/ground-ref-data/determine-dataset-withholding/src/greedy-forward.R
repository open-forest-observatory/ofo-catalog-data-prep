# ==============================================================================
# PROPORTIONAL STRATIFIED SAMPLING: GREEDY FORWARD SELECTION
# ==============================================================================
# Purpose: Select ~20% of drone footprint groups to withhold for validation
#          such that withheld plots match the catalog's distribution
#
# Method: Greedy forward selection - iteratively add groups that bring the
#         selection closest to matching the target distribution
#
# Input: Single data frame with plots and their group_id assignments
#
# Author: [Your Name]
# Date: 2025-11-04
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

# Greedy algorithm parameters
N_RUNS <- 4                  # Number of independent runs
TOP_K_CANDIDATES <- 5        # Consider top K groups when selecting which to add
STOCHASTIC_TEMP <- 3         # Temperature for probability weighting

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
    return(0)
  }
  
  distances <- map_dbl(names(factorial_catalog), function(combo_name) {
    
    catalog_dist <- factorial_catalog[[combo_name]]$distribution
    selected_dist <- factorial_selected[[combo_name]]$distribution
    
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
  
  single_var_distance <- overall_distribution_distance(dist_selected, dist_catalog)
  
  if (is.null(factorial_selected) || is.null(factorial_catalog)) {
    return(single_var_distance)
  }
  
  factorial_distance <- factorial_distribution_distance(factorial_selected, factorial_catalog)
  
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
  
  selected_plots <- plots_df |>
    filter(group_id %in% selected_groups)
  
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

# ==============================================================================
# ALGORITHM: GREEDY FORWARD SELECTION
# ==============================================================================

#' Run one iteration of greedy forward selection
#'
#' @param plots_df Data frame of plots (with bin columns and group_id)
#' @param catalog_distribution Distribution of full catalog
#' @param catalog_factorial Factorial distribution of catalog
#' @param total_groups Total number of groups
#' @param total_plots Total number of plots
#' @param total_trees Total number of trees
#' @param run_id Run identifier for random seed
#' @param required_groups Vector of group_ids that must be included
#' @param previous_selected_plots Data frame of plots already selected in previous tiers (optional)
#' @return Data frame of states at each iteration
greedy_forward_selection <- function(plots_df, catalog_distribution, catalog_factorial,
                                    total_groups, total_plots, total_trees, run_id,
                                    required_groups = c(),
                                    previous_selected_plots = NULL) {
  
  # Set run-specific seed
  set.seed(42 + run_id)
  
  # Check if we're filling gaps from previous selections
  filling_gaps <- !is.null(previous_selected_plots)
  
  # Initialize: start with required groups
  selected_groups <- required_groups
  all_groups <- unique(plots_df$group_id)
  available_groups <- setdiff(all_groups, required_groups)
  
  # Pre-calculate group statistics for faster lookup
  group_stats <- plots_df |>
    group_by(group_id) |>
    summarise(
      n_plots = n(),
      n_trees = sum(n_trees, na.rm = TRUE),
      .groups = "drop"
    )
  
  group_n_plots <- setNames(group_stats$n_plots, group_stats$group_id)
  group_n_trees <- setNames(group_stats$n_trees, group_stats$group_id)
  
  # Storage for state history
  state_history <- list()
  iteration <- 1
  
  # Calculate and save initial state
  metrics <- calculate_metrics(selected_groups, plots_df,
                               total_groups, total_plots, total_trees)
  
  if (length(selected_groups) > 0) {
    selected_plots <- plots_df |> filter(group_id %in% selected_groups)
    
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
  } else {
    dist_distance <- Inf  # No groups selected yet
  }
  
  state_history[[iteration]] <- list(
    iteration = iteration,
    selected_groups = selected_groups,
    metrics = metrics,
    dist_distance = dist_distance
  )
  
  # Greedy forward selection loop
  while (metrics$pct_groups < MAX_GROUPS_PCT && 
         metrics$pct_plots < MAX_PCT && 
         metrics$pct_trees < MAX_PCT &&
         length(available_groups) > 0) {
    
    # Pre-allocate candidate_scores
    n_available <- length(available_groups)
    candidate_scores <- tibble(
      group_id = rep(NA_real_, n_available),
      dist_distance = rep(NA_real_, n_available),
      n_plots = rep(NA_integer_, n_available),
      test_pct_plots = rep(NA_real_, n_available),
      test_pct_trees = rep(NA_real_, n_available)
    )
    
    valid_idx <- 0
    
    for (candidate_group in available_groups) {
      
      # Fast lookup of group stats
      candidate_n_plots <- group_n_plots[as.character(candidate_group)]
      candidate_n_trees <- group_n_trees[as.character(candidate_group)]
      
      # Quick check: Would adding this group violate maximum constraints?
      test_n_groups <- length(selected_groups) + 1
      test_n_plots <- metrics$n_plots + candidate_n_plots
      test_n_trees <- metrics$n_trees + candidate_n_trees
      test_pct_groups <- (test_n_groups / total_groups) * 100
      test_pct_plots <- (test_n_plots / total_plots) * 100
      test_pct_trees <- (test_n_trees / total_trees) * 100
      
      # Skip if would exceed maximum constraints
      if (test_pct_groups > MAX_GROUPS_PCT ||
          test_pct_plots > MAX_PCT || 
          test_pct_trees > MAX_PCT) {
        next
      }
      
      # Calculate distribution distance if we ADD this group
      test_groups <- c(selected_groups, candidate_group)
      test_plots <- plots_df |> filter(group_id %in% test_groups)
      
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
      
      # Store candidate
      valid_idx <- valid_idx + 1
      candidate_scores$group_id[valid_idx] <- candidate_group
      candidate_scores$dist_distance[valid_idx] <- test_dist_distance
      candidate_scores$n_plots[valid_idx] <- candidate_n_plots
      candidate_scores$test_pct_plots[valid_idx] <- test_pct_plots
      candidate_scores$test_pct_trees[valid_idx] <- test_pct_trees
    }
    
    # Remove unused rows
    if (valid_idx == 0) {
      break  # No valid additions possible
    }
    candidate_scores <- candidate_scores |> slice(1:valid_idx)
    
    # Sort by distribution distance (ascending - lower is better)
    candidate_scores <- candidate_scores |>
      arrange(dist_distance)
    
    # Stochastic selection from top K candidates
    n_candidates <- min(TOP_K_CANDIDATES, nrow(candidate_scores))
    top_candidates <- candidate_scores |>
      slice(1:n_candidates)
    
    # Calculate selection probabilities
    weights <- exp(-((1:n_candidates) - 1) / STOCHASTIC_TEMP)
    weights <- weights / sum(weights)
    
    # Select group to add
    selected_idx <- sample(1:n_candidates, size = 1, prob = weights)
    group_to_add <- top_candidates$group_id[selected_idx]
    
    # Add the selected group
    selected_groups <- c(selected_groups, group_to_add)
    available_groups <- setdiff(available_groups, group_to_add)
    
    # Update metrics and save state
    iteration <- iteration + 1
    metrics <- calculate_metrics(selected_groups, plots_df,
                                total_groups, total_plots, total_trees)
    selected_plots <- plots_df |> filter(group_id %in% selected_groups)
    
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
  
  if (is.null(reference_df)) {
    reference_df <- plots_df
  }
  
  factorial_dists <- list()
  
  for (combo in FACTORIAL_COMBINATIONS) {
    var1 <- combo$var1
    var2 <- combo$var2
    
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
        pct_plots_selected = pmin(pct_plots_selected_raw, 40),
        pct_trees_selected = pmin(pct_trees_selected_raw, 40),
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
        pct_plots_selected = pmin(pct_plots_selected_raw, 40),
        pct_trees_selected = pmin(pct_trees_selected_raw, 40),
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

#' Main function to select withheld groups via greedy forward selection
#'
#' @param plots_df Data frame with columns: plot_id, group_id, attributes, n_trees
#' @param required_groups Vector of group_ids that must be in the withheld set
#' @param reference_distribution Optional: Pre-calculated distribution to match
#' @param reference_factorial Optional: Pre-calculated factorial distribution to match
#' @param reference_plots_df Optional: Full reference dataset for binning consistency
#' @param previous_selected_plots Optional: Plots already selected in previous tiers
#' @return List containing withheld and training plots with diagnostics
select_withheld_groups <- function(plots_df, required_groups = c(),
                                  reference_distribution = NULL,
                                  reference_factorial = NULL,
                                  reference_plots_df = NULL,
                                  previous_selected_plots = NULL) {
  
  cat("Starting proportional stratified group selection (GREEDY FORWARD)...\n")
  cat("Method: Greedy forward selection (add groups that best match distribution)\n")
  
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
  
  # Prepare quantile bins
  bin_info <- prepare_quantile_bins(binning_df, n_quantiles)
  plots_df_binned <- assign_quantile_bins(plots_df, bin_info)
  
  # Calculate catalog distribution
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
    tier_factorial <- calculate_factorial_distribution(plots_df_binned, reference_df = reference_plots_df)
  } else {
    tier_distribution <- NULL
    tier_factorial <- NULL
  }
  
  cat("  Catalog distribution prepared\n\n")
  
  # ============================================================================
  # STEP 1.5: Handle required groups
  # ============================================================================
  
  required_groups <- intersect(required_groups, unique(plots_df$group_id))
  
  if (length(required_groups) > 0) {
    cat(sprintf("Step 1.5: Checking required groups (%d groups)...\n", 
                length(required_groups)))
    
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
      
      cat(sprintf("  ✓ Required groups already in acceptable range [%d%%, %d%%]\n",
                  MIN_PCT, MAX_PCT))
      cat("  Using required groups as final selection\n\n")
      
      required_plots <- plots_df_binned |> filter(group_id %in% required_groups)
      required_dist <- calculate_distribution(required_plots)
      dist_distance <- overall_distribution_distance(required_dist, catalog_distribution)
      
      training_groups <- setdiff(unique(plots_df$group_id), required_groups)
      training_plots <- plots_df_binned |> filter(group_id %in% training_groups)
      
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
      
      cat(sprintf("  ⚠ Warning: Required groups exceed MAX_PCT (%.1f%% plots, %.1f%% trees)\n",
                  required_metrics$pct_plots, required_metrics$pct_trees))
      cat(sprintf("  Accepting as-is since these are required from higher hierarchy levels\n\n"))
      
      required_plots <- plots_df_binned |> filter(group_id %in% required_groups)
      required_dist <- calculate_distribution(required_plots)
      dist_distance <- overall_distribution_distance(required_dist, catalog_distribution)
      
      training_groups <- setdiff(unique(plots_df$group_id), required_groups)
      training_plots <- plots_df_binned |> filter(group_id %in% training_groups)
      
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
      cat(sprintf("  Required groups below target. Will proceed with greedy forward selection\n"))
      cat(sprintf("  (Required groups locked in, algorithm will add more)\n\n"))
    }
  }
  
  # ============================================================================
  # STEP 2: Run greedy forward selection
  # ============================================================================
  
  cat(sprintf("Step 2: Running %d independent greedy forward selection runs...\n", N_RUNS))
  
  # Set up parallel processing
  if (PARALLELIZE) {
    n_cores <- availableCores()
    cat(sprintf("  Using parallel processing with %d cores\n", n_cores))
    plan(multisession, workers = n_cores)
  } else {
    cat("  Running sequentially (PARALLELIZE = FALSE)\n")
    plan(sequential)
  }
  
  # Run greedy forward selection in parallel
  all_states <- future_map(1:N_RUNS, function(run) {
    
    states <- greedy_forward_selection(
      plots_df = plots_df_binned,
      catalog_distribution = catalog_distribution,
      catalog_factorial = catalog_factorial_ref,
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
  
  # ============================================================================
  # STEP 3: Find best solution
  # ============================================================================
  
  cat("Step 3: Selecting best solution...\n")
  
  cat(sprintf("  Found %d valid solutions\n", nrow(valid_states)))
  
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
    arrange(dist_distance, target_distance) |>
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
  
  withheld_plots <- plots_df_binned |>
    filter(group_id %in% withheld_groups)
  
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
  
  # Calculate factorial distribution
  if (use_reference_distribution) {
    withheld_factorial <- calculate_factorial_distribution(withheld_plots, reference_df = reference_plots_df)
  } else {
    withheld_factorial <- calculate_factorial_distribution(withheld_plots, reference_df = plots_df_binned)
  }
  
  # Generate diagnostics
  if (use_reference_distribution) {
    cat("  Generating dual diagnostics (reference + tier-specific)...\n")
    
    diagnostics$distribution_comparison_reference <- create_distribution_comparison(
      catalog_distribution,
      withheld_distribution
    )
    
    diagnostics$distribution_comparison_tier <- create_distribution_comparison(
      tier_distribution,
      withheld_distribution
    )
    
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
    
    diagnostics$factorial_distributions <- list(
      reference_catalog = catalog_factorial_ref,
      tier_catalog = tier_factorial,
      withheld = withheld_factorial
    )
    
    diagnostics$factorial_plots <- create_dual_factorial_plots(
      catalog_factorial_ref,
      tier_factorial,
      withheld_factorial
    )
    
    diagnostics$target_distance <- best_solution$target_distance
    diagnostics$distribution_distance_reference <- best_solution$dist_distance
    diagnostics$distribution_distance_tier <- overall_distribution_distance(withheld_distribution, tier_distribution)
    
  } else {
    cat("  Generating standard diagnostics (tier-specific only)...\n")
    
    diagnostics$distribution_comparison <- create_distribution_comparison(
      catalog_distribution,
      withheld_distribution
    )
    
    diagnostics$summary_stats <- create_summary_stats_comparison(
      plots_df,
      withheld_plots
    )
    
    catalog_factorial <- calculate_factorial_distribution(plots_df_binned)
    diagnostics$factorial_distributions <- list(
      catalog = catalog_factorial,
      withheld = withheld_factorial
    )
    diagnostics$factorial_plots <- create_factorial_plots(
      catalog_factorial,
      withheld_factorial
    )
    
    diagnostics$target_distance <- best_solution$target_distance
    diagnostics$distribution_distance <- best_solution$dist_distance
  }
  
  cat("  Diagnostics complete\n\n")
  
  # ============================================================================
  # Return results
  # ============================================================================
  
  results <- list(
    withheld_group_ids = withheld_groups,
    training_group_ids = training_groups,
    withheld_plots = withheld_plots,
    training_plots = training_plots,
    diagnostics = diagnostics,
    all_states = all_states_df,
    best_solution = best_solution,
    config = list(
      algorithm = "greedy_forward",
      n_quantiles = n_quantiles,
      target_pct = TARGET_PCT,
      min_pct = MIN_PCT,
      max_pct = MAX_PCT,
      n_runs = N_RUNS,
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
  cat("PROPORTIONAL STRATIFIED SAMPLING REPORT (GREEDY FORWARD)\n")
  cat("=" |> str_pad(width = 80, side = "both", pad = "="), "\n\n")
  
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
  
  cat(sprintf("Target distance from %.0f%%: %.2f\n", 
              results$config$target_pct,
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
    cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
    for (var in CONTINUOUS_VARS) {
      cat(sprintf("\n%s:\n", toupper(var)))
      print(results$diagnostics$distribution_comparison_tier[[var]] |> 
            select(bin, 
                   catalog_n_plots, catalog_pct_plots, 
                   selected_n_plots, selected_pct_plots,
                   diff_pct_plots), 
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
  cat(sprintf("  Runs: %d\n", results$config$n_runs))
  cat(sprintf("  Reference distribution: %s\n", 
              ifelse(use_reference, "Yes (cross-tier)", "No (within-tier)")))
  cat("\n")
}

#' Save factorial distribution plots to files
#'
#' @param results Output from select_withheld_groups()
#' @param output_dir Directory to save plots
#' @param width Plot width in inches
#' @param height Plot height in inches
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
#' @param all_plots_df Data frame containing ALL plots across all tiers
#' @return List with reference_distribution, reference_factorial, and binned data
prepare_reference_distributions <- function(all_plots_df) {
  
  cat("Preparing reference distributions from full dataset...\n")
  cat(sprintf("  Total plots: %d\n", nrow(all_plots_df)))
  
  n_quantiles <- calculate_optimal_n_quantiles(nrow(all_plots_df))
  cat(sprintf("  Using %d quantiles for continuous variables\n", n_quantiles))
  
  bin_info <- prepare_quantile_bins(all_plots_df, n_quantiles)
  all_plots_binned <- assign_quantile_bins(all_plots_df, bin_info)
  
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
# EXAMPLE USAGE
# ==============================================================================

# # Standard within-tier stratification
# plots <- read_csv("path/to/plots.csv")
# results <- select_withheld_groups(plots_df = plots)
# print_selection_report(results)
# write_csv(results$withheld_plots, "withheld_plots.csv")
# save_factorial_plots(results, output_dir = "factorial_plots")