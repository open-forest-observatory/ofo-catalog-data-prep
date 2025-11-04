# ==============================================================================
# SIMPLIFIED GREEDY FORWARD SELECTION FOR STRATIFIED SAMPLING
# ==============================================================================
# Purpose: Select groups to withhold for validation matching catalog distribution
# Method: Greedy forward selection with count-based distance metric
# ==============================================================================

library(tidyverse)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

set.seed(42)

# Continuous variables to stratify
CONTINUOUS_VARS <- c("ppt", "trees_per_ha", "mean_ba_live", "area_ha")

# Categorical variables to stratify
CATEGORICAL_VARS <- c("ecoregion", "sp_comp_group", "project_name", "pairing_tier")

# Binning parameters
TARGET_PLOTS_PER_BIN <- 10
MAX_BINS <- 5
N_BINS_FACTORIAL <- 3

# Selection parameters
MIN_PCT <- 12
MAX_PCT <- 22
MIN_GROUPS_PCT <- 15
MAX_GROUPS_PCT <- 25
FACTORIAL_WEIGHT <- 0.5

# Factorial combinations to track
FACTORIAL_COMBINATIONS <- list(
  list(var1 = "mean_ba_live", var2 = "sp_comp_group"),
  list(var1 = "trees_per_ha", var2 = "sp_comp_group"),
  list(var1 = "ppt", var2 = "sp_comp_group"),
  list(var1 = "mean_ba_live", var2 = "ecoregion"),
  list(var1 = "trees_per_ha", var2 = "ecoregion"),
  list(var1 = "ppt", var2 = "ecoregion"),
  list(var1 = "trees_per_ha", var2 = "mean_ba_live"),
  list(var1 = "pairing_tier", var2 = "sp_comp_group")
)

# ==============================================================================
# BINNING FUNCTIONS
# ==============================================================================

calculate_optimal_n_quantiles <- function(n_plots_total) {
  max_quantiles <- floor(n_plots_total / TARGET_PLOTS_PER_BIN)
  n_quantiles <- min(max_quantiles, MAX_BINS)
  max(2, n_quantiles)
}

prepare_quantile_bins <- function(plots_df, n_quantiles) {
  bin_info <- list()
  for (var in CONTINUOUS_VARS) {
    probs <- seq(0, 1, length.out = n_quantiles + 1)
    breaks <- unique(quantile(plots_df[[var]], probs = probs, na.rm = TRUE))
    labels <- if (length(breaks) > 2) paste0("Q", 1:(length(breaks) - 1)) else c("Low", "High")
    bin_info[[var]] <- list(breaks = breaks, labels = labels)
  }
  bin_info
}

assign_quantile_bins <- function(plots_df, bin_info) {
  result <- plots_df
  for (var in CONTINUOUS_VARS) {
    bin_col <- paste0(var, "_bin")
    result[[bin_col]] <- cut(result[[var]], 
                             breaks = bin_info[[var]]$breaks,
                             labels = bin_info[[var]]$labels,
                             include.lowest = TRUE)
  }
  result
}

# ==============================================================================
# DISTRIBUTION CALCULATION
# ==============================================================================

calculate_distribution <- function(plots_df) {
  distribution <- list()
  
  # Continuous variables
  for (var in CONTINUOUS_VARS) {
    bin_col <- paste0(var, "_bin")
    dist <- plots_df |>
      filter(!is.na(.data[[bin_col]])) |>
      group_by(.data[[bin_col]]) |>
      summarise(n_plots = n(), n_trees = sum(n_trees, na.rm = TRUE), .groups = "drop") |>
      rename(bin = !!bin_col)
    distribution[[var]] <- dist
  }
  
  # Categorical variables
  for (var in CATEGORICAL_VARS) {
    dist <- plots_df |>
      filter(!is.na(.data[[var]])) |>
      group_by(.data[[var]]) |>
      summarise(n_plots = n(), n_trees = sum(n_trees, na.rm = TRUE), .groups = "drop") |>
      rename(category = !!var)
    distribution[[var]] <- dist
  }
  
  distribution
}

calculate_factorial_distribution <- function(plots_df) {
  factorial_dists <- list()
  
  for (combo in FACTORIAL_COMBINATIONS) {
    var1 <- combo$var1
    var2 <- combo$var2
    plots_working <- plots_df
    
    # Create factorial bins for continuous variables
    for (var in c(var1, var2)) {
      if (var %in% CONTINUOUS_VARS) {
        col_name <- paste0(var, "_factorial_bin")
        probs <- seq(0, 1, length.out = N_BINS_FACTORIAL + 1)
        breaks <- unique(quantile(plots_df[[var]], probs = probs, na.rm = TRUE))
        if (length(breaks) > 1) {
          plots_working[[col_name]] <- cut(plots_working[[var]], breaks = breaks, 
                                          include.lowest = TRUE, dig.lab = 3)
        } else {
          plots_working[[col_name]] <- as.factor(paste0(var, "_all"))
        }
      }
    }
    
    col1 <- if (var1 %in% CONTINUOUS_VARS) paste0(var1, "_factorial_bin") else var1
    col2 <- if (var2 %in% CONTINUOUS_VARS) paste0(var2, "_factorial_bin") else var2
    
    dist <- plots_working |>
      filter(!is.na(.data[[col1]]) & !is.na(.data[[col2]])) |>
      group_by(.data[[col1]], .data[[col2]]) |>
      summarise(n_plots = n(), n_trees = sum(n_trees, na.rm = TRUE), .groups = "drop") |>
      rename(var1_value = !!col1, var2_value = !!col2)
    
    factorial_dists[[paste(var1, var2, sep = "_x_")]] <- list(distribution = dist)
  }
  
  factorial_dists
}

# ==============================================================================
# DISTANCE METRICS
# ==============================================================================

distribution_distance_counts <- function(dist_selected, dist_catalog, var, target_n_plots) {
  key_col <- if (var %in% CONTINUOUS_VARS) "bin" else "category"
  
  catalog_with_targets <- dist_catalog |>
    select(all_of(key_col), catalog_n_plots = n_plots) |>
    mutate(target_count = (catalog_n_plots / sum(catalog_n_plots)) * target_n_plots)
  
  comparison <- catalog_with_targets |>
    full_join(dist_selected |> select(all_of(key_col), selected_n_plots = n_plots), by = key_col) |>
    mutate(target_count = replace_na(target_count, 0),
           selected_n_plots = replace_na(selected_n_plots, 0),
           abs_diff = abs(target_count - selected_n_plots))
  
  sum(comparison$abs_diff)
}

overall_distribution_distance_counts <- function(dist_selected, dist_catalog, target_n_plots) {
  all_vars <- c(CONTINUOUS_VARS, CATEGORICAL_VARS)
  distances <- map_dbl(all_vars, ~distribution_distance_counts(
    dist_selected[[.x]], dist_catalog[[.x]], .x, target_n_plots))
  mean(distances) / target_n_plots
}

factorial_distribution_distance_counts <- function(factorial_selected, factorial_catalog, target_n_plots) {
  if (length(FACTORIAL_COMBINATIONS) == 0) return(0)
  
  distances <- map_dbl(names(factorial_catalog), function(combo_name) {
    catalog_dist <- factorial_catalog[[combo_name]]$distribution
    selected_dist <- factorial_selected[[combo_name]]$distribution
    
    comparison <- catalog_dist |>
      rename(catalog_n_plots = n_plots) |>
      mutate(target_count = (catalog_n_plots / sum(catalog_n_plots)) * target_n_plots) |>
      full_join(selected_dist |> rename(selected_n_plots = n_plots), 
                by = c("var1_value", "var2_value")) |>
      mutate(target_count = replace_na(target_count, 0),
             selected_n_plots = replace_na(selected_n_plots, 0),
             abs_diff = abs(target_count - selected_n_plots))
    
    sum(comparison$abs_diff)
  })
  
  mean(distances) / target_n_plots
}

combined_distribution_distance_counts <- function(dist_selected, dist_catalog,
                                                  factorial_selected, 
                                                  factorial_catalog,
                                                  target_n_plots) {
  single_var_distance <- overall_distribution_distance_counts(dist_selected, dist_catalog, target_n_plots)
  factorial_distance <- factorial_distribution_distance_counts(factorial_selected, factorial_catalog, target_n_plots)
  (1 - FACTORIAL_WEIGHT) * single_var_distance + FACTORIAL_WEIGHT * factorial_distance
}

# ==============================================================================
# GREEDY FORWARD SELECTION
# ==============================================================================

greedy_forward_selection <- function(plots_df, catalog_distribution, catalog_factorial,
                                    total_groups, total_plots, total_trees,
                                    target_n_plots) {
  
  selected_groups <- c()
  available_groups <- unique(plots_df$group_id)
  
  # Pre-calculate group statistics
  group_stats <- plots_df |>
    group_by(group_id) |>
    summarise(n_plots = n(), n_trees = sum(n_trees, na.rm = TRUE), .groups = "drop")
  
  iteration <- 0
  
  while (length(available_groups) > 0) {
    
    iteration <- iteration + 1
    
    # Calculate current metrics
    n_groups <- length(selected_groups)
    if (n_groups > 0) {
      selected_plots <- plots_df |> filter(group_id %in% selected_groups)
      n_plots <- nrow(selected_plots)
      n_trees <- sum(selected_plots$n_trees, na.rm = TRUE)
    } else {
      n_plots <- 0
      n_trees <- 0
    }
    
    pct_groups <- (n_groups / total_groups) * 100
    pct_plots <- (n_plots / total_plots) * 100
    pct_trees <- (n_trees / total_trees) * 100
    
    # Progress indicator
    if (iteration %% 1 == 0 || pct_plots >= MIN_PCT) {
      cat(sprintf("\r  Iteration %d: %d groups (%.1f%%), %d plots (%.1f%%), %d trees (%.1f%%)   ",
                  iteration, n_groups, pct_groups, n_plots, pct_plots, n_trees, pct_trees))
      flush.console()
    }
    
    # Check if we've reached max constraints
    if (pct_groups >= MAX_GROUPS_PCT || pct_plots >= MAX_PCT || pct_trees >= MAX_PCT) {
      cat("\n")
      break
    }
    
    # Evaluate all candidates
    best_group <- NULL
    best_distance <- Inf
    
    for (candidate_group in available_groups) {
      
      # Quick constraint check
      candidate_stats <- group_stats |> filter(group_id == candidate_group)
      test_n_groups <- n_groups + 1
      test_n_plots <- n_plots + candidate_stats$n_plots
      test_n_trees <- n_trees + candidate_stats$n_trees
      test_pct_groups <- (test_n_groups / total_groups) * 100
      test_pct_plots <- (test_n_plots / total_plots) * 100
      test_pct_trees <- (test_n_trees / total_trees) * 100
      
      if (test_pct_groups > MAX_GROUPS_PCT || test_pct_plots > MAX_PCT || test_pct_trees > MAX_PCT) {
        next
      }
      
      # Calculate distance with this candidate
      test_groups <- c(selected_groups, candidate_group)
      test_plots <- plots_df |> filter(group_id %in% test_groups)
      test_dist <- calculate_distribution(test_plots)
      test_factorial <- calculate_factorial_distribution(test_plots)
      
      test_distance <- combined_distribution_distance_counts(
        test_dist, catalog_distribution, test_factorial, catalog_factorial, target_n_plots)
      
      if (test_distance < best_distance) {
        best_distance <- test_distance
        best_group <- candidate_group
      }
    }
    
    if (is.null(best_group)) {
      cat("\n")
      break  # No valid additions possible
    }
    
    # Add the best group
    selected_groups <- c(selected_groups, best_group)
    available_groups <- setdiff(available_groups, best_group)
  }
  
  # Clear progress line if we exited normally
  if (!is.null(best_group)) {
    cat("\n")
  }
  
  selected_groups
}

# ==============================================================================
# MAIN FUNCTION
# ==============================================================================

select_withheld_groups <- function(plots_df) {
  
  cat("Starting greedy forward selection...\n")
  
  # Validate input data
  required_cols <- c("group_id", "n_trees", CONTINUOUS_VARS, CATEGORICAL_VARS)
  missing_cols <- setdiff(required_cols, names(plots_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Calculate totals
  total_groups <- n_distinct(plots_df$group_id)
  total_plots <- nrow(plots_df)
  total_trees <- sum(plots_df$n_trees, na.rm = TRUE)
  
  # Create bins
  n_quantiles <- calculate_optimal_n_quantiles(total_plots)
  bin_info <- prepare_quantile_bins(plots_df, n_quantiles)
  plots_df_binned <- assign_quantile_bins(plots_df, bin_info)
  
  # Check for and warn about plots with missing values in continuous variables
  na_counts <- sapply(CONTINUOUS_VARS, function(var) {
    bin_col <- paste0(var, "_bin")
    sum(is.na(plots_df_binned[[bin_col]]))
  })
  
  total_excluded <- sum(na_counts > 0)
  if (total_excluded > 0) {
    cat("\n⚠ Warning: Plots with missing values will be excluded from selection:\n")
    for (var in names(na_counts[na_counts > 0])) {
      cat(sprintf("  - %s: %d plots with NA values\n", var, na_counts[var]))
    }
    # Count unique plots with any NA
    any_na <- apply(plots_df_binned[, paste0(CONTINUOUS_VARS, "_bin"), drop = FALSE], 1, 
                    function(x) any(is.na(x)))
    n_excluded <- sum(any_na)
    cat(sprintf("  Total plots excluded: %d (%.1f%%)\n\n", n_excluded, 100 * n_excluded / total_plots))
  }
  
  # Calculate catalog distribution
  catalog_distribution <- calculate_distribution(plots_df_binned)
  catalog_factorial <- calculate_factorial_distribution(plots_df_binned)
  
  # Calculate target
  target_n_plots <- round(total_plots * ((MIN_PCT + MAX_PCT) / 2) / 100)
  
  cat(sprintf("Total: %d groups, %d plots, %d trees\n", total_groups, total_plots, total_trees))
  cat(sprintf("Target: %d plots\n", target_n_plots))
  
  # Run greedy selection
  cat("Running greedy selection...\n")
  withheld_groups <- greedy_forward_selection(
    plots_df_binned, catalog_distribution, catalog_factorial,
    total_groups, total_plots, total_trees, target_n_plots)
  
  # Extract final sets
  withheld_plots <- plots_df_binned |> filter(group_id %in% withheld_groups)
  training_groups <- setdiff(unique(plots_df$group_id), withheld_groups)
  training_plots <- plots_df_binned |> filter(group_id %in% training_groups)
  
  # Calculate final metrics
  final_pct_groups <- (length(withheld_groups) / total_groups) * 100
  final_pct_plots <- (nrow(withheld_plots) / total_plots) * 100
  final_pct_trees <- (sum(withheld_plots$n_trees) / total_trees) * 100
  
  cat(sprintf("\nFinal selection:\n"))
  cat(sprintf("  Withheld: %d groups (%.1f%%), %d plots (%.1f%%), %d trees (%.1f%%)\n",
              length(withheld_groups), final_pct_groups,
              nrow(withheld_plots), final_pct_plots,
              sum(withheld_plots$n_trees), final_pct_trees))
  
  # Calculate distributions for reporting
  withheld_distribution <- calculate_distribution(withheld_plots)
  withheld_factorial <- calculate_factorial_distribution(withheld_plots)
  
  # Return results
  list(
    withheld_group_ids = withheld_groups,
    training_group_ids = training_groups,
    withheld_plots = withheld_plots,
    training_plots = training_plots,
    catalog_distribution = catalog_distribution,
    catalog_factorial = catalog_factorial,
    withheld_distribution = withheld_distribution,
    withheld_factorial = withheld_factorial,
    config = list(
      n_quantiles = n_quantiles,
      target_n_plots = target_n_plots
    )
  )
}

# ==============================================================================
# REPORTING FUNCTIONS
# ==============================================================================

#' Print selection report with count-based metrics
#'
#' @param results Output from select_withheld_groups()
print_selection_report <- function(results) {
  
  cat("=" |> str_pad(width = 80, side = "both", pad = "="), "\n")
  cat("PROPORTIONAL STRATIFIED SAMPLING REPORT\n")
  cat("=" |> str_pad(width = 80, side = "both", pad = "="), "\n\n")
  
  cat("METHOD: Greedy forward selection with count-based distance metric\n\n")
  
  # Summary metrics
  withheld_plots <- results$withheld_plots
  total_plots <- nrow(withheld_plots) + nrow(results$training_plots)
  total_trees <- sum(withheld_plots$n_trees, na.rm = TRUE) + sum(results$training_plots$n_trees, na.rm = TRUE)
  total_groups <- length(results$withheld_group_ids) + length(results$training_group_ids)
  
  target_pct <- ((MIN_PCT + MAX_PCT) / 2)
  
  cat("SUMMARY METRICS\n")
  cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
  
  summary_df <- tibble(
    metric = c("Groups", "Plots", "Trees"),
    n_withheld = c(length(results$withheld_group_ids),
                   nrow(withheld_plots),
                   sum(withheld_plots$n_trees, na.rm = TRUE)),
    n_total = c(total_groups, total_plots, total_trees),
    pct_withheld = (n_withheld / n_total) * 100,
    distance_from_target = abs(pct_withheld - target_pct)
  )
  
  print(summary_df, n = Inf)
  cat("\n")
  
  cat(sprintf("Target: %d plots (%.0f%% of catalog)\n",
              results$config$target_n_plots, target_pct))
  cat(sprintf("Target distance from %.0f%%: %.2f\n\n", 
              target_pct, sum(summary_df$distance_from_target[2:3])))
  
  # Use pre-calculated distributions from results
  catalog_dist <- results$catalog_distribution
  withheld_dist <- results$withheld_distribution
  
  # Distribution comparisons for continuous variables
  cat("CONTINUOUS VARIABLE DISTRIBUTIONS\n")
  cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
  
  for (var in CONTINUOUS_VARS) {
    cat(sprintf("\n%s:\n", toupper(var)))
    
    comparison <- catalog_dist[[var]] |>
      rename(catalog_n_plots = n_plots, catalog_n_trees = n_trees) |>
      mutate(target_n_plots = (catalog_n_plots / sum(catalog_n_plots)) * results$config$target_n_plots) |>
      full_join(withheld_dist[[var]] |> rename(selected_n_plots = n_plots, selected_n_trees = n_trees),
                by = "bin") |>
      mutate(selected_n_plots = replace_na(selected_n_plots, 0),
             selected_n_trees = replace_na(selected_n_trees, 0),
             abs_diff_plots = abs(target_n_plots - selected_n_plots)) |>
      select(bin, catalog_n_plots, target_n_plots, selected_n_plots, 
             catalog_n_trees, selected_n_trees, abs_diff_plots)
    
    print(comparison, n = Inf)
    cat(sprintf("  Total absolute difference (plots): %.2f\n", sum(comparison$abs_diff_plots)))
  }
  cat("\n")
  
  # Distribution comparisons for categorical variables
  cat("CATEGORICAL VARIABLE DISTRIBUTIONS\n")
  cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
  
  for (var in CATEGORICAL_VARS) {
    cat(sprintf("\n%s:\n", toupper(var)))
    
    comparison <- catalog_dist[[var]] |>
      rename(catalog_n_plots = n_plots, catalog_n_trees = n_trees) |>
      mutate(target_n_plots = (catalog_n_plots / sum(catalog_n_plots)) * results$config$target_n_plots) |>
      full_join(withheld_dist[[var]] |> rename(selected_n_plots = n_plots, selected_n_trees = n_trees),
                by = "category") |>
      mutate(selected_n_plots = replace_na(selected_n_plots, 0),
             selected_n_trees = replace_na(selected_n_trees, 0),
             abs_diff_plots = abs(target_n_plots - selected_n_plots)) |>
      select(category, catalog_n_plots, target_n_plots, selected_n_plots,
             catalog_n_trees, selected_n_trees, abs_diff_plots) |>
      arrange(desc(catalog_n_plots))
    
    print(comparison, n = Inf)
    cat(sprintf("  Total absolute difference (plots): %.2f\n", sum(comparison$abs_diff_plots)))
  }
  cat("\n")
  
  # Summary statistics comparison
  cat("SUMMARY STATISTICS COMPARISON\n")
  cat("-" |> str_pad(width = 80, side = "both", pad = "-"), "\n")
  
  catalog_plots <- bind_rows(withheld_plots, results$training_plots)
  
  stats <- map_dfr(CONTINUOUS_VARS, function(var) {
    tibble(
      variable = var,
      catalog_mean = mean(catalog_plots[[var]], na.rm = TRUE),
      selected_mean = mean(withheld_plots[[var]], na.rm = TRUE),
      mean_pct_diff = ((selected_mean - catalog_mean) / catalog_mean) * 100,
      catalog_sd = sd(catalog_plots[[var]], na.rm = TRUE),
      selected_sd = sd(withheld_plots[[var]], na.rm = TRUE),
      sd_pct_diff = ((selected_sd - catalog_sd) / catalog_sd) * 100
    )
  })
  
  print(stats, n = Inf)
  cat("\n")
  
  cat("=" |> str_pad(width = 80, side = "both", pad = "="), "\n\n")
  
  cat("Withheld group IDs:\n")
  cat(paste(results$withheld_group_ids, collapse = ", "), "\n")
}

#' Create factorial distribution heatmaps
#'
#' @param results Output from select_withheld_groups()
#' @return List of ggplot objects
create_factorial_plots <- function(results) {
  
  # Use pre-calculated factorial distributions from results
  catalog_factorial <- results$catalog_factorial
  withheld_factorial <- results$withheld_factorial
  
  plots <- list()
  
  for (combo_name in names(catalog_factorial)) {
    combo_info <- catalog_factorial[[combo_name]]
    var1 <- strsplit(combo_name, "_x_")[[1]][1]
    var2 <- strsplit(combo_name, "_x_")[[1]][2]
    
    combined <- catalog_factorial[[combo_name]]$distribution |>
      rename(catalog_n_plots = n_plots, catalog_n_trees = n_trees) |>
      full_join(withheld_factorial[[combo_name]]$distribution |>
                  rename(selected_n_plots = n_plots, selected_n_trees = n_trees),
                by = c("var1_value", "var2_value")) |>
      mutate(
        catalog_n_plots = replace_na(catalog_n_plots, 0),
        catalog_n_trees = replace_na(catalog_n_trees, 0),
        selected_n_plots = replace_na(selected_n_plots, 0),
        selected_n_trees = replace_na(selected_n_trees, 0),
        pct_plots_selected = ifelse(catalog_n_plots > 0, 
                                     100 * selected_n_plots / catalog_n_plots, 0),
        pct_trees_selected = ifelse(catalog_n_trees > 0,
                                     100 * selected_n_trees / catalog_n_trees, 0)
      )
    
    # Heatmap for plots
    p_plots <- ggplot(combined, aes(x = var1_value, y = var2_value)) +
      geom_tile(aes(fill = pct_plots_selected), color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("C:%d\nS:%d\n%.0f%%", 
                                     catalog_n_plots, selected_n_plots, pct_plots_selected)),
                size = 3, lineheight = 0.8) +
      scale_fill_gradient2(low = "#d73027", mid = "#fee08b", high = "#1a9850",
                          midpoint = 20, limits = c(0, 100), name = "% Selected") +
      labs(title = paste("Factorial Distribution:", var1, "×", var2, "(Plots)"),
           subtitle = "Colors show % of catalog selected in each bin",
           x = var1, y = var2) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold"))
    
    # Heatmap for trees
    p_trees <- ggplot(combined, aes(x = var1_value, y = var2_value)) +
      geom_tile(aes(fill = pct_trees_selected), color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("C:%d\nS:%d\n%.0f%%", 
                                     catalog_n_trees, selected_n_trees, pct_trees_selected)),
                size = 3, lineheight = 0.8) +
      scale_fill_gradient2(low = "#d73027", mid = "#fee08b", high = "#1a9850",
                          midpoint = 20, limits = c(0, 100), name = "% Selected") +
      labs(title = paste("Factorial Distribution:", var1, "×", var2, "(Trees)"),
           subtitle = "Colors show % of catalog selected in each bin",
           x = var1, y = var2) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold"))
    
    plots[[paste0(combo_name, "_plots")]] <- p_plots
    plots[[paste0(combo_name, "_trees")]] <- p_trees
  }
  
  plots
}
