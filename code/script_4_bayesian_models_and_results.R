#!/usr/bin/env Rscript

# =============================================================================
# SCRIPT 4: BAYESIAN MODELS AND RESULTS
# =============================================================================
#
# This script reproduces the final Bayesian results for sequence structure and
# acoustic structure. It starts from the two model-ready tables created from
# the distance matrices in Scripts 2 and 3.
#
# The response is the standardized distance between two recording sessions.
# For every acoustic or sequence metric, the reported stage effect is:
#
#       expected distance after pairing - expected distance before pairing
#
# Negative values therefore indicate convergence and positive values indicate
# divergence. Results are averaged equally over the three study pairs and,
# where a combined result is shown, equally over the four metrics.
#
# Reproducing the reported results does not run Stan. By default, this script
# loads the saved brms model objects and posterior draws supplied with the
# repository, recalculates the reported summaries, and recreates the figures.
# Model sampling occurs only when the script is called with --refit:
#
#   Rscript code/script_4_bayesian_models_and_results.R
#   Rscript code/script_4_bayesian_models_and_results.R --refit
#
# Refitting takes substantially longer and requires a working Stan toolchain.
# =============================================================================

required_packages <- c(
  "tidyverse",
  "brms",
  "tidybayes",
  "bayesplot",
  "posterior",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the required R packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(tidybayes)
  library(bayesplot)
  library(posterior)
  library(patchwork)
})

# -----------------------------------------------------------------------------
# Paths and run settings
# -----------------------------------------------------------------------------

get_script_directory <- function() {
  full_arguments <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", full_arguments, value = TRUE)

  if (length(file_argument)) {
    return(dirname(normalizePath(sub("^--file=", "", file_argument[[1]]))))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    active_file <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(active_file)) {
      return(dirname(normalizePath(active_file)))
    }
  }

  normalizePath(getwd())
}

script_directory <- get_script_directory()
project_directory <- normalizePath(
  file.path(script_directory, ".."),
  mustWork = TRUE
)

arguments <- commandArgs(trailingOnly = TRUE)
allowed_arguments <- "--refit"
unknown_arguments <- setdiff(arguments, allowed_arguments)
if (length(unknown_arguments)) {
  stop(
    "Unknown argument: ",
    paste(unknown_arguments, collapse = ", "),
    ". Use no argument for cached reproduction or --refit to sample the models."
  )
}

# Expensive model fitting is deliberately disabled in a normal run.
REFIT_MODELS <- FALSE
if ("--refit" %in% arguments) {
  REFIT_MODELS <- TRUE
}

model_data_directory <- file.path(project_directory, "data", "models")
model_cache_directory <- file.path(model_data_directory, "cached_models")
posterior_cache_directory <- file.path(
  model_data_directory,
  "posterior_draws"
)

result_directory <- file.path(
  project_directory,
  "results",
  "bayesian_models"
)
table_directory <- file.path(result_directory, "tables")
figure_directory <- file.path(result_directory, "figures")
diagnostic_directory <- file.path(result_directory, "diagnostics")
refitted_model_directory <- file.path(result_directory, "refitted_models")

walk(
  c(table_directory, figure_directory, diagnostic_directory),
  ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE)
)
if (REFIT_MODELS) {
  dir.create(
    refitted_model_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

sequence_data_path <- file.path(
  model_data_directory,
  "sequence_session_distances.csv"
)
acoustic_data_path <- file.path(
  model_data_directory,
  "acoustic_session_distances.csv"
)
sequence_model_path <- file.path(
  model_cache_directory,
  "sequence_model.rds"
)
acoustic_model_path <- file.path(
  model_cache_directory,
  "acoustic_model.rds"
)
pair_draw_path <- file.path(
  posterior_cache_directory,
  "pair_specific_stage_effect_draws.rds"
)

required_files <- c(
  sequence_data_path,
  acoustic_data_path,
  sequence_model_path,
  acoustic_model_path
)
if (!REFIT_MODELS) {
  required_files <- c(required_files, pair_draw_path)
}

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop(
    "The public analysis folder is incomplete. These required files were not found:\n  - ",
    paste(missing_files, collapse = "\n  - ")
  )
}

# -----------------------------------------------------------------------------
# Model inputs and specification
# -----------------------------------------------------------------------------

SEQUENCE_METRICS <- c(
  "transition_matrix",
  "bigram",
  "repeat_A_len",
  "local_alignment"
)

ACOUSTIC_METRICS <- c(
  "T_PC_pairwise",
  "M_PC_pairwise",
  "spectro_DTW",
  "chatter_cosine"
)

metric_labels <- c(
  "transition_matrix" = "Transition probabilities",
  "bigram" = "Bigram distribution",
  "repeat_A_len" = "Repeat distribution of phee",
  "local_alignment" = "Local alignment",
  "T_PC_pairwise" = "Spectro-temporal parameters (STP)",
  "M_PC_pairwise" = "Mel-frequency cepstral coefficients (MFCC)",
  "spectro_DTW" = "Dynamic time warping (DTW)",
  "chatter_cosine" = "Variational autoencoder (VAE)"
)

panel_metric_labels <- c(
  "transition_matrix" = "Transition\nprobabilities",
  "bigram" = "Bigram\ndistribution",
  "repeat_A_len" = "Repeat distribution\nof phee",
  "local_alignment" = "Local alignment",
  "T_PC_pairwise" = "Spectro-temporal\nparameters (STP)",
  "M_PC_pairwise" = "Mel-frequency cepstral\ncoefficients (MFCC)",
  "spectro_DTW" = "Dynamic time\nwarping (DTW)",
  "chatter_cosine" = "Variational autoencoder\n(VAE)"
)

PARTNER_COLOUR <- "#2c728e"
NONPARTNER_COLOUR <- "#b8de29"

require_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

prepare_model_data <- function(data, metric_order, label) {
  required <- c(
    "pair", "context", "stage", "session_1_id", "session_2_id",
    "comparison_id", "metric", "distance", "z_distance"
  )
  require_columns(data, required, label)

  if (anyNA(data[required])) {
    stop(label, " contains missing values in required columns")
  }
  if (!setequal(unique(data$stage), c("before", "after"))) {
    stop(label, " must contain before and after stages")
  }
  if (!setequal(unique(data$context), c("Partner", "Non-partner"))) {
    stop(label, " must contain Partner and Non-partner contexts")
  }
  if (!setequal(unique(data$metric), metric_order)) {
    stop(label, " does not contain the expected four metrics")
  }
  if (any(data$session_1_id == data$session_2_id)) {
    stop(label, " contains a comparison of a session with itself")
  }
  if (anyDuplicated(data[c("comparison_id", "metric")])) {
    stop(label, " contains duplicated comparison-by-metric rows")
  }

  metrics_per_comparison <- data %>%
    distinct(comparison_id, metric) %>%
    count(comparison_id, name = "number_of_metrics")
  if (any(metrics_per_comparison$number_of_metrics != length(metric_order))) {
    stop(label, " has a comparison without all four metrics")
  }

  scaling_check <- data %>%
    group_by(metric) %>%
    summarise(
      mean_z = mean(z_distance),
      population_sd_z = sd(z_distance) * sqrt((n() - 1) / n()),
      .groups = "drop"
    )
  if (any(abs(scaling_check$mean_z) > 1e-8) ||
      any(abs(scaling_check$population_sd_z - 1) > 1e-8)) {
    stop(label, " is not standardized within metric as expected")
  }

  session_levels <- union(
    as.character(data$session_1_id),
    as.character(data$session_2_id)
  )

  data %>%
    mutate(
      stage = factor(stage, levels = c("before", "after")),
      context = factor(context, levels = c("Non-partner", "Partner")),
      metric = factor(metric, levels = metric_order),
      pair = factor(pair),
      comparison_id = factor(comparison_id),
      session_1_id = factor(session_1_id, levels = session_levels),
      session_2_id = factor(session_2_id, levels = session_levels)
    )
}

sequence_data <- read_csv(sequence_data_path, show_col_types = FALSE) %>%
  prepare_model_data(SEQUENCE_METRICS, "Sequence model data")
acoustic_data <- read_csv(acoustic_data_path, show_col_types = FALSE) %>%
  prepare_model_data(ACOUSTIC_METRICS, "Acoustic model data")

if (n_distinct(sequence_data$pair) != 3L ||
    n_distinct(acoustic_data$pair) != 3L) {
  stop("Both analyses must contain the three experimental pairs")
}

final_formula <- brms::bf(
  z_distance ~ stage * context * metric +
    (1 + stage * context || pair) +
    (1 | mm(session_1_id, session_2_id))
)

model_priors <- c(
  brms::prior(normal(0, 1), class = "Intercept"),
  brms::prior(normal(0, 1), class = "b"),
  brms::prior(exponential(1), class = "sd"),
  brms::prior(exponential(1), class = "sigma")
)

# -----------------------------------------------------------------------------
# Load the saved models, or explicitly refit them
# -----------------------------------------------------------------------------

message("Loading the saved sequence and acoustic model objects ...")
sequence_model <- readRDS(sequence_model_path)
acoustic_model <- readRDS(acoustic_model_path)

fit_model <- function(data, output_stem) {
  message("Sampling ", output_stem, " ...")
  detected_cores <- parallel::detectCores(logical = FALSE)
  if (is.na(detected_cores)) {
    detected_cores <- 1L
  }
  brms::brm(
    formula = final_formula,
    data = data,
    family = gaussian(),
    prior = model_priors,
    chains = 4,
    cores = min(4L, detected_cores),
    iter = 4000,
    warmup = 2000,
    seed = 123,
    control = list(adapt_delta = 0.995, max_treedepth = 15),
    backend = "rstan",
    save_pars = brms::save_pars(all = TRUE),
    file = file.path(refitted_model_directory, output_stem),
    file_refit = "always",
    refresh = 100
  )
}

extract_pair_stage_draws <- function(model, data, analysis_label) {
  prediction_grid <- expand_grid(
    pair = levels(data$pair),
    context = levels(data$context),
    metric = levels(data$metric),
    stage = levels(data$stage)
  ) %>%
    mutate(
      pair = factor(pair, levels = levels(data$pair)),
      context = factor(context, levels = levels(data$context)),
      metric = factor(metric, levels = levels(data$metric)),
      stage = factor(stage, levels = levels(data$stage)),
      grid_row = as.character(row_number())
    )

  # Pair effects are retained, while the multi-membership session intercept is
  # omitted. The result is a pair-specific expectation rather than a prediction
  # for a particular recording session.
  expected_distance <- brms::posterior_epred(
    model,
    newdata = prediction_grid,
    re_formula = ~(1 + stage * context || pair)
  )

  expected_wide <- as_tibble(as.data.frame(expected_distance))
  names(expected_wide) <- as.character(seq_len(ncol(expected_wide)))

  expected_wide %>%
    mutate(.draw = row_number()) %>%
    pivot_longer(
      cols = -.draw,
      names_to = "grid_row",
      values_to = "expected_z_distance"
    ) %>%
    left_join(prediction_grid, by = "grid_row") %>%
    group_by(.draw, pair, context, metric) %>%
    summarise(
      stage_effect =
        expected_z_distance[stage == "after"] -
        expected_z_distance[stage == "before"],
      .groups = "drop"
    ) %>%
    mutate(
      analysis = analysis_label,
      pair = as.character(pair),
      context = as.character(context),
      metric = as.character(metric),
      metric_label = recode(metric, !!!metric_labels)
    ) %>%
    select(
      .draw,
      pair,
      analysis,
      context,
      metric,
      metric_label,
      stage_effect
    )
}

if (REFIT_MODELS) {
  message(
    "--refit was supplied. The two Bayesian models will now be sampled; ",
    "this may take several hours."
  )
  sequence_model <- fit_model(sequence_data, "sequence_model_refitted")
  acoustic_model <- fit_model(acoustic_data, "acoustic_model_refitted")

  pair_stage_draws <- bind_rows(
    extract_pair_stage_draws(
      sequence_model,
      sequence_data,
      "Sequence structure"
    ),
    extract_pair_stage_draws(
      acoustic_model,
      acoustic_data,
      "Call structure"
    )
  )
} else {
  message("Loading the saved pair-specific posterior stage effects ...")
  pair_stage_draws <- readRDS(pair_draw_path)
}

require_columns(
  pair_stage_draws,
  c(
    ".draw", "pair", "analysis", "context", "metric", "metric_label",
    "stage_effect"
  ),
  "Posterior draw cache"
)

expected_analyses <- c("Sequence structure", "Call structure")
if (!setequal(unique(pair_stage_draws$analysis), expected_analyses) ||
    !setequal(unique(pair_stage_draws$context), c("Partner", "Non-partner"))) {
  stop("The posterior draw cache does not match the final analyses")
}

# -----------------------------------------------------------------------------
# Posterior summaries
# -----------------------------------------------------------------------------

summarise_posterior <- function(draws, value_column, grouping_columns) {
  draws %>%
    group_by(across(all_of(grouping_columns))) %>%
    summarise(
      estimate_median = median(.data[[value_column]]),
      estimate_mean = mean(.data[[value_column]]),
      lower_95 = quantile(.data[[value_column]], 0.025),
      upper_95 = quantile(.data[[value_column]], 0.975),
      Pr_effect_lt_0 = mean(.data[[value_column]] < 0),
      Pr_effect_gt_0 = mean(.data[[value_column]] > 0),
      posterior_two_sided_tail_probability =
        2 * min(Pr_effect_lt_0, Pr_effect_gt_0),
      .groups = "drop"
    )
}

# Every pair receives equal weight, irrespective of the number of session
# comparisons available for that pair.
stage_draws <- pair_stage_draws %>%
  group_by(.draw, analysis, context, metric, metric_label) %>%
  summarise(stage_effect = mean(stage_effect), .groups = "drop")

saveRDS(
  stage_draws,
  file.path(
    table_directory,
    "pair_marginal_stage_effect_draws_reproduced.rds"
  )
)

model_data_long <- bind_rows(
  sequence_data %>% mutate(analysis = "Sequence structure"),
  acoustic_data %>% mutate(analysis = "Call structure")
) %>%
  mutate(
    across(
      c(pair, stage, context, metric, comparison_id, session_1_id, session_2_id),
      as.character
    )
  )

data_counts <- model_data_long %>%
  group_by(analysis, context, metric) %>%
  summarise(
    N_rows = n(),
    N_comparisons = n_distinct(comparison_id),
    N_sessions = n_distinct(c(session_1_id, session_2_id)),
    .groups = "drop"
  )

stage_effect_summary <- summarise_posterior(
  stage_draws,
  "stage_effect",
  c("analysis", "context", "metric", "metric_label")
) %>%
  left_join(data_counts, by = c("analysis", "context", "metric")) %>%
  mutate(
    estimand_definition = "expected after - expected before",
    direction = "negative = convergence; positive = divergence"
  )

write_csv(
  stage_effect_summary,
  file.path(
    table_directory,
    "stage_effects_by_context_metric_reproduced.csv"
  )
)

average_stage_draws <- stage_draws %>%
  group_by(.draw, analysis, context) %>%
  summarise(stage_effect = mean(stage_effect), .groups = "drop")

average_stage_summary <- summarise_posterior(
  average_stage_draws,
  "stage_effect",
  c("analysis", "context")
) %>%
  mutate(
    estimand_definition = "expected after - expected before",
    direction = "negative = convergence; positive = divergence",
    pair_weighting = "equal over three observed pairs",
    metric_weighting = "equal over four metrics"
  )

write_csv(
  average_stage_summary,
  file.path(
    table_directory,
    "average_stage_effects_by_context_analysis_reproduced.csv"
  )
)
saveRDS(
  average_stage_draws,
  file.path(
    table_directory,
    "pair_marginal_combined_stage_effect_draws_reproduced.rds"
  )
)

metric_context_contrast_draws <- stage_draws %>%
  select(.draw, analysis, metric, metric_label, context, stage_effect) %>%
  pivot_wider(names_from = context, values_from = stage_effect) %>%
  mutate(context_contrast = Partner - .data[["Non-partner"]])

metric_context_contrast_summary <- summarise_posterior(
  metric_context_contrast_draws,
  "context_contrast",
  c("analysis", "metric", "metric_label")
) %>%
  mutate(
    estimand_definition =
      "(after - before) Partner - (after - before) Non-partner",
    direction = "positive = the Non-partner stage change is more negative"
  )

write_csv(
  metric_context_contrast_summary,
  file.path(
    table_directory,
    "partner_minus_nonpartner_stage_effects_by_metric_reproduced.csv"
  )
)

average_context_contrast_draws <- average_stage_draws %>%
  pivot_wider(names_from = context, values_from = stage_effect) %>%
  mutate(context_contrast = Partner - .data[["Non-partner"]])

average_context_contrast_summary <- summarise_posterior(
  average_context_contrast_draws,
  "context_contrast",
  "analysis"
) %>%
  mutate(
    estimand_definition =
      "(after - before) Partner - (after - before) Non-partner",
    direction = "positive = the Non-partner stage change is more negative",
    pair_weighting = "equal over three observed pairs",
    metric_weighting = "equal over four metrics",
    contains_model_extrapolation = TRUE
  )

write_csv(
  average_context_contrast_summary,
  file.path(
    table_directory,
    "average_partner_minus_nonpartner_stage_effects_reproduced.csv"
  )
)

# Identify the pair-by-context cells for which both stages were observed. This
# supports the common-support sensitivity analysis reported with the models.
pair_context_support <- model_data_long %>%
  distinct(analysis, pair, context, stage, comparison_id) %>%
  count(analysis, pair, context, stage, name = "session_comparisons") %>%
  complete(
    analysis,
    pair,
    context = c("Partner", "Non-partner"),
    stage = c("before", "after"),
    fill = list(session_comparisons = 0L)
  ) %>%
  group_by(analysis, pair, context) %>%
  summarise(
    before_comparisons = sum(session_comparisons[stage == "before"]),
    after_comparisons = sum(session_comparisons[stage == "after"]),
    has_before_and_after =
      before_comparisons > 0 & after_comparisons > 0,
    .groups = "drop"
  )

write_csv(
  pair_context_support,
  file.path(table_directory, "pair_context_stage_support_reproduced.csv")
)

common_support_pairs <- pair_context_support %>%
  group_by(analysis, pair) %>%
  summarise(
    supported_contexts = sum(has_before_and_after),
    .groups = "drop"
  ) %>%
  filter(supported_contexts == 2L)

common_support_metric_draws <- pair_stage_draws %>%
  semi_join(common_support_pairs, by = c("analysis", "pair")) %>%
  select(.draw, pair, analysis, metric, metric_label, context, stage_effect) %>%
  pivot_wider(names_from = context, values_from = stage_effect) %>%
  mutate(context_contrast = Partner - .data[["Non-partner"]]) %>%
  group_by(.draw, analysis, metric, metric_label) %>%
  summarise(
    context_contrast = mean(context_contrast),
    pairs_with_direct_support = n_distinct(pair),
    .groups = "drop"
  )

common_support_metric_summary <- summarise_posterior(
  common_support_metric_draws,
  "context_contrast",
  c("analysis", "metric", "metric_label", "pairs_with_direct_support")
) %>%
  mutate(
    estimand_definition = paste(
      "(after - before) Partner - (after - before) Non-partner;",
      "positive values mean the Non-partner stage change is more negative"
    )
  )

write_csv(
  common_support_metric_summary,
  file.path(
    table_directory,
    "common_support_partner_minus_nonpartner_by_metric_reproduced.csv"
  )
)

common_support_combined_draws <- common_support_metric_draws %>%
  group_by(.draw, analysis, pairs_with_direct_support) %>%
  summarise(context_contrast = mean(context_contrast), .groups = "drop")

common_support_combined_summary <- summarise_posterior(
  common_support_combined_draws,
  "context_contrast",
  c("analysis", "pairs_with_direct_support")
) %>%
  mutate(
    estimand_definition = paste(
      "(after - before) Partner - (after - before) Non-partner;",
      "positive values mean the Non-partner stage change is more negative"
    )
  )

write_csv(
  common_support_combined_summary,
  file.path(
    table_directory,
    "common_support_combined_partner_minus_nonpartner_reproduced.csv"
  )
)

direct_support_counts <- pair_context_support %>%
  group_by(analysis, context) %>%
  summarise(
    directly_supported_pairs = sum(has_before_and_after),
    .groups = "drop"
  )
common_support_pair_counts <- common_support_pairs %>%
  count(analysis, name = "directly_supported_pairs")

primary_results <- bind_rows(
  average_stage_summary %>%
    left_join(direct_support_counts, by = c("analysis", "context")) %>%
    transmute(
      analysis,
      estimand = paste(context, "stage effect"),
      estimand_definition,
      direction,
      pairs_in_summary = 3L,
      directly_supported_pairs,
      metrics_in_summary = 4L,
      contains_model_extrapolation = directly_supported_pairs < pairs_in_summary,
      estimate_median,
      lower_95,
      upper_95,
      posterior_probability_for_stated_direction = Pr_effect_lt_0,
      reported_direction = "effect < 0 (convergence)"
    ),
  average_context_contrast_summary %>%
    left_join(common_support_pair_counts, by = "analysis") %>%
    transmute(
      analysis,
      estimand = "Primary stage-by-context interaction",
      estimand_definition,
      direction,
      pairs_in_summary = 3L,
      directly_supported_pairs,
      metrics_in_summary = 4L,
      contains_model_extrapolation,
      estimate_median,
      lower_95,
      upper_95,
      posterior_probability_for_stated_direction = Pr_effect_gt_0,
      reported_direction = "interaction > 0"
    ),
  common_support_combined_summary %>%
    transmute(
      analysis,
      estimand = "Common-support stage-by-context interaction",
      estimand_definition,
      direction = "positive = the Non-partner stage change is more negative",
      pairs_in_summary = pairs_with_direct_support,
      directly_supported_pairs = pairs_with_direct_support,
      metrics_in_summary = 4L,
      contains_model_extrapolation = FALSE,
      estimate_median,
      lower_95,
      upper_95,
      posterior_probability_for_stated_direction = Pr_effect_gt_0,
      reported_direction = "interaction > 0"
    )
)

write_csv(
  primary_results,
  file.path(table_directory, "manuscript_primary_results_reproduced.csv")
)

# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------

wrap_metric_label <- function(x, width = 17) {
  stringr::str_wrap(x, width = width)
}

make_metric_stage_plot <- function(
    draws,
    metric_order,
    analysis_label,
    x_limits) {
  wrapped_labels <- wrap_metric_label(unname(metric_labels[metric_order]))

  plot_draws <- draws %>%
    filter(analysis == analysis_label) %>%
    mutate(
      context = factor(context, levels = c("Partner", "Non-partner")),
      metric = factor(metric, levels = metric_order),
      metric_label_wrapped = factor(
        wrap_metric_label(metric_label),
        levels = wrapped_labels
      )
    )

  plot_summary <- plot_draws %>%
    group_by(context, metric_label_wrapped) %>%
    summarise(
      estimate = median(stage_effect),
      lower = quantile(stage_effect, 0.025),
      upper = quantile(stage_effect, 0.975),
      .groups = "drop"
    )

  ggplot(plot_draws, aes(x = stage_effect, y = 0)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    tidybayes::stat_halfeye(
      aes(fill = context),
      slab_alpha = 0.6,
      .width = 0,
      point_interval = NULL
    ) +
    geom_segment(
      data = plot_summary,
      aes(x = lower, xend = upper, y = 0, yend = 0, colour = context),
      linewidth = 0.9,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = plot_summary,
      aes(x = estimate, y = 0, colour = context),
      size = 2.1,
      inherit.aes = FALSE
    ) +
    facet_grid(context ~ metric_label_wrapped, switch = "y") +
    scale_y_continuous(NULL, breaks = NULL) +
    scale_fill_manual(
      values = c(
        "Partner" = PARTNER_COLOUR,
        "Non-partner" = NONPARTNER_COLOUR
      ),
      drop = FALSE
    ) +
    scale_colour_manual(
      values = c(
        "Partner" = PARTNER_COLOUR,
        "Non-partner" = NONPARTNER_COLOUR
      ),
      drop = FALSE
    ) +
    coord_cartesian(xlim = x_limits) +
    labs(
      title = analysis_label,
      x = "Posterior stage effect: after - before (metric-specific SD)",
      y = NULL,
      fill = NULL,
      colour = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(size = 10, lineheight = 0.95),
      panel.spacing = grid::unit(1, "lines"),
      legend.position = "bottom"
    )
}

shared_limit <- as.numeric(
  quantile(abs(stage_draws$stage_effect), 0.995, names = FALSE)
) * 1.05
shared_limits <- c(-shared_limit, shared_limit)

credible_symmetric_limits <- function(
    summary_data,
    padding = 1.10,
    rounding_step = 0.5,
    minimum_limit = 1.0) {
  interval_extent <- max(
    abs(summary_data$lower_95),
    abs(summary_data$upper_95),
    na.rm = TRUE
  )
  limit <- ceiling(interval_extent * padding / rounding_step) * rounding_step
  c(-max(limit, minimum_limit), max(limit, minimum_limit))
}

combined_limits <- credible_symmetric_limits(average_stage_summary)

sequence_metric_figure <- make_metric_stage_plot(
  stage_draws,
  SEQUENCE_METRICS,
  "Sequence structure",
  shared_limits
)
acoustic_metric_figure <- make_metric_stage_plot(
  stage_draws,
  ACOUSTIC_METRICS,
  "Call structure",
  shared_limits
)

ggsave(
  file.path(
    figure_directory,
    "Fig_sequence_stage_by_context_metric_reproduced.png"
  ),
  sequence_metric_figure,
  width = 12,
  height = 4.8,
  dpi = 300
)
ggsave(
  file.path(
    figure_directory,
    "Fig_acoustic_stage_by_context_metric_reproduced.png"
  ),
  acoustic_metric_figure,
  width = 12,
  height = 4.8,
  dpi = 300
)

average_plot_draws <- average_stage_draws %>%
  mutate(
    context = factor(context, levels = c("Partner", "Non-partner")),
    analysis = factor(
      analysis,
      levels = c("Sequence structure", "Call structure")
    ),
    label = paste(context, analysis, sep = "\n")
  )

average_plot_summary <- average_plot_draws %>%
  group_by(context, analysis, label) %>%
  summarise(
    estimate = median(stage_effect),
    lower = quantile(stage_effect, 0.025),
    upper = quantile(stage_effect, 0.975),
    .groups = "drop"
  )

average_label_order <- c(
  "Partner\nSequence structure",
  "Non-partner\nSequence structure",
  "Partner\nCall structure",
  "Non-partner\nCall structure"
)

average_figure <- ggplot(
  average_plot_draws,
  aes(x = stage_effect, y = factor(label, levels = average_label_order))
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  tidybayes::stat_halfeye(
    aes(fill = context),
    slab_alpha = 0.55,
    .width = 0,
    point_interval = NULL
  ) +
  geom_segment(
    data = average_plot_summary,
    aes(
      x = lower,
      xend = upper,
      y = factor(label, levels = average_label_order),
      yend = factor(label, levels = average_label_order),
      colour = context
    ),
    linewidth = 0.9,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = average_plot_summary,
    aes(
      x = estimate,
      y = factor(label, levels = average_label_order),
      colour = context
    ),
    shape = 16,
    size = 2.2,
    inherit.aes = FALSE
  ) +
  scale_y_discrete(limits = rev(average_label_order)) +
  scale_fill_manual(
    values = c(
      "Partner" = PARTNER_COLOUR,
      "Non-partner" = NONPARTNER_COLOUR
    )
  ) +
  scale_colour_manual(
    values = c(
      "Partner" = PARTNER_COLOUR,
      "Non-partner" = NONPARTNER_COLOUR
    )
  ) +
  coord_cartesian(xlim = shared_limits) +
  labs(
    title = "Equal-weight average across metrics",
    x = "Posterior stage effect: after - before",
    y = NULL,
    fill = NULL,
    colour = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(axis.text.y = element_text(hjust = 0.5), legend.position = "bottom")

ggsave(
  file.path(figure_directory, "Fig_average_stage_posterior_reproduced.png"),
  average_figure,
  width = 8,
  height = 4.8,
  dpi = 300
)

make_combined_analysis_plot <- function(
    draws,
    combined_draws,
    metric_order,
    analysis_label,
    show_x_title,
    metric_x_limits,
    combined_x_limits,
    text_scale = 1) {
  pad_title <- function(x, number_of_lines = 2L) {
    vapply(x, function(label) {
      current_lines <- stringr::str_count(label, fixed("\n")) + 1L
      if (current_lines < number_of_lines) {
        paste0(label, str_dup("\n ", number_of_lines - current_lines))
      } else {
        label
      }
    }, character(1))
  }

  metric_titles <- pad_title(unname(panel_metric_labels[metric_order]))
  combined_title <- pad_title("Combined result")
  row_levels <- c("Partner", "Non-partner")

  metric_plot_draws <- draws %>%
    filter(analysis == analysis_label, metric %in% metric_order) %>%
    transmute(
      .draw,
      analysis,
      context,
      stage_effect,
      panel_label = pad_title(
        unname(panel_metric_labels[as.character(metric)])
      )
    ) %>%
    mutate(
      context = factor(context, levels = row_levels),
      row_label = factor(context, levels = row_levels),
      panel_label = factor(panel_label, levels = metric_titles)
    )

  combined_plot_draws <- combined_draws %>%
    filter(analysis == analysis_label) %>%
    transmute(
      .draw,
      analysis,
      context,
      stage_effect,
      panel_label = combined_title
    ) %>%
    mutate(
      context = factor(context, levels = row_levels),
      row_label = factor(context, levels = row_levels),
      panel_label = factor(panel_label, levels = combined_title)
    )

  summarise_for_plot <- function(plot_draws) {
    plot_draws %>%
      group_by(analysis, context, row_label, panel_label) %>%
      summarise(
        estimate = median(stage_effect),
        lower = quantile(stage_effect, 0.025),
        upper = quantile(stage_effect, 0.975),
        .groups = "drop"
      )
  }

  metric_plot_summary <- summarise_for_plot(metric_plot_draws)
  combined_plot_summary <- summarise_for_plot(combined_plot_draws)

  direction_zones <- expand_grid(
    row_label = factor(row_levels, levels = row_levels),
    direction = c("Convergence", "Divergence")
  ) %>%
    mutate(
      panel_label = factor(combined_title, levels = combined_title),
      xmin = if_else(direction == "Convergence", combined_x_limits[[1]], 0),
      xmax = if_else(direction == "Convergence", 0, combined_x_limits[[2]])
    )

  direction_labels <- tibble(
    row_label = factor(row_levels[[1]], levels = row_levels),
    panel_label = factor(combined_title, levels = combined_title),
    direction = c("Convergence", "Divergence"),
    x = c(
      combined_x_limits[[1]] + diff(combined_x_limits) * 0.02,
      combined_x_limits[[2]] - diff(combined_x_limits) * 0.02
    ),
    hjust = c(0, 1)
  )

  make_distribution_panel <- function(
      plot_draws,
      plot_summary,
      shade_directions,
      show_row_labels,
      panel_x_limits) {
    plot <- ggplot(plot_draws, aes(x = stage_effect, y = 0))

    if (shade_directions) {
      plot <- plot +
        geom_rect(
          data = filter(direction_zones, direction == "Convergence"),
          aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
          fill = "#dceff2",
          alpha = 0.58,
          inherit.aes = FALSE
        ) +
        geom_rect(
          data = filter(direction_zones, direction == "Divergence"),
          aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
          fill = "#f6e5d8",
          alpha = 0.58,
          inherit.aes = FALSE
        )
    }

    plot <- plot +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.55,
        colour = "grey50"
      ) +
      tidybayes::stat_halfeye(
        aes(fill = context),
        slab_alpha = 0.66,
        .width = 0,
        point_interval = NULL
      ) +
      geom_segment(
        data = plot_summary,
        aes(x = lower, xend = upper, y = 0, yend = 0, colour = context),
        linewidth = 1.05,
        inherit.aes = FALSE
      ) +
      geom_point(
        data = plot_summary,
        aes(x = estimate, y = 0, colour = context),
        shape = 16,
        size = 2.8,
        inherit.aes = FALSE
      )

    if (shade_directions) {
      plot <- plot +
        geom_text(
          data = filter(direction_labels, direction == "Convergence"),
          aes(x = x, y = Inf, label = direction, hjust = hjust),
          colour = "#326d78",
          fontface = "bold",
          size = 3.2 * text_scale,
          vjust = 1.15,
          inherit.aes = FALSE
        ) +
        geom_text(
          data = filter(direction_labels, direction == "Divergence"),
          aes(x = x, y = Inf, label = direction, hjust = hjust),
          colour = "#955133",
          fontface = "bold",
          size = 3.2 * text_scale,
          vjust = 1.15,
          inherit.aes = FALSE
        )
    }

    plot +
      facet_grid(row_label ~ panel_label, switch = "y") +
      scale_y_continuous(breaks = NULL) +
      scale_x_continuous(breaks = function(x) pretty(x, n = 3)) +
      scale_fill_manual(
        values = c(
          "Partner" = PARTNER_COLOUR,
          "Non-partner" = NONPARTNER_COLOUR
        ),
        name = NULL,
        drop = FALSE
      ) +
      scale_colour_manual(
        values = c(
          "Partner" = PARTNER_COLOUR,
          "Non-partner" = NONPARTNER_COLOUR
        ),
        guide = "none",
        drop = FALSE
      ) +
      coord_cartesian(xlim = panel_x_limits) +
      labs(
        x = if (show_x_title) {
          "Posterior stage effect (after - before)"
        } else {
          NULL
        },
        y = if (show_row_labels) analysis_label else NULL
      ) +
      guides(fill = guide_legend(order = 1, override.aes = list(alpha = 0.8))) +
      theme_classic(base_size = 16 * text_scale) +
      theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(
          size = if (shade_directions) {
            13.2 * text_scale
          } else {
            9.8 * text_scale
          },
          lineheight = 0.94
        ),
        strip.text.y.left = if (show_row_labels) {
          element_text(
            size = 13.2 * text_scale,
            angle = 0,
            hjust = 0.5,
            lineheight = 0.95
          )
        } else {
          element_blank()
        },
        panel.spacing.x = grid::unit(
          if (shade_directions) 0.65 else 0.35,
          "lines"
        ),
        panel.spacing.y = grid::unit(0.48, "lines"),
        axis.title.x = element_text(
          size = 16 * text_scale,
          margin = margin(t = 8)
        ),
        axis.title.y = element_text(
          size = 15 * text_scale,
          angle = 90,
          margin = margin(r = 9)
        ),
        axis.text.x = element_text(size = 12.5 * text_scale),
        axis.line = element_line(linewidth = 0.55),
        axis.ticks = element_line(linewidth = 0.55),
        legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(size = 13.2 * text_scale),
        legend.key.height = grid::unit(1.15, "lines"),
        legend.spacing.y = grid::unit(0.35, "lines"),
        plot.margin = margin(5, 4, 3, 4)
      )
  }

  combined_plot <- make_distribution_panel(
    combined_plot_draws,
    combined_plot_summary,
    shade_directions = TRUE,
    show_row_labels = TRUE,
    panel_x_limits = combined_x_limits
  )

  metric_plot <- make_distribution_panel(
    metric_plot_draws,
    metric_plot_summary,
    shade_directions = FALSE,
    show_row_labels = FALSE,
    panel_x_limits = metric_x_limits
  )

  (combined_plot | metric_plot) +
    patchwork::plot_layout(
      widths = c(1.45, 4.6),
      guides = "auto",
      axis_titles = "collect_x"
    )
}

make_full_results_figure <- function(
    stage_draw_data,
    combined_draw_data,
    metric_x_limits,
    combined_x_limits,
    text_scale = 1.45) {
  sequence_panel <- make_combined_analysis_plot(
    stage_draw_data,
    combined_draw_data,
    SEQUENCE_METRICS,
    "Sequence structure",
    show_x_title = FALSE,
    metric_x_limits = metric_x_limits,
    combined_x_limits = combined_x_limits,
    text_scale = text_scale
  )

  acoustic_panel <- make_combined_analysis_plot(
    stage_draw_data,
    combined_draw_data,
    ACOUSTIC_METRICS,
    "Call structure",
    show_x_title = TRUE,
    metric_x_limits = metric_x_limits,
    combined_x_limits = combined_x_limits,
    text_scale = text_scale
  )

  analysis_grid <- (sequence_panel / acoustic_panel) +
    patchwork::plot_layout(
      heights = c(1, 1),
      guides = "collect",
      axis_titles = "collect_x"
    ) &
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.justification = c(0, 1),
      legend.box = "vertical",
      legend.box.just = "top",
      legend.text.align = 0,
      legend.margin = margin(5, 4, 5, 8)
    )

  model_axis_label <- ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = "Fitted models",
      angle = 90,
      size = 5.7 * text_scale
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))

  (model_axis_label | analysis_grid) +
    patchwork::plot_layout(widths = c(0.045, 1))
}

combined_figure <- make_full_results_figure(
  stage_draws,
  average_stage_draws,
  metric_x_limits = shared_limits,
  combined_x_limits = combined_limits
)

ggsave(
  file.path(
    figure_directory,
    "Fig_combined_interaction_panel_reproduced.png"
  ),
  combined_figure,
  width = 16.5,
  height = 8.5,
  dpi = 300
)

# -----------------------------------------------------------------------------
# Diagnostics from the saved posterior samples
# -----------------------------------------------------------------------------

finite_statistic <- function(values, summary_function) {
  values <- values[is.finite(values)]
  if (!length(values)) {
    return(NA_real_)
  }
  match.fun(summary_function)(values)
}

model_diagnostics <- function(model, model_label) {
  draw_summary <- posterior::summarise_draws(
    posterior::as_draws_array(model)
  )
  sampler_diagnostics <- brms::nuts_params(model)

  tibble(
    model = model_label,
    max_rhat = finite_statistic(draw_summary$rhat, "max"),
    min_ess_bulk = finite_statistic(draw_summary$ess_bulk, "min"),
    min_ess_tail = finite_statistic(draw_summary$ess_tail, "min"),
    divergences = sum(
      sampler_diagnostics$Parameter == "divergent__" &
        sampler_diagnostics$Value == 1,
      na.rm = TRUE
    ),
    max_treedepth_hits = sum(
      sampler_diagnostics$Parameter == "treedepth__" &
        sampler_diagnostics$Value >= 15,
      na.rm = TRUE
    )
  ) %>%
    mutate(
      diagnostic_pass =
        max_rhat <= 1.01 &
        min_ess_bulk >= 400 &
        min_ess_tail >= 400 &
        divergences == 0 &
        max_treedepth_hits == 0
    )
}

diagnostic_summary <- bind_rows(
  model_diagnostics(sequence_model, "sequence_model"),
  model_diagnostics(acoustic_model, "acoustic_model")
)

write_csv(
  diagnostic_summary,
  file.path(diagnostic_directory, "diagnostics_summary_reproduced.csv")
)

save_diagnostic_plots <- function(model, output_stem) {
  set.seed(123)
  predictive_check <- brms::pp_check(model, ndraws = 200)
  ggsave(
    file.path(
      diagnostic_directory,
      paste0(output_stem, "_pp_check_reproduced.png")
    ),
    predictive_check,
    width = 9,
    height = 6,
    dpi = 180
  )

  draw_array <- posterior::as_draws_array(model)
  parameter_names <- posterior::variables(draw_array)
  trace_parameters <- unique(c(
    grep("^b_", parameter_names, value = TRUE),
    grep("^sd_", parameter_names, value = TRUE),
    "sigma"
  ))
  trace_parameters <- intersect(trace_parameters, parameter_names)

  trace_plot <- bayesplot::mcmc_trace(
    draw_array,
    pars = trace_parameters
  )
  ggsave(
    file.path(
      diagnostic_directory,
      paste0(output_stem, "_trace_reproduced.png")
    ),
    trace_plot,
    width = 14,
    height = max(8, length(trace_parameters) * 0.55),
    dpi = 160,
    limitsize = FALSE
  )
}

save_diagnostic_plots(sequence_model, "sequence_model")
save_diagnostic_plots(acoustic_model, "acoustic_model")

if (!all(diagnostic_summary$diagnostic_pass)) {
  warning(
    "At least one saved model did not pass the diagnostic thresholds. See ",
    file.path(diagnostic_directory, "diagnostics_summary_reproduced.csv")
  )
}

message("\nScript 4 completed successfully.")
message("Model fitting requested: ", REFIT_MODELS)
message(
  "Primary result table: ",
  file.path(table_directory, "manuscript_primary_results_reproduced.csv")
)
message(
  "Primary figure: ",
  file.path(
    figure_directory,
    "Fig_combined_interaction_panel_reproduced.png"
  )
)
