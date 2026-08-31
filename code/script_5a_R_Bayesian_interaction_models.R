#!/usr/bin/env Rscript

# ============================================================
# SCRIPT 5a: BAYESIAN INTERACTION MODELS
# ============================================================
#
# This script is the final step of the analysis. It starts from the two
# model-ready tables made during Script 4a:
#
#   glm data/glm_data_seq_interaction.csv
#   glm data/glm_data_call_interaction.csv
#
# Those tables are treated as precomputed inputs in this public repository.
# The default command validates them and uses saved brms models if they are
# available. It NEVER starts Stan sampling.
#
# Safe default (no sampling):
#   Rscript code/script_5a_R_Bayesian_interaction_models.R
#
# Fit a model only when its saved model is missing:
#   Rscript code/script_5a_R_Bayesian_interaction_models.R --fit
#
# Deliberately refit both models from scratch:
#   Rscript code/script_5a_R_Bayesian_interaction_models.R --refit
#
# The model asks whether the before-to-after change in distance differs
# between Partner and Non-partner comparisons, and whether that pattern
# differs among distance metrics.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(tidybayes)
  library(patchwork)
})


# 1. Paths and command-line options ---------------------------------------

get_script_dir <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)

  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))))
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

args <- commandArgs(trailingOnly = TRUE)
allowed_args <- c("--fit", "--refit")
unknown_args <- setdiff(args, allowed_args)

if (length(unknown_args)) {
  stop(
    "Unknown option(s): ", paste(unknown_args, collapse = ", "), "\n",
    "Use no option for the safe cached-model mode, --fit to fit missing ",
    "models, or --refit to refit both models."
  )
}

fit_missing <- "--fit" %in% args
force_refit <- "--refit" %in% args

if (fit_missing && force_refit) {
  stop("Choose either --fit or --refit, not both.")
}

script_dir <- get_script_dir()
repo_dir <- normalizePath(file.path(script_dir, ".."))
data_dir <- file.path(repo_dir, "glm data")

sequence_csv <- file.path(data_dir, "glm_data_seq_interaction.csv")
call_csv <- file.path(data_dir, "glm_data_call_interaction.csv")

output_dir <- file.path(data_dir, "brms_interaction_outputs")
model_dir <- file.path(output_dir, "models")
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")

sequence_model_name <- "sequence_interaction_model"
call_model_name <- "call_structure_interaction_model"


# 2. Metric names and model specification --------------------------------

sequence_metrics <- c(
  "transition_matrix",
  "bigram",
  "repeat_A_len",
  "local_alignment"
)

call_metrics <- c(
  "T_PC_pairwise",
  "M_PC_pairwise",
  "spectro_DTW"
)

metric_labels <- c(
  "transition_matrix" = "Transition probabilities",
  "bigram" = "Bigram distribution",
  "repeat_A_len" = "Phee-repeat distribution",
  "local_alignment" = "Local alignment",
  "T_PC_pairwise" = "Spectro-temporal parameters (STP)",
  "M_PC_pairwise" = "Mel-frequency cepstral coefficients (MFCC)",
  "spectro_DTW" = "Dynamic time warping (DTW)"
)

partner_colour <- "#2c728e"
nonpartner_colour <- "#b8de29"

# stage, context, and metric are factors in the data. With "before" and
# "Non-partner" as their reference levels, the interaction coefficients
# describe how the before-to-after change varies across the design.
model_formula <- bf(
  z_distance ~ stage * context * metric +
    (1 | pair) +
    (1 | mm(session_1_id, session_2_id))
)

model_priors <- c(
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)

model_control <- list(adapt_delta = 0.995, max_treedepth = 15)
available_cores <- parallel::detectCores(logical = FALSE)
if (is.na(available_cores)) available_cores <- 1L
model_cores <- max(1L, min(4L, available_cores))


# 3. Read and check the precomputed Script 4a tables ----------------------

require_columns <- function(data, columns, label) {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns)) {
    stop(
      label, " is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

read_model_table <- function(path, expected_metrics, label) {
  if (!file.exists(path)) {
    stop(
      "Missing precomputed Script 4a input:\n  ", path, "\n\n",
      "Script 5a does not rebuild the pairwise distance tables. Add the ",
      "precomputed file before running this script."
    )
  }

  data <- read_csv(path, show_col_types = FALSE)
  required_columns <- c(
    "pair", "context", "stage", "session_1_id", "session_2_id",
    "comparison_id", "metric", "distance", "z_distance"
  )
  require_columns(data, required_columns, label)

  if (anyNA(data[required_columns])) {
    missing_counts <- colSums(is.na(data[required_columns]))
    details <- missing_counts[missing_counts > 0]
    stop(
      label, " contains missing required values:\n",
      paste(names(details), details, sep = " = ", collapse = "\n")
    )
  }

  if (!setequal(unique(data$stage), c("before", "after"))) {
    stop(label, " must contain both before and after stages.")
  }
  if (!setequal(unique(data$context), c("Non-partner", "Partner"))) {
    stop(label, " must contain Partner and Non-partner contexts.")
  }
  if (!setequal(unique(data$metric), expected_metrics)) {
    stop(
      label, " has unexpected metrics. Found: ",
      paste(sort(unique(data$metric)), collapse = ", ")
    )
  }
  if (any(data$session_1_id == data$session_2_id)) {
    stop(label, " contains a self-comparison.")
  }
  if (anyDuplicated(data[c("comparison_id", "metric")])) {
    stop(label, " contains duplicate comparison_id x metric rows.")
  }

  metrics_per_comparison <- data %>%
    distinct(comparison_id, metric) %>%
    count(comparison_id, name = "n_metrics")
  if (any(metrics_per_comparison$n_metrics != length(expected_metrics))) {
    stop(
      label, " has comparisons without exactly ",
      length(expected_metrics), " metric rows."
    )
  }

  # Script 4a standardizes distances separately for each metric. The
  # population SD is used here because that is how Script 4a calculated z.
  scaling_check <- data %>%
    group_by(metric) %>%
    summarise(
      z_mean = mean(z_distance),
      z_sd = sd(z_distance) * sqrt((n() - 1) / n()),
      .groups = "drop"
    )
  if (any(abs(scaling_check$z_mean) > 1e-6) ||
      any(abs(scaling_check$z_sd - 1) > 1e-6)) {
    stop(label, " has an unexpected z_distance standardization.")
  }

  # Both columns in the multi-membership term must use the same factor levels.
  session_levels <- union(
    as.character(data$session_1_id),
    as.character(data$session_2_id)
  )

  data %>%
    mutate(
      stage = factor(stage, levels = c("before", "after")),
      context = factor(context, levels = c("Non-partner", "Partner")),
      metric = factor(metric, levels = expected_metrics),
      pair = factor(pair),
      comparison_id = factor(comparison_id),
      session_1_id = factor(session_1_id, levels = session_levels),
      session_2_id = factor(session_2_id, levels = session_levels)
    )
}

cat("Reading precomputed Script 4a inputs...\n")
sequence_data <- read_model_table(
  sequence_csv,
  sequence_metrics,
  "Sequence table"
)
call_data <- read_model_table(
  call_csv,
  call_metrics,
  "Call-structure table"
)

cat("  Sequence rows:      ", nrow(sequence_data), "\n", sep = "")
cat("  Call-structure rows: ", nrow(call_data), "\n", sep = "")


# 4. Load cached models, or fit only with explicit permission -------------

model_path <- function(model_name) {
  file.path(model_dir, paste0(model_name, ".rds"))
}

same_values <- function(cached, current) {
  if (is.numeric(current)) {
    return(isTRUE(all.equal(
      as.numeric(cached), as.numeric(current),
      tolerance = 1e-10, check.attributes = FALSE
    )))
  }
  identical(as.character(cached), as.character(current))
}

validate_cached_model <- function(fit, current_data, label) {
  if (!inherits(fit, "brmsfit")) {
    stop(label, " is not a brmsfit object.")
  }
  if (is.null(fit$data)) {
    stop(label, " does not contain the data used for fitting.")
  }

  comparison_columns <- c(
    "z_distance", "stage", "context", "metric", "pair",
    "session_1_id", "session_2_id"
  )
  require_columns(fit$data, comparison_columns, paste(label, "cached data"))

  if (nrow(fit$data) != nrow(current_data)) {
    stop(
      label, " was fitted to ", nrow(fit$data), " rows, but the current ",
      "input has ", nrow(current_data), " rows. Use --refit deliberately ",
      "if the input table has changed."
    )
  }

  changed_columns <- comparison_columns[!vapply(
    comparison_columns,
    function(column) same_values(fit$data[[column]], current_data[[column]]),
    logical(1)
  )]
  if (length(changed_columns)) {
    stop(
      label, " does not match the current input in: ",
      paste(changed_columns, collapse = ", "),
      ". Use --refit deliberately if the input table has changed."
    )
  }

  cached_formula <- gsub(
    "\\s+", "",
    paste(deparse(fit$formula$formula), collapse = "")
  )
  expected_formula <- gsub(
    "\\s+", "",
    paste(deparse(model_formula$formula), collapse = "")
  )
  if (!identical(cached_formula, expected_formula)) {
    stop(
      label, " was fitted with a different formula. ",
      "Use --refit to replace it."
    )
  }

  invisible(TRUE)
}

load_cached_model <- function(model_name, data, label) {
  path <- model_path(model_name)
  if (!file.exists(path)) return(NULL)

  cat("Loading cached ", label, ":\n  ", path, "\n", sep = "")
  fit <- tryCatch(
    readRDS(path),
    error = function(error) {
      stop("Could not read ", path, ": ", conditionMessage(error))
    }
  )
  validate_cached_model(fit, data, label)
  fit
}

fit_model <- function(data, model_name, label, refit = FALSE) {
  dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)
  cat("\nFitting ", label, "...\n", sep = "")

  fit <- brm(
    formula = model_formula,
    data = data,
    family = gaussian(),
    prior = model_priors,
    chains = 4,
    cores = model_cores,
    iter = 4000,
    warmup = 2000,
    seed = 123,
    control = model_control,
    backend = "rstan",
    save_pars = save_pars(all = TRUE),
    file = file.path(model_dir, model_name),
    file_refit = if (refit) "always" else "never",
    refresh = 100
  )
  validate_cached_model(fit, data, label)
  fit
}

get_model <- function(data, model_name, label) {
  if (force_refit) {
    return(fit_model(data, model_name, label, refit = TRUE))
  }

  cached <- load_cached_model(model_name, data, label)
  if (!is.null(cached)) return(cached)

  if (fit_missing) {
    return(fit_model(data, model_name, label, refit = FALSE))
  }

  NULL
}

sequence_fit <- get_model(
  sequence_data,
  sequence_model_name,
  "sequence model"
)
call_fit <- get_model(
  call_data,
  call_model_name,
  "call-structure model"
)

if (is.null(sequence_fit) || is.null(call_fit)) {
  missing_models <- c(
    if (is.null(sequence_fit)) model_path(sequence_model_name),
    if (is.null(call_fit)) model_path(call_model_name)
  )

  cat("\nInput validation succeeded, but the following cached model(s) are missing:\n")
  cat(paste0("  - ", missing_models, collapse = "\n"), "\n")
  cat(
    "\nNo sampling was started. Use --fit only if you intend to fit the ",
    "missing model(s).\n",
    sep = ""
  )
  quit(save = "no", status = 0)
}


# 5. Calculate before-to-after posterior contrasts ------------------------

extract_stage_effects <- function(fit, data, analysis_label) {
  prediction_grid <- expand_grid(
    context = levels(data$context),
    metric = levels(data$metric),
    stage = levels(data$stage)
  ) %>%
    mutate(
      context = factor(context, levels = levels(data$context)),
      metric = factor(metric, levels = levels(data$metric)),
      stage = factor(stage, levels = levels(data$stage)),
      grid_row = as.character(row_number())
    )

  # re_formula = NA gives the population-level comparison. Pair and session
  # deviations are not assigned to the hypothetical rows in this grid.
  expected <- posterior_epred(
    fit,
    newdata = prediction_grid,
    re_formula = NA
  )

  expected_table <- as_tibble(as.data.frame(expected))
  names(expected_table) <- as.character(seq_len(ncol(expected_table)))

  expected_table %>%
    mutate(.draw = row_number()) %>%
    pivot_longer(
      cols = -.draw,
      names_to = "grid_row",
      values_to = "expected_z_distance"
    ) %>%
    left_join(prediction_grid, by = "grid_row") %>%
    group_by(.draw, context, metric) %>%
    summarise(
      stage_effect =
        expected_z_distance[stage == "after"] -
        expected_z_distance[stage == "before"],
      .groups = "drop"
    ) %>%
    mutate(
      analysis = analysis_label,
      context = as.character(context),
      metric = as.character(metric),
      metric_label = recode(metric, !!!metric_labels)
    )
}

stage_draws <- bind_rows(
  extract_stage_effects(sequence_fit, sequence_data, "Sequence structure"),
  extract_stage_effects(call_fit, call_data, "Call structure")
)

summarise_draws <- function(data, value, groups) {
  data %>%
    group_by(across(all_of(groups))) %>%
    summarise(
      estimate_median = median(.data[[value]]),
      estimate_mean = mean(.data[[value]]),
      lower_95 = quantile(.data[[value]], 0.025),
      upper_95 = quantile(.data[[value]], 0.975),
      Pr_effect_lt_0 = mean(.data[[value]] < 0),
      p_two_sided = 2 * min(Pr_effect_lt_0, 1 - Pr_effect_lt_0),
      .groups = "drop"
    )
}

stage_summary <- summarise_draws(
  stage_draws,
  "stage_effect",
  c("analysis", "context", "metric", "metric_label")
)

# This is the interaction in a directly interpretable form:
# (Partner after-before) - (Non-partner after-before).
context_contrast_draws <- stage_draws %>%
  select(.draw, analysis, metric, metric_label, context, stage_effect) %>%
  pivot_wider(names_from = context, values_from = stage_effect) %>%
  mutate(context_contrast = Partner - `Non-partner`)

context_contrast_summary <- summarise_draws(
  context_contrast_draws,
  "context_contrast",
  c("analysis", "metric", "metric_label")
)

# An equal-weight mean gives every metric the same influence. It is a summary
# of the posterior contrasts, not a separately fitted pooled model.
average_stage_draws <- stage_draws %>%
  group_by(.draw, analysis, context) %>%
  summarise(stage_effect = mean(stage_effect), .groups = "drop")

average_stage_summary <- summarise_draws(
  average_stage_draws,
  "stage_effect",
  c("analysis", "context")
)

dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(stage_summary, file.path(table_dir, "stage_effects_by_context_metric.csv"))
write_csv(
  context_contrast_summary,
  file.path(table_dir, "partner_minus_nonpartner_by_metric.csv")
)
write_csv(
  average_stage_summary,
  file.path(table_dir, "average_stage_effects.csv")
)


# 6. Plot the posterior contrasts ----------------------------------------

make_metric_plot <- function(draws, metrics, analysis_label, x_limits) {
  label_order <- stringr::str_wrap(unname(metric_labels[metrics]), width = 18)

  plot_draws <- draws %>%
    filter(analysis == analysis_label) %>%
    mutate(
      context = factor(context, levels = c("Partner", "Non-partner")),
      metric_label = factor(
        stringr::str_wrap(metric_label, width = 18),
        levels = label_order
      )
    )

  plot_summary <- plot_draws %>%
    group_by(context, metric_label) %>%
    summarise(
      estimate = median(stage_effect),
      lower = quantile(stage_effect, 0.025),
      upper = quantile(stage_effect, 0.975),
      .groups = "drop"
    )

  ggplot(plot_draws, aes(x = stage_effect, y = 0)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    stat_halfeye(
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
    facet_grid(context ~ metric_label, switch = "y") +
    scale_y_continuous(NULL, breaks = NULL) +
    scale_fill_manual(
      values = c("Partner" = partner_colour, "Non-partner" = nonpartner_colour)
    ) +
    scale_colour_manual(
      values = c("Partner" = partner_colour, "Non-partner" = nonpartner_colour)
    ) +
    coord_cartesian(xlim = x_limits) +
    labs(
      title = analysis_label,
      x = "Posterior stage effect: after - before (metric-specific SD)",
      fill = NULL,
      colour = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(size = 10),
      legend.position = "bottom"
    )
}

shared_limit <- as.numeric(
  quantile(abs(stage_draws$stage_effect), 0.995, names = FALSE)
) * 1.05
shared_limits <- c(-shared_limit, shared_limit)

sequence_plot <- make_metric_plot(
  stage_draws,
  sequence_metrics,
  "Sequence structure",
  shared_limits
)
call_plot <- make_metric_plot(
  stage_draws,
  call_metrics,
  "Call structure",
  shared_limits
)

average_plot_data <- average_stage_draws %>%
  mutate(
    context = factor(context, levels = c("Partner", "Non-partner")),
    analysis = factor(
      analysis,
      levels = c("Sequence structure", "Call structure")
    ),
    label = paste(context, analysis, sep = "\n")
  )

average_plot_summary <- average_plot_data %>%
  group_by(context, analysis, label) %>%
  summarise(
    estimate = median(stage_effect),
    lower = quantile(stage_effect, 0.025),
    upper = quantile(stage_effect, 0.975),
    .groups = "drop"
  )

average_order <- c(
  "Partner\nSequence structure",
  "Non-partner\nSequence structure",
  "Partner\nCall structure",
  "Non-partner\nCall structure"
)

average_plot <- ggplot(
  average_plot_data,
  aes(x = stage_effect, y = factor(label, levels = average_order))
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  stat_halfeye(
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
      y = factor(label, levels = average_order),
      yend = factor(label, levels = average_order),
      colour = context
    ),
    linewidth = 0.9,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = average_plot_summary,
    aes(
      x = estimate,
      y = factor(label, levels = average_order),
      colour = context,
      shape = analysis
    ),
    size = 2.2,
    inherit.aes = FALSE
  ) +
  scale_y_discrete(limits = rev(average_order)) +
  scale_fill_manual(
    values = c("Partner" = partner_colour, "Non-partner" = nonpartner_colour)
  ) +
  scale_colour_manual(
    values = c("Partner" = partner_colour, "Non-partner" = nonpartner_colour)
  ) +
  coord_cartesian(xlim = shared_limits) +
  labs(
    title = "Equal-weight average across metrics",
    x = "Posterior stage effect: after - before",
    y = NULL,
    fill = NULL,
    colour = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom")

combined_plot <- (
  average_plot + theme(legend.position = "none")
) | (
  (sequence_plot + theme(legend.position = "none")) /
    (call_plot + theme(legend.position = "none"))
)
combined_plot <- combined_plot +
  plot_layout(widths = c(0.42, 1)) +
  plot_annotation(tag_levels = "a")

dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(
  file.path(figure_dir, "sequence_stage_by_context_metric.png"),
  sequence_plot,
  width = 12,
  height = 4.8,
  dpi = 300
)
ggsave(
  file.path(figure_dir, "call_stage_by_context_metric.png"),
  call_plot,
  width = 12,
  height = 4.8,
  dpi = 300
)
ggsave(
  file.path(figure_dir, "combined_interaction_panel.png"),
  combined_plot,
  width = 17,
  height = 9,
  dpi = 300
)

cat("\nScript 5a finished.\n")
cat("Tables:  ", table_dir, "\n", sep = "")
cat("Figures: ", figure_dir, "\n", sep = "")
