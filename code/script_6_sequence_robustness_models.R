#!/usr/bin/env Rscript

# ============================================================
# SCRIPT 6: NON-PARTNER SEQUENCE ROBUSTNESS MODELS
#
# This reviewer-directed sensitivity analysis is intentionally restricted to
# the outcome/context for which the primary Bayesian analysis (Script 4)
# found evidence of convergence:
# Non-partner sequence structure.
#
# Model 1 asks whether repertoire distances are greater across stages than
# within stages. It compares three cells separately:
#   * Before vs Before sessions
#   * After vs After sessions
#   * Before vs After sessions
# The reported contrast is:
#   Across - 0.5 * (Within-before + Within-after)
# Positive values mean that the across-stage distance exceeds the mean of the
# within-Before and within-After distances.
#
# Model 2 asks whether the Before-to-After decline is greater for the three
# bonded opposite-sex dyads than for the six nonbonded opposite-sex dyads.
# The reported contrast is:
#   (After - Before)_bonded - (After - Before)_nonbonded
# Negative values mean greater convergence among bonded dyads.
#
# Both models use the exact frozen 107-repertoire manifest and four 107 x 107
# sequence-distance matrices used by the primary analysis. Reused sessions
# enter through multiple-membership effects. Posterior summaries give metrics
# equal weight and Model 2 gives dyads equal weight rather than weighting them
# by Cartesian row counts.
#
# Usage:
#   Rscript code/script_6_sequence_robustness_models.R --prep-only
#   Rscript code/script_6_sequence_robustness_models.R --validate-only
#   Rscript code/script_6_sequence_robustness_models.R --refit
#
# Development fits must use both:
#   MARMOSET_QUICK_FIT=1
#   MARMOSET_6_OUTPUT_DIR=/an/explicit/noncanonical/path
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(tidybayes)
  library(posterior)
  library(patchwork)
  library(reticulate)
})

# -------------------------
# Paths and run settings
# -------------------------
get_script_dir <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg)) {
    script_path <- sub("^--file=", "", file_arg[[1]])
    script_path <- gsub("~\\+~", " ", script_path)
    return(dirname(normalizePath(script_path)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path)) {
    return(dirname(normalizePath(rstudioapi::getActiveDocumentContext()$path)))
  }
  normalizePath(getwd())
}

script_dir <- get_script_dir()
root_override <- Sys.getenv("MARMOSET_PROJECT_ROOT", unset = "")
root <- if (nzchar(root_override)) {
  normalizePath(root_override, mustWork = TRUE)
} else {
  normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
prep_only <- "--prep-only" %in% args
validate_only <- "--validate-only" %in% args
force_refit <- "--refit" %in% args
unknown_args <- setdiff(args, c("--prep-only", "--validate-only", "--refit"))
if (length(unknown_args)) {
  stop("Unknown command-line arguments: ", paste(unknown_args, collapse = ", "))
}
if (sum(c(prep_only, validate_only, force_refit)) > 1L) {
  stop("Choose only one of --prep-only, --validate-only, and --refit")
}

quick_fit <- identical(Sys.getenv("MARMOSET_QUICK_FIT", unset = "0"), "1")
if (!prep_only && !validate_only && !force_refit && !quick_fit) {
  stop(
    "No action requested. Use --validate-only to check the analysis or ",
    "--refit to sample both models."
  )
}
run_mode <- case_when(
  prep_only ~ "prep_only",
  validate_only ~ "validate_only",
  quick_fit ~ "quick_fit",
  TRUE ~ "final_fit"
)

manifest_path <- file.path(
  root, "data", "sequence", "session_order_107.csv"
)
matrix_dir <- file.path(
  root, "data", "sequence", "distances"
)
matrix_paths <- c(
  transition_matrix = file.path(
    matrix_dir, "transition_probability_107.npy"
  ),
  bigram = file.path(matrix_dir, "bigram_107.npy"),
  repeat_A_len = file.path(matrix_dir, "phee_repeat_107.npy"),
  local_alignment = file.path(matrix_dir, "local_alignment_107.npy")
)

canonical_out_dir <- normalizePath(
  file.path(root, "results", "sequence_robustness"),
  mustWork = FALSE
)
output_override <- Sys.getenv("MARMOSET_6_OUTPUT_DIR", unset = "")
if (quick_fit && !nzchar(output_override)) {
  stop(
    "MARMOSET_QUICK_FIT=1 requires an explicit MARMOSET_6_OUTPUT_DIR; ",
    "development runs may not overwrite the canonical results."
  )
}
out_dir <- if (nzchar(output_override)) {
  normalizePath(output_override, mustWork = FALSE)
} else if (prep_only) {
  paste0(canonical_out_dir, "_prep_only")
} else if (validate_only) {
  paste0(canonical_out_dir, "_validate_only")
} else {
  canonical_out_dir
}
if (run_mode != "final_fit" && identical(out_dir, canonical_out_dir)) {
  stop("A non-final run may not write to the canonical results directory")
}

model_input_dir <- file.path(out_dir, "model_inputs")
model_dir <- file.path(out_dir, "models")
table_dir <- file.path(out_dir, "tables")
figure_dir <- file.path(out_dir, "figures")
diagnostic_dir <- file.path(out_dir, "diagnostics")
walk(
  c(out_dir, model_input_dir, model_dir, table_dir, figure_dir, diagnostic_dir),
  ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE)
)

# -------------------------
# Frozen inputs and constants
# -------------------------
SEQ_METRICS <- c(
  "transition_matrix", "bigram", "repeat_A_len", "local_alignment"
)
METRIC_LABELS <- c(
  transition_matrix = "Transition probabilities",
  bigram = "Bigram distribution",
  repeat_A_len = "Phee-repeat distribution",
  local_alignment = "Local alignment"
)
COMPARISON_TYPES <- c("within_before", "within_after", "across_stage")

BONDED_PAIRS <- tribble(
  ~pair_id, ~male_id, ~female_id,
  "Tabor-Lola", "Tabor", "Lola",
  "Odin-Nougatti", "Odin", "Nougatti",
  "Wuschel-Olympia", "Wuschel", "Olympia"
)

EXPECTED_HASHES <- c(
  manifest = "b986a7ec0ff783b3a1344ed275334c22c32c24dc9f3a367a73b2b33580503daa",
  transition_matrix = "50472b5436c0a1aa9d71a964366913a866b026e69165306c99bccce0fc6d903f",
  bigram = "8ff6648e707d404b01c5e1793f8d531014711f3b58e903742f92fab329ec6a44",
  repeat_A_len = "4f915ef1ef0bf0164459999879415c5166f236af4f664b40304e67fabecabde7",
  local_alignment = "0203611c811af4e0a1630d7e94359c2b82a3e13868af02a9ca35fa26e5adc202"
)

priors_gaussian <- c(
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)

model1_formula <- bf(
  z_distance ~ comparison_type * metric +
    (1 | focal_id) +
    (1 | conspecific_id) +
    (1 | mm(session_1_id, session_2_id))
)

model2_formula <- bf(
  z_distance ~ stage * bonded_status * metric +
    (1 + stage || dyad_id) +
    (1 | mm(individual_1, individual_2)) +
    (1 | mm(session_1_id, session_2_id))
)

brms_control <- list(adapt_delta = 0.995, max_treedepth = 15)
model_chains <- if (quick_fit) 2L else 4L
model_iterations <- if (quick_fit) 600L else 4000L
model_warmup <- if (quick_fit) 300L else 2000L
available_cores <- parallel::detectCores(logical = FALSE)
if (is.na(available_cores) || available_cores < 1L) available_cores <- 4L
model_cores <- max(1L, min(model_chains, available_cores))

# -------------------------
# General helpers
# -------------------------
require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

require_files <- function(paths, label) {
  missing <- unname(paths[!file.exists(paths)])
  if (length(missing)) {
    stop(label, " is missing:\n", paste(missing, collapse = "\n"))
  }
}

sha256_file <- function(path) digest::digest(file = path, algo = "sha256")

project_relative <- function(path) {
  absolute <- normalizePath(path, mustWork = FALSE)
  root_prefix <- paste0(normalizePath(root), .Platform$file.sep)
  if (startsWith(absolute, root_prefix)) {
    substring(absolute, nchar(root_prefix) + 1L)
  } else {
    absolute
  }
}

population_sd <- function(x) sqrt(mean((x - mean(x))^2))

format_formula <- function(formula) {
  gsub("[[:space:]]+", " ", paste(deparse(formula$formula), collapse = " "))
}

clean_id_component <- function(value) {
  if (length(value) != 1L || is.na(value)) {
    stop("Identifier components may not be missing")
  }
  text <- if (is.numeric(value) && abs(value - round(value)) < 1e-8) {
    format(round(value), scientific = FALSE, trim = TRUE)
  } else {
    trimws(as.character(value))
  }
  if (!nzchar(text)) stop("Identifier components may not be empty")
  gsub("__", "_", text, fixed = TRUE)
}

make_session_ids <- function(focal, receiver, stage, session_number) {
  mapply(
    function(f, r, s, n) {
      paste(
        clean_id_component(f), clean_id_component(r),
        clean_id_component(s), clean_id_component(n), sep = "__"
      )
    },
    focal, receiver, stage, session_number,
    USE.NAMES = FALSE
  )
}

make_comparison_ids <- function(session_1_id, session_2_id) {
  map2_chr(
    session_1_id, session_2_id,
    ~ paste(sort(c(.x, .y)), collapse = "||")
  )
}

canonical_dyad_id <- function(individual_1, individual_2) {
  map2_chr(
    individual_1, individual_2,
    ~ paste(sort(c(.x, .y)), collapse = "--")
  )
}

configure_project_python <- function() {
  configured_python <- Sys.getenv("RETICULATE_PYTHON", unset = "")
  if (nzchar(configured_python) && !file.exists(configured_python)) {
    stop("RETICULATE_PYTHON does not exist: ", configured_python)
  }

  candidates <- c(
    configured_python,
    file.path(root, ".venv", "bin", "python"),
    file.path(root, ".venv", "Scripts", "python.exe")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (length(candidates)) {
    reticulate::use_python(candidates[[1]], required = TRUE)
  }
  if (!reticulate::py_module_available("numpy")) {
    stop(
      "Python package 'numpy' is required. Create the documented .venv, ",
      "or set RETICULATE_PYTHON to a Python interpreter with NumPy installed."
    )
  }
  reticulate::py_config()$python
}

numpy_module <- NULL
read_numpy_matrix <- function(path) {
  if (is.null(numpy_module)) {
    numpy_module <<- reticulate::import("numpy", convert = TRUE)
  }
  numpy_module$load(path, allow_pickle = FALSE)
}

validate_distance_matrix <- function(matrix, label) {
  if (!identical(dim(matrix), c(107L, 107L))) {
    stop(label, " must be 107 x 107")
  }
  if (any(!is.finite(matrix))) stop(label, " contains non-finite values")
  if (max(abs(matrix - t(matrix))) > 1e-6) stop(label, " is not symmetric")
  if (max(abs(diag(matrix))) > 1e-6) stop(label, " diagonal is not zero")
  invisible(TRUE)
}

add_metric_distances <- function(comparisons) {
  distances <- imap_dfr(matrix_paths, function(path, metric_name) {
    matrix <- read_numpy_matrix(path)
    validate_distance_matrix(matrix, metric_name)
    values <- matrix[cbind(
      comparisons$source_index_1 + 1L,
      comparisons$source_index_2 + 1L
    )]
    tibble(
      comparison_id = comparisons$comparison_id,
      metric = metric_name,
      distance = as.numeric(values)
    )
  })

  comparisons %>%
    left_join(distances, by = "comparison_id") %>%
    group_by(metric) %>%
    mutate(
      metric_mean = mean(distance),
      metric_population_sd = population_sd(distance),
      z_distance = (distance - metric_mean) / metric_population_sd
    ) %>%
    ungroup()
}

validate_metric_rows <- function(data, label) {
  if (anyDuplicated(data[c("comparison_id", "metric")])) {
    stop(label, " has duplicate comparison x metric rows")
  }
  per_comparison <- data %>%
    distinct(comparison_id, metric) %>%
    count(comparison_id, name = "n_metrics")
  if (any(per_comparison$n_metrics != length(SEQ_METRICS))) {
    stop(label, " does not have four metrics for every comparison")
  }
  scaling <- data %>%
    group_by(metric) %>%
    summarise(
      mean = mean(z_distance),
      population_sd = population_sd(z_distance),
      .groups = "drop"
    )
  if (any(abs(scaling$mean) > 1e-8) ||
      any(abs(scaling$population_sd - 1) > 1e-8)) {
    stop(label, " metric standardization failed")
  }
  invisible(TRUE)
}

factor_mm_columns <- function(data, first_column, second_column) {
  levels_union <- union(
    as.character(data[[first_column]]), as.character(data[[second_column]])
  )
  data[[first_column]] <- factor(data[[first_column]], levels = levels_union)
  data[[second_column]] <- factor(data[[second_column]], levels = levels_union)
  data
}

summarise_posterior <- function(data, value_column, grouping_columns) {
  data %>%
    group_by(across(all_of(grouping_columns))) %>%
    summarise(
      estimate_median = median(.data[[value_column]]),
      estimate_mean = mean(.data[[value_column]]),
      lower_95 = quantile(.data[[value_column]], 0.025),
      upper_95 = quantile(.data[[value_column]], 0.975),
      Pr_lt_0 = mean(.data[[value_column]] < 0),
      Pr_gt_0 = mean(.data[[value_column]] > 0),
      .groups = "drop"
    )
}

# -------------------------
# Verify and load frozen inputs
# -------------------------
artifact_paths <- c(manifest = manifest_path, matrix_paths)
require_files(artifact_paths, "Script 6 frozen inputs")
observed_hashes <- map_chr(artifact_paths, sha256_file)
if (any(observed_hashes != EXPECTED_HASHES[names(observed_hashes)])) {
  mismatch <- tibble(
    artifact = names(observed_hashes),
    expected = EXPECTED_HASHES[names(observed_hashes)],
    observed = observed_hashes
  ) %>% filter(expected != observed)
  stop(
    "Frozen input hash validation failed:\n",
    paste(capture.output(print(mismatch, n = Inf)), collapse = "\n")
  )
}
write_csv(
  tibble(
    artifact = names(artifact_paths),
    project_relative_path = map_chr(artifact_paths, project_relative),
    sha256 = observed_hashes
  ),
  file.path(model_input_dir, "frozen_input_hashes.csv")
)

invisible(configure_project_python())

manifest <- read_csv(
  manifest_path, show_col_types = FALSE, name_repair = "minimal"
) %>% arrange(group_id)
require_columns(
  manifest,
  c(
    "group_id", "focal ID", "conspecific_ID", "stage", "session_number",
    "paired_status", "n_seq", "pair_id", "sex"
  ),
  "107-repertoire manifest"
)
if (nrow(manifest) != 107L ||
    !identical(as.integer(manifest$group_id), 0:106)) {
  stop("Manifest must contain ordered group_id values 0:106")
}
if (anyDuplicated(manifest[c(
  "focal ID", "conspecific_ID", "stage", "session_number"
)])) {
  stop("Manifest contains duplicate focal x receiver x stage x session keys")
}

sessions <- manifest %>%
  transmute(
    focal_id = .data[["focal ID"]],
    conspecific_id = conspecific_ID,
    sex = tolower(sex),
    stage = tolower(stage),
    context = recode(
      paired_status,
      "partner" = "Partner",
      "non-partner" = "Non-partner",
      .default = NA_character_
    ),
    session_number = session_number,
    session_id = make_session_ids(
      .data[["focal ID"]], conspecific_ID, stage, session_number
    ),
    source_index = as.integer(group_id),
    n_sequences = as.integer(n_seq),
    source_pair = pair_id
  )
if (anyNA(sessions$context)) stop("Manifest has an unknown paired_status")
if (anyDuplicated(sessions$session_id)) stop("Session IDs are not unique")

sex_table <- sessions %>% distinct(focal_id, sex)
if (nrow(sex_table) != 6L ||
    sum(sex_table$sex == "male") != 3L ||
    sum(sex_table$sex == "female") != 3L) {
  stop("Expected exactly six focal animals: three male and three female")
}

nonpartner_sessions <- sessions %>% filter(context == "Non-partner")
if (nrow(nonpartner_sessions) != 74L) {
  stop("Expected 74 eligible Non-partner sequence repertoires; observed ",
       nrow(nonpartner_sessions))
}

# -------------------------
# Model 1: within- versus across-stage repertoire distances
# -------------------------
common_receiver_strata <- nonpartner_sessions %>%
  distinct(focal_id, conspecific_id, stage) %>%
  count(focal_id, conspecific_id, name = "n_stages") %>%
  filter(n_stages == 2L) %>%
  select(focal_id, conspecific_id)

model1_sessions <- nonpartner_sessions %>%
  inner_join(common_receiver_strata, by = c("focal_id", "conspecific_id"))

build_within_individual_comparisons <- function(data) {
  chunks <- list()
  chunk_index <- 1L
  strata <- data %>% distinct(focal_id, conspecific_id)
  for (row_index in seq_len(nrow(strata))) {
    focal_value <- strata$focal_id[[row_index]]
    receiver_value <- strata$conspecific_id[[row_index]]
    block <- data %>%
      filter(
        focal_id == focal_value,
        conspecific_id == receiver_value
      ) %>%
      arrange(source_index)
    if (nrow(block) < 2L) next
    indices <- t(combn(seq_len(nrow(block)), 2L))
    left <- block[indices[, 1], ]
    right <- block[indices[, 2], ]
    chunks[[chunk_index]] <- tibble(
      focal_id = focal_value,
      conspecific_id = receiver_value,
      stage_1 = left$stage,
      stage_2 = right$stage,
      comparison_type = case_when(
        stage_1 != stage_2 ~ "across_stage",
        stage_1 == "before" ~ "within_before",
        TRUE ~ "within_after"
      ),
      session_1_id = left$session_id,
      session_2_id = right$session_id,
      source_index_1 = left$source_index,
      source_index_2 = right$source_index,
      session_number_1 = left$session_number,
      session_number_2 = right$session_number,
      n_sequences_1 = left$n_sequences,
      n_sequences_2 = right$n_sequences
    ) %>%
      mutate(
        comparison_id = make_comparison_ids(session_1_id, session_2_id)
      )
    chunk_index <- chunk_index + 1L
  }
  bind_rows(chunks)
}

model1_comparisons <- build_within_individual_comparisons(model1_sessions)
if (nrow(model1_comparisons) != 208L) {
  stop("Expected 208 Model 1 session comparisons; observed ",
       nrow(model1_comparisons))
}
expected_model1_types <- c(
  within_before = 37L, within_after = 59L, across_stage = 112L
)
observed_model1_types <- table(model1_comparisons$comparison_type)
if (any(observed_model1_types[names(expected_model1_types)] !=
        expected_model1_types)) {
  stop("Model 1 comparison-type counts changed")
}
if (anyDuplicated(model1_comparisons$comparison_id)) {
  stop("Model 1 contains duplicate session comparisons")
}

model1_unfactored <- add_metric_distances(model1_comparisons)
validate_metric_rows(model1_unfactored, "Model 1")
if (nrow(model1_unfactored) != 832L) {
  stop("Expected 832 Model 1 comparison x metric rows")
}

model1_data <- model1_unfactored %>%
  mutate(
    comparison_type = factor(comparison_type, levels = COMPARISON_TYPES),
    metric = factor(metric, levels = SEQ_METRICS),
    focal_id = factor(focal_id),
    conspecific_id = factor(conspecific_id),
    comparison_id = factor(comparison_id)
  ) %>%
  factor_mm_columns("session_1_id", "session_2_id")

# -------------------------
# Model 2: bonded versus nonbonded opposite-sex dyads
# -------------------------
is_bonded_dyad <- function(male_id, female_id) {
  map2_lgl(
    male_id, female_id,
    ~ any(BONDED_PAIRS$male_id == .x & BONDED_PAIRS$female_id == .y)
  )
}

build_opposite_sex_comparisons <- function(data) {
  males <- sex_table %>% filter(sex == "male") %>% pull(focal_id) %>% sort()
  females <- sex_table %>% filter(sex == "female") %>% pull(focal_id) %>% sort()
  chunks <- list()
  chunk_index <- 1L
  for (male_value in males) {
    for (female_value in females) {
      for (stage_value in c("before", "after")) {
        left <- data %>%
          filter(focal_id == male_value, stage == stage_value) %>%
          arrange(conspecific_id, session_number)
        right <- data %>%
          filter(focal_id == female_value, stage == stage_value) %>%
          arrange(conspecific_id, session_number)
        if (!nrow(left) || !nrow(right)) next
        indices <- expand_grid(
          left_row = seq_len(nrow(left)), right_row = seq_len(nrow(right))
        )
        chunks[[chunk_index]] <- tibble(
          dyad_id = canonical_dyad_id(male_value, female_value),
          bonded_status = if_else(
            is_bonded_dyad(male_value, female_value),
            "bonded", "nonbonded"
          ),
          stage = stage_value,
          individual_1 = male_value,
          individual_2 = female_value,
          male_receiver = left$conspecific_id[indices$left_row],
          female_receiver = right$conspecific_id[indices$right_row],
          session_1_id = left$session_id[indices$left_row],
          session_2_id = right$session_id[indices$right_row],
          source_index_1 = left$source_index[indices$left_row],
          source_index_2 = right$source_index[indices$right_row],
          session_number_1 = left$session_number[indices$left_row],
          session_number_2 = right$session_number[indices$right_row],
          n_sequences_1 = left$n_sequences[indices$left_row],
          n_sequences_2 = right$n_sequences[indices$right_row]
        ) %>%
          mutate(
            comparison_id = make_comparison_ids(session_1_id, session_2_id),
            direct_candidate_interaction =
              male_receiver == individual_2 |
              female_receiver == individual_1
          )
        chunk_index <- chunk_index + 1L
      }
    }
  }
  bind_rows(chunks)
}

model2_comparisons <- build_opposite_sex_comparisons(nonpartner_sessions)
if (nrow(model2_comparisons) != 689L) {
  stop("Expected 689 Model 2 session comparisons; observed ",
       nrow(model2_comparisons))
}
if (n_distinct(model2_comparisons$dyad_id) != 9L ||
    n_distinct(model2_comparisons$dyad_id[
      model2_comparisons$bonded_status == "bonded"
    ]) != 3L ||
    n_distinct(model2_comparisons$dyad_id[
      model2_comparisons$bonded_status == "nonbonded"
    ]) != 6L) {
  stop("Model 2 must contain 3 bonded and 6 nonbonded opposite-sex dyads")
}
expected_model2_counts <- tribble(
  ~bonded_status, ~stage, ~expected,
  "bonded", "before", 98L,
  "bonded", "after", 171L,
  "nonbonded", "before", 154L,
  "nonbonded", "after", 266L
)
observed_model2_counts <- model2_comparisons %>%
  count(bonded_status, stage, name = "observed") %>%
  left_join(expected_model2_counts,
            by = c("bonded_status", "stage"))
if (any(observed_model2_counts$observed != observed_model2_counts$expected)) {
  stop("Model 2 bond-status x stage counts changed")
}
if (anyDuplicated(model2_comparisons$comparison_id)) {
  stop("Model 2 contains duplicate session comparisons")
}

model2_unfactored <- add_metric_distances(model2_comparisons)
validate_metric_rows(model2_unfactored, "Model 2")
if (nrow(model2_unfactored) != 2756L) {
  stop("Expected 2,756 Model 2 comparison x metric rows")
}

individual_levels <- sort(unique(c(
  model2_unfactored$individual_1, model2_unfactored$individual_2
)))
session_levels_model2 <- sort(unique(c(
  model2_unfactored$session_1_id, model2_unfactored$session_2_id
)))
model2_data <- model2_unfactored %>%
  mutate(
    stage = factor(stage, levels = c("before", "after")),
    bonded_status = factor(
      bonded_status, levels = c("nonbonded", "bonded")
    ),
    metric = factor(metric, levels = SEQ_METRICS),
    dyad_id = factor(dyad_id),
    individual_1 = factor(individual_1, levels = individual_levels),
    individual_2 = factor(individual_2, levels = individual_levels),
    male_receiver = factor(male_receiver),
    female_receiver = factor(female_receiver),
    session_1_id = factor(session_1_id, levels = session_levels_model2),
    session_2_id = factor(session_2_id, levels = session_levels_model2),
    comparison_id = factor(comparison_id)
  )

# -------------------------
# Save prepared inputs and audit tables
# -------------------------
model1_input_path <- file.path(model_input_dir, "model1_individual_shift_data.csv")
model2_input_path <- file.path(model_input_dir, "model2_pair_specificity_data.csv")
write_csv(model1_unfactored, model1_input_path)
write_csv(model2_unfactored, model2_input_path)

write_csv(
  model1_comparisons %>%
    count(focal_id, conspecific_id, comparison_type, name = "N_comparisons"),
  file.path(model_input_dir, "model1_support_counts.csv")
)
write_csv(
  model2_comparisons %>%
    count(
      dyad_id, bonded_status, stage, male_receiver, female_receiver,
      direct_candidate_interaction,
      name = "N_comparisons"
    ),
  file.path(model_input_dir, "model2_support_counts.csv")
)

scaling_summary <- bind_rows(
  model1_unfactored %>% mutate(model = "Model 1"),
  model2_unfactored %>% mutate(model = "Model 2")
) %>%
  group_by(model, metric) %>%
  summarise(
    N_rows = n(),
    raw_mean = first(metric_mean),
    raw_population_sd = first(metric_population_sd),
    z_mean = mean(z_distance),
    z_population_sd = population_sd(z_distance),
    .groups = "drop"
  )
write_csv(scaling_summary, file.path(model_input_dir, "scaling_summary.csv"))

model_specification <- list(
  run_mode = run_mode,
  scope = "Non-partner sequence structure only",
  generated_by = "code/script_6_sequence_robustness_models.R",
  family = paste(
    "gaussian(identity) with a common residual SD, matching the primary",
    "Bayesian analysis (Script 4)"
  ),
  model_1_formula = format_formula(model1_formula),
  model_2_formula = format_formula(model2_formula),
  model_1_primary_estimand =
    paste(
      "Across-stage - 0.5 * (Within-before + Within-after); positive means that",
      "the across-stage distance exceeds the mean within-stage distance"
    ),
  model_2_primary_estimand =
    "(After - Before)_bonded - (After - Before)_nonbonded; negative supports additional bonded-dyad convergence",
  scaling = paste(
    "Population z score within metric after constructing each model's full",
    "comparison universe; never within stage, individual, dyad, or bond status."
  ),
  sampling = list(
    chains = model_chains,
    iterations = model_iterations,
    warmup = model_warmup,
    seed = 123L,
    adapt_delta = brms_control$adapt_delta,
    max_treedepth = brms_control$max_treedepth
  ),
  input_files = map_chr(artifact_paths, project_relative),
  input_sha256 = as.list(observed_hashes)
)
jsonlite::write_json(
  model_specification,
  file.path(out_dir, "analysis_provenance.json"),
  pretty = TRUE, auto_unbox = TRUE
)

if (prep_only) {
  cat("Script 6 preparation completed successfully.\n")
  quit(save = "no", status = 0)
}

if (validate_only) {
  cat("Building Stan data for Script 6 Model 1...\n")
  invisible(make_standata(
    model1_formula, data = model1_data, family = gaussian(),
    prior = priors_gaussian, backend = "rstan"
  ))
  cat("Building Stan data for Script 6 Model 2...\n")
  invisible(make_standata(
    model2_formula, data = model2_data, family = gaussian(),
    prior = priors_gaussian, backend = "rstan"
  ))
  cat("Script 6 validation completed successfully.\n")
  quit(save = "no", status = 0)
}

# -------------------------
# Fit both Bayesian sensitivity models
# -------------------------
input_fingerprints <- c(
  model1 = sha256_file(model1_input_path),
  model2 = sha256_file(model2_input_path)
)

fit_model <- function(formula, data, model_name, input_fingerprint) {
  formula_text <- format_formula(formula)
  fingerprint <- digest::digest(
    list(
      input_sha256 = input_fingerprint,
      formula = formula_text,
      priors = capture.output(print(priors_gaussian)),
      sampling = model_specification$sampling,
      frozen_inputs = observed_hashes
    ),
    algo = "sha256"
  )
  cache_name <- paste0(model_name, "_", substr(fingerprint, 1, 12))
  cache_path <- file.path(model_dir, cache_name)
  cat("\n--- Fitting/loading ", cache_name, " ---\n", sep = "")
  fit <- brm(
    formula = formula,
    data = data,
    family = gaussian(),
    prior = priors_gaussian,
    chains = model_chains,
    cores = model_cores,
    iter = model_iterations,
    warmup = model_warmup,
    seed = 123,
    control = brms_control,
    backend = "rstan",
    save_pars = save_pars(all = TRUE),
    file = cache_path,
    file_refit = if (force_refit) "always" else "on_change",
    refresh = 100
  )
  list(fit = fit, fingerprint = fingerprint, cache_path = paste0(cache_path, ".rds"))
}

run_marker <- list(
  status = "in_progress",
  run_mode = run_mode,
  manuscript_output_eligible = FALSE,
  started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
jsonlite::write_json(
  run_marker, file.path(out_dir, "analysis_completion.json"),
  pretty = TRUE, auto_unbox = TRUE
)
jsonlite::write_json(
  run_marker, file.path(out_dir, "manuscript_readiness.json"),
  pretty = TRUE, auto_unbox = TRUE
)

model1_result <- fit_model(
  model1_formula, model1_data, "model1_individual_stage_shift",
  input_fingerprints[["model1"]]
)
model2_result <- fit_model(
  model2_formula, model2_data, "model2_pair_specificity",
  input_fingerprints[["model2"]]
)
fit_model1 <- model1_result$fit
fit_model2 <- model2_result$fit

# -------------------------
# Model 1 posterior estimands
# -------------------------
model1_strata <- model1_comparisons %>%
  distinct(focal_id, conspecific_id)
model1_grid <- crossing(
  model1_strata,
  comparison_type = COMPARISON_TYPES,
  metric = SEQ_METRICS
) %>%
  mutate(
    focal_id = factor(focal_id, levels = levels(model1_data$focal_id)),
    conspecific_id = factor(
      conspecific_id, levels = levels(model1_data$conspecific_id)
    ),
    comparison_type = factor(
      comparison_type, levels = levels(model1_data$comparison_type)
    ),
    metric = factor(metric, levels = levels(model1_data$metric)),
    session_1_id = factor(
      levels(model1_data$session_1_id)[1],
      levels = levels(model1_data$session_1_id)
    ),
    session_2_id = factor(
      levels(model1_data$session_2_id)[1],
      levels = levels(model1_data$session_2_id)
    ),
    grid_row = as.character(row_number())
  )

model1_expected <- posterior_epred(
  fit_model1,
  newdata = model1_grid,
  re_formula = ~(1 | focal_id) + (1 | conspecific_id)
)
model1_expected_long <- as_tibble(as.data.frame(model1_expected))
names(model1_expected_long) <- as.character(seq_len(ncol(model1_expected_long)))
model1_expected_long <- model1_expected_long %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(-.draw, names_to = "grid_row", values_to = "expected") %>%
  left_join(
    model1_grid %>%
      transmute(
        grid_row, focal_id = as.character(focal_id),
        conspecific_id = as.character(conspecific_id),
        comparison_type = as.character(comparison_type),
        metric = as.character(metric)
      ),
    by = "grid_row"
  )

model1_individual_draws <- model1_expected_long %>%
  group_by(.draw, focal_id, metric, comparison_type) %>%
  summarise(expected = mean(expected), .groups = "drop") %>%
  pivot_wider(names_from = comparison_type, values_from = expected) %>%
  mutate(
    shift_excess =
      across_stage - 0.5 * (within_before + within_after),
    across_minus_before = across_stage - within_before,
    across_minus_after = across_stage - within_after,
    after_minus_before_dispersion = within_after - within_before
  )
saveRDS(
  model1_individual_draws,
  file.path(table_dir, "model1_individual_contrast_draws.rds")
)

model1_metric_draws <- model1_individual_draws %>%
  group_by(.draw, metric) %>%
  summarise(shift_excess = mean(shift_excess), .groups = "drop") %>%
  mutate(scope = "Population contrast across eligible common-receiver sessions")
model1_combined_draws <- model1_metric_draws %>%
  group_by(.draw, scope) %>%
  summarise(shift_excess = mean(shift_excess), .groups = "drop") %>%
  mutate(metric = "Combined")
model1_draws <- bind_rows(model1_metric_draws, model1_combined_draws)
saveRDS(model1_draws, file.path(table_dir, "model1_primary_draws.rds"))
model1_summary <- summarise_posterior(
  model1_draws, "shift_excess", c("scope", "metric")
) %>%
  mutate(
    estimand = "Across - 0.5 * (Within-before + Within-after)",
    direction = paste(
      "positive = across-stage distance exceeds the mean of within-Before",
      "and within-After variation"
    )
  )
write_csv(
  model1_summary,
  file.path(table_dir, "model1_individual_stage_shift_results.csv")
)

model1_component_draws_metric <- model1_individual_draws %>%
  group_by(.draw, metric) %>%
  summarise(
    across_minus_before = mean(across_minus_before),
    across_minus_after = mean(across_minus_after),
    after_minus_before_dispersion = mean(after_minus_before_dispersion),
    .groups = "drop"
  )
model1_component_draws_combined <- model1_component_draws_metric %>%
  group_by(.draw) %>%
  summarise(
    across_minus_before = mean(across_minus_before),
    across_minus_after = mean(across_minus_after),
    after_minus_before_dispersion = mean(after_minus_before_dispersion),
    .groups = "drop"
  ) %>%
  mutate(metric = "Combined")
model1_component_draws <- bind_rows(
  model1_component_draws_metric,
  model1_component_draws_combined
) %>%
  pivot_longer(
    cols = c(
      across_minus_before,
      across_minus_after,
      after_minus_before_dispersion
    ),
    names_to = "contrast",
    values_to = "contrast_value"
  )
write_csv(
  summarise_posterior(
    model1_component_draws,
    "contrast_value",
    c("metric", "contrast")
  ) %>%
    mutate(
      contrast_definition = recode(
        contrast,
        across_minus_before = "Before-After minus Before-Before",
        across_minus_after = "Before-After minus After-After",
        after_minus_before_dispersion = "After-After minus Before-Before"
      )
    ),
  file.path(table_dir, "model1_component_contrasts.csv")
)

# -------------------------
# Model 2 posterior estimands
# -------------------------
dyad_specs <- model2_comparisons %>%
  distinct(dyad_id, bonded_status, individual_1, individual_2)
model2_grid <- crossing(
  dyad_specs,
  stage = c("before", "after"),
  metric = SEQ_METRICS
) %>%
  mutate(
    dyad_id = factor(dyad_id, levels = levels(model2_data$dyad_id)),
    bonded_status = factor(
      bonded_status, levels = levels(model2_data$bonded_status)
    ),
    individual_1 = factor(individual_1, levels = individual_levels),
    individual_2 = factor(individual_2, levels = individual_levels),
    stage = factor(stage, levels = levels(model2_data$stage)),
    metric = factor(metric, levels = levels(model2_data$metric)),
    session_1_id = factor(
      levels(model2_data$session_1_id)[1],
      levels = levels(model2_data$session_1_id)
    ),
    session_2_id = factor(
      levels(model2_data$session_2_id)[1],
      levels = levels(model2_data$session_2_id)
    ),
    grid_row = as.character(row_number())
  )

if (nrow(model2_grid) != 72L) {
  stop("Expected a 72-row dyad x stage x metric Model 2 prediction grid")
}

model2_expected <- posterior_epred(
  fit_model2,
  newdata = model2_grid,
  re_formula =
    ~(1 + stage || dyad_id) +
     (1 | mm(individual_1, individual_2))
)
model2_expected_long <- as_tibble(as.data.frame(model2_expected))
names(model2_expected_long) <- as.character(seq_len(ncol(model2_expected_long)))
model2_expected_long <- model2_expected_long %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(-.draw, names_to = "grid_row", values_to = "expected") %>%
  left_join(
    model2_grid %>%
      transmute(
        grid_row,
        dyad_id = as.character(dyad_id),
        bonded_status = as.character(bonded_status),
        stage = as.character(stage),
        metric = as.character(metric)
      ),
    by = "grid_row"
  )

model2_dyad_draws <- model2_expected_long %>%
  group_by(.draw, dyad_id, bonded_status, stage, metric) %>%
  summarise(expected = mean(expected), .groups = "drop") %>%
  pivot_wider(names_from = stage, values_from = expected) %>%
  mutate(stage_effect = after - before)
saveRDS(
  model2_dyad_draws,
  file.path(table_dir, "model2_dyad_stage_effect_draws.rds")
)

model2_group_metric_draws <- model2_dyad_draws %>%
  group_by(.draw, bonded_status, metric) %>%
  summarise(stage_effect = mean(stage_effect), .groups = "drop")
model2_group_combined_draws <- model2_group_metric_draws %>%
  group_by(.draw, bonded_status) %>%
  summarise(stage_effect = mean(stage_effect), .groups = "drop") %>%
  mutate(metric = "Combined")
model2_group_draws <- bind_rows(
  model2_group_metric_draws,
  model2_group_combined_draws
)
write_csv(
  summarise_posterior(
    model2_group_draws, "stage_effect", c("bonded_status", "metric")
  ) %>%
    mutate(
      estimand = "Expected After - Before",
      direction = "negative = convergence"
    ),
  file.path(table_dir, "model2_group_stage_effects_by_metric.csv")
)

model2_specificity_draws_metric <- model2_group_metric_draws %>%
  select(.draw, bonded_status, metric, stage_effect) %>%
  pivot_wider(names_from = bonded_status, values_from = stage_effect) %>%
  mutate(specificity_contrast = bonded - nonbonded)
model2_specificity_draws_combined <- model2_specificity_draws_metric %>%
  group_by(.draw) %>%
  summarise(
    specificity_contrast = mean(specificity_contrast),
    .groups = "drop"
  ) %>%
  mutate(metric = "Combined")
model2_specificity_draws <- bind_rows(
  model2_specificity_draws_metric %>%
    select(.draw, metric, specificity_contrast),
  model2_specificity_draws_combined
)
saveRDS(
  model2_specificity_draws,
  file.path(table_dir, "model2_pair_specificity_draws.rds")
)
model2_summary <- summarise_posterior(
  model2_specificity_draws,
  "specificity_contrast",
  "metric"
) %>%
  mutate(
    estimand = "(After - Before)_bonded - (After - Before)_nonbonded",
    direction = "negative = greater convergence among bonded dyads"
  )
write_csv(
  model2_summary,
  file.path(table_dir, "model2_pair_specificity_results.csv")
)

# -------------------------
# Diagnostics and figures
# -------------------------
fit_diagnostics <- function(fit, model_name) {
  draws <- posterior::as_draws_array(fit)
  summaries <- posterior::summarise_draws(
    draws, posterior::default_convergence_measures()
  ) %>% filter(variable != "lp__")
  sampler <- brms::nuts_params(fit)
  tibble(
    model = model_name,
    max_Rhat = max(summaries$rhat, na.rm = TRUE),
    min_bulk_ESS = min(summaries$ess_bulk, na.rm = TRUE),
    min_tail_ESS = min(summaries$ess_tail, na.rm = TRUE),
    divergences = sum(
      sampler$Parameter == "divergent__" & sampler$Value > 0
    ),
    max_treedepth_hits = sum(
      sampler$Parameter == "treedepth__" &
      sampler$Value >= brms_control$max_treedepth
    )
  )
}

diagnostics <- bind_rows(
  fit_diagnostics(fit_model1, "Model 1: stage-variation contrast"),
  fit_diagnostics(fit_model2, "Model 2: pair specificity")
)
write_csv(diagnostics, file.path(diagnostic_dir, "diagnostics_summary.csv"))

save_ppc <- function(fit, filename) {
  plot <- pp_check(fit, type = "dens_overlay", ndraws = 100)
  ggsave(
    file.path(diagnostic_dir, filename), plot,
    width = 7, height = 5, dpi = 300
  )
}
save_ppc(fit_model1, "model1_posterior_predictive_density.png")
save_ppc(fit_model2, "model2_posterior_predictive_density.png")

metric_plot_order <- c(
  "Combined",
  "Transition probabilities",
  "Bigram distribution",
  "Phee-repeat distribution",
  "Local alignment"
)
posterior_fill <- "#b8de29"
posterior_line <- "#86a900"

prepare_posterior_plot_data <- function(data, value_column) {
  data %>%
    transmute(
      .draw,
      metric_label = recode(
        as.character(metric), !!!METRIC_LABELS, Combined = "Combined"
      ),
      posterior_value = .data[[value_column]]
    ) %>%
    mutate(
      metric_label = factor(
        metric_label,
        levels = rev(metric_plot_order)
      )
    )
}

make_posterior_panel <- function(plot_data, title, x_label) {
  plot_summary <- plot_data %>%
    group_by(metric_label) %>%
    summarise(
      estimate = median(posterior_value),
      lower = quantile(posterior_value, 0.025),
      upper = quantile(posterior_value, 0.975),
      .groups = "drop"
    )

  ggplot(plot_data, aes(x = posterior_value, y = metric_label)) +
    geom_vline(
      xintercept = 0, linetype = "dashed", colour = "grey50",
      linewidth = 0.45
    ) +
    stat_halfeye(
      fill = posterior_fill,
      colour = NA,
      slab_alpha = 0.72,
      .width = 0,
      point_interval = NULL,
      normalize = "groups",
      height = 0.72
    ) +
    geom_segment(
      data = plot_summary,
      aes(
        x = lower, xend = upper,
        y = metric_label, yend = metric_label
      ),
      colour = posterior_line,
      linewidth = 1.05,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = plot_summary,
      aes(x = estimate, y = metric_label),
      colour = posterior_line,
      size = 2.2,
      inherit.aes = FALSE
    ) +
    labs(title = title, x = x_label, y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      plot.margin = margin(8, 10, 8, 8)
    )
}

model1_plot_data <- prepare_posterior_plot_data(
  model1_draws, "shift_excess"
)
model1_plot <- make_posterior_panel(
  model1_plot_data,
  "Across- versus within-stage variation",
  "Across-stage - mean within-stage distance\n(metric-specific SD units)"
)

model2_plot_data <- prepare_posterior_plot_data(
  model2_specificity_draws, "specificity_contrast"
)
model2_plot <- make_posterior_panel(
  model2_plot_data,
  "Bonded versus nonbonded change",
  "Bonded - nonbonded stage change\n(metric-specific SD units)"
)

combined_robustness_plot <-
  (model1_plot | model2_plot) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))

# Retain panel-level outputs for checking and write one canonical two-panel
# supplementary figure for the manuscript.
ggsave(
  file.path(figure_dir, "model1_individual_stage_shift.png"),
  model1_plot, width = 6.2, height = 5.2, dpi = 300
)
ggsave(
  file.path(figure_dir, "model1_individual_stage_shift.pdf"),
  model1_plot, width = 6.2, height = 5.2
)
ggsave(
  file.path(figure_dir, "model2_pair_specificity.png"),
  model2_plot, width = 6.2, height = 5.2, dpi = 300
)
ggsave(
  file.path(figure_dir, "model2_pair_specificity.pdf"),
  model2_plot, width = 6.2, height = 5.2
)
ggsave(
  file.path(figure_dir, "Fig_S9_sequence_robustness_posteriors.png"),
  combined_robustness_plot, width = 12.4, height = 5.2, dpi = 300
)
ggsave(
  file.path(figure_dir, "Fig_S9_sequence_robustness_posteriors.pdf"),
  combined_robustness_plot, width = 12.4, height = 5.2
)

diagnostic_gate <- diagnostics %>%
  summarise(
    passed = all(
      max_Rhat <= 1.01,
      min_bulk_ESS >= 400,
      min_tail_ESS >= 400,
      divergences == 0,
      max_treedepth_hits == 0
    )
  ) %>% pull(passed)
manuscript_ready <- !quick_fit && isTRUE(diagnostic_gate)

completion <- list(
  status = "complete",
  run_mode = run_mode,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  manuscript_output_eligible = manuscript_ready,
  diagnostic_gate_passed = isTRUE(diagnostic_gate),
  model_1_fingerprint = model1_result$fingerprint,
  model_2_fingerprint = model2_result$fingerprint,
  output_directory = project_relative(out_dir),
  primary_tables = c(
    project_relative(file.path(
      table_dir, "model1_individual_stage_shift_results.csv"
    )),
    project_relative(file.path(
      table_dir, "model2_pair_specificity_results.csv"
    ))
  ),
  primary_figures = c(
    project_relative(file.path(
      figure_dir, "Fig_S9_sequence_robustness_posteriors.png"
    )),
    project_relative(file.path(
      figure_dir, "Fig_S9_sequence_robustness_posteriors.pdf"
    ))
  )
)
jsonlite::write_json(
  completion, file.path(out_dir, "analysis_completion.json"),
  pretty = TRUE, auto_unbox = TRUE
)
jsonlite::write_json(
  completion, file.path(out_dir, "manuscript_readiness.json"),
  pretty = TRUE, auto_unbox = TRUE
)

cat("\nScript 6 completed.\n")
cat("Output directory: ", out_dir, "\n", sep = "")
cat("Manuscript output eligible: ", manuscript_ready, "\n", sep = "")
