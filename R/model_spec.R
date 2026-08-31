# Model specification and input validation for the two confirmatory analyses.
#
# This file only defines functions and constants. Sourcing it never loads a
# fitted model or starts Stan.

SEQUENCE_METRICS <- c(
  "transition_probability",
  "bigram",
  "phee_repeat",
  "local_alignment"
)

CALL_METRICS <- c("stp", "mfcc", "dtw", "vae")

MODEL_METRICS <- list(
  sequence = SEQUENCE_METRICS,
  call = CALL_METRICS
)

METRIC_LABELS <- c(
  transition_probability = "Transition probabilities",
  bigram = "Bigram distribution",
  phee_repeat = "Repeat distribution of phee",
  local_alignment = "Local alignment",
  stp = "Spectro-temporal parameters (STP)",
  mfcc = "Mel-frequency cepstral coefficients (MFCC)",
  dtw = "Dynamic time warping (DTW)",
  vae = "Variational autoencoder (VAE)"
)

# Aliases make the model layer compatible with both the historical tables and
# the concise names used by the rebuilt preprocessing workflow. All exported
# posterior products use the canonical names above.
METRIC_ALIASES <- c(
  transition_matrix = "transition_probability",
  transition_probabilities = "transition_probability",
  transition_probability = "transition_probability",
  bigram = "bigram",
  bigram_distribution = "bigram",
  repeat_A_len = "phee_repeat",
  phee_repeat = "phee_repeat",
  repeat_distribution = "phee_repeat",
  repeat_distribution_of_phee = "phee_repeat",
  local_alignment = "local_alignment",
  T_PC_pairwise = "stp",
  STP = "stp",
  stp = "stp",
  M_PC_pairwise = "mfcc",
  MFCC = "mfcc",
  mfcc = "mfcc",
  spectro_DTW = "dtw",
  DTW = "dtw",
  dtw = "dtw",
  VAE_pairwise = "vae",
  VAE = "vae",
  vae = "vae"
)

EXPECTED_ANALYSIS_COUNTS <- list(
  sequence = c(rows = 1324L, comparisons = 331L, metrics = 4L),
  call = c(rows = 1724L, comparisons = 431L, metrics = 4L)
)

model_formula <- function() {
  stats::as.formula(
    paste0(
      "z_distance ~ stage * context * metric + ",
      "(1 + stage * context || pair) + ",
      "(1 | mm(session_1_id, session_2_id))"
    )
  )
}

pair_prediction_formula <- function() {
  # Pair effects are retained when obtaining manuscript contrasts. The
  # session multi-membership effect is deliberately marginalized out.
  stats::as.formula("~ (1 + stage * context || pair)")
}

normalise_analysis_name <- function(analysis) {
  if (length(analysis) != 1L || is.na(analysis)) {
    stop("analysis must be one of: sequence, call", call. = FALSE)
  }
  key <- tolower(trimws(as.character(analysis)))
  aliases <- c(
    sequence = "sequence",
    "sequence structure" = "sequence",
    call = "call",
    spectral = "call",
    acoustic = "call",
    "call structure" = "call"
  )
  if (!key %in% names(aliases)) {
    stop("analysis must be one of: sequence, call", call. = FALSE)
  }
  unname(aliases[[key]])
}

analysis_label <- function(analysis) {
  switch(
    normalise_analysis_name(analysis),
    sequence = "Sequence structure",
    call = "Call structure"
  )
}

canonicalise_metric <- function(metric) {
  metric <- as.character(metric)
  unknown <- sort(unique(metric[!metric %in% names(METRIC_ALIASES)]))
  if (length(unknown)) {
    stop(
      "Unknown metric code(s): ", paste(unknown, collapse = ", "),
      ". See METRIC_ALIASES in R/model_spec.R.",
      call. = FALSE
    )
  }
  unname(METRIC_ALIASES[metric])
}

require_columns <- function(data, columns, label = "Model data") {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(
      label, " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

normalise_stage <- function(stage) {
  value <- tolower(trimws(as.character(stage)))
  value[value %in% c("pre", "baseline")] <- "before"
  value[value %in% c("post")] <- "after"
  value
}

normalise_context <- function(context) {
  value <- tolower(trimws(as.character(context)))
  value <- gsub("_", "-", value, fixed = TRUE)
  value[value %in% c("partner", "paired")] <- "Partner"
  value[value %in% c(
    "non-partner", "nonpartner", "non-paired", "stranger", "unpaired"
  )] <-
    "Non-partner"
  value
}

population_sd <- function(x) {
  sqrt(mean((x - mean(x))^2))
}

copy_column_alias <- function(data, canonical, alias) {
  if (!canonical %in% names(data) && alias %in% names(data)) {
    data[[canonical]] <- data[[alias]]
  } else if (all(c(canonical, alias) %in% names(data))) {
    canonical_value <- as.character(data[[canonical]])
    alias_value <- as.character(data[[alias]])
    if (anyNA(canonical_value) || anyNA(alias_value) ||
        !identical(canonical_value, alias_value)) {
      stop(
        "Columns ", canonical, " and ", alias, " disagree",
        call. = FALSE
      )
    }
  }
  data
}

centered_stage_contrasts <- function() {
  matrix(
    c(-0.5, 0.5),
    nrow = 2L,
    dimnames = list(c("before", "after"), "after_minus_before")
  )
}

validate_standardisation <- function(data, tolerance = 1e-6) {
  metrics <- levels(data$metric)
  problems <- vapply(metrics, function(metric_name) {
    values <- data$z_distance[data$metric == metric_name]
    mean_bad <- abs(mean(values)) > tolerance
    sd_bad <- abs(population_sd(values) - 1) > tolerance
    mean_bad || sd_bad
  }, logical(1))
  if (any(problems)) {
    stop(
      "z_distance must have mean 0 and population SD 1 within each metric; ",
      "failed for: ", paste(metrics[problems], collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

complete_support_pairs <- function(data) {
  require_columns(data, c("pair", "stage", "context"))
  observed <- unique(data.frame(
    pair = as.character(data$pair),
    stage = as.character(data$stage),
    context = as.character(data$context),
    stringsAsFactors = FALSE
  ))
  required_cells <- expand.grid(
    stage = c("before", "after"),
    context = c("Non-partner", "Partner"),
    stringsAsFactors = FALSE
  )
  pairs <- sort(unique(observed$pair))
  pairs[vapply(pairs, function(pair_name) {
    pair_cells <- observed[observed$pair == pair_name, c("stage", "context")]
    all(paste(required_cells$stage, required_cells$context) %in%
          paste(pair_cells$stage, pair_cells$context))
  }, logical(1))]
}

prepare_model_data <- function(
    data,
    analysis,
    check_standardisation = TRUE,
    expected_pairs = 3L,
    expected_common_support_pairs = 2L) {
  analysis <- normalise_analysis_name(analysis)
  label <- analysis_label(analysis)
  data <- copy_column_alias(data, "pair", "pair_id")
  data <- copy_column_alias(data, "session_1_id", "repertoire_a_id")
  data <- copy_column_alias(data, "session_2_id", "repertoire_b_id")
  if (!"z_distance" %in% names(data) &&
      "standardized_distance" %in% names(data)) {
    data$z_distance <- data$standardized_distance
  } else if (all(c("z_distance", "standardized_distance") %in% names(data))) {
    z_value <- suppressWarnings(as.numeric(as.character(data$z_distance)))
    standard_value <- suppressWarnings(as.numeric(
      as.character(data$standardized_distance)
    ))
    if (anyNA(z_value) || anyNA(standard_value) ||
        any(abs(z_value - standard_value) > 1e-12)) {
      stop(
        label,
        " z_distance and standardized_distance columns disagree",
        call. = FALSE
      )
    }
  }
  required <- c(
    "pair", "context", "stage", "session_1_id", "session_2_id",
    "comparison_id", "metric", "z_distance"
  )
  require_columns(data, required, paste(label, "model data"))

  if (!nrow(data)) {
    stop(label, " model data has no rows", call. = FALSE)
  }
  if (anyNA(data[required])) {
    missing_counts <- colSums(is.na(data[required]))
    details <- paste(
      names(missing_counts)[missing_counts > 0],
      missing_counts[missing_counts > 0],
      sep = "=",
      collapse = ", "
    )
    stop(label, " model data contains missing values: ", details, call. = FALSE)
  }

  data$stage <- normalise_stage(data$stage)
  data$context <- normalise_context(data$context)
  data$metric <- canonicalise_metric(data$metric)
  data$z_distance <- suppressWarnings(as.numeric(as.character(data$z_distance)))

  if (any(!is.finite(data$z_distance))) {
    stop(label, " z_distance contains non-finite values", call. = FALSE)
  }
  if (!setequal(unique(data$stage), c("before", "after"))) {
    stop(label, " model data must contain before and after stages", call. = FALSE)
  }
  expected_stage_centered <- ifelse(data$stage == "after", 0.5, -0.5)
  if ("stage_centered" %in% names(data)) {
    observed_stage_centered <- suppressWarnings(as.numeric(
      as.character(data$stage_centered)
    ))
    if (anyNA(observed_stage_centered) ||
        any(abs(observed_stage_centered - expected_stage_centered) > 1e-12)) {
      stop(
        label,
        " stage_centered must code before=-0.5 and after=+0.5",
        call. = FALSE
      )
    }
  } else {
    data$stage_centered <- expected_stage_centered
  }
  if (!setequal(unique(data$context), c("Non-partner", "Partner"))) {
    stop(
      label,
      " model data must contain Partner and Non-partner contexts",
      call. = FALSE
    )
  }

  expected_metrics <- MODEL_METRICS[[analysis]]
  if (!setequal(unique(data$metric), expected_metrics)) {
    stop(
      label, " metric mismatch. Expected: ",
      paste(expected_metrics, collapse = ", "),
      "; observed: ", paste(sort(unique(data$metric)), collapse = ", "),
      call. = FALSE
    )
  }
  if (any(as.character(data$session_1_id) == as.character(data$session_2_id))) {
    stop(label, " model data contains a self-comparison", call. = FALSE)
  }

  duplicate_key <- paste(data$comparison_id, data$metric, sep = "\r")
  if (anyDuplicated(duplicate_key)) {
    stop(
      label, " model data contains duplicate comparison_id x metric rows",
      call. = FALSE
    )
  }

  metric_count <- table(as.character(data$comparison_id))
  if (any(metric_count != length(expected_metrics))) {
    stop(
      label, " model data has comparisons without exactly ",
      length(expected_metrics), " metrics",
      call. = FALSE
    )
  }

  pair_names <- sort(unique(as.character(data$pair)))
  if (!is.null(expected_pairs) && length(pair_names) != expected_pairs) {
    stop(
      label, " model data must contain exactly ", expected_pairs,
      " pairs; observed ", length(pair_names),
      call. = FALSE
    )
  }

  session_levels <- sort(unique(c(
    as.character(data$session_1_id),
    as.character(data$session_2_id)
  )))
  data$stage <- factor(data$stage, levels = c("before", "after"))
  contrasts(data$stage) <- centered_stage_contrasts()
  data$context <- factor(data$context, levels = c("Non-partner", "Partner"))
  data$metric <- factor(data$metric, levels = expected_metrics)
  data$pair <- factor(as.character(data$pair), levels = pair_names)
  data$comparison_id <- factor(as.character(data$comparison_id))
  data$session_1_id <- factor(as.character(data$session_1_id), levels = session_levels)
  data$session_2_id <- factor(as.character(data$session_2_id), levels = session_levels)

  common_pairs <- complete_support_pairs(data)
  if (!is.null(expected_common_support_pairs) &&
      length(common_pairs) != expected_common_support_pairs) {
    stop(
      label, " model data must contain exactly ",
      expected_common_support_pairs,
      " pairs with complete stage-by-context support; observed ",
      length(common_pairs),
      call. = FALSE
    )
  }

  if (check_standardisation) {
    validate_standardisation(data)
  }
  data
}

read_model_data <- function(
    path,
    analysis,
    check_standardisation = TRUE,
    expected_pairs = 3L,
    expected_common_support_pairs = 2L) {
  if (!file.exists(path)) {
    stop(
      "Missing ", normalise_analysis_name(analysis), " model input: ", path,
      ". Run the model-table preprocessing step first.",
      call. = FALSE
    )
  }
  data <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  prepare_model_data(
    data,
    analysis = analysis,
    check_standardisation = check_standardisation,
    expected_pairs = expected_pairs,
    expected_common_support_pairs = expected_common_support_pairs
  )
}

validate_expected_counts <- function(data, analysis) {
  analysis <- normalise_analysis_name(analysis)
  expected <- EXPECTED_ANALYSIS_COUNTS[[analysis]]
  observed_rows <- nrow(data)
  observed_comparisons <- length(unique(as.character(data$comparison_id)))
  observed_metrics <- length(unique(as.character(data$metric)))
  observed <- c(
    rows = observed_rows,
    comparisons = observed_comparisons,
    metrics = observed_metrics
  )
  if (!identical(as.integer(observed), as.integer(expected))) {
    stop(
      analysis_label(analysis), " manuscript-count mismatch. Expected ",
      paste(names(expected), expected, sep = "=", collapse = ", "),
      "; observed ",
      paste(names(observed), observed, sep = "=", collapse = ", "),
      call. = FALSE
    )
  }
  invisible(observed)
}
