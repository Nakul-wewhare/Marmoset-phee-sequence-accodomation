repo_root <- normalizePath(file.path(getwd()), mustWork = TRUE)
if (!file.exists(file.path(repo_root, "R", "model_spec.R"))) {
  repo_root <- normalizePath(file.path(getwd(), "..", ".."), mustWork = TRUE)
}

source(file.path(repo_root, "R", "model_spec.R"))
source(file.path(repo_root, "R", "model_engine.R"))
source(file.path(repo_root, "R", "posterior_results.R"))
source(file.path(repo_root, "R", "run_models.R"))

synthetic_model_data <- function(analysis = "sequence") {
  analysis <- normalise_analysis_name(analysis)
  metrics <- MODEL_METRICS[[analysis]]
  grid <- expand.grid(
    pair = c("Pair A", "Pair B", "Pair C"),
    context = c("Non-partner", "Partner"),
    stage = c("before", "after"),
    stringsAsFactors = FALSE
  )
  # Pair C is the manuscript's unsupported Partner-after cell.
  grid <- grid[!(
    grid$pair == "Pair C" &
      grid$context == "Partner" &
      grid$stage == "after"
  ), ]
  grid$comparison_id <- paste(
    grid$pair,
    grid$context,
    grid$stage,
    sep = "__"
  )
  grid$session_1_id <- paste0(grid$comparison_id, "__s1")
  grid$session_2_id <- paste0(grid$comparison_id, "__s2")
  data <- merge(
    grid,
    data.frame(metric = metrics, stringsAsFactors = FALSE),
    all = TRUE
  )
  data <- data[order(data$comparison_id, match(data$metric, metrics)), ]
  data$z_distance <- rep(seq_along(metrics), length.out = nrow(data))
  rownames(data) <- NULL
  data
}

synthetic_stage_draws <- function(analysis = "sequence", n_draws = 2L) {
  analysis <- normalise_analysis_name(analysis)
  metrics <- MODEL_METRICS[[analysis]]
  grid <- expand.grid(
    .draw = seq_len(n_draws),
    pair = c("Pair A", "Pair B", "Pair C"),
    context = c("Non-partner", "Partner"),
    metric = metrics,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pair_value <- c(`Pair A` = 0, `Pair B` = 2, `Pair C` = 100)
  metric_value <- stats::setNames(seq_along(metrics) - 1, metrics)
  grid$delta <-
    unname(pair_value[grid$pair]) +
    ifelse(grid$context == "Partner", 10, 0) +
    unname(metric_value[grid$metric]) +
    (grid$.draw - 1)
  grid$analysis <- analysis_label(analysis)
  grid$metric_label <- unname(METRIC_LABELS[grid$metric])
  grid$expected_before <- 0
  grid$expected_after <- grid$delta
  grid$common_support_pair <- grid$pair %in% c("Pair A", "Pair B")
  grid[required_draw_columns()]
}

synthetic_sampling_metadata <- function(n_draws = 2L, seed = 42L) {
  model_sampling_metadata(
    seed = seed,
    chains = 1L,
    iter = as.integer(n_draws) + 1L,
    warmup = 1L,
    backend = "rstan"
  )
}

synthetic_public_model_data <- function(analysis = "sequence") {
  data <- synthetic_model_data(analysis)
  data$z_distance <- seq_len(nrow(data))
  data$standardized_distance <- ave(
    data$z_distance,
    data$metric,
    FUN = function(value) {
      (value - mean(value)) / sqrt(mean((value - mean(value))^2))
    }
  )
  data$stage_centered <- ifelse(data$stage == "after", 0.5, -0.5)
  data$pair_id <- data$pair
  data$repertoire_a_id <- data$session_1_id
  data$repertoire_b_id <- data$session_2_id
  data[c("pair", "session_1_id", "session_2_id", "z_distance")] <- NULL
  data
}
