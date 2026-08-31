# Explicit fitting and posterior-cache functions.
#
# None of these functions run automatically. In particular, fit_model() is
# only called by the command-line driver after an explicit --refit flag.

model_priors <- function() {
  c(
    brms::prior("normal(0, 1)", class = "Intercept"),
    brms::prior("normal(0, 1)", class = "b"),
    brms::prior("exponential(1)", class = "sd"),
    brms::prior("exponential(1)", class = "sigma")
  )
}

default_model_cores <- function() {
  detected <- parallel::detectCores(logical = FALSE)
  if (is.na(detected) || detected < 1L) detected <- 1L
  min(4L, as.integer(detected))
}

validate_brms_model <- function(data) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop(
      "Package 'brms' is required for Stan-data validation. ",
      "Install the locked R environment before using --validate-only.",
      call. = FALSE
    )
  }
  invisible(brms::make_standata(
    formula = brms::bf(model_formula()),
    data = data,
    family = stats::gaussian(),
    prior = model_priors(),
    backend = "rstan"
  ))
}

fit_model <- function(
    data,
    analysis,
    chains = 4L,
    iter = 4000L,
    warmup = 2000L,
    cores = default_model_cores(),
    seed = 123L,
    backend = "rstan",
    refresh = 100L) {
  analysis <- normalise_analysis_name(analysis)
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop(
      "Package 'brms' is required for --refit. ",
      "The default cached-report path does not require brms.",
      call. = FALSE
    )
  }
  if (iter <= warmup) {
    stop("iter must be greater than warmup", call. = FALSE)
  }
  message(
    "Explicit --refit requested: fitting the ", analysis_label(analysis),
    " model. This is the expensive path."
  )
  brms::brm(
    formula = brms::bf(model_formula()),
    data = data,
    family = stats::gaussian(),
    prior = model_priors(),
    chains = as.integer(chains),
    iter = as.integer(iter),
    warmup = as.integer(warmup),
    cores = as.integer(cores),
    seed = as.integer(seed),
    backend = backend,
    control = list(adapt_delta = 0.995, max_treedepth = 15L),
    save_pars = brms::save_pars(all = TRUE),
    refresh = as.integer(refresh)
  )
}

prediction_grid <- function(data) {
  cells <- expand.grid(
    pair = levels(data$pair),
    context = levels(data$context),
    metric = levels(data$metric),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  before <- cells
  after <- cells
  before$stage <- "before"
  after$stage <- "after"
  grid <- rbind(before, after)

  grid$stage <- factor(grid$stage, levels = levels(data$stage))
  contrasts(grid$stage) <- contrasts(data$stage)
  grid$context <- factor(grid$context, levels = levels(data$context))
  grid$metric <- factor(grid$metric, levels = levels(data$metric))
  grid$pair <- factor(grid$pair, levels = levels(data$pair))

  # These columns are not used because re_formula excludes the session term,
  # but supplying valid levels keeps prediction robust across brms versions.
  grid$session_1_id <- factor(
    levels(data$session_1_id)[[1]],
    levels = levels(data$session_1_id)
  )
  grid$session_2_id <- factor(
    levels(data$session_2_id)[[1]],
    levels = levels(data$session_2_id)
  )
  list(cells = cells, grid = grid)
}

extract_pair_stage_draws <- function(fit, data, analysis) {
  analysis <- normalise_analysis_name(analysis)
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required to extract draws from a fit", call. = FALSE)
  }
  grids <- prediction_grid(data)
  n_cells <- nrow(grids$cells)

  expected <- brms::posterior_epred(
    fit,
    newdata = grids$grid,
    re_formula = pair_prediction_formula(),
    allow_new_levels = FALSE
  )
  if (!is.matrix(expected) || ncol(expected) != 2L * n_cells) {
    stop(
      "Unexpected posterior_epred shape while constructing pair contrasts",
      call. = FALSE
    )
  }

  expected_before <- expected[, seq_len(n_cells), drop = FALSE]
  expected_after <- expected[, n_cells + seq_len(n_cells), drop = FALSE]
  delta <- expected_after - expected_before
  common_pairs <- complete_support_pairs(data)

  draws <- data.frame(
    .draw = rep(seq_len(nrow(delta)), each = n_cells),
    analysis = analysis_label(analysis),
    pair = rep(grids$cells$pair, times = nrow(delta)),
    context = rep(grids$cells$context, times = nrow(delta)),
    metric = rep(grids$cells$metric, times = nrow(delta)),
    metric_label = unname(METRIC_LABELS[
      rep(grids$cells$metric, times = nrow(delta))
    ]),
    expected_before = as.vector(t(expected_before)),
    expected_after = as.vector(t(expected_after)),
    delta = as.vector(t(delta)),
    stringsAsFactors = FALSE
  )
  draws$common_support_pair <- draws$pair %in% common_pairs
  validate_stage_draws(draws, analysis = analysis)
}

required_draw_columns <- function() {
  c(
    ".draw", "analysis", "pair", "context", "metric", "metric_label",
    "expected_before", "expected_after", "delta", "common_support_pair"
  )
}

coerce_logical_strict <- function(x, column = "common_support_pair") {
  if (is.logical(x)) return(x)
  value <- tolower(trimws(as.character(x)))
  result <- rep(NA, length(value))
  result[value %in% c("true", "t", "1")] <- TRUE
  result[value %in% c("false", "f", "0")] <- FALSE
  if (anyNA(result)) {
    stop(column, " must contain only TRUE or FALSE", call. = FALSE)
  }
  result
}

validate_stage_draws <- function(draws, analysis = NULL) {
  require_columns(draws, required_draw_columns(), "Posterior stage-draw cache")
  if (!nrow(draws)) {
    stop("Posterior stage-draw cache has no rows", call. = FALSE)
  }

  draws$.draw <- suppressWarnings(as.integer(draws$.draw))
  draws$expected_before <- suppressWarnings(as.numeric(draws$expected_before))
  draws$expected_after <- suppressWarnings(as.numeric(draws$expected_after))
  draws$delta <- suppressWarnings(as.numeric(draws$delta))
  draws$common_support_pair <- coerce_logical_strict(draws$common_support_pair)
  draws$metric <- canonicalise_metric(draws$metric)
  draws$metric_label <- unname(METRIC_LABELS[draws$metric])

  if (anyNA(draws$.draw) || any(draws$.draw < 1L)) {
    stop("Posterior stage-draw cache has invalid .draw identifiers", call. = FALSE)
  }
  numeric_columns <- c("expected_before", "expected_after", "delta")
  if (any(!is.finite(as.matrix(draws[numeric_columns])))) {
    stop("Posterior stage-draw cache contains non-finite values", call. = FALSE)
  }
  if (any(abs(
    draws$delta - (draws$expected_after - draws$expected_before)
  ) > 1e-10)) {
    stop(
      "Posterior stage-draw cache has delta values inconsistent with ",
      "expected_after - expected_before",
      call. = FALSE
    )
  }
  if (!setequal(unique(draws$context), c("Non-partner", "Partner"))) {
    stop(
      "Posterior stage-draw cache must contain both contexts",
      call. = FALSE
    )
  }

  analyses <- unique(as.character(draws$analysis))
  if (length(analyses) != 1L) {
    stop("Each posterior draw cache must contain one analysis", call. = FALSE)
  }
  inferred_analysis <- normalise_analysis_name(analyses)
  if (!is.null(analysis) &&
      inferred_analysis != normalise_analysis_name(analysis)) {
    stop(
      "Posterior cache analysis does not match requested analysis: ", analysis,
      call. = FALSE
    )
  }
  draws$analysis <- analysis_label(inferred_analysis)
  expected_metrics <- MODEL_METRICS[[inferred_analysis]]
  if (!setequal(unique(draws$metric), expected_metrics)) {
    stop(
      "Posterior cache metric mismatch for ", inferred_analysis,
      ". Expected: ", paste(expected_metrics, collapse = ", "),
      call. = FALSE
    )
  }

  key <- paste(
    draws$.draw, draws$pair, draws$context, draws$metric,
    sep = "\r"
  )
  if (anyDuplicated(key)) {
    stop(
      "Posterior cache contains duplicate draw x pair x context x metric rows",
      call. = FALSE
    )
  }
  expected_rows_per_draw <-
    length(unique(draws$pair)) * 2L * length(expected_metrics)
  rows_per_draw <- table(draws$.draw)
  if (any(rows_per_draw != expected_rows_per_draw)) {
    stop("Posterior cache has incomplete prediction grids", call. = FALSE)
  }
  support_by_pair <- tapply(
    draws$common_support_pair,
    draws$pair,
    function(value) length(unique(value)) == 1L && unique(value)
  )
  if (anyNA(support_by_pair) || !any(support_by_pair)) {
    stop(
      "Posterior cache has no consistently identified common-support pairs",
      call. = FALSE
    )
  }
  if (length(support_by_pair) != 3L || sum(support_by_pair) != 2L) {
    stop(
      "Current posterior caches must contain three pairs, exactly two of ",
      "which have common stage-by-context support",
      call. = FALSE
    )
  }
  draws
}

cache_formula_text <- function(formula) {
  gsub(
    "[[:space:]]+",
    " ",
    trimws(paste(deparse(formula, width.cutoff = 500L), collapse = " "))
  )
}

require_cache_provenance_packages <- function() {
  missing <- c("digest", "jsonlite")[!vapply(
    c("digest", "jsonlite"),
    requireNamespace,
    logical(1),
    quietly = TRUE
  )]
  if (length(missing)) {
    stop(
      "Packages required for posterior-cache provenance are missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

cache_sha256_file <- function(path) {
  require_cache_provenance_packages()
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("Cannot checksum missing provenance input: ", path, call. = FALSE)
  }
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

cache_sha256_object <- function(value) {
  require_cache_provenance_packages()
  encoded <- jsonlite::toJSON(
    value,
    auto_unbox = FALSE,
    dataframe = "columns",
    digits = NA,
    na = "string",
    null = "null",
    POSIXt = "ISO8601",
    pretty = FALSE
  )
  digest::digest(encoded, algo = "sha256", serialize = FALSE)
}

canonical_model_table <- function(data, analysis) {
  analysis <- normalise_analysis_name(analysis)
  columns <- c(
    "comparison_id", "pair", "context", "stage", "metric", "z_distance",
    "session_1_id", "session_2_id"
  )
  require_columns(data, columns, paste(analysis_label(analysis), "model table"))
  canonical <- data[columns]
  character_columns <- setdiff(columns, "z_distance")
  canonical[character_columns] <- lapply(
    canonical[character_columns],
    as.character
  )
  canonical$z_distance <- as.numeric(canonical$z_distance)
  if (anyNA(canonical) || any(!is.finite(canonical$z_distance))) {
    stop("Canonical model table contains missing or non-finite values", call. = FALSE)
  }
  canonical
}

canonical_model_table_sha256 <- function(data, analysis) {
  cache_sha256_object(canonical_model_table(data, analysis))
}

canonical_stage_draws_sha256 <- function(draws, analysis) {
  draws <- validate_stage_draws(draws, analysis)
  canonical <- draws[required_draw_columns()]
  canonical$.draw <- as.integer(canonical$.draw)
  canonical$expected_before <- as.numeric(canonical$expected_before)
  canonical$expected_after <- as.numeric(canonical$expected_after)
  canonical$delta <- as.numeric(canonical$delta)
  canonical$common_support_pair <- as.logical(canonical$common_support_pair)
  cache_sha256_object(canonical)
}

package_version_or_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) return("not-installed")
  as.character(utils::packageVersion(package))
}

model_software_metadata <- function(backend = "rstan") {
  backend <- as.character(backend)
  list(
    R = R.version.string,
    platform = R.version$platform,
    brms = package_version_or_missing("brms"),
    posterior = package_version_or_missing("posterior"),
    backend = backend,
    backend_package_version = package_version_or_missing(backend)
  )
}

model_sampling_metadata <- function(
    seed,
    chains,
    iter,
    warmup,
    backend,
    adapt_delta = 0.995,
    max_treedepth = 15L) {
  list(
    seed = as.integer(seed),
    chains = as.integer(chains),
    iter = as.integer(iter),
    warmup = as.integer(warmup),
    backend = as.character(backend),
    control = list(
      adapt_delta = as.numeric(adapt_delta),
      max_treedepth = as.integer(max_treedepth)
    )
  )
}

default_sampling_metadata <- function() {
  model_sampling_metadata(
    seed = 123L,
    chains = 4L,
    iter = 4000L,
    warmup = 2000L,
    backend = "rstan"
  )
}

RELEASE_SAMPLING_CONTRACT <- c(
  chains = 4L,
  iter = 4000L,
  warmup = 2000L,
  post_warmup_draws = 8000L
)

portable_draw_metadata_path <- function(path) {
  paste0(path, ".metadata.json")
}

normalise_cache_metadata <- function(metadata) {
  if (!is.list(metadata)) {
    stop("Posterior-cache metadata must be a JSON/R list", call. = FALSE)
  }
  required_top <- c(
    "schema_version", "created_utc", "analysis", "formula",
    "pair_prediction_formula", "metrics", "pairs", "common_support_pairs",
    "model_table", "sampling", "software", "draws", "provenance_sha256"
  )
  missing <- setdiff(required_top, names(metadata))
  if (length(missing)) {
    stop(
      "Posterior-cache metadata is missing field(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  required_nested <- list(
    model_table = c(
      "path", "file_sha256", "content_sha256", "bytes", "n_rows", "columns"
    ),
    sampling = c("seed", "chains", "iter", "warmup", "backend", "control"),
    software = c(
      "R", "platform", "brms", "posterior", "backend",
      "backend_package_version"
    ),
    draws = c(
      "n_draws", "n_rows", "content_sha256", "portable_sha256"
    )
  )
  for (section in names(required_nested)) {
    if (!is.list(metadata[[section]])) {
      stop("Posterior-cache metadata section is not an object: ", section,
           call. = FALSE)
    }
    section_missing <- setdiff(required_nested[[section]], names(metadata[[section]]))
    if (length(section_missing)) {
      stop(
        "Posterior-cache metadata ", section, " section is missing field(s): ",
        paste(section_missing, collapse = ", "),
        call. = FALSE
      )
    }
  }
  if (!is.list(metadata$sampling$control) ||
      !all(c("adapt_delta", "max_treedepth") %in%
           names(metadata$sampling$control))) {
    stop(
      "Posterior-cache sampling.control must record adapt_delta and ",
      "max_treedepth",
      call. = FALSE
    )
  }

  metadata$schema_version <- as.integer(metadata$schema_version)
  metadata$created_utc <- as.character(metadata$created_utc)
  metadata$analysis <- as.character(metadata$analysis)
  metadata$formula <- as.character(metadata$formula)
  metadata$pair_prediction_formula <- as.character(metadata$pair_prediction_formula)
  metadata$metrics <- as.character(unlist(metadata$metrics, use.names = FALSE))
  metadata$pairs <- as.character(unlist(metadata$pairs, use.names = FALSE))
  metadata$common_support_pairs <- as.character(unlist(
    metadata$common_support_pairs,
    use.names = FALSE
  ))
  metadata$provenance_sha256 <- as.character(metadata$provenance_sha256)

  for (field in c(
    "path", "file_sha256", "content_sha256"
  )) metadata$model_table[[field]] <- as.character(metadata$model_table[[field]])
  metadata$model_table$bytes <- as.numeric(metadata$model_table$bytes)
  metadata$model_table$n_rows <- as.integer(metadata$model_table$n_rows)
  metadata$model_table$columns <- as.character(unlist(
    metadata$model_table$columns,
    use.names = FALSE
  ))

  for (field in c("seed", "chains", "iter", "warmup")) {
    metadata$sampling[[field]] <- as.integer(metadata$sampling[[field]])
  }
  metadata$sampling$backend <- as.character(metadata$sampling$backend)
  metadata$sampling$control$adapt_delta <- as.numeric(
    metadata$sampling$control$adapt_delta
  )
  metadata$sampling$control$max_treedepth <- as.integer(
    metadata$sampling$control$max_treedepth
  )

  for (field in c(
    "R", "platform", "brms", "posterior", "backend",
    "backend_package_version"
  )) metadata$software[[field]] <- as.character(metadata$software[[field]])

  metadata$draws$n_draws <- as.integer(metadata$draws$n_draws)
  metadata$draws$n_rows <- as.integer(metadata$draws$n_rows)
  metadata$draws$content_sha256 <- as.character(metadata$draws$content_sha256)
  metadata$draws$portable_sha256 <- as.character(metadata$draws$portable_sha256)
  metadata
}

cache_metadata_fingerprint <- function(metadata) {
  metadata <- normalise_cache_metadata(metadata)
  metadata$provenance_sha256 <- NULL
  cache_sha256_object(metadata)
}

stage_draw_cache_metadata <- function(
    draws,
    data,
    analysis,
    model_table_path,
    sampling,
    software = model_software_metadata(sampling$backend),
    portable_sha256 = "") {
  analysis <- normalise_analysis_name(analysis)
  if (!file.exists(model_table_path)) {
    stop(
      "The canonical model table is required to write posterior provenance: ",
      model_table_path,
      call. = FALSE
    )
  }
  canonical_table <- canonical_model_table(data, analysis)
  metadata <- list(
    schema_version = 2L,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    analysis = analysis,
    formula = cache_formula_text(model_formula()),
    pair_prediction_formula = cache_formula_text(pair_prediction_formula()),
    metrics = MODEL_METRICS[[analysis]],
    pairs = levels(data$pair),
    common_support_pairs = complete_support_pairs(data),
    model_table = list(
      path = basename(model_table_path),
      file_sha256 = cache_sha256_file(model_table_path),
      content_sha256 = cache_sha256_object(canonical_table),
      bytes = unname(file.info(model_table_path)$size),
      n_rows = nrow(canonical_table),
      columns = names(canonical_table)
    ),
    sampling = sampling,
    software = software,
    draws = list(
      n_draws = length(unique(draws$.draw)),
      n_rows = nrow(draws),
      content_sha256 = canonical_stage_draws_sha256(draws, analysis),
      portable_sha256 = portable_sha256
    ),
    provenance_sha256 = ""
  )
  metadata <- normalise_cache_metadata(metadata)
  metadata$provenance_sha256 <- cache_metadata_fingerprint(metadata)
  metadata
}

write_gzip_csv <- function(data, path) {
  connection <- gzfile(path, open = "wt")
  tryCatch(
    utils::write.csv(data, connection, row.names = FALSE, na = ""),
    finally = close(connection)
  )
  invisible(path)
}

write_cache_metadata_json <- function(metadata, path) {
  require_cache_provenance_packages()
  metadata <- normalise_cache_metadata(metadata)
  metadata$provenance_sha256 <- cache_metadata_fingerprint(metadata)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(
    metadata,
    path,
    pretty = TRUE,
    auto_unbox = TRUE,
    digits = NA,
    null = "null"
  )
  invisible(metadata)
}

write_stage_draw_cache <- function(
    draws,
    data,
    analysis,
    rds_path,
    portable_path,
    model_table_path,
    sampling,
    software = model_software_metadata(sampling$backend)) {
  analysis <- normalise_analysis_name(analysis)
  draws <- validate_stage_draws(draws, analysis)
  dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(portable_path), recursive = TRUE, showWarnings = FALSE)

  write_gzip_csv(draws, portable_path)
  portable_sha256 <- cache_sha256_file(portable_path)
  metadata <- stage_draw_cache_metadata(
    draws,
    data,
    analysis,
    model_table_path = model_table_path,
    sampling = sampling,
    software = software,
    portable_sha256 = portable_sha256
  )
  payload <- list(
    metadata = metadata,
    draws = draws
  )
  saveRDS(payload, rds_path, compress = "xz")
  write_cache_metadata_json(metadata, portable_draw_metadata_path(portable_path))
  invisible(metadata)
}

write_portable_stage_draws <- function(draws, path, analysis = NULL) {
  metadata <- attr(draws, "cache_metadata", exact = TRUE)
  draws <- validate_stage_draws(draws, analysis)
  if (is.null(metadata)) {
    stop(
      "Refusing to write portable posterior draws without bound cache metadata",
      call. = FALSE
    )
  }
  metadata <- normalise_cache_metadata(metadata)
  expected_content <- canonical_stage_draws_sha256(draws, analysis)
  if (!identical(metadata$draws$content_sha256, expected_content)) {
    stop("Posterior draws no longer match their bound metadata", call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) {
    stop(
      "Missing portable posterior draws: ", path,
      ". Report mode will not reconstruct a release cache implicitly.",
      call. = FALSE
    )
  }
  observed_sha256 <- cache_sha256_file(path)
  if (!identical(metadata$draws$portable_sha256, observed_sha256)) {
    stop(
      "Existing portable posterior draws fail their SHA-256 check: ", path,
      call. = FALSE
    )
  }
  metadata$provenance_sha256 <- cache_metadata_fingerprint(metadata)
  write_cache_metadata_json(metadata, portable_draw_metadata_path(path))
  path
}

validate_stage_draw_cache_metadata <- function(
    metadata,
    data,
    analysis,
    model_table_path,
    draws = NULL,
    portable_path = NULL,
    enforce_release_sampling = TRUE) {
  metadata <- normalise_cache_metadata(metadata)
  if (!identical(metadata$schema_version, 2L)) {
    stop("Unsupported posterior cache schema; expected version 2", call. = FALSE)
  }
  expected_fingerprint <- cache_metadata_fingerprint(metadata)
  if (!identical(metadata$provenance_sha256, expected_fingerprint)) {
    stop(
      "Posterior-cache provenance fingerprint is invalid or was edited",
      call. = FALSE
    )
  }
  requested_analysis <- normalise_analysis_name(analysis)
  expected_fields <- list(
    analysis = requested_analysis,
    formula = cache_formula_text(model_formula()),
    pair_prediction_formula = cache_formula_text(pair_prediction_formula()),
    metrics = MODEL_METRICS[[requested_analysis]],
    pairs = levels(data$pair),
    common_support_pairs = complete_support_pairs(data)
  )
  mismatched <- names(expected_fields)[!vapply(
    names(expected_fields),
    function(field) identical(metadata[[field]], expected_fields[[field]]),
    logical(1)
  )]
  if (length(mismatched)) {
    stop(
      "Posterior cache metadata is incompatible with the current model ",
      "specification (", paste(mismatched, collapse = ", "), ")",
      call. = FALSE
    )
  }
  if (!file.exists(model_table_path)) {
    stop("Missing canonical model table: ", model_table_path, call. = FALSE)
  }
  canonical <- canonical_model_table(data, analysis)
  table_checks <- c(
    file_sha256 = identical(
      metadata$model_table$file_sha256,
      cache_sha256_file(model_table_path)
    ),
    content_sha256 = identical(
      metadata$model_table$content_sha256,
      cache_sha256_object(canonical)
    ),
    bytes = identical(
      metadata$model_table$bytes,
      as.numeric(unname(file.info(model_table_path)$size))
    ),
    n_rows = identical(metadata$model_table$n_rows, as.integer(nrow(canonical))),
    columns = identical(metadata$model_table$columns, names(canonical))
  )
  if (!all(table_checks)) {
    stop(
      "Posterior cache was generated from a different or stale canonical ",
      "model table (failed: ",
      paste(names(table_checks)[!table_checks], collapse = ", "),
      ")",
      call. = FALSE
    )
  }
  sampling <- metadata$sampling
  sampling_valid <-
    all(vapply(sampling[c("seed", "chains", "iter", "warmup")], function(x) {
      length(x) == 1L && is.finite(x) && x >= 1L
    }, logical(1))) &&
    sampling$iter > sampling$warmup &&
    sampling$backend %in% c("rstan", "cmdstanr") &&
    identical(sampling$control$adapt_delta, 0.995) &&
    identical(sampling$control$max_treedepth, 15L)
  expected_draws <- sampling$chains * (sampling$iter - sampling$warmup)
  sampling_valid <- sampling_valid &&
    identical(metadata$draws$n_draws, as.integer(expected_draws))
  if (enforce_release_sampling) {
    release_values <- c(
      chains = sampling$chains,
      iter = sampling$iter,
      warmup = sampling$warmup,
      post_warmup_draws = metadata$draws$n_draws
    )
    sampling_valid <- sampling_valid && identical(
      release_values,
      RELEASE_SAMPLING_CONTRACT
    )
  }
  software_fields <- c(
    "R", "platform", "brms", "posterior", "backend",
    "backend_package_version"
  )
  software_valid <- all(vapply(
    metadata$software[software_fields],
    function(x) length(x) == 1L && !is.na(x) && nzchar(x) && x != "not-installed",
    logical(1)
  )) && identical(metadata$software$backend, sampling$backend)
  if (!sampling_valid || !software_valid) {
    stop(
      "Posterior cache lacks valid exact sampler/software provenance",
      call. = FALSE
    )
  }
  if (!is.null(draws)) {
    draws <- validate_stage_draws(draws, analysis)
    draw_checks <- c(
      n_rows = identical(metadata$draws$n_rows, as.integer(nrow(draws))),
      n_draws = identical(
        metadata$draws$n_draws,
        as.integer(length(unique(draws$.draw)))
      ),
      content_sha256 = identical(
        metadata$draws$content_sha256,
        canonical_stage_draws_sha256(draws, analysis)
      )
    )
    if (!all(draw_checks)) {
      stop(
        "Posterior draw payload fails its bound metadata (failed: ",
        paste(names(draw_checks)[!draw_checks], collapse = ", "),
        ")",
        call. = FALSE
      )
    }
  }
  if (!is.null(portable_path)) {
    if (!file.exists(portable_path) || !identical(
      metadata$draws$portable_sha256,
      cache_sha256_file(portable_path)
    )) {
      stop(
        "Portable posterior draws fail their bound SHA-256 check: ",
        portable_path,
        call. = FALSE
      )
    }
  }
  metadata
}

validate_portable_sidecar_against_metadata <- function(portable_path, metadata) {
  metadata_path <- portable_draw_metadata_path(portable_path)
  if (!file.exists(metadata_path)) {
    stop(
      "Portable posterior draws require adjacent provenance metadata: ",
      metadata_path,
      call. = FALSE
    )
  }
  require_cache_provenance_packages()
  portable_metadata <- normalise_cache_metadata(jsonlite::read_json(
    metadata_path,
    simplifyVector = TRUE
  ))
  if (!identical(
    portable_metadata$provenance_sha256,
    cache_metadata_fingerprint(portable_metadata)
  )) {
    stop(
      "Portable posterior-cache provenance fingerprint is invalid: ",
      metadata_path,
      call. = FALSE
    )
  }
  if (!identical(portable_metadata, normalise_cache_metadata(metadata))) {
    stop(
      "RDS and portable posterior-cache metadata disagree: ", metadata_path,
      call. = FALSE
    )
  }
  invisible(TRUE)
}

read_rds_draw_cache <- function(
    path,
    analysis,
    data,
    model_table_path,
    portable_path = NULL,
    enforce_release_sampling = TRUE) {
  payload <- readRDS(path)
  if (!is.list(payload) || is.null(payload$metadata) || is.null(payload$draws)) {
    stop(
      "Current posterior cache RDS must contain metadata and draws: ", path,
      ". Legacy fitted-model RDS files are not accepted as posterior caches.",
      call. = FALSE
    )
  }
  if (enforce_release_sampling &&
      (is.null(portable_path) || !file.exists(portable_path))) {
    stop(
      "A release posterior cache requires both RDS and portable draw files",
      call. = FALSE
    )
  }
  metadata <- validate_stage_draw_cache_metadata(
    payload$metadata,
    data,
    analysis,
    model_table_path,
    draws = payload$draws,
    portable_path = if (!is.null(portable_path) && file.exists(portable_path)) {
      portable_path
    } else NULL,
    enforce_release_sampling = enforce_release_sampling
  )
  if (!is.null(portable_path) && file.exists(portable_path)) {
    validate_portable_sidecar_against_metadata(portable_path, metadata)
  }
  draws <- validate_stage_draws(payload$draws, analysis)
  attr(draws, "cache_metadata") <- metadata
  draws
}

read_portable_draw_cache <- function(
    path,
    analysis,
    data,
    model_table_path,
    enforce_release_sampling = TRUE) {
  metadata_path <- portable_draw_metadata_path(path)
  if (!file.exists(metadata_path)) {
    stop(
      "Portable posterior draws require adjacent provenance metadata: ",
      metadata_path,
      call. = FALSE
    )
  }
  require_cache_provenance_packages()
  metadata <- jsonlite::read_json(metadata_path, simplifyVector = TRUE)
  draws <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  metadata <- validate_stage_draw_cache_metadata(
    metadata,
    data,
    analysis,
    model_table_path,
    draws = draws,
    portable_path = path,
    enforce_release_sampling = enforce_release_sampling
  )
  draws <- validate_stage_draws(draws, analysis)
  attr(draws, "cache_metadata") <- metadata
  draws
}

load_cached_stage_draws <- function(
    rds_path,
    portable_path,
    analysis,
    data = NULL,
    model_table_path = NULL,
    enforce_release_sampling = TRUE) {
  if (file.exists(rds_path)) {
    if (is.null(data) || is.null(model_table_path)) {
      stop(
        "Current model data and its canonical CSV path are required to ",
        "validate posterior provenance",
        call. = FALSE
      )
    }
    return(read_rds_draw_cache(
      rds_path,
      analysis,
      data,
      model_table_path,
      portable_path,
      enforce_release_sampling = enforce_release_sampling
    ))
  }
  if (file.exists(portable_path)) {
    if (is.null(data) || is.null(model_table_path)) {
      stop(
        "Current model data and its canonical CSV path are required to ",
        "validate posterior provenance",
        call. = FALSE
      )
    }
    message("RDS cache absent; using portable posterior draws: ", portable_path)
    return(read_portable_draw_cache(
      portable_path,
      analysis,
      data,
      model_table_path,
      enforce_release_sampling = enforce_release_sampling
    ))
  }
  stop(
    "Missing current posterior-draw cache for ",
    analysis_label(analysis), ". Checked:\n  - ", rds_path,
    "\n  - ", portable_path,
    "\nRun Rscript R/run_models.R --refit only if an intentional, expensive ",
    "model refit is required. Cached/report-only mode never fits silently.",
    call. = FALSE
  )
}

validate_cache_against_data <- function(draws, data, analysis) {
  analysis <- normalise_analysis_name(analysis)
  draws <- validate_stage_draws(draws, analysis)
  data_pairs <- sort(levels(data$pair))
  draw_pairs <- sort(unique(draws$pair))
  if (!identical(data_pairs, draw_pairs)) {
    stop(
      "Posterior cache pair IDs do not match the current ", analysis,
      " model table",
      call. = FALSE
    )
  }
  expected_common <- sort(complete_support_pairs(data))
  cached_common <- sort(unique(draws$pair[draws$common_support_pair]))
  if (!identical(expected_common, cached_common)) {
    stop(
      "Posterior cache common-support pairs do not match the current ",
      analysis, " model table",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

finite_stat <- function(values, fun) {
  values <- values[is.finite(values)]
  if (!length(values)) return(NA_real_)
  match.fun(fun)(values)
}

collect_fit_diagnostics <- function(fit, model_name, max_treedepth = 15L) {
  if (!requireNamespace("posterior", quietly = TRUE) ||
      !requireNamespace("brms", quietly = TRUE)) {
    stop("Packages 'posterior' and 'brms' are required for diagnostics")
  }
  draw_summary <- posterior::summarise_draws(
    posterior::as_draws_array(fit)
  )
  sampler <- brms::nuts_params(fit)
  data.frame(
    model = model_name,
    max_rhat = finite_stat(draw_summary$rhat, "max"),
    min_ess_bulk = finite_stat(draw_summary$ess_bulk, "min"),
    min_ess_tail = finite_stat(draw_summary$ess_tail, "min"),
    divergences = sum(
      sampler$Parameter == "divergent__" & sampler$Value == 1,
      na.rm = TRUE
    ),
    max_treedepth_hits = sum(
      sampler$Parameter == "treedepth__" &
        sampler$Value >= max_treedepth,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}

save_diagnostic_plots <- function(fit, output_dir, model_name) {
  if (!requireNamespace("bayesplot", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("posterior", quietly = TRUE) ||
      !requireNamespace("brms", quietly = TRUE)) {
    stop(
      "brms, posterior, bayesplot, and ggplot2 are required for ",
      "diagnostic plots",
      call. = FALSE
    )
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ppc_path <- file.path(output_dir, paste0(model_name, "_pp_check.png"))
  trace_path <- file.path(output_dir, paste0(model_name, "_trace.png"))
  pp <- brms::pp_check(fit, ndraws = 200)
  ggplot2::ggsave(
    ppc_path,
    pp,
    width = 9,
    height = 6,
    dpi = 180
  )

  draw_array <- posterior::as_draws_array(fit)
  parameters <- posterior::variables(draw_array)
  trace_parameters <- unique(c(
    grep("^b_", parameters, value = TRUE),
    grep("^sd_", parameters, value = TRUE),
    "sigma"
  ))
  trace_parameters <- intersect(trace_parameters, parameters)
  trace <- bayesplot::mcmc_trace(draw_array, pars = trace_parameters)
  ggplot2::ggsave(
    trace_path,
    trace,
    width = 14,
    height = max(8, length(trace_parameters) * 0.45),
    dpi = 160,
    limitsize = FALSE
  )
  if (!all(file.exists(c(ppc_path, trace_path)))) {
    stop("Required refit diagnostic plots were not written", call. = FALSE)
  }
  invisible(c(ppc = ppc_path, trace = trace_path))
}
