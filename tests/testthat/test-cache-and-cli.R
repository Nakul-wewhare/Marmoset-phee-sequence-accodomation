testthat::test_that("missing posterior caches never trigger a fit", {
  location <- tempfile("missing-cache-")
  testthat::expect_error(
    load_cached_stage_draws(
      paste0(location, ".rds"),
      paste0(location, ".csv.gz"),
      "sequence"
    ),
    "never fits silently"
  )
})

testthat::test_that("portable posterior draws are a valid report-only cache", {
  directory <- tempfile("draw-cache-")
  dir.create(directory)
  model_table_path <- file.path(directory, "model_sequence.csv")
  utils::write.csv(
    synthetic_model_data("sequence"),
    model_table_path,
    row.names = FALSE
  )
  data <- read_model_data(
    model_table_path,
    "sequence",
    check_standardisation = FALSE
  )
  draws <- synthetic_stage_draws("sequence")
  rds_path <- file.path(directory, "sequence.rds")
  portable_path <- file.path(directory, "sequence.csv.gz")
  write_stage_draw_cache(
    draws,
    data,
    "sequence",
    rds_path,
    portable_path,
    model_table_path = model_table_path,
    sampling = synthetic_sampling_metadata()
  )
  testthat::expect_error(
    load_cached_stage_draws(
      rds_path,
      portable_path,
      "sequence",
      data = data,
      model_table_path = model_table_path
    ),
    "sampler/software provenance"
  )
  rds_loaded <- load_cached_stage_draws(
    rds_path,
    portable_path,
    "sequence",
    data = data,
    model_table_path = model_table_path,
    enforce_release_sampling = FALSE
  )
  testthat::expect_equal(rds_loaded$delta, draws$delta)
  unlink(rds_path)
  loaded <- load_cached_stage_draws(
    rds_path,
    portable_path,
    "sequence",
    data = data,
    model_table_path = model_table_path,
    enforce_release_sampling = FALSE
  )
  testthat::expect_equal(loaded$delta, draws$delta)
  testthat::expect_true(file.exists(portable_draw_metadata_path(portable_path)))
  testthat::expect_silent(
    validate_cache_against_data(loaded, data, "sequence")
  )
})

testthat::test_that("portable cache metadata binds the exact model table", {
  directory <- tempfile("bound-cache-")
  dir.create(directory)
  model_table_path <- file.path(directory, "model_sequence.csv")
  utils::write.csv(
    synthetic_model_data("sequence"),
    model_table_path,
    row.names = FALSE
  )
  data <- read_model_data(
    model_table_path,
    "sequence",
    check_standardisation = FALSE
  )
  rds_path <- file.path(directory, "sequence.rds")
  portable_path <- file.path(directory, "sequence.csv.gz")
  metadata <- write_stage_draw_cache(
    synthetic_stage_draws("sequence"),
    data,
    "sequence",
    rds_path,
    portable_path,
    model_table_path = model_table_path,
    sampling = synthetic_sampling_metadata(seed = 42L)
  )
  testthat::expect_identical(metadata$schema_version, 2L)
  testthat::expect_identical(metadata$sampling$seed, 42L)
  testthat::expect_identical(metadata$formula, cache_formula_text(model_formula()))
  testthat::expect_identical(metadata$metrics, SEQUENCE_METRICS)
  testthat::expect_match(metadata$model_table$file_sha256, "^[0-9a-f]{64}$")
  testthat::expect_match(metadata$provenance_sha256, "^[0-9a-f]{64}$")
  testthat::expect_false(metadata$software$brms == "not-installed")

  unlink(rds_path)
  write("", model_table_path, append = TRUE)
  testthat::expect_error(
    load_cached_stage_draws(
      rds_path,
      portable_path,
      "sequence",
      data = data,
      model_table_path = model_table_path,
      enforce_release_sampling = FALSE
    ),
    "different or stale canonical model table"
  )
})

testthat::test_that("cache validation rejects semantic input and metadata drift", {
  directory <- tempfile("mismatch-cache-")
  dir.create(directory)
  model_table_path <- file.path(directory, "model_sequence.csv")
  utils::write.csv(
    synthetic_model_data("sequence"),
    model_table_path,
    row.names = FALSE
  )
  data <- read_model_data(
    model_table_path,
    "sequence",
    check_standardisation = FALSE
  )
  rds_path <- file.path(directory, "sequence.rds")
  portable_path <- file.path(directory, "sequence.csv.gz")
  write_stage_draw_cache(
    synthetic_stage_draws("sequence"),
    data,
    "sequence",
    rds_path,
    portable_path,
    model_table_path = model_table_path,
    sampling = synthetic_sampling_metadata()
  )

  mismatched_data <- data
  mismatched_data$z_distance[[1]] <- mismatched_data$z_distance[[1]] + 0.25
  testthat::expect_error(
    load_cached_stage_draws(
      rds_path,
      portable_path,
      "sequence",
      data = mismatched_data,
      model_table_path = model_table_path,
      enforce_release_sampling = FALSE
    ),
    "different or stale canonical model table"
  )

  unlink(rds_path)
  metadata_path <- portable_draw_metadata_path(portable_path)
  metadata <- jsonlite::read_json(metadata_path, simplifyVector = TRUE)
  metadata$formula <- "z_distance ~ 1"
  jsonlite::write_json(metadata, metadata_path, auto_unbox = TRUE, pretty = TRUE)
  testthat::expect_error(
    load_cached_stage_draws(
      rds_path,
      portable_path,
      "sequence",
      data = data,
      model_table_path = model_table_path,
      enforce_release_sampling = FALSE
    ),
    "provenance fingerprint is invalid"
  )
})

testthat::test_that("portable draw bytes are cryptographically bound", {
  directory <- tempfile("draw-tamper-")
  dir.create(directory)
  model_table_path <- file.path(directory, "model_call.csv")
  utils::write.csv(
    synthetic_model_data("call"),
    model_table_path,
    row.names = FALSE
  )
  data <- read_model_data(
    model_table_path,
    "call",
    check_standardisation = FALSE
  )
  rds_path <- file.path(directory, "call.rds")
  portable_path <- file.path(directory, "call.csv.gz")
  write_stage_draw_cache(
    synthetic_stage_draws("call"),
    data,
    "call",
    rds_path,
    portable_path,
    model_table_path = model_table_path,
    sampling = synthetic_sampling_metadata()
  )
  unlink(rds_path)
  tampered <- utils::read.csv(portable_path, stringsAsFactors = FALSE)
  tampered$expected_after[[1]] <- tampered$expected_after[[1]] + 1
  tampered$delta[[1]] <- tampered$delta[[1]] + 1
  write_gzip_csv(tampered, portable_path)
  testthat::expect_error(
    load_cached_stage_draws(
      rds_path,
      portable_path,
      "call",
      data = data,
      model_table_path = model_table_path,
      enforce_release_sampling = FALSE
    ),
    "bound metadata"
  )
})

testthat::test_that("legacy model RDS files are not accepted as draw caches", {
  directory <- tempfile("legacy-cache-")
  dir.create(directory)
  model_table_path <- file.path(directory, "model_call.csv")
  utils::write.csv(synthetic_model_data("call"), model_table_path, row.names = FALSE)
  data <- read_model_data(model_table_path, "call", check_standardisation = FALSE)
  rds_path <- tempfile(fileext = ".rds")
  saveRDS(list(old_brms_fit = TRUE), rds_path)
  testthat::expect_error(
    load_cached_stage_draws(
      rds_path,
      tempfile(),
      "call",
      data = data,
      model_table_path = model_table_path,
      enforce_release_sampling = FALSE
    ),
    "Legacy fitted-model RDS files are not accepted"
  )
})

testthat::test_that("CLI defaults to cached reporting and gates expensive work", {
  config <- parse_model_cli(character(), project_root = tempdir())
  testthat::expect_false(config$refit)
  testthat::expect_false(config$validate_only)
  testthat::expect_match(
    config$sequence_data,
    "data/derived/model_sequence.csv$"
  )
  testthat::expect_match(config$manifest, "results/model_manifest.json$")
  testthat::expect_identical(
    formals(save_stage_effect_figure)$basename,
    "figure_3_regenerated"
  )

  refit <- parse_model_cli(c("--refit"), project_root = tempdir())
  testthat::expect_true(refit$refit)
  testthat::expect_true(refit$diagnostics)
  testthat::expect_identical(refit$chains, 4L)
  testthat::expect_error(
    parse_model_cli(
      c("--refit", "--chains=2", "--iter=20", "--warmup=10"),
      project_root = tempdir()
    ),
    "exactly 8,000"
  )
  testthat::expect_error(
    parse_model_cli(c("--diagnostics"), project_root = tempdir()),
    "only valid with --refit"
  )
  testthat::expect_error(
    parse_model_cli(
      c("--refit", "--validate-only"),
      project_root = tempdir()
    ),
    "Choose only one"
  )
})

testthat::test_that("diagnostics use the canonical notebook schema", {
  diagnostics <- data.frame(
    analysis = c("Sequence structure", "Call structure"),
    max_rhat = c(1.0042, 1.0041),
    min_bulk_ess = c(1905, 1177),
    min_tail_ess = c(2751, 2049),
    n_divergent = c(0, 0),
    n_max_treedepth = c(0, 0),
    ppc_status = c("passed", "passed"),
    ppc_path = c("sequence_ppc.png", "call_ppc.png"),
    trace_path = c("sequence_trace.png", "call_trace.png"),
    diagnostic_status = c("pass", "pass"),
    stringsAsFactors = FALSE
  )
  validated <- validate_diagnostics_table(diagnostics)
  testthat::expect_identical(validated$min_bulk_ess, c(1905, 1177))

  diagnostics$min_tail_ess[[2]] <- 999
  testthat::expect_error(
    validate_diagnostics_table(diagnostics),
    "ESS >= 1000"
  )
})

testthat::test_that("cached mode renders reports end to end without fitting", {
  project <- tempfile("cached-report-")
  directories <- c(
    "data/derived",
    "data/cache/posterior_draws",
    "results/posterior_draws",
    "results/tables"
  )
  for (directory in directories) {
    dir.create(file.path(project, directory), recursive = TRUE)
  }

  inputs <- list(
    sequence = synthetic_public_model_data("sequence"),
    call = synthetic_public_model_data("call")
  )
  for (analysis in names(inputs)) {
    input_path <- file.path(
      project,
      "data",
      "derived",
      paste0("model_", analysis, ".csv")
    )
    utils::write.csv(inputs[[analysis]], input_path, row.names = FALSE)
    prepared <- read_model_data(input_path, analysis)
    write_stage_draw_cache(
      synthetic_stage_draws(analysis, n_draws = 10L),
      prepared,
      analysis,
      file.path(
        project,
        "data",
        "cache",
        "posterior_draws",
        paste0(analysis, "_stage_effect_draws.rds")
      ),
      file.path(
        project,
        "results",
        "posterior_draws",
        paste0(analysis, "_stage_effect_draws.csv.gz")
      ),
      model_table_path = input_path,
      sampling = synthetic_sampling_metadata(n_draws = 10L)
    )
  }

  diagnostic_paths <- c(
    "results/figures/main/diagnostics/sequence/sequence_pp_check.png",
    "results/figures/main/diagnostics/call/call_pp_check.png",
    "results/figures/main/diagnostics/sequence/sequence_trace.png",
    "results/figures/main/diagnostics/call/call_trace.png"
  )
  for (path in diagnostic_paths) {
    dir.create(dirname(file.path(project, path)), recursive = TRUE, showWarnings = FALSE)
    file.create(file.path(project, path))
  }

  diagnostics <- data.frame(
    analysis = c("Sequence structure", "Call structure"),
    max_rhat = c(1.0042, 1.0041),
    min_bulk_ess = c(1905, 1177),
    min_tail_ess = c(2751, 2049),
    n_divergent = c(0, 0),
    n_max_treedepth = c(0, 0),
    ppc_status = c("passed", "passed"),
    ppc_path = diagnostic_paths[1:2],
    trace_path = diagnostic_paths[3:4],
    diagnostic_status = c("pass", "pass"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    diagnostics,
    file.path(
      project,
      "data",
      "cache",
      "posterior_draws",
      "model_diagnostics.csv"
    ),
    row.names = FALSE
  )

  config <- parse_model_cli(
    c("--skip-manuscript-count-check", "--skip-figures"),
    project_root = project
  )
  config$enforce_release_sampling <- FALSE
  result <- run_model_pipeline(config)
  testthat::expect_identical(result$mode, "cached-report")
  testthat::expect_true(file.exists(file.path(
    project,
    "results",
    "model_manifest.json"
  )))
  testthat::expect_true(file.exists(file.path(
    project,
    "results",
    "tables",
    "psi_combined.csv"
  )))
})
