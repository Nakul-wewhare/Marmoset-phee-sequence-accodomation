#!/usr/bin/env Rscript

# Cached posterior reporting is the default. The Stan models are fitted only
# when this script receives the explicit --refit flag.

current_file <- function() {
  frames <- sys.frames()
  if (length(frames)) {
    for (index in rev(seq_along(frames))) {
      path <- frames[[index]]$ofile
      if (!is.null(path) && nzchar(path)) return(normalizePath(path))
    }
  }
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]])))
  }
  normalizePath("R/run_models.R")
}

.model_r_dir <- dirname(current_file())
source(file.path(.model_r_dir, "model_spec.R"))
source(file.path(.model_r_dir, "model_engine.R"))
source(file.path(.model_r_dir, "posterior_results.R"))
rm(.model_r_dir)

MIN_MODEL_ESS <- 1000

model_cli_defaults <- function(project_root = dirname(dirname(current_file()))) {
  project_root <- normalizePath(project_root, mustWork = FALSE)
  cache_dir <- file.path(project_root, "data", "cache", "posterior_draws")
  draws_dir <- file.path(project_root, "results", "posterior_draws")
  table_dir <- file.path(project_root, "results", "tables")
  figure_dir <- file.path(project_root, "results", "figures", "main")
  list(
    project_root = project_root,
    sequence_data = file.path(project_root, "data", "derived", "model_sequence.csv"),
    call_data = file.path(project_root, "data", "derived", "model_call.csv"),
    cache_dir = cache_dir,
    draws_dir = draws_dir,
    table_dir = table_dir,
    figure_dir = figure_dir,
    manifest = file.path(project_root, "results", "model_manifest.json"),
    diagnostics_cache = file.path(cache_dir, "model_diagnostics.csv"),
    refit = FALSE,
    validate_only = FALSE,
    validate_cache_only = FALSE,
    diagnostics = FALSE,
    enforce_release_sampling = TRUE,
    skip_figures = FALSE,
    check_counts = TRUE,
    chains = 4L,
    iter = 4000L,
    warmup = 2000L,
    cores = default_model_cores(),
    seed = 123L,
    backend = "rstan",
    help = FALSE
  )
}

model_cli_usage <- function() {
  paste(
    "Usage: Rscript R/run_models.R [options]",
    "",
    "Default: validate current model tables and cached posterior draws, then",
    "write manuscript summaries and Figure 3. It never fits a model.",
    "",
    "Modes:",
    "  --refit                 Explicitly fit both expensive brms models",
    "  --validate-only         Validate input data and Stan data; do not fit",
    "  --validate-cache-only   Validate inputs, portable/RDS draws, diagnostics",
    "  --report-only           Explicit alias for the safe default mode",
    "",
    "Paths (use --name=PATH):",
    "  --sequence-data, --call-data, --cache-dir, --draws-dir, --table-dir,",
    "  --figure-dir, --manifest",
    "",
    "Refit settings:",
    "  --chains, --iter, --warmup, --cores, --seed, --backend",
    "  Every --refit saves mandatory posterior-predictive and trace plots.",
    "  --diagnostics is accepted with --refit for compatibility but is redundant.",
    "",
    "Other:",
    "  --skip-figures",
    "  --skip-manuscript-count-check",
    "  --help",
    sep = "\n"
  )
}

parse_model_cli <- function(
    args = commandArgs(trailingOnly = TRUE),
    project_root = dirname(dirname(current_file()))) {
  config <- model_cli_defaults(project_root)
  flag_map <- c(
    "--refit" = "refit",
    "--validate-only" = "validate_only",
    "--validate-cache-only" = "validate_cache_only",
    "--diagnostics" = "diagnostics",
    "--skip-figures" = "skip_figures",
    "--help" = "help"
  )
  path_map <- c(
    "sequence-data" = "sequence_data",
    "call-data" = "call_data",
    "cache-dir" = "cache_dir",
    "draws-dir" = "draws_dir",
    "table-dir" = "table_dir",
    "figure-dir" = "figure_dir",
    "manifest" = "manifest"
  )
  integer_map <- c(
    chains = "chains",
    iter = "iter",
    warmup = "warmup",
    cores = "cores",
    seed = "seed"
  )

  for (argument in args) {
    if (argument == "--report-only") next
    if (argument == "--skip-manuscript-count-check") {
      config$check_counts <- FALSE
      next
    }
    if (argument %in% names(flag_map)) {
      config[[flag_map[[argument]]]] <- TRUE
      next
    }
    if (!startsWith(argument, "--") || !grepl("=", argument, fixed = TRUE)) {
      stop("Unknown command-line argument: ", argument, call. = FALSE)
    }
    pieces <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1]]
    name <- pieces[[1]]
    value <- paste(pieces[-1], collapse = "=")
    if (!nzchar(value)) stop("Missing value for --", name, call. = FALSE)
    if (name %in% names(path_map)) {
      config[[path_map[[name]]]] <- value
    } else if (name %in% names(integer_map)) {
      parsed <- suppressWarnings(as.integer(value))
      if (is.na(parsed) || parsed < 1L) {
        stop("--", name, " must be a positive integer", call. = FALSE)
      }
      config[[integer_map[[name]]]] <- parsed
    } else if (name == "backend") {
      if (!value %in% c("rstan", "cmdstanr")) {
        stop("--backend must be rstan or cmdstanr", call. = FALSE)
      }
      config$backend <- value
    } else {
      stop("Unknown command-line option: --", name, call. = FALSE)
    }
  }

  active_modes <- sum(c(
    refit = config$refit,
    validate_only = config$validate_only,
    validate_cache_only = config$validate_cache_only
  ))
  if (active_modes > 1L) {
    stop(
      "Choose only one of --refit, --validate-only, or --validate-cache-only",
      call. = FALSE
    )
  }
  if (config$diagnostics && !config$refit) {
    stop("--diagnostics is only valid with --refit", call. = FALSE)
  }
  if (config$refit) config$diagnostics <- TRUE
  if (config$refit && config$enforce_release_sampling) {
    requested <- c(
      chains = config$chains,
      iter = config$iter,
      warmup = config$warmup
    )
    canonical <- RELEASE_SAMPLING_CONTRACT[c("chains", "iter", "warmup")]
    if (!identical(requested, canonical)) {
      stop(
        "Release --refit requires --chains=4, --iter=4000, and ",
        "--warmup=2000 so the cache contains exactly 8,000 post-warmup ",
        "draws",
        call. = FALSE
      )
    }
  }
  if (config$iter <= config$warmup) {
    stop("--iter must be greater than --warmup", call. = FALSE)
  }

  # Recompute dependent defaults when the user changes a parent directory.
  config$diagnostics_cache <- file.path(config$cache_dir, "model_diagnostics.csv")
  config
}

draw_paths <- function(config, analysis) {
  analysis <- normalise_analysis_name(analysis)
  list(
    rds = file.path(
      config$cache_dir,
      paste0(analysis, "_stage_effect_draws.rds")
    ),
    portable = file.path(
      config$draws_dir,
      paste0(analysis, "_stage_effect_draws.csv.gz")
    ),
    model = file.path(
      config$cache_dir,
      paste0(analysis, "_interaction_brmsfit.rds")
    )
  )
}

make_diagnostics_row <- function(
    fit,
    analysis,
    ppc_path,
    trace_path) {
  raw <- collect_fit_diagnostics(
    fit,
    paste0(normalise_analysis_name(analysis), "_interaction_model")
  )
  status <- if (
    raw$max_rhat <= 1.01 &&
      raw$min_ess_bulk >= MIN_MODEL_ESS &&
      raw$min_ess_tail >= MIN_MODEL_ESS &&
      raw$divergences == 0 &&
      raw$max_treedepth_hits == 0 &&
      file.exists(ppc_path) &&
      file.exists(trace_path)
  ) "pass" else "review"
  result <- data.frame(
    analysis = analysis_label(analysis),
    max_rhat = raw$max_rhat,
    min_bulk_ess = raw$min_ess_bulk,
    min_tail_ess = raw$min_ess_tail,
    n_divergent = raw$divergences,
    n_max_treedepth = raw$max_treedepth_hits,
    ppc_status = "generated",
    ppc_path = ppc_path,
    trace_path = trace_path,
    diagnostic_status = status,
    stringsAsFactors = FALSE
  )
  if (status != "pass") {
    stop(
      analysis_label(analysis),
      " refit failed the release diagnostic contract: R-hat <= 1.01, ",
      "bulk/tail ESS >= ", MIN_MODEL_ESS,
      ", zero divergences/tree-depth hits, and both PPC/trace plots are ",
      "mandatory",
      call. = FALSE
    )
  }
  result
}

project_relative_path <- function(path, project_root) {
  root <- normalizePath(project_root, mustWork = TRUE)
  normalized <- normalizePath(path, mustWork = TRUE)
  prefix <- paste0(root, .Platform$file.sep)
  if (startsWith(normalized, prefix)) {
    return(substring(normalized, nchar(prefix) + 1L))
  }
  normalized
}

validate_refit_diagnostic_dependencies <- function() {
  required <- c("brms", "posterior", "bayesplot", "ggplot2", "digest", "jsonlite")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Cannot start --refit because mandatory fitting/diagnostic/provenance ",
      "packages are missing: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

fit_and_cache_analysis <- function(data, analysis, config) {
  paths <- draw_paths(config, analysis)
  analysis <- normalise_analysis_name(analysis)
  model_table_path <- if (analysis == "sequence") {
    config$sequence_data
  } else {
    config$call_data
  }
  fit <- fit_model(
    data,
    analysis = analysis,
    chains = config$chains,
    iter = config$iter,
    warmup = config$warmup,
    cores = config$cores,
    seed = config$seed,
    backend = config$backend
  )
  draws <- extract_pair_stage_draws(fit, data, analysis)
  diagnostic_dir <- file.path(config$figure_dir, "diagnostics", analysis)
  diagnostic_paths <- save_diagnostic_plots(
    fit,
    diagnostic_dir,
    paste0(analysis, "_interaction_model")
  )
  diagnostics <- make_diagnostics_row(
    fit,
    analysis,
    ppc_path = diagnostic_paths[["ppc"]],
    trace_path = diagnostic_paths[["trace"]]
  )
  diagnostics$ppc_path <- project_relative_path(
    diagnostics$ppc_path,
    config$project_root
  )
  diagnostics$trace_path <- project_relative_path(
    diagnostics$trace_path,
    config$project_root
  )

  dir.create(config$cache_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(fit, paths$model, compress = FALSE)
  sampling <- model_sampling_metadata(
    seed = config$seed,
    chains = config$chains,
    iter = config$iter,
    warmup = config$warmup,
    backend = config$backend
  )
  metadata <- write_stage_draw_cache(
    draws,
    data,
    analysis,
    rds_path = paths$rds,
    portable_path = paths$portable,
    model_table_path = model_table_path,
    sampling = sampling,
    software = model_software_metadata(config$backend)
  )
  attr(draws, "cache_metadata") <- metadata
  list(
    draws = draws,
    diagnostics = diagnostics,
    paths = paths
  )
}

read_and_validate_inputs <- function(config) {
  sequence_data <- read_model_data(config$sequence_data, "sequence")
  call_data <- read_model_data(config$call_data, "call")
  if (config$check_counts) {
    validate_expected_counts(sequence_data, "sequence")
    validate_expected_counts(call_data, "call")
  }
  list(sequence = sequence_data, call = call_data)
}

run_model_pipeline <- function(config) {
  inputs <- read_and_validate_inputs(config)

  if (config$validate_only) {
    message("Validation-only mode: constructing Stan data without sampling.")
    validate_brms_model(inputs$sequence)
    validate_brms_model(inputs$call)
    message("Validation successful. No models were fitted and no caches changed.")
    return(invisible(list(mode = "validate-only")))
  }

  if (config$refit) {
    validate_refit_diagnostic_dependencies()
    sequence <- fit_and_cache_analysis(inputs$sequence, "sequence", config)
    call <- fit_and_cache_analysis(inputs$call, "call", config)
    diagnostics <- rbind(sequence$diagnostics, call$diagnostics)
    diagnostics <- validate_diagnostics_table(diagnostics)
    dir.create(config$cache_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(
      diagnostics,
      config$diagnostics_cache,
      row.names = FALSE
    )
    sequence_draws <- sequence$draws
    call_draws <- call$draws
    mode <- "explicit-refit"
  } else {
    message(
      "Cached/report-only mode: no Stan sampling will run. ",
      "Use --refit explicitly to replace current caches."
    )
    sequence_paths <- draw_paths(config, "sequence")
    call_paths <- draw_paths(config, "call")
    sequence_draws <- load_cached_stage_draws(
      sequence_paths$rds,
      sequence_paths$portable,
      "sequence",
      data = inputs$sequence,
      model_table_path = config$sequence_data,
      enforce_release_sampling = config$enforce_release_sampling
    )
    call_draws <- load_cached_stage_draws(
      call_paths$rds,
      call_paths$portable,
      "call",
      data = inputs$call,
      model_table_path = config$call_data,
      enforce_release_sampling = config$enforce_release_sampling
    )
    diagnostics <- load_cached_diagnostics(
      config$diagnostics_cache,
      file.path(config$table_dir, "model_diagnostics.csv"),
      project_root = config$project_root
    )
    mode <- "cached-report"
  }

  validate_cache_against_data(sequence_draws, inputs$sequence, "sequence")
  validate_cache_against_data(call_draws, inputs$call, "call")

  if (config$validate_cache_only) {
    message(
      "Cache validation successful. No models were fitted and no outputs changed."
    )
    return(invisible(list(mode = "validate-cache-only")))
  }

  sequence_portable <- write_portable_stage_draws(
    sequence_draws,
    draw_paths(config, "sequence")$portable,
    "sequence"
  )
  call_portable <- write_portable_stage_draws(
    call_draws,
    draw_paths(config, "call")$portable,
    "call"
  )

  products <- combine_posterior_products(
    build_posterior_products(sequence_draws),
    build_posterior_products(call_draws)
  )
  summaries <- summarise_posterior_products(products)
  counts <- rbind(
    model_data_counts(inputs$sequence, "sequence"),
    model_data_counts(inputs$call, "call")
  )
  summaries <- add_manuscript_counts(summaries, counts)
  table_paths <- write_summary_tables(summaries, config$table_dir)
  diagnostics_path <- write_report_diagnostics(
    diagnostics,
    file.path(config$table_dir, "model_diagnostics.csv")
  )

  figure_paths <- character()
  if (!config$skip_figures) {
    figure_paths <- save_stage_effect_figure(products, config$figure_dir)
  }

  output_paths <- c(
    sequence_portable,
    call_portable,
    portable_draw_metadata_path(sequence_portable),
    portable_draw_metadata_path(call_portable),
    table_paths,
    diagnostics = diagnostics_path,
    diagnostic_artifact_paths(diagnostics, config$project_root),
    figure_paths
  )
  manifest <- write_results_manifest(
    output_paths,
    config$manifest,
    config$project_root,
    mode
  )
  message("Model reporting complete. Manifest: ", manifest)
  invisible(list(
    mode = mode,
    draws = list(sequence = sequence_draws, call = call_draws),
    products = products,
    summaries = summaries,
    diagnostics = diagnostics,
    outputs = c(output_paths, manifest = manifest)
  ))
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  config <- parse_model_cli(args)
  if (config$help) {
    cat(model_cli_usage(), "\n")
    return(invisible(0L))
  }
  run_model_pipeline(config)
  invisible(0L)
}

if (sys.nframe() == 0L) {
  main()
}
