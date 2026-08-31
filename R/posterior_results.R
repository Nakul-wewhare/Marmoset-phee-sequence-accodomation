# Posterior summaries, manuscript tables, figures, and provenance manifests.
# These helpers operate on portable pair-level contrast draws and never fit a
# model.

support_sets <- function(draws) {
  common_pairs <- sort(unique(draws$pair[draws$common_support_pair]))
  list(
    all_pairs = sort(unique(draws$pair)),
    common_support_A_B = common_pairs
  )
}

collapse_pair_equal <- function(draws, pairs, support) {
  selected <- draws[draws$pair %in% pairs, , drop = FALSE]
  if (!nrow(selected)) {
    stop("No posterior rows remain for support set: ", support, call. = FALSE)
  }
  group_columns <- c(
    ".draw", "analysis", "context", "metric", "metric_label"
  )
  collapsed <- stats::aggregate(
    selected[c("expected_before", "expected_after", "delta")],
    by = selected[group_columns],
    FUN = mean
  )
  collapsed$support <- support
  collapsed$n_pairs <- length(pairs)
  collapsed$pairs_included <- paste(sort(pairs), collapse = ";")
  collapsed[c(
    ".draw", "analysis", "support", "context", "metric", "metric_label",
    "expected_before", "expected_after", "delta", "n_pairs",
    "pairs_included"
  )]
}

pair_equal_metric_draws <- function(draws) {
  draws <- validate_stage_draws(draws)
  sets <- support_sets(draws)
  do.call(rbind, lapply(names(sets), function(name) {
    collapse_pair_equal(draws, sets[[name]], name)
  }))
}

combine_metrics_equal_weight <- function(metric_draws) {
  required <- c(
    ".draw", "analysis", "support", "context", "metric", "delta",
    "n_pairs", "pairs_included"
  )
  require_columns(metric_draws, required, "Pair-equal metric draws")
  group_columns <- c(
    ".draw", "analysis", "support", "context", "n_pairs",
    "pairs_included"
  )
  combined <- stats::aggregate(
    metric_draws[c("expected_before", "expected_after", "delta")],
    by = metric_draws[group_columns],
    FUN = mean
  )
  combined$metric <- "combined"
  combined$metric_label <- "Combined (equal metric weight)"
  combined$n_metrics <- 4L
  combined[c(
    ".draw", "analysis", "support", "context", "metric", "metric_label",
    "expected_before", "expected_after", "delta", "n_pairs",
    "pairs_included", "n_metrics"
  )]
}

partner_minus_nonpartner <- function(stage_draws) {
  required <- c(
    ".draw", "analysis", "support", "context", "metric", "metric_label",
    "delta", "n_pairs", "pairs_included"
  )
  require_columns(stage_draws, required, "Stage-effect draws")
  partner <- stage_draws[
    stage_draws$context == "Partner",
    setdiff(names(stage_draws), "context"),
    drop = FALSE
  ]
  nonpartner <- stage_draws[
    stage_draws$context == "Non-partner",
    setdiff(names(stage_draws), "context"),
    drop = FALSE
  ]
  join_columns <- c(
    ".draw", "analysis", "support", "metric", "metric_label",
    "n_pairs", "pairs_included"
  )
  joined <- merge(
    partner[c(join_columns, "delta")],
    nonpartner[c(join_columns, "delta")],
    by = join_columns,
    suffixes = c("_partner", "_nonpartner"),
    all = FALSE,
    sort = FALSE
  )
  joined$psi <- joined$delta_partner - joined$delta_nonpartner
  joined[c(join_columns, "delta_partner", "delta_nonpartner", "psi")]
}

build_posterior_products <- function(draws) {
  metric_stage <- pair_equal_metric_draws(draws)
  combined_stage <- combine_metrics_equal_weight(metric_stage)
  combined_psi <- partner_minus_nonpartner(combined_stage)
  combined_psi$n_metrics <- 4L
  list(
    metric_stage = metric_stage,
    combined_stage = combined_stage,
    metric_psi = partner_minus_nonpartner(metric_stage),
    combined_psi = combined_psi
  )
}

combine_posterior_products <- function(...) {
  products <- list(...)
  names_to_bind <- names(products[[1]])
  stats::setNames(lapply(names_to_bind, function(name) {
    do.call(rbind, lapply(products, `[[`, name))
  }), names_to_bind)
}

summarise_effect_draws <- function(data, value_column, grouping_columns) {
  require_columns(
    data,
    c(value_column, grouping_columns),
    "Posterior effect draws"
  )
  key <- interaction(data[grouping_columns], drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(data)), key)
  rows <- lapply(groups, function(index) {
    values <- data[[value_column]][index]
    first <- data[index[[1]], grouping_columns, drop = FALSE]
    first$estimate_median <- stats::median(values)
    first$estimate_mean <- mean(values)
    first$lower_95 <- unname(stats::quantile(values, 0.025))
    first$upper_95 <- unname(stats::quantile(values, 0.975))
    first$Pr_effect_lt_0 <- mean(values < 0)
    first$Pr_effect_gt_0 <- mean(values > 0)
    first$p_two_sided <- 2 * min(
      first$Pr_effect_lt_0,
      first$Pr_effect_gt_0
    )
    first$n_draws <- length(unique(data$.draw[index]))
    first
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

summarise_posterior_products <- function(products) {
  list(
    stage_effects_by_context_metric = summarise_effect_draws(
      products$metric_stage,
      "delta",
      c(
        "analysis", "support", "context", "metric", "metric_label",
        "n_pairs", "pairs_included"
      )
    ),
    stage_effects_combined = summarise_effect_draws(
      products$combined_stage,
      "delta",
      c(
        "analysis", "support", "context", "metric", "metric_label",
        "n_pairs", "pairs_included", "n_metrics"
      )
    ),
    psi_by_metric = summarise_effect_draws(
      products$metric_psi,
      "psi",
      c(
        "analysis", "support", "metric", "metric_label", "n_pairs",
        "pairs_included"
      )
    ),
    psi_combined = summarise_effect_draws(
      products$combined_psi,
      "psi",
      c(
        "analysis", "support", "metric", "metric_label", "n_pairs",
        "pairs_included", "n_metrics"
      )
    )
  )
}

model_data_counts <- function(data, analysis) {
  analysis <- normalise_analysis_name(analysis)
  sets <- list(
    all_pairs = sort(levels(data$pair)),
    common_support_A_B = sort(complete_support_pairs(data))
  )
  do.call(rbind, lapply(names(sets), function(support_name) {
    selected <- data[data$pair %in% sets[[support_name]], , drop = FALSE]
    context_rows <- do.call(rbind, lapply(levels(data$context), function(context) {
      context_data <- selected[selected$context == context, , drop = FALSE]
      before <- context_data[context_data$stage == "before", , drop = FALSE]
      after <- context_data[context_data$stage == "after", , drop = FALSE]
      data.frame(
        analysis = analysis_label(analysis),
        support = support_name,
        context = context,
        N_comparisons = length(unique(as.character(context_data$comparison_id))),
        N_before = length(unique(as.character(before$comparison_id))),
        N_after = length(unique(as.character(after$comparison_id))),
        stringsAsFactors = FALSE
      )
    }))
    all_contexts <- data.frame(
      analysis = analysis_label(analysis),
      support = support_name,
      context = "All contexts",
      N_comparisons = length(unique(as.character(selected$comparison_id))),
      N_before = length(unique(as.character(
        selected$comparison_id[selected$stage == "before"]
      ))),
      N_after = length(unique(as.character(
        selected$comparison_id[selected$stage == "after"]
      ))),
      stringsAsFactors = FALSE
    )
    rbind(context_rows, all_contexts)
  }))
}

add_manuscript_counts <- function(summaries, counts) {
  stage_keys <- c("analysis", "support", "context")
  summaries$stage_effects_by_context_metric <- merge(
    summaries$stage_effects_by_context_metric,
    counts[counts$context != "All contexts", ],
    by = stage_keys,
    all.x = TRUE,
    sort = FALSE
  )
  summaries$stage_effects_combined <- merge(
    summaries$stage_effects_combined,
    counts[counts$context != "All contexts", ],
    by = stage_keys,
    all.x = TRUE,
    sort = FALSE
  )
  psi_counts <- counts[counts$context == "All contexts", ]
  psi_counts$context <- NULL
  summaries$psi_by_metric <- merge(
    summaries$psi_by_metric,
    psi_counts,
    by = c("analysis", "support"),
    all.x = TRUE,
    sort = FALSE
  )
  summaries$psi_combined <- merge(
    summaries$psi_combined,
    psi_counts,
    by = c("analysis", "support"),
    all.x = TRUE,
    sort = FALSE
  )
  summaries
}

write_summary_tables <- function(summaries, table_dir) {
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(table_dir, paste0(names(summaries), ".csv"))
  for (index in seq_along(summaries)) {
    utils::write.csv(summaries[[index]], paths[[index]], row.names = FALSE)
  }
  stats::setNames(paths, names(summaries))
}

figure_3_draw_data <- function(products, support = "all_pairs") {
  metric <- products$metric_stage[
    products$metric_stage$support == support,
    ,
    drop = FALSE
  ]
  combined <- products$combined_stage[
    products$combined_stage$support == support,
    names(metric),
    drop = FALSE
  ]
  plot_data <- rbind(combined, metric)
  expected_cells <- 2L * 2L * 5L
  observed_cells <- unique(plot_data[c("analysis", "context", "metric")])
  if (nrow(observed_cells) != expected_cells) {
    stop(
      "Figure 3 requires two analyses x two context rows x five columns ",
      "(combined plus four metrics)",
      call. = FALSE
    )
  }
  plot_data$context <- factor(
    plot_data$context,
    levels = c("Partner", "Non-partner")
  )
  plot_data$analysis <- factor(
    plot_data$analysis,
    levels = c("Sequence structure", "Call structure")
  )
  plot_data
}

figure_3_analysis_plot <- function(plot_data, analysis) {
  selected <- plot_data[
    as.character(plot_data$analysis) == analysis,
    ,
    drop = FALSE
  ]
  preferred <- c(
    "Combined (equal metric weight)",
    unname(METRIC_LABELS[MODEL_METRICS[[normalise_analysis_name(analysis)]]])
  )
  selected$metric_label <- factor(selected$metric_label, levels = preferred)
  summary <- summarise_effect_draws(
    selected,
    "delta",
    c("metric_label", "context")
  )
  combined_background <- expand.grid(
    metric_label = factor("Combined (equal metric weight)", levels = preferred),
    context = factor(c("Partner", "Non-partner"), levels = levels(selected$context)),
    side = c("Convergence", "Divergence"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  combined_background$xmin <- ifelse(
    combined_background$side == "Convergence",
    -Inf,
    0
  )
  combined_background$xmax <- ifelse(
    combined_background$side == "Convergence",
    0,
    Inf
  )
  combined_background$fill_colour <- ifelse(
    combined_background$side == "Convergence",
    "#e8f3f5",
    "#f8eee8"
  )
  background_labels <- data.frame(
    metric_label = factor(
      rep("Combined (equal metric weight)", 2L),
      levels = preferred
    ),
    context = factor(rep("Partner", 2L), levels = levels(selected$context)),
    x = c(-Inf, Inf),
    y = Inf,
    label = c("Convergence", "Divergence"),
    hjust = c(-0.1, 1.1),
    text_colour = c("#2b6f79", "#9a5637"),
    stringsAsFactors = FALSE
  )
  colours <- c(Partner = "#2c728e", `Non-partner` = "#b8de29")

  ggplot2::ggplot(
    selected,
    ggplot2::aes(x = delta, colour = context, fill = context)
  ) +
    ggplot2::geom_rect(
      data = combined_background,
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = -Inf,
        ymax = Inf,
        fill = fill_colour
      ),
      inherit.aes = FALSE,
      alpha = 0.85,
      colour = NA
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey55",
      linewidth = 0.4
    ) +
    ggplot2::geom_density(alpha = 0.22, linewidth = 0.65) +
    ggplot2::geom_text(
      data = background_labels,
      ggplot2::aes(
        x = x,
        y = y,
        label = label,
        hjust = hjust,
        colour = text_colour
      ),
      inherit.aes = FALSE,
      vjust = 1.2,
      fontface = "bold",
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = summary,
      ggplot2::aes(
        x = lower_95,
        xend = upper_95,
        y = 0,
        yend = 0,
        colour = context
      ),
      inherit.aes = FALSE,
      linewidth = 0.9
    ) +
    ggplot2::geom_point(
      data = summary,
      ggplot2::aes(x = estimate_median, y = 0, colour = context),
      inherit.aes = FALSE,
      size = 1.8
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(context),
      cols = ggplot2::vars(metric_label),
      scales = "free",
      switch = "y"
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        colours,
        `#2b6f79` = "#2b6f79",
        `#9a5637` = "#9a5637"
      ),
      breaks = names(colours),
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = c(colours, `#e8f3f5` = "#e8f3f5", `#f8eee8` = "#f8eee8"),
      breaks = names(colours),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = analysis,
      x = NULL,
      y = NULL,
      colour = NULL,
      fill = NULL
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 9),
      strip.placement = "outside",
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.spacing.x = grid::unit(0.25, "lines"),
      legend.position = "right",
      plot.title = ggplot2::element_text(size = 13, face = "bold")
    )
}

make_stage_effect_figure <- function(products, support = "all_pairs") {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "Packages 'ggplot2' and 'patchwork' are required to draw Figure 3",
      call. = FALSE
    )
  }
  plot_data <- figure_3_draw_data(products, support = support)
  plots <- lapply(
    c("Sequence structure", "Call structure"),
    function(analysis) figure_3_analysis_plot(plot_data, analysis)
  )
  patchwork::wrap_plots(plots, ncol = 1L, guides = "collect") +
    patchwork::plot_annotation(
      caption = "Posterior stage effect (after - before; metric-specific SD)"
    ) &
    ggplot2::theme(legend.position = "right")
}

save_stage_effect_figure <- function(
    products,
    figure_dir,
    basename = "figure_3_regenerated") {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  figure <- make_stage_effect_figure(products)
  png_path <- file.path(figure_dir, paste0(basename, ".png"))
  svg_path <- file.path(figure_dir, paste0(basename, ".svg"))
  ggplot2::ggsave(png_path, figure, width = 15, height = 7.8, dpi = 300)
  ggplot2::ggsave(svg_path, figure, width = 15, height = 7.8)
  c(png = png_path, svg = svg_path)
}

resolve_project_artifact_path <- function(path, project_root) {
  if (!nzchar(path)) return(path)
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) return(path)
  file.path(project_root, path)
}

diagnostic_artifact_paths <- function(diagnostics, project_root) {
  paths <- unlist(diagnostics[c("ppc_path", "trace_path")], use.names = FALSE)
  unique(vapply(
    paths,
    resolve_project_artifact_path,
    character(1),
    project_root = project_root
  ))
}

validate_diagnostics_table <- function(diagnostics, project_root = NULL) {
  required <- c(
    "analysis", "max_rhat", "min_bulk_ess", "min_tail_ess",
    "n_divergent", "n_max_treedepth", "ppc_status", "ppc_path",
    "trace_path", "diagnostic_status"
  )
  require_columns(diagnostics, required, "Model diagnostics cache")
  expected <- c("Sequence structure", "Call structure")
  if (!setequal(diagnostics$analysis, expected)) {
    stop(
      "Model diagnostics cache must contain sequence and call rows",
      call. = FALSE
    )
  }
  numeric <- c(
    "max_rhat", "min_bulk_ess", "min_tail_ess", "n_divergent",
    "n_max_treedepth"
  )
  diagnostics[numeric] <- lapply(diagnostics[numeric], as.numeric)
  if (any(!is.finite(as.matrix(diagnostics[numeric])))) {
    stop("Model diagnostics cache has non-finite values", call. = FALSE)
  }
  complete_ppc <- tolower(trimws(as.character(diagnostics$ppc_status))) %in%
    c("pass", "passed", "ok", "complete", "completed", "generated")
  nonempty_paths <-
    nzchar(trimws(as.character(diagnostics$ppc_path))) &
    nzchar(trimws(as.character(diagnostics$trace_path)))
  passing <-
    diagnostics$max_rhat <= 1.01 &
    diagnostics$min_bulk_ess >= 1000 &
    diagnostics$min_tail_ess >= 1000 &
    diagnostics$n_divergent == 0 &
    diagnostics$n_max_treedepth == 0 &
    complete_ppc &
    nonempty_paths &
    tolower(trimws(as.character(diagnostics$diagnostic_status))) == "pass"
  if (!all(passing)) {
    stop(
      "Model diagnostics do not satisfy the release contract: R-hat <= 1.01, ",
      "bulk/tail ESS >= 1000, zero divergences/tree-depth hits, completed ",
      "posterior-predictive checks, trace plots, and diagnostic_status=pass ",
      "are mandatory",
      call. = FALSE
    )
  }
  if (!is.null(project_root)) {
    paths <- diagnostic_artifact_paths(diagnostics, project_root)
    missing_paths <- paths[!file.exists(paths)]
    if (length(missing_paths)) {
      stop(
        "Model diagnostics cache references missing required plot(s): ",
        paste(missing_paths, collapse = ", "),
        call. = FALSE
      )
    }
  }
  diagnostics
}

read_cached_diagnostics <- function(path, project_root = NULL) {
  if (!file.exists(path)) {
    stop(
      "Missing current model-diagnostics cache: ", path,
      ". Cached/report-only mode cannot reconstruct R-hat, ESS, or sampler ",
      "diagnostics from portable contrast draws. Use the diagnostics produced ",
      "by the current explicit --refit workflow.",
      call. = FALSE
    )
  }
  validate_diagnostics_table(
    utils::read.csv(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    project_root = project_root
  )
}

load_cached_diagnostics <- function(
    cache_path,
    portable_path,
    project_root = NULL) {
  if (file.exists(cache_path)) {
    return(read_cached_diagnostics(cache_path, project_root = project_root))
  }
  if (file.exists(portable_path)) {
    message(
      "Diagnostics cache absent; using portable diagnostics table: ",
      portable_path
    )
    return(read_cached_diagnostics(portable_path, project_root = project_root))
  }
  stop(
    "Missing current model diagnostics. Checked:\n  - ", cache_path,
    "\n  - ", portable_path,
    "\nCached/report-only mode never infers sampler diagnostics from contrast ",
    "draws and never refits silently.",
    call. = FALSE
  )
}

write_report_diagnostics <- function(diagnostics, path) {
  diagnostics <- validate_diagnostics_table(diagnostics)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(diagnostics, path, row.names = FALSE)
  path
}

sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required to write the results manifest")
  }
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

write_results_manifest <- function(files, manifest_path, project_root, mode) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to write the results manifest")
  }
  files <- unique(unname(files[file.exists(files)]))
  root <- normalizePath(project_root, mustWork = TRUE)
  normalized <- normalizePath(files, mustWork = TRUE)
  inside_project <- startsWith(normalized, paste0(root, .Platform$file.sep))
  if (any(!inside_project)) {
    stop(
      "A release model manifest can contain only project-relative outputs; ",
      "outside project_root: ", paste(normalized[!inside_project], collapse = ", "),
      call. = FALSE
    )
  }
  relative <- substring(normalized, nchar(root) + 2L)
  entries <- lapply(seq_along(files), function(index) {
    info <- file.info(files[[index]])
    list(
      path = relative[[index]],
      bytes = unname(info$size),
      sha256 = sha256_file(files[[index]])
    )
  })
  manifest <- list(
    schema_version = 1L,
    generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    mode = mode,
    model_formula = paste(
      deparse(model_formula(), width.cutoff = 500L),
      collapse = " "
    ),
    contrast_definition = paste(
      "Delta is pair-equal after minus before; combined effects are",
      "equal-weight means of four metrics; Psi is Partner Delta minus",
      "Non-partner Delta. Session multi-membership deviations are excluded",
      "from predictions."
    ),
    files = entries
  )
  dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(
    manifest,
    manifest_path,
    pretty = TRUE,
    auto_unbox = TRUE,
    digits = NA
  )
  manifest_path
}
