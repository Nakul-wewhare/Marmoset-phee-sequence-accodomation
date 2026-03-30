# ============================================================
# FINAL BAYESIAN ANALYSIS SCRIPT (brms) — manuscript-ready outputs
# Outcome: z_distance (SD units)
# Models: NO random slopes anywhere
#
# Outputs (for manuscript):
#   Tables:
#     - pooled_results.csv  (β median/mean, 95% CrI, Pr(β<0), two-sided p, N)
#     - per_metric_results.csv (same; with manuscript metric names + metric_code)
#     - diagnostics_summary_ALL.csv (Rhat, ESS, divergences, Pareto-k, LOO elpd)
#   Figures (aesthetic consistent with your Bayesian plots):
#     - Fig_main_pooled_posterior.png/svg
#     - Fig_grid_sequence_metrics_fullposterior.png/svg
#     - Fig_grid_spectral_metrics_fullposterior.png/svg
#   Optional auxiliary outputs when requested with environment variables:
#     - per-model diagnostics (pp_check, trace plots, loo files)
#     - full fitted model objects
#
# Notes:
#   - stage is coded as centered dummy: before=-0.5, after=+0.5
#     => beta == (after - before) on z_distance (SD units).
#   - Metric labels in plots/tables use manuscript names, while metric_code retains CSV codes.
# ============================================================

# 0) Working directory (supports RStudio and command-line Rscript execution)
script_path <- NULL

if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::hasFun("getActiveDocumentContext")) {
  script_path <- tryCatch(
    rstudioapi::getActiveDocumentContext()$path,
    error = function(e) NULL
  )
}

if (is.null(script_path) || !nzchar(script_path)) {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1) {
    script_path <- sub("^--file=", "", file_arg)
  }
}

if (!is.null(script_path) && nzchar(script_path)) {
  setwd(dirname(normalizePath(script_path)))
  message("Working directory set to script folder: ", getwd())
}

library(tidyverse)
library(brms)
library(tidybayes)
library(bayesplot)
library(posterior)

# -------------------------
# User settings
# -------------------------
root    <- normalizePath("..")
csv_dir <- file.path(root, "glm data")

distance_col <- "z_distance"

# Metric codes (from your CSVs) — ordered to match manuscript flow
SEQ_METRICS  <- c("transition_matrix", "bigram", "repeat_A_len", "local_alignment")
SPEC_METRICS <- c("T_PC_pairwise", "M_PC_pairwise", "spectro_DTW")

# Manuscript-friendly metric labels (codes -> names)
metric_label_map <- c(
  # sequence metrics
  "transition_matrix" = "Transition probabilities",
  "bigram"            = "Bigram distribution",
  "repeat_A_len"      = "Repeat distribution of phee",
  "local_alignment"   = "Local alignment",
  
  # spectral metrics
  "T_PC_pairwise"     = "Spectro-temporal parameters (STP)",
  "M_PC_pairwise"     = "Mel-frequency cepstral coefficients (MFCC)",
  "spectro_DTW"       = "Dynamic time warping (DTW)"
)

# Colours & shapes (match your established aesthetic)
PARTNER_COL    <- "#2c728e"
NONPARTNER_COL <- "#b8de29"

# Output dirs
out_dir    <- file.path(csv_dir, "brms_manuscript_outputs-newsession")
fig_dir    <- file.path(out_dir, "figures")
diag_dir   <- file.path(out_dir, "diagnostics")
model_dir  <- file.path(out_dir, "models")
table_dir  <- file.path(out_dir, "tables")
save_auxiliary_outputs <- identical(Sys.getenv("MANUSCRIPT_REPO_SAVE_AUXILIARY"), "1")
save_model_objects <- identical(Sys.getenv("MANUSCRIPT_REPO_SAVE_MODELS"), "1")

dir.create(out_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

if (save_auxiliary_outputs) {
  dir.create(diag_dir, showWarnings = FALSE, recursive = TRUE)
}

if (save_model_objects) {
  dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)
}

# -------------------------
# Data loader
# -------------------------
prep_df <- function(fname) {
  df <- read_csv(file.path(csv_dir, fname), show_col_types = FALSE)
  if (!distance_col %in% names(df)) stop(sprintf("Column '%s' missing in %s", distance_col, fname))
  
  df %>%
    mutate(
      stage = factor(stage, levels = c("before", "after")),
      stage_after = if_else(stage == "after", 0.5, -0.5),
      
      metric = factor(metric),
      pair = factor(pair),
      
      # --- extract focal/partner IDs from "pair" assuming focal-first ---
      focal_id   = sub("-.*$", "", as.character(pair)),
      partner_id = sub("^.*-", "", as.character(pair)),
      
      # --- make session IDs unique by individual + stage + session number ---
      session_focal_id   = paste(focal_id,   stage, session_focal,   sep = "_"),
      session_partner_id = paste(partner_id, stage, session_partner, sep = "_"),
      
      session_focal_id   = factor(session_focal_id),
      session_partner_id = factor(session_partner_id)
    ) %>%
    rename(y = all_of(distance_col))
}

# Add metric_code (CSV code) and metric_label (manuscript label)
add_metric_labels <- function(df) {
  df %>%
    mutate(
      metric_code  = as.character(metric),
      metric_label = dplyr::recode(metric_code, !!!metric_label_map)
    )
}

cat("\nLoading data...\n")
spec_partner    <- prep_df("glm_data_spec_paired.csv")      %>% mutate(dataset_tag = "Partner\nCall structure")       %>% add_metric_labels()
spec_nonpartner <- prep_df("glm_data_spec_non_paired.csv")  %>% mutate(dataset_tag = "Non-partner\nCall structure")   %>% add_metric_labels()
seq_partner     <- prep_df("glm_data_seq_partner.csv")      %>% mutate(dataset_tag = "Partner\nSequence structure")   %>% add_metric_labels()
seq_nonpartner  <- prep_df("glm_data_seq_non_partner.csv")  %>% mutate(dataset_tag = "Non-partner\nSequence structure") %>% add_metric_labels()

all_dat <- bind_rows(spec_partner, spec_nonpartner, seq_partner, seq_nonpartner) %>%
  mutate(
    dataset_tag = factor(
      dataset_tag,
      levels = c(
        "Non-partner\nSequence structure",
        "Non-partner\nCall structure",
        "Partner\nSequence structure",
        "Partner\nCall structure"
      )
    )
  )

cat("\nMetric codes present:\n")
print(sort(unique(all_dat$metric_code)))
cat("\nMetric labels present:\n")
print(sort(unique(all_dat$metric_label)))

if (identical(Sys.getenv("MANUSCRIPT_REPO_VALIDATE_ONLY"), "1")) {
  message("Validation mode: data loaded successfully; exiting before model fitting.")
  quit(save = "no", status = 0)
}

# -------------------------
# brms settings
# -------------------------
brms_ctrl <- list(adapt_delta = 0.995, max_treedepth = 15)

priors_gaussian <- c(
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)

# -------------------------
# Utility: posterior summaries for beta
# -------------------------
summarise_beta <- function(beta_draws) {
  p_lt0 <- mean(beta_draws < 0)
  p_2s  <- 2 * min(p_lt0, 1 - p_lt0)  # two-sided posterior tail area
  tibble(
    beta_median = median(beta_draws),
    beta_mean   = mean(beta_draws),
    beta_lo95   = as.numeric(quantile(beta_draws, 0.025)),
    beta_hi95   = as.numeric(quantile(beta_draws, 0.975)),
    Pr_beta_lt0 = p_lt0,
    p_two_sided = p_2s
  )
}

# -------------------------
# Utility: save standard plots per model
# -------------------------
save_model_plots <- function(fit, model_name) {
  if (!save_auxiliary_outputs) {
    return(invisible(FALSE))
  }
  
  png(file.path(diag_dir, paste0(model_name, "_pp_check.png")),
      width = 1400, height = 1000, res = 150)
  print(pp_check(fit, ndraws = 200))
  dev.off()
  
  png(file.path(diag_dir, paste0(model_name, "_trace_stage_sigma.png")),
      width = 1600, height = 1000, res = 150)
  bayesplot::mcmc_trace(as_draws_array(fit), pars = c("b_stage_after", "sigma"))
  dev.off()
  
  invisible(TRUE)
}

# -------------------------
# Utility: diagnostics summary per brms fit (Rhat/ESS/divergences/LOO)
# -------------------------
extract_diagnostics <- function(fit, model_name) {
  
  draws_sum <- posterior::summarise_draws(as_draws_array(fit))
  
  max_rhat <- suppressWarnings(max(draws_sum$rhat, na.rm = TRUE))
  min_ess_bulk <- suppressWarnings(min(draws_sum$ess_bulk, na.rm = TRUE))
  min_ess_tail <- suppressWarnings(min(draws_sum$ess_tail, na.rm = TRUE))
  
  # divergences
  sp <- fit$fit@sim$sampler_diagnostics
  divs <- 0
  if (!is.null(sp)) {
    divs <- sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1)))
  }
  
  # LOO (moment match only if needed)
  loo_res <- loo(fit)
  k <- loo_res$diagnostics$pareto_k
  max_k <- suppressWarnings(max(k, na.rm = TRUE))
  
  used_mm <- FALSE
  if (any(is.finite(k) & k > 0.7)) {
    used_mm <- TRUE
    loo_res <- loo(fit, moment_match = TRUE)
    k <- loo_res$diagnostics$pareto_k
    max_k <- suppressWarnings(max(k, na.rm = TRUE))
  }
  
  if (save_auxiliary_outputs) {
    saveRDS(loo_res, file.path(diag_dir, paste0(model_name, "_loo.rds")))
    write_csv(as.data.frame(loo_res$pointwise),
              file.path(diag_dir, paste0(model_name, "_loo_pointwise.csv")))
  }
  
  tibble(
    model = model_name,
    max_rhat = max_rhat,
    min_ess_bulk = min_ess_bulk,
    min_ess_tail = min_ess_tail,
    divergences = divs,
    max_pareto_k = max_k,
    loo_elpd = loo_res$estimates["elpd_loo", "Estimate"],
    loo_se   = loo_res$estimates["elpd_loo", "SE"],
    loo_used_moment_match = used_mm
  )
}

# ============================================================
# 1) POOLED MODELS (4 datasets)
# y ~ stage + (1|metric) + (1|pair) + (1|session_focal) + (1|session_partner)
# ============================================================
fit_pooled <- function(data, model_name) {
  fml <- y ~ stage_after +
    (1 | metric) +
    (1 | pair) +
    (1 | session_focal_id) +
    (1 | session_partner_id)
  
  cat("\n--- Fitting pooled:", model_name, "---\n")
  fit <- brm(
    formula = fml,
    data = data,
    family = gaussian(),
    prior = priors_gaussian,
    chains = 4, cores = 4,
    iter = 4000, warmup = 2000,
    seed = 123,
    control = brms_ctrl,
    backend = "rstan",
    save_pars = save_pars(all = TRUE)
  )
  
  if (save_model_objects) {
    saveRDS(fit, file.path(model_dir, paste0(model_name, ".rds")))
  }
  save_model_plots(fit, model_name)
  fit
}

fit_seq_nonpar   <- fit_pooled(seq_nonpartner,  "pooled_seq_nonpartner")
fit_seq_partner  <- fit_pooled(seq_partner,     "pooled_seq_partner")
fit_spec_nonpar  <- fit_pooled(spec_nonpartner, "pooled_spec_nonpartner")
fit_spec_partner <- fit_pooled(spec_partner,    "pooled_spec_partner")

# Pooled results table (manuscript-ready)
pooled_results <- tibble(
  dataset_tag = c("Non-partner\nSequence structure",
                  "Partner\nSequence structure",
                  "Non-partner\nCall structure",
                  "Partner\nCall structure"),
  model = c("pooled_seq_nonpartner",
            "pooled_seq_partner",
            "pooled_spec_nonpartner",
            "pooled_spec_partner"),
  N = c(nrow(seq_nonpartner), nrow(seq_partner), nrow(spec_nonpartner), nrow(spec_partner))
) %>%
  mutate(
    pairing  = if_else(grepl("^Partner", dataset_tag), "Partner", "Non-partner"),
    modality = if_else(grepl("Sequence", dataset_tag), "Sequence structure", "Call structure")
  )

pooled_betas <- list(
  pooled_seq_nonpartner = as_draws_df(fit_seq_nonpar)$b_stage_after,
  pooled_seq_partner    = as_draws_df(fit_seq_partner)$b_stage_after,
  pooled_spec_nonpartner= as_draws_df(fit_spec_nonpar)$b_stage_after,
  pooled_spec_partner   = as_draws_df(fit_spec_partner)$b_stage_after
)

pooled_summ <- imap_dfr(pooled_betas, ~ summarise_beta(.x) %>% mutate(model = .y))

pooled_results <- pooled_results %>% left_join(pooled_summ, by = "model")

write_csv(pooled_results, file.path(table_dir, "pooled_results.csv"))

diag_pooled <- bind_rows(
  extract_diagnostics(fit_seq_nonpar,   "pooled_seq_nonpartner"),
  extract_diagnostics(fit_seq_partner,  "pooled_seq_partner"),
  extract_diagnostics(fit_spec_nonpar,  "pooled_spec_nonpartner"),
  extract_diagnostics(fit_spec_partner, "pooled_spec_partner")
)
write_csv(diag_pooled, file.path(table_dir, "diagnostics_summary_pooled.csv"))

# ============================================================
# 2) FIGURE: pooled posterior (same aesthetic)
# ============================================================
get_stage_draws <- function(fit, pairing, modality) {
  as_draws_df(fit) %>%
    transmute(
      pairing  = pairing,
      modality = modality,
      b_stage  = b_stage_after,
      ylab = paste0(pairing, "\n", modality)
    )
}

post_stage_pooled <- bind_rows(
  get_stage_draws(fit_seq_partner,  "Partner",     "Sequence structure"),
  get_stage_draws(fit_seq_nonpar,   "Non-partner", "Sequence structure"),
  get_stage_draws(fit_spec_partner, "Partner",     "Call structure"),
  get_stage_draws(fit_spec_nonpar,  "Non-partner", "Call structure")
) %>%
  mutate(
    ylab = factor(
      ylab,
      levels = c(
        "Partner\nSequence structure",
        "Non-partner\nSequence structure",
        "Partner\nCall structure",
        "Non-partner\nCall structure"
      )
    ),
    pairing  = factor(pairing,  levels = c("Partner","Non-partner")),
    modality = factor(modality, levels = c("Sequence structure","Call structure"))
  )

post_summ_pooled <- post_stage_pooled %>%
  group_by(pairing, modality, ylab) %>%
  summarise(
    estimate = median(b_stage),
    lo = quantile(b_stage, 0.025),
    hi = quantile(b_stage, 0.975),
    .groups = "drop"
  )

fig_main <- ggplot(post_stage_pooled, aes(x = b_stage, y = ylab)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  stat_halfeye(
    aes(fill = pairing),
    slab_alpha = 0.55,
    .width = 0,
    point_interval = NULL
  ) +
  geom_segment(
    data = post_summ_pooled,
    aes(x = lo, xend = hi, y = ylab, yend = ylab, colour = pairing),
    linewidth = 0.9,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = post_summ_pooled,
    aes(x = estimate, y = ylab, colour = pairing, shape = modality),
    size = 2.2,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c("Partner" = PARTNER_COL,
                               "Non-partner" = NONPARTNER_COL)) +
  scale_color_manual(values = c("Partner" = PARTNER_COL,
                                "Non-partner" = NONPARTNER_COL)) +
  scale_shape_manual(values = c("Sequence structure" = 16,
                                "Call structure"     = 17)) +
  labs(
    x = "Posterior for stage coefficient (β)",
    y = "Fitted models",
    fill = NULL, colour = NULL, shape = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_blank(),
    axis.text.y = element_text(hjust = 0.5)
  )
# Desired top-to-bottom order
ylab_order <- c(
  "Partner\nSequence structure",
  "Non-partner\nSequence structure",
  "Partner\nCall structure",
  "Non-partner\nCall structure"
)

post_stage_pooled <- post_stage_pooled %>%
  mutate(ylab = factor(ylab, levels = ylab_order))

# IMPORTANT: reverse so first level appears at TOP in ggplot
fig_main <- fig_main + scale_y_discrete(limits = rev(ylab_order))
print(fig_main)

ggsave(file.path(fig_dir, "Fig_main_pooled_posterior.png"),
       fig_main, dpi = 300, width = 8, height = 4.5)
ggsave(file.path(fig_dir, "Fig_main_pooled_posterior.svg"),
       fig_main, width = 8, height = 4.5)

# ============================================================
# 3) PER-METRIC MODELS (robustness)
# y ~ stage + (1|pair) + (1|session_focal) + (1|session_partner)
# ============================================================
fit_per_metric <- function(dat_sub, model_name) {
  fml <- y ~ stage_after +
    (1 | pair) +
    (1 | session_focal_id) +
    (1 | session_partner_id)
  
  cat("\n--- Fitting per-metric:", model_name, "---\n")
  fit <- brm(
    formula = fml,
    data = dat_sub,
    family = gaussian(),
    prior = priors_gaussian,
    chains = 4, cores = 4,
    iter = 4000, warmup = 2000,
    seed = 123,
    control = brms_ctrl,
    backend = "rstan",
    save_pars = save_pars(all = TRUE)
  )
  
  if (save_model_objects) {
    saveRDS(fit, file.path(model_dir, paste0(model_name, ".rds")))
  }
  save_model_plots(fit, model_name)
  fit
}

run_per_metric_set <- function(dat, dataset_tag, keep_metric_codes) {
  
  dat2 <- dat %>%
    filter(metric_code %in% keep_metric_codes) %>%
    droplevels()
  
  out_rows  <- list()
  out_draws <- list()
  out_diag  <- list()
  
  for (m in keep_metric_codes) {
    dat_m <- dat2 %>% filter(metric_code == m)
    if (nrow(dat_m) == 0) next
    
    m_label <- metric_label_map[[m]]
    
    model_name <- paste0(
      "metric_", gsub("\\s+","_", tolower(gsub("[^A-Za-z0-9_]", "_", dataset_tag))),
      "__", m
    )
    
    fit_m <- fit_per_metric(dat_m, model_name)
    
    beta_draws <- as_draws_df(fit_m)$b_stage_after
    
    summ <- summarise_beta(beta_draws) %>%
      mutate(
        dataset_tag  = dataset_tag,
        metric_code  = m,
        metric       = m_label,
        model        = model_name,
        N            = nrow(dat_m),
        pairing      = if_else(grepl("^Partner", dataset_tag), "Partner", "Non-partner"),
        modality     = if_else(grepl("Sequence", dataset_tag), "Sequence structure", "Call structure")
      )
    
    out_rows[[m]] <- summ
    
    out_draws[[m]] <- tibble(
      dataset_tag  = dataset_tag,
      metric_code  = m,
      metric_label = m_label,
      pairing      = if_else(grepl("^Partner", dataset_tag), "Partner", "Non-partner"),
      modality     = if_else(grepl("Sequence", dataset_tag), "Sequence structure", "Call structure"),
      b_stage      = beta_draws
    )
    
    out_diag[[m]] <- extract_diagnostics(fit_m, model_name)
  }
  
  list(
    results = bind_rows(out_rows),
    draws   = bind_rows(out_draws),
    diag    = bind_rows(out_diag)
  )
}

pm_seq_nonpar   <- run_per_metric_set(seq_nonpartner,  "Non-partner\nSequence structure", SEQ_METRICS)
pm_seq_partner  <- run_per_metric_set(seq_partner,     "Partner\nSequence structure",     SEQ_METRICS)
pm_spec_nonpar  <- run_per_metric_set(spec_nonpartner, "Non-partner\nCall structure",     SPEC_METRICS)
pm_spec_partner <- run_per_metric_set(spec_partner,    "Partner\nCall structure",         SPEC_METRICS)

per_metric_results <- bind_rows(
  pm_seq_nonpar$results,
  pm_seq_partner$results,
  pm_spec_nonpar$results,
  pm_spec_partner$results
) %>%
  mutate(
    pairing  = factor(pairing, levels = c("Partner","Non-partner")),
    modality = factor(modality, levels = c("Sequence structure","Call structure")),
    metric   = factor(metric, levels = unname(metric_label_map[c(SEQ_METRICS, SPEC_METRICS)])),
    metric_code = factor(metric_code, levels = c(SEQ_METRICS, SPEC_METRICS))
  )

write_csv(per_metric_results, file.path(table_dir, "per_metric_results.csv"))

per_metric_diag <- bind_rows(
  pm_seq_nonpar$diag,
  pm_seq_partner$diag,
  pm_spec_nonpar$diag,
  pm_spec_partner$diag
)
write_csv(per_metric_diag, file.path(table_dir, "diagnostics_summary_per_metric.csv"))

per_metric_draws <- bind_rows(
  pm_seq_nonpar$draws,
  pm_seq_partner$draws,
  pm_spec_nonpar$draws,
  pm_spec_partner$draws
) %>%
  mutate(
    pairing  = factor(pairing, levels = c("Partner","Non-partner")),
    modality = factor(modality, levels = c("Sequence structure","Call structure")),
    metric_label = factor(metric_label, levels = unname(metric_label_map[c(SEQ_METRICS, SPEC_METRICS)])),
    metric_code = factor(metric_code, levels = c(SEQ_METRICS, SPEC_METRICS))
  )

# ============================================================
# 4) GRID FIGURES (full posterior per cell; manuscript metric names)
# ============================================================
library(stringr)

wrap_lines <- function(x, width = 16, max_lines = 3) {
  w <- str_wrap(x, width = width)
  parts <- strsplit(w, "\n")
  vapply(parts, function(p) {
    if (length(p) <= max_lines) return(paste(p, collapse = "\n"))
    paste(c(p[1:(max_lines-1)], paste(p[max_lines:length(p)], collapse = " ")), collapse = "\n")
  }, character(1))
}

make_metric_grid_bayes_fullposterior_v2 <- function(draws, keep_metric_codes, row_levels,
                                                    wrap_width = 16, wrap_lines_n = 3) {
  labels <- unname(metric_label_map[keep_metric_codes])
  labels_wrapped <- wrap_lines(labels, width = wrap_width, max_lines = wrap_lines_n)
  
  dd <- draws %>%
    filter(as.character(metric_code) %in% keep_metric_codes) %>%
    mutate(
      dataset_tag = factor(as.character(dataset_tag), levels = row_levels),
      metric_label = factor(as.character(metric_label), levels = labels),
      metric_label_wrapped = factor(
        wrap_lines(as.character(metric_label), width = wrap_width, max_lines = wrap_lines_n),
        levels = labels_wrapped
      )
    )
  
  # One summary per panel (NO pairing/modality in grouping => no duplicates)
  summ <- dd %>%
    group_by(dataset_tag, metric_label_wrapped) %>%
    summarise(
      pairing  = if_else(grepl("^Partner", first(dataset_tag)), "Partner", "Non-partner"),
      modality = if_else(grepl("Sequence", first(dataset_tag)), "Sequence structure", "Call structure"),
      estimate = median(b_stage),
      lo = quantile(b_stage, 0.025),
      hi = quantile(b_stage, 0.975),
      .groups = "drop"
    ) %>%
    mutate(
      pairing  = factor(pairing, levels = c("Partner","Non-partner")),
      modality = factor(modality, levels = c("Sequence structure","Call structure"))
    )
  
  p <- ggplot(dd, aes(x = b_stage, y = 0)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    stat_halfeye(
      aes(fill = if_else(grepl("^Partner", dataset_tag), "Partner", "Non-partner")),
      slab_alpha = 0.6,
      .width = 0,
      point_interval = NULL
    ) +
    geom_segment(
      data = summ,
      aes(x = lo, xend = hi, y = 0, yend = 0,
          colour = pairing),
      linewidth = 0.9,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = summ,
      aes(x = estimate, y = 0, colour = pairing, shape = modality),
      size = 2.0,
      inherit.aes = FALSE
    ) +
    facet_grid(dataset_tag ~ metric_label_wrapped, switch = "y") +
    scale_y_continuous(NULL, breaks = NULL) +
    scale_fill_manual(values = c("Partner" = PARTNER_COL,
                                 "Non-partner" = NONPARTNER_COL)) +
    scale_color_manual(values = c("Partner" = PARTNER_COL,
                                  "Non-partner" = NONPARTNER_COL)) +
    scale_shape_manual(values = c("Sequence structure" = 16,
                                  "Call structure"     = 17)) +
    labs(
      x = "Posterior for stage coefficient (β)",
      y = NULL,
      fill = NULL, colour = NULL, shape = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      # REMOVE the black boxes around facet strips
      strip.background = element_blank(),
      strip.placement  = "outside",
      strip.text       = element_text(size = 10, lineheight = 0.95),
      panel.spacing    = unit(1.0, "lines")
      
    )
  theme(
    strip.text.y.left = element_blank(),
    strip.text.y.right = element_blank(),
    strip.background.y = element_blank()
  )
  
  p
  p <- p + theme(
    strip.text.y.left = element_blank(),
    strip.text.y.right = element_blank()
  )
}
seq_rows <- c("Partner\nSequence structure", "Non-partner\nSequence structure")
fig_seq_grid <- make_metric_grid_bayes_fullposterior_v2(per_metric_draws, SEQ_METRICS, seq_rows,
                                                        wrap_width = 16, wrap_lines_n = 2)

spec_rows <- c("Partner\nCall structure", "Non-partner\nCall structure")
fig_spec_grid <- make_metric_grid_bayes_fullposterior_v2(per_metric_draws, SPEC_METRICS, spec_rows,
                                                         wrap_width = 14, wrap_lines_n = 3)
ggsave(file.path(fig_dir, "Fig_grid_sequence_metrics_fullposterior.png"),
       fig_seq_grid, width = 12, height = 4.5, dpi = 300)
ggsave(file.path(fig_dir, "Fig_grid_spectral_metrics_fullposterior.png"),
       fig_spec_grid, width = 12, height = 4.5, dpi = 300)
# ============================================================
# 5) One combined diagnostics table (pooled + per-metric)
# ============================================================
diagnostics_all <- bind_rows(
  diag_pooled %>% mutate(model_set = "pooled"),
  per_metric_diag %>% mutate(model_set = "per_metric")
)

write_csv(diagnostics_all, file.path(table_dir, "diagnostics_summary_ALL.csv"))

cat("\n============================================================\n")
cat("DONE. Manuscript outputs written to:\n", out_dir, "\n")
cat("Tables:\n - ", file.path(table_dir, "pooled_results.csv"), "\n",
    " - ", file.path(table_dir, "per_metric_results.csv"), "\n",
    " - ", file.path(table_dir, "diagnostics_summary_ALL.csv"), "\n", sep="")
cat("Figures:\n - ", file.path(fig_dir, "Fig_main_pooled_posterior.png"), "\n",
    " - ", file.path(fig_dir, "Fig_grid_sequence_metrics_fullposterior.png"), "\n",
    " - ", file.path(fig_dir, "Fig_grid_spectral_metrics_fullposterior.png"), "\n", sep="")
if (save_auxiliary_outputs) {
  cat("Diagnostics folder:\n - ", diag_dir, "\n", sep="")
}
if (save_model_objects) {
  cat("Model objects folder:\n - ", model_dir, "\n", sep="")
}
cat("============================================================\n")



library(patchwork)

# 1) remove legends from each plot (collect later)
fig_main_noleg <- fig_main + theme(legend.position = "none")
fig_seq_noleg  <- fig_seq_grid + theme(legend.position = "none")
fig_spec_noleg <- fig_spec_grid + theme(legend.position = "none")

# 2) tighten margins so rows align better
tight_margins <- theme(plot.margin = margin(6, 6, 6, 6))
fig_main_noleg <- fig_main_noleg + tight_margins
fig_seq_noleg  <- fig_seq_noleg  + tight_margins
fig_spec_noleg <- fig_spec_noleg + tight_margins

# 3) right column with equal total height (2 + 2 = 4 units)
right_col <- fig_seq_noleg / fig_spec_noleg +
  plot_layout(heights = c(2, 2))

# 4) combine with left column narrower
panel <- (fig_main_noleg | right_col) +
  plot_layout(widths = c(0.4, 1), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

# 5) bigger text everywhere
big_text_theme <- theme(
  text = element_text(size = 16),
  axis.title = element_text(size = 18),
  axis.text  = element_text(size = 14),
  strip.text = element_text(size = 14),
  legend.text = element_text(size = 14)
)

panel <- panel & big_text_theme

ggsave(file.path(fig_dir, "Fig_combined_panel.png"),
       panel, width = 16, height = 8.5, dpi = 300)
ggsave(file.path(fig_dir, "Fig_combined_panel.svg"),
       panel, width = 16, height = 8.5)

print(panel)
