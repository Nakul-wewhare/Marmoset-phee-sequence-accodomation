testthat::test_that("pair and metric averaging are explicitly equal weight", {
  draws <- synthetic_stage_draws("sequence")
  products <- build_posterior_products(draws)

  all_metric <- products$metric_stage[
    products$metric_stage$.draw == 1 &
      products$metric_stage$support == "all_pairs" &
      products$metric_stage$context == "Non-partner" &
      products$metric_stage$metric == SEQUENCE_METRICS[[1]],
  ]
  testthat::expect_equal(all_metric$delta, mean(c(0, 2, 100)))
  testthat::expect_equal(all_metric$n_pairs, 3)

  common_metric <- products$metric_stage[
    products$metric_stage$.draw == 1 &
      products$metric_stage$support == "common_support_A_B" &
      products$metric_stage$context == "Non-partner" &
      products$metric_stage$metric == SEQUENCE_METRICS[[1]],
  ]
  testthat::expect_equal(common_metric$delta, 1)
  testthat::expect_equal(common_metric$n_pairs, 2)
  testthat::expect_identical(
    common_metric$pairs_included,
    "Pair A;Pair B"
  )

  common_combined <- products$combined_stage[
    products$combined_stage$.draw == 1 &
      products$combined_stage$support == "common_support_A_B" &
      products$combined_stage$context == "Non-partner",
  ]
  testthat::expect_equal(common_combined$delta, mean(1 + 0:3))
  testthat::expect_equal(common_combined$n_metrics, 4)
})

testthat::test_that("Psi is Partner Delta minus Non-partner Delta", {
  products <- build_posterior_products(synthetic_stage_draws("call"))
  testthat::expect_true(all(products$metric_psi$psi == 10))
  testthat::expect_true(all(products$combined_psi$psi == 10))
})

testthat::test_that("posterior summaries retain both support analyses", {
  products <- build_posterior_products(synthetic_stage_draws("sequence"))
  summaries <- summarise_posterior_products(products)
  testthat::expect_setequal(
    summaries$psi_combined$support,
    c("all_pairs", "common_support_A_B")
  )
  testthat::expect_true(all(summaries$psi_combined$n_draws == 2L))
  testthat::expect_true(all(c(
    "estimate_median", "lower_95", "upper_95", "Pr_effect_gt_0"
  ) %in% names(summaries$psi_combined)))
})

testthat::test_that("Figure 3 uses separate context rows and five columns", {
  products <- combine_posterior_products(
    build_posterior_products(synthetic_stage_draws("sequence")),
    build_posterior_products(synthetic_stage_draws("call"))
  )
  plot_data <- figure_3_draw_data(products)
  cells <- unique(plot_data[c("analysis", "context", "metric")])
  testthat::expect_equal(nrow(cells), 20L)
  testthat::expect_setequal(
    as.character(unique(plot_data$context)),
    c("Partner", "Non-partner")
  )
  counts <- table(cells$analysis, cells$context)
  testthat::expect_true(all(counts == 5L))
  testthat::expect_false("psi" %in% names(plot_data))
})

testthat::test_that("draw-cache validation rejects inconsistent deltas", {
  draws <- synthetic_stage_draws("sequence")
  draws$delta[[1]] <- draws$delta[[1]] + 1
  testthat::expect_error(
    validate_stage_draws(draws, "sequence"),
    "inconsistent"
  )
})
