testthat::test_that("the model formula matches the manuscript", {
  formula_text <- paste(deparse(model_formula(), width.cutoff = 500L), collapse = " ")
  formula_text <- gsub("[[:space:]]+", " ", formula_text)
  testthat::expect_identical(
    formula_text,
    paste0(
      "z_distance ~ stage * context * metric + ",
      "(1 + stage * context || pair) + ",
      "(1 | mm(session_1_id, session_2_id))"
    )
  )
  pair_text <- gsub(
    "[[:space:]]+",
    " ",
    paste(deparse(pair_prediction_formula()), collapse = " ")
  )
  testthat::expect_identical(
    pair_text,
    "~(1 + stage * context || pair)"
  )
})

testthat::test_that("both analyses use four canonical metrics", {
  testthat::expect_length(SEQUENCE_METRICS, 4L)
  testthat::expect_length(CALL_METRICS, 4L)
  testthat::expect_true("vae" %in% CALL_METRICS)
  testthat::expect_identical(
    canonicalise_metric(c("T_PC_pairwise", "spectro_DTW", "VAE")),
    c("stp", "dtw", "vae")
  )
})

testthat::test_that("model data validation identifies common support", {
  data <- prepare_model_data(
    synthetic_model_data("sequence"),
    "sequence",
    check_standardisation = FALSE
  )
  testthat::expect_identical(
    complete_support_pairs(data),
    c("Pair A", "Pair B")
  )
  testthat::expect_identical(levels(data$stage), c("before", "after"))
  testthat::expect_equal(
    unname(contrasts(data$stage)),
    matrix(c(-0.5, 0.5), nrow = 2L)
  )
  testthat::expect_identical(
    levels(data$context),
    c("Non-partner", "Partner")
  )
  testthat::expect_identical(levels(data$metric), SEQUENCE_METRICS)
  testthat::expect_identical(
    levels(data$session_1_id),
    levels(data$session_2_id)
  )
})

testthat::test_that("public model-table identifier names map to the formula", {
  data <- synthetic_model_data("call")
  data$pair_id <- data$pair
  data$repertoire_a_id <- data$session_1_id
  data$repertoire_b_id <- data$session_2_id
  data$standardized_distance <- data$z_distance
  data$stage_centered <- ifelse(data$stage == "after", 0.5, -0.5)
  data[c("pair", "session_1_id", "session_2_id", "z_distance")] <- NULL
  prepared <- prepare_model_data(
    data,
    "call",
    check_standardisation = FALSE
  )
  testthat::expect_true(all(c(
    "pair", "session_1_id", "session_2_id", "z_distance"
  ) %in% names(prepared)))
})

testthat::test_that("standardized_distance is accepted as the input alias", {
  data <- synthetic_model_data("sequence")
  data$standardized_distance <- data$z_distance
  data$z_distance <- NULL
  prepared <- prepare_model_data(
    data,
    "sequence",
    check_standardisation = FALSE
  )
  testthat::expect_true("z_distance" %in% names(prepared))
  testthat::expect_equal(prepared$z_distance, data$standardized_distance)
})

testthat::test_that("validation rejects incomplete metric grids", {
  data <- synthetic_model_data("call")
  data <- data[-1, ]
  testthat::expect_error(
    prepare_model_data(
      data,
      "call",
      check_standardisation = FALSE
    ),
    "without exactly 4 metrics"
  )
})

testthat::test_that("missing model tables fail with an actionable message", {
  testthat::expect_error(
    read_model_data(tempfile(), "sequence"),
    "Missing sequence model input"
  )
})
