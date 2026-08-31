testthat::test_that("release model manifests reject outside-project outputs", {
  testthat::skip_if_not_installed("digest")
  testthat::skip_if_not_installed("jsonlite")
  project <- tempfile("manifest-root-")
  dir.create(project)
  outside <- tempfile("outside-result-")
  writeLines("result", outside)
  manifest <- file.path(project, "results", "model_manifest.json")

  testthat::expect_error(
    write_results_manifest(outside, manifest, project, "cached-report"),
    "project-relative outputs"
  )
})
