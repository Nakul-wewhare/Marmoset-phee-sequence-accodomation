# Install the R packages used by Scripts 4 and 6.

if (identical(getOption("repos")[["CRAN"]], "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

required_packages <- c(
  "tidyverse",
  "brms",
  "tidybayes",
  "bayesplot",
  "posterior",
  "patchwork",
  "reticulate",
  "digest",
  "jsonlite",
  "rstan",
  "rstudioapi"
)

packages_to_install <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(packages_to_install)) {
  install.packages(packages_to_install)
} else {
  message("All required R packages are already installed.")
}
