cran_packages <- c(
  "bayesplot",
  "broom.mixed",
  "brms",
  "digest",
  "emmeans",
  "lme4",
  "lmerTest",
  "loo",
  "jsonlite",
  "patchwork",
  "posterior",
  "rstan",
  "rstudioapi",
  "stringr",
  "testthat",
  "tidybayes",
  "tidyverse"
)

installed <- rownames(installed.packages())
to_install <- setdiff(cran_packages, installed)

if (length(to_install)) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
} else {
  message("All listed R packages are already installed.")
}

message("Note: the final Bayesian script uses brms with backend = 'rstan'.")
message("If model fitting fails, check that your local Stan toolchain is configured correctly.")
