cran_packages <- c(
  "brms",
  "patchwork",
  "rstan",
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

message("Script 5a can validate inputs and load cached models without sampling.")
message("The --fit and --refit options use brms with backend = 'rstan'.")
message("Those options require a working local Stan/C++ toolchain.")
