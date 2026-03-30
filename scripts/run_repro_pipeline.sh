#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR/code"

run_notebook() {
  local notebook="$1"
  local tmp_dir
  tmp_dir="$(mktemp -d)"
  jupyter nbconvert --to notebook --execute --output-dir "$tmp_dir" "$notebook" >/dev/null
  rm -rf "$tmp_dir"
}

echo "[1/4] Recomputing spectral distances from bundled processed inputs"
run_notebook "script_2_spectral_measurements.ipynb"

echo "[2/4] Recomputing sequence distances"
run_notebook "script_3_seq_measurements.ipynb"

echo "[3/4] Rebuilding model-ready analysis tables"
run_notebook "script_4_GLM preprocessing.ipynb"

echo "[4/4] Rerunning Bayesian models"
Rscript "script_5_R_Bayesian_models.R"

echo "Done. Final outputs are available in:"
echo "$ROOT_DIR/glm data/brms_manuscript_outputs-newsession/"
