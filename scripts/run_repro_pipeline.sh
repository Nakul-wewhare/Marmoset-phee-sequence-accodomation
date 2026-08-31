#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON:-python}"
MODE="cached"
STRICT=0
ASSEMBLE_CALL_MODEL=0
RECOMPUTE_SEQUENCE_ALIGNMENT=0
REFIT_MODELS=0

usage() {
  cat <<'EOF'
Usage: bash scripts/run_repro_pipeline.sh [options]

Safe default:
  --mode cached       Validate every artifact present; never recompute or fit.
  --strict            Also require every current-analysis release artifact.

Explicit cache-only write:
  --assemble-call-model
                      Build and register both call tables from released caches.

Explicit full-mode options:
  --mode full
  --recompute-sequence-alignment
                      Rebuild local alignment and the sequence model table.
  --refit-models      Run both expensive brms models after inputs validate.

Audio-dependent DTW and VAE rebuild drivers are deliberately separate from
this runner. Use `marmoset-repro preflight-audio --audio-dir PATH` before one.

Environment:
  PYTHON=/path/to/python may select the interpreter. Install the package first
  with `python -m pip install -e .`.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)
      [[ $# -ge 2 ]] || { echo "--mode requires cached or full" >&2; exit 2; }
      MODE="$2"
      shift 2
      ;;
    --mode=*)
      MODE="${1#*=}"
      shift
      ;;
    --strict)
      STRICT=1
      shift
      ;;
    --assemble-call-model)
      ASSEMBLE_CALL_MODEL=1
      shift
      ;;
    --recompute-sequence-alignment)
      RECOMPUTE_SEQUENCE_ALIGNMENT=1
      shift
      ;;
    --refit-models)
      REFIT_MODELS=1
      shift
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ "$MODE" != "cached" && "$MODE" != "full" ]]; then
  echo "--mode must be cached or full" >&2
  exit 2
fi
if [[ "$MODE" == "cached" && ( $RECOMPUTE_SEQUENCE_ALIGNMENT -eq 1 || $REFIT_MODELS -eq 1 ) ]]; then
  echo "Recompute/refit flags require --mode full" >&2
  exit 2
fi

cd "$ROOT_DIR"
if ! "$PYTHON_BIN" -c 'import marmoset_convergence' 2>/dev/null; then
  echo "The local package is not installed for $PYTHON_BIN." >&2
  echo "Run: $PYTHON_BIN -m pip install -e ." >&2
  exit 2
fi

if [[ $ASSEMBLE_CALL_MODEL -eq 1 ]]; then
  echo "Explicit cache-only component: assembling the four-metric call model table"
  "$PYTHON_BIN" scripts/assemble_call_model_table.py --root "$ROOT_DIR"
fi

if [[ "$MODE" == "full" ]]; then
  if [[ $RECOMPUTE_SEQUENCE_ALIGNMENT -eq 1 ]]; then
    echo "Explicit component: rebuilding the complete sequence local-alignment cache"
    "$PYTHON_BIN" scripts/prepare_cached_artifacts.py \
      --root "$ROOT_DIR" --recompute-sequence-alignment
  fi

  if [[ $REFIT_MODELS -eq 1 ]]; then
    echo "Explicit component: fitting two 4-chain brms models"
    Rscript R/run_models.R --refit
  fi

  if [[ $RECOMPUTE_SEQUENCE_ALIGNMENT -eq 0 && $REFIT_MODELS -eq 0 ]]; then
    echo "Full mode selected, but no component flag was supplied; validating only."
  fi
fi

VALIDATE_ARGS=(--root "$ROOT_DIR")
if [[ $STRICT -eq 1 ]]; then
  VALIDATE_ARGS+=(--strict)
fi
"$PYTHON_BIN" scripts/validate_release.py "${VALIDATE_ARGS[@]}"

if [[ -f data/derived/model_sequence.csv && -f data/derived/model_call.csv \
      && -f results/posterior_draws/sequence_stage_effect_draws.csv.gz \
      && -f results/posterior_draws/call_stage_effect_draws.csv.gz ]]; then
  echo "Validating installed posterior caches without writing reports or sampling"
  Rscript R/run_models.R --validate-cache-only
else
  echo "Posterior-cache validation skipped: the manifest identifies the missing release assets."
fi

echo "Validation complete. No unrequested costly step was run."
