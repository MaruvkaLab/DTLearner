#!/usr/bin/env bash
################################################################################
# run_all.sh
# Runs every analysis script in the order their inputs require.
# Shell equivalent of run_all.R.
#
# Put this file in the same folder as the other scripts (WORK_DIR, next to
# input_files/) and run it:
#
#     bash run_all.sh
#   or, after  chmod +x run_all.sh :
#     ./run_all.sh
#
# Options:
#     --keep-going     carry on after a script fails and report every failure
#                      at the end (default: stop at the first failure)
#     --from SCRIPT    start at SCRIPT (e.g. --from plot_figure_3.R), skipping
#                      the ones before it
#
# ORDER
#   1. plot_figure_1.R           writes figure_1/sup_tab2.xlsx, which
#                                plot_figure_2.R reads for panel 2f
#   2. plot_figure_2.R           modelling + Figure 2
#   3. predict_hcc.R             -> hcc_pred/         (Figure 3)
#   4. predict_lgg.R             -> lgg_pred/         (Figure 3)
#   5. predict_hcc_LOOCV.R       -> loocv_hcc/        (Supp. Figure 3)
#   6. predict_tcga.R            -> tcga_pred/        (Figure 4)
#   7. predict_drug_cohort.R     -> drug_cohort_pred/ (Figure 5)
#   8. predict_ccle.R            -> ccle_pred/        (Figure 5)
#   9. plot_figure_3.R
#  10. plot_figure_4.R
#  11. plot_figure_5.R
#
# plot_figure_2_from_cache.R is not run: it only redraws Figure 2 from the
# files plot_figure_2.R has already written.
#
# Each script runs in its own R session (Rscript), and everything it prints
# goes to logs/<script>.log. Requires Rscript on the PATH (set RSCRIPT to its
# full path otherwise, e.g.  RSCRIPT=/usr/local/bin/Rscript bash run_all.sh).
################################################################################

set -u

# ── Settings ──────────────────────────────────────────────────────────────────

# Scripts to run, in order. Comment out a line to skip that script.
SCRIPTS=(
  "plot_figure_1.R"
  "plot_figure_2.R"
  "predict_hcc.R"
  "predict_lgg.R"
  "predict_hcc_LOOCV.R"
  "predict_tcga.R"
  "predict_drug_cohort.R"
  "predict_ccle.R"
  "plot_figure_3.R"
  "plot_figure_4.R"
  "plot_figure_5.R"
)

STOP_ON_ERROR=1
START_FROM=""
RSCRIPT="${RSCRIPT:-Rscript}"

# ── Command-line options ──────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --keep-going) STOP_ON_ERROR=0; shift ;;
    --from)       START_FROM="${2:-}"; shift 2 ;;
    -h|--help)    sed -n '2,20p' "$0"; exit 0 ;;
    *)            echo "Unknown option: $1 (see --help)" >&2; exit 2 ;;
  esac
done

# ── WORK_DIR = the folder this file is in ─────────────────────────────────────
WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "WORK_DIR: $WORK_DIR"

if ! command -v "$RSCRIPT" >/dev/null 2>&1; then
  echo "Rscript not found. Add R to your PATH or set RSCRIPT=/path/to/Rscript." >&2
  exit 1
fi

# Apply --from: drop the scripts before START_FROM
if [[ -n "$START_FROM" ]]; then
  found=0; kept=()
  for s in "${SCRIPTS[@]}"; do
    [[ "$s" == "$START_FROM" ]] && found=1
    [[ $found -eq 1 ]] && kept+=("$s")
  done
  if [[ $found -eq 0 ]]; then
    echo "--from: '$START_FROM' is not in the script list." >&2; exit 2
  fi
  SCRIPTS=("${kept[@]}")
fi

# ── Check that every script is there before starting ─────────────────────────
missing=()
for s in "${SCRIPTS[@]}"; do
  [[ -f "$WORK_DIR/$s" ]] || missing+=("$s")
done
if [[ ${#missing[@]} -gt 0 ]]; then
  echo "Script(s) not found in $WORK_DIR: ${missing[*]}" >&2
  exit 1
fi

LOG_DIR="$WORK_DIR/logs"
mkdir -p "$LOG_DIR"
SUMMARY="$LOG_DIR/run_all_summary.csv"
echo "script,status,minutes" > "$SUMMARY"

# ── Run ───────────────────────────────────────────────────────────────────────
cd "$WORK_DIR" || exit 1
n=${#SCRIPTS[@]}
statuses=()
failed=0
t_all=$(date +%s)

for i in "${!SCRIPTS[@]}"; do
  s="${SCRIPTS[$i]}"
  log="$LOG_DIR/${s%.R}.log"
  echo ""
  echo "[$((i + 1))/$n] $s  (log: logs/$(basename "$log"))"

  t0=$(date +%s)
  "$RSCRIPT" "$WORK_DIR/$s" > "$log" 2>&1
  code=$?
  mins=$(awk -v s="$(( $(date +%s) - t0 ))" 'BEGIN { printf "%.1f", s / 60 }')

  if [[ $code -eq 0 ]]; then
    status="OK"
  else
    status="FAILED (exit $code)"
    failed=1
  fi
  statuses[$i]="$status"
  echo "       $status in $mins min"
  echo "$s,$status,$mins" >> "$SUMMARY"

  if [[ $code -ne 0 ]]; then
    echo "       Last lines of the log:"
    tail -n 15 "$log" | sed 's/^/         /'
    if [[ $STOP_ON_ERROR -eq 1 ]]; then
      for j in $(seq $((i + 1)) $((n - 1))); do
        statuses[$j]="not run"
        echo "${SCRIPTS[$j]},not run," >> "$SUMMARY"
      done
      break
    fi
  fi
done

# ── Summary ───────────────────────────────────────────────────────────────────
total=$(awk -v s="$(( $(date +%s) - t_all ))" 'BEGIN { printf "%.1f", s / 60 }')
echo ""
echo "================ SUMMARY ($total min total) ================"
for i in "${!SCRIPTS[@]}"; do
  printf "  %-26s %s\n" "${SCRIPTS[$i]}" "${statuses[$i]:-not run}"
done

if [[ $failed -ne 0 ]]; then
  echo ""
  echo "Not every script completed — see the summary above and the logs in $LOG_DIR" >&2
  exit 1
fi
echo ""
echo "All scripts completed."
