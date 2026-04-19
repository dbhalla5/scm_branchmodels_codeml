#!/bin/bash
set -euo pipefail

# ------------------------------------------------------------------
# Usage:
#   bash run_branch_models_subset_parallel.sh CYP4A5 Primary.Lifestyle 20
#
# Args:
#   $1 = Gene name
#   $2 = Trait name
#   $3 = (optional) number of SCM trees to sample (default = 20)
# ------------------------------------------------------------------

GENE="$1"
TRAIT="$2"
NSAMPLE="${3:-20}"    # number of SCM trees to run
MAX_JOBS=20           # limit concurrent codeml runs

BASE="/path/to/SCM"
PHY="$BASE/sub_trees/${GENE}.phy"
TREE_DIR="$BASE/${GENE}/${TRAIT}/trees_paml_trimmed"
OUT_DIR="$BASE/branchmodel_files/${GENE}/${TRAIT}"

mkdir -p "$OUT_DIR"

# ---- Sanity checks ----
[[ -f "$PHY" ]] || { echo "[ERROR] Phy file not found: $PHY" >&2; exit 1; }
[[ -d "$TREE_DIR" ]] || { echo "[ERROR] Tree directory not found: $TREE_DIR" >&2; exit 1; }

# ---- Select subset of trees ----
echo "[INFO] Sampling up to $NSAMPLE SCM trees from $TREE_DIR"
cd "$TREE_DIR"
ALL_TREES=(*_paml.nh)
NTOTAL=${#ALL_TREES[@]}
(( NTOTAL > 0 )) || { echo "[ERROR] No *_paml.nh trees found in $TREE_DIR" >&2; exit 1; }

# Pick random subset reproducibly if desired
shuf -n "$NSAMPLE" < <(printf "%s\n" "${ALL_TREES[@]}") > "$OUT_DIR/selected_trees.txt"
echo "[INFO] Selected $(wc -l < "$OUT_DIR/selected_trees.txt") trees (saved to selected_trees.txt)"

# ---- Helper: limit concurrent jobs ----
running_jobs() { jobs -rp | wc -l; }

# ---- Main loop ----
for TREE in $(cat "$OUT_DIR/selected_trees.txt"); do
  TREE_PATH="$TREE_DIR/$TREE"
  BASE_NAME=$(basename "$TREE" .nh)
  RUN_DIR="$OUT_DIR/$BASE_NAME"
  mkdir -p "$RUN_DIR"

  echo "[INFO] Preparing codeml configs for $BASE_NAME"

  # --- ALT model ---
  cat > "$RUN_DIR/alt.ctl" <<EOF
      seqfile = $PHY
      treefile = $TREE_PATH
      outfile  = $RUN_DIR/${BASE_NAME}_alt.txt
      noisy = 3
      verbose = 1
      runmode = 0
      seqtype = 1
      CodonFreq = 2
      model = 2
      fix_omega = 0
      omega = 0.2
      cleandata = 1
EOF

  # --- NULL model ---
  cat > "$RUN_DIR/null.ctl" <<EOF
      seqfile = $PHY
      treefile = $TREE_PATH
      outfile  = $RUN_DIR/${BASE_NAME}_null.txt
      noisy = 3
      verbose = 1
      runmode = 0
      seqtype = 1
      CodonFreq = 2
      model = 2
      fix_omega = 1
      omega = 1
      cleandata = 1
EOF

  # --- Launch both models sequentially ---
  (
    cd "$RUN_DIR"
    echo "[RUN] $(basename "$RUN_DIR")"
    codeml alt.ctl > alt.log 2>&1
    codeml null.ctl > null.log 2>&1
  ) &

  # --- Parallel throttle ---
  while (( $(running_jobs) >= MAX_JOBS )); do
    sleep 15
  done
done

wait
echo "[INFO] All codeml branch-model runs completed successfully."

