#!/bin/bash

# ============================================================
# CYP4V2 branch-site CODEML batch runner
#
# Runs 10 SCM workflows concurrently.
#
# Within each SCM:
#       NULL -> ALT (omega=1) -> ALT (omega=2)
#
# Each CODEML command is executed from its own model directory.
#
# Existing CODEML output files are moved to a timestamped
# backup directory before a new run.
#
# Control files (*.ctl) and input trees/alignment files are
# NEVER removed or modified.
#
# Usage:
#       bash run_branch_site_batch.sh 1
#       bash run_branch_site_batch.sh 2
#       ...
#       bash run_branch_site_batch.sh 5
#
# Batch 1 = SCM_001-010
# Batch 2 = SCM_011-020
# Batch 3 = SCM_021-030
# Batch 4 = SCM_031-040
# Batch 5 = SCM_041-050
# ============================================================

set -u


# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

BASE="/path/to/SCM/branch_site_model/CYP4V2/Migration"    ####### change this accordingly. (this is for "CYP4V2" set)
RUNS="${BASE}/runs"


# ------------------------------------------------------------
# Check argument
# ------------------------------------------------------------

if [ "$#" -ne 1 ]; then
    echo "Usage: bash $0 <batch>"
    echo
    echo "  1 = SCM_001-010"
    echo "  2 = SCM_011-020"
    echo "  3 = SCM_021-030"
    echo "  4 = SCM_031-040"
    echo "  5 = SCM_041-050"
    exit 1
fi


case "$1" in

    1)
        START=1
        END=10
        ;;

    2)
        START=11
        END=20
        ;;

    3)
        START=21
        END=30
        ;;

    4)
        START=31
        END=40
        ;;

    5)
        START=41
        END=50
        ;;

    *)
        echo "[ERROR] Batch must be 1, 2, 3, 4, or 5."
        exit 1
        ;;

esac


# ------------------------------------------------------------
# Find CODEML
# ------------------------------------------------------------

CODEML=$(command -v codeml || true)

if [ -z "$CODEML" ]; then
    echo "[ERROR] codeml was not found in PATH."
    echo "[ERROR] Activate the PAML environment first."
    exit 1
fi


# ------------------------------------------------------------
# Check base directory
# ------------------------------------------------------------

if [ ! -d "$RUNS" ]; then
    echo "[ERROR] Run directory does not exist:"
    echo "        $RUNS"
    exit 1
fi


echo "============================================================"
echo " CYP4V2 branch-site CODEML"                                     ####### change this accordingly. (this is for "CYP4V2" set)
echo " Batch $1"
echo " SCM range: $(printf '%03d' "$START")-$(printf '%03d' "$END")"
echo " Concurrent SCM workflows: 10"
echo "============================================================"
echo
echo "[INFO] CODEML: $CODEML"
echo "[INFO] Runs:   $RUNS"
echo


# ------------------------------------------------------------
# Check all required directories and CTL files BEFORE
# launching anything.
# ------------------------------------------------------------

echo "[INFO] Checking batch..."

for ((i=START; i<=END; i++)); do

    SCM=$(printf "%03d" "$i")
    SCM_DIR="${RUNS}/SCM_${SCM}"

    if [ ! -d "$SCM_DIR" ]; then
        echo "[ERROR] Missing directory:"
        echo "        $SCM_DIR"
        exit 1
    fi

    for MODEL in null alt omega_2; do

        MODEL_DIR="${SCM_DIR}/${MODEL}"

        if [ ! -d "$MODEL_DIR" ]; then
            echo "[ERROR] Missing directory:"
            echo "        $MODEL_DIR"
            exit 1
        fi

    done

    if [ ! -f "${SCM_DIR}/null/null.ctl" ]; then
        echo "[ERROR] Missing:"
        echo "        ${SCM_DIR}/null/null.ctl"
        exit 1
    fi

    if [ ! -f "${SCM_DIR}/alt/alt.ctl" ]; then
        echo "[ERROR] Missing:"
        echo "        ${SCM_DIR}/alt/alt.ctl"
        exit 1
    fi

    if [ ! -f "${SCM_DIR}/omega_2/omega_2.ctl" ]; then
        echo "[ERROR] Missing:"
        echo "        ${SCM_DIR}/omega_2/omega_2.ctl"
        exit 1
    fi

done

echo "[OK] All SCM directories and CTL files found."
echo


# ------------------------------------------------------------
# Function to back up existing CODEML output
# ------------------------------------------------------------

backup_existing_output() {

    MODEL_DIR="$1"

    # CODEML-generated files seen in the CYP4V2 runs.
    # Do NOT include *.ctl or input trees.
    OUTPUT_FILES=(
        "2NG.dN"
        "2NG.dS"
        "2NG.t"
        "2ML.dN"
        "2ML.dS"
        "2ML.t"
        "4fold.nuc"
        "lnf"
        "rst"
        "rst1"
        "rub"
        "mlc"
    )

    BACKUP_DIR="${MODEL_DIR}/previous_run_$(date '+%Y%m%d_%H%M%S')"

    FOUND=0

    for FILE in "${OUTPUT_FILES[@]}"; do

        if [ -e "${MODEL_DIR}/${FILE}" ]; then
            FOUND=1
            break
        fi

    done

    if [ "$FOUND" -eq 0 ]; then
        return 0
    fi

    mkdir -p "$BACKUP_DIR"

    for FILE in "${OUTPUT_FILES[@]}"; do

        if [ -e "${MODEL_DIR}/${FILE}" ]; then
            mv "${MODEL_DIR}/${FILE}" "${BACKUP_DIR}/"
        fi

    done

    echo "[BACKUP] ${MODEL_DIR}"
    echo "         Existing CODEML output moved to:"
    echo "         ${BACKUP_DIR}"

}


# ------------------------------------------------------------
# Run one SCM
# ------------------------------------------------------------

run_scm() {

    SCM="$1"
    SCM_DIR="${RUNS}/SCM_${SCM}"

    LOG="${SCM_DIR}/batch.log"

    # Start a fresh log for this workflow.
    {
        echo "============================================================"
        echo "SCM_${SCM}"
        echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "============================================================"
        echo
    } > "$LOG"


    # ========================================================
    # NULL
    # ========================================================

    MODEL_DIR="${SCM_DIR}/null"

    echo "[START] SCM_${SCM} NULL" >> "$LOG"
    echo "Directory: ${MODEL_DIR}" >> "$LOG"
    echo "Time: $(date '+%Y-%m-%d %H:%M:%S')" >> "$LOG"

    backup_existing_output "$MODEL_DIR" >> "$LOG" 2>&1

    (
        cd "$MODEL_DIR" || exit 1
        "$CODEML" null.ctl
    ) >> "$LOG" 2>&1

    STATUS=$?

    echo "CODEML exit status: ${STATUS}" >> "$LOG"

    if [ "$STATUS" -eq 0 ]; then
        echo "[DONE] SCM_${SCM} NULL" >> "$LOG"
    else
        echo "[WARNING] SCM_${SCM} NULL returned exit status ${STATUS}" >> "$LOG"
        echo "[INFO] Continuing to ALT omega=1." >> "$LOG"
    fi

    echo >> "$LOG"


    # ========================================================
    # ALT omega = 1
    # ========================================================

    MODEL_DIR="${SCM_DIR}/alt"

    echo "[START] SCM_${SCM} ALT omega=1" >> "$LOG"
    echo "Directory: ${MODEL_DIR}" >> "$LOG"
    echo "Time: $(date '+%Y-%m-%d %H:%M:%S')" >> "$LOG"

    backup_existing_output "$MODEL_DIR" >> "$LOG" 2>&1

    (
        cd "$MODEL_DIR" || exit 1
        "$CODEML" alt.ctl
    ) >> "$LOG" 2>&1

    STATUS=$?

    echo "CODEML exit status: ${STATUS}" >> "$LOG"

    if [ "$STATUS" -eq 0 ]; then
        echo "[DONE] SCM_${SCM} ALT omega=1" >> "$LOG"
    else
        echo "[WARNING] SCM_${SCM} ALT omega=1 returned exit status ${STATUS}" >> "$LOG"
        echo "[INFO] Continuing to ALT omega=2." >> "$LOG"
    fi

    echo >> "$LOG"


    # ========================================================
    # ALT omega = 2
    # ========================================================

    MODEL_DIR="${SCM_DIR}/omega_2"

    echo "[START] SCM_${SCM} ALT omega=2" >> "$LOG"
    echo "Directory: ${MODEL_DIR}" >> "$LOG"
    echo "Time: $(date '+%Y-%m-%d %H:%M:%S')" >> "$LOG"

    backup_existing_output "$MODEL_DIR" >> "$LOG" 2>&1

    (
        cd "$MODEL_DIR" || exit 1
        "$CODEML" omega_2.ctl
    ) >> "$LOG" 2>&1

    STATUS=$?

    echo "CODEML exit status: ${STATUS}" >> "$LOG"

    if [ "$STATUS" -eq 0 ]; then
        echo "[DONE] SCM_${SCM} ALT omega=2" >> "$LOG"
    else
        echo "[WARNING] SCM_${SCM} ALT omega=2 returned exit status ${STATUS}" >> "$LOG"
    fi

    echo >> "$LOG"


    # ========================================================
    # Finish SCM
    # ========================================================

    {
        echo "============================================================"
        echo "SCM_${SCM} workflow finished"
        echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "============================================================"
    } >> "$LOG"

}


# ------------------------------------------------------------
# Launch the 10 SCM workflows
# ------------------------------------------------------------

PIDS=()
SCMS=()

for ((i=START; i<=END; i++)); do

    SCM=$(printf "%03d" "$i")

    run_scm "$SCM" &

    PID=$!

    PIDS+=("$PID")
    SCMS+=("$SCM")

    echo "[LAUNCHED] SCM_${SCM} (PID ${PID})"

done


echo
echo "[INFO] All 10 SCM workflows launched."
echo "[INFO] Waiting for completion..."
echo


# ------------------------------------------------------------
# Wait for all SCM workflows
# ------------------------------------------------------------

FAILED=0

for ((j=0; j<${#PIDS[@]}; j++)); do

    PID="${PIDS[$j]}"
    SCM="${SCMS[$j]}"

    if wait "$PID"; then
        echo "[COMPLETE] SCM_${SCM}"
    else
        echo "[FAILED] SCM_${SCM} workflow"
        FAILED=1
    fi

done


echo

# ------------------------------------------------------------
# Final status
# ------------------------------------------------------------

if [ "$FAILED" -ne 0 ]; then

    echo "============================================================"
    echo "Batch $1 finished with workflow failures."
    echo "Check:"
    echo "${RUNS}/SCM_###/batch.log"
    echo "============================================================"

    exit 1

fi


echo "============================================================"
echo "Batch $1 completed."
echo "SCM range: $(printf '%03d' "$START")-$(printf '%03d' "$END")"
echo "============================================================"
