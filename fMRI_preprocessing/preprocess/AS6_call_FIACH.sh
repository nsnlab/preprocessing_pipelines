#!/bin/bash
# ============================================================
# AS6_call_FIACH.sh
# Runs the packaged MATLAB FIACH (fiach.m) headless via AS6_fiach_pipeline.m
#
# Inputs (exported by wrapper_function_preprocess.sh):
#   SUB, DIR, WORK_DIR, RUN_NAME, fiach_file, MOTION_PAR, BOLD_JSON
#   FIACH_MODE (task|rs), FIACH_GMM (R|ModifiedICH), FIACH_MASK (user|auto)
#   RS_INCLUDE_MOTION (true|false), SPM_DIR
#
# FIACH code: fiach.m and fiach_bold_contrast.m sit in this folder (preprocess/)
#
# Output used downstream (unchanged name):
#   ${WORK_DIR}/FIACH/rclean_${fiach_file}.nii
#   ${WORK_DIR}/FIACH/noisyPCs_${fiach_file}.tsv   (6 PCs, for AS7)
# ============================================================

sDIR="$(dirname $(realpath $0))/"

for var in SUB fiach_file WORK_DIR DIR RUN_NAME MOTION_PAR; do
    if [ -z "${!var}" ]; then echo "ERROR: ${var} not defined."; exit 1; fi
done

FIACH_MODE=${FIACH_MODE:-task}
FIACH_GMM=${FIACH_GMM:-R}
FIACH_MASK=${FIACH_MASK:-user}
RS_INCLUDE_MOTION=${RS_INCLUDE_MOTION:-true}

if [[ "$FIACH_MODE" != "task" && "$FIACH_MODE" != "rs" ]]; then
    echo "ERROR: FIACH_MODE must be 'task' or 'rs' (got '${FIACH_MODE}')"; exit 1
fi

# SPM: by default use the SPM already on your MATLAB path (saved path / startup.m).
# Set SPM_DIR only to force a specific copy:  SPM_DIR=/path/to/spm12 ./wrapper_function_preprocess.sh
SPM_DIR="${SPM_DIR:-}"
if [ -n "${SPM_DIR}" ]; then
    if [ ! -f "${SPM_DIR}/spm.m" ]; then
        echo "ERROR: SPM not found at ${SPM_DIR} (no spm.m)."; exit 1
    fi
    echo "  SPM path: ${SPM_DIR} (forced)"
else
    echo "  SPM path: from MATLAB's own path"
fi

# --- Resolve JSON sidecar ---
if [ -n "${BOLD_JSON}" ] && [ -f "${BOLD_JSON}" ]; then
    JSON="${BOLD_JSON}"
else
    echo "WARNING: BOLD_JSON not set, falling back to glob (may be ambiguous)"
    JSON=$(ls "${DIR}${SUB}/func/"*bold.json 2>/dev/null | head -n 1)
    if [ -z "${JSON}" ]; then
        echo "ERROR: No bold.json sidecar found for ${SUB}"
        exit 1
    fi
fi

echo "Using JSON: ${JSON}"
TR=$(jq -r '.RepetitionTime' "${JSON}")
TE=$(jq -r '.EchoTime'       "${JSON}")
if [ "${TR}" = "null" ] || [ -z "${TR}" ]; then echo "ERROR: No RepetitionTime in ${JSON}"; exit 1; fi
if [ "${TE}" = "null" ] || [ -z "${TE}" ]; then echo "ERROR: No EchoTime in ${JSON}";       exit 1; fi
# fiach_bold_contrast needs TE as whole milliseconds (BIDS stores seconds)
TE_MS=$(awk -v te="${TE}" 'BEGIN { if (te < 1) te = te * 1000; printf "%d", te + 0.5 }')
echo "TR: ${TR}s  TE: ${TE}s (${TE_MS} ms)"

B0=7
FIACH_DIR="${WORK_DIR}/FIACH"

INPUT="${FIACH_DIR}/${fiach_file}.nii"
if [ ! -f "${INPUT}" ]; then echo "ERROR: FIACH input not found: ${INPUT}"; exit 1; fi
if [ ! -f "${MOTION_PAR}" ]; then echo "ERROR: Motion parameters not found: ${MOTION_PAR}"; exit 1; fi

if [ "${FIACH_MASK}" = "user" ]; then
    MASK="${FIACH_DIR}/rfBrainMask_${RUN_NAME}.nii"
    if [ ! -f "${MASK}" ]; then echo "ERROR: Brain mask not found: ${MASK}"; exit 1; fi
else
    MASK=""
fi

echo "FIACH mode: ${FIACH_MODE} | GMM: ${FIACH_GMM} | mask: ${FIACH_MASK}"

# Statements go straight to matlab -batch (no temp .m file), so run names
# with hyphens such as task-face stay as they are everywhere.
# cfg.fiachRoot = preprocess/ (this folder), where fiach.m and fiach_bold_contrast.m live.
MATLAB_CMD="cfg = struct();"
MATLAB_CMD+=" cfg.fiachRoot = '${sDIR}';"
MATLAB_CMD+=" cfg.spmRoot = '${SPM_DIR}';"
MATLAB_CMD+=" cfg.fiachDir = '${FIACH_DIR}';"
MATLAB_CMD+=" cfg.input = '${INPUT}';"
MATLAB_CMD+=" cfg.mask = '${MASK}';"
MATLAB_CMD+=" cfg.motionPar = '${MOTION_PAR}';"
MATLAB_CMD+=" cfg.TR = ${TR};"
MATLAB_CMD+=" cfg.TEms = ${TE_MS};"
MATLAB_CMD+=" cfg.B0 = ${B0};"
MATLAB_CMD+=" cfg.mode = '${FIACH_MODE}';"
MATLAB_CMD+=" cfg.gmmMethod = '${FIACH_GMM}';"
MATLAB_CMD+=" cfg.posteriorThreshold = 0.5;"
MATLAB_CMD+=" cfg.contrastProfile = 'Conservative';"
MATLAB_CMD+=" cfg.nMads = 1.96;"
MATLAB_CMD+=" cfg.freq = 128;"
MATLAB_CMD+=" cfg.maxGap = 1;"
MATLAB_CMD+=" cfg.rsIncludeMotion = ${RS_INCLUDE_MOTION};"
MATLAB_CMD+=" run('${sDIR}AS6_fiach_pipeline.m');"

LOG_FILE="${FIACH_DIR}/fiach_log_${RUN_NAME}.txt"
echo "Running FIACH via matlab -batch"
matlab -nosplash -batch "${MATLAB_CMD}" 2>&1 | tee "${LOG_FILE}"
MATLAB_EXIT=${PIPESTATUS[0]}

if [ ${MATLAB_EXIT} -ne 0 ]; then
    echo "ERROR: MATLAB exited with code ${MATLAB_EXIT} during FIACH"
    exit 1
fi