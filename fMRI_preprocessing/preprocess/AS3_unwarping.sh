#!/bin/bash
# ============================================================
# AS3_unwarping.sh
# GRE fieldmap-based EPI unwarping using FSL.
# Expects env vars: SUB, RUN_NAME, mc_file, BOLD_JSON,
#                   WORK_DIR, FMAP_DIR
# ============================================================

if [ -z "${SUB}" ];       then echo "ERROR: SUB not defined.";       exit 1; fi
if [ -z "${mc_file}" ];   then echo "ERROR: mc_file not defined.";   exit 1; fi
if [ -z "${WORK_DIR}" ];  then echo "ERROR: WORK_DIR not defined.";  exit 1; fi
if [ -z "${FMAP_DIR}" ] || [ ! -d "${FMAP_DIR}" ]; then
    echo "ERROR: FMAP_DIR not set or does not exist: ${FMAP_DIR}"; exit 1
fi

echo "  Using fmap directory: ${FMAP_DIR}"

# Quick check that the B0 and EPI are in the same space
echo '----------- fieldmap and EPI need to be in the same space-----------'
fslinfo .../fmap/fieldmap_rads.nii.gz | grep -E "^dim[1-3]|^pixdim[1-3]"
fslinfo .../mc/<mc_file>_mcf.nii.gz   | grep -E "^dim[1-3]|^pixdim[1-3]"




MC_DIR="${WORK_DIR}/mc/"
FMAP_OUT="${WORK_DIR}/fmap/"
FIACH_DIR="${WORK_DIR}/FIACH/"
mkdir -p "${FMAP_OUT}"

# --- Find fieldmap files ---
MAGN=$(ls "${FMAP_DIR}"*_e1.nii.gz 2>/dev/null | head -n 1)
PHASE=$(ls "${FMAP_DIR}"*_e2_ph.nii.gz 2>/dev/null | head -n 1)
if [ -z "${MAGN}" ];  then echo "ERROR: No magnitude (*_e1.nii.gz) in ${FMAP_DIR}";     exit 1; fi
if [ -z "${PHASE}" ]; then echo "ERROR: No phase diff (*_e2_ph.nii.gz) in ${FMAP_DIR}"; exit 1; fi
echo "  Magnitude : ${MAGN}"
echo "  Phase diff: ${PHASE}"

# --- Read dTE from JSON ---
dTE=""
E1_JSON=$(ls "${FMAP_DIR}"*_e1.json 2>/dev/null | head -n 1)
E2_JSON=$(ls "${FMAP_DIR}"*_e2.json 2>/dev/null | head -n 1)

if [ -n "${E1_JSON}" ] && [ -n "${E2_JSON}" ]; then
    TE1=$(jq -r '.EchoTime // empty' "${E1_JSON}" 2>/dev/null)
    TE2=$(jq -r '.EchoTime // empty' "${E2_JSON}" 2>/dev/null)
    if [ -n "${TE1}" ] && [ -n "${TE2}" ]; then
        dTE=$(echo "scale=6; (${TE2} - ${TE1}) * 1000" | bc)
        echo "  dTE from e1/e2 JSONs: ${TE1}s / ${TE2}s → ${dTE} ms"
    else
        echo "  WARNING: Could not read EchoTime from e1/e2 JSONs"
    fi
else
    echo "  WARNING: Could not find e1/e2 JSON files in ${FMAP_DIR}"
fi

if [ -z "${dTE}" ]; then
    echo "ERROR: could not read dTE from the e1/e2 fieldmap JSONs in ${FMAP_DIR}"; exit 1
fi

# --- Brain extract magnitude ---
echo "  Brain extracting magnitude..."
MAGN_BRAIN="${FMAP_OUT}magnitude_brain"
bet "${MAGN}" "${MAGN_BRAIN}" -f 0.5 -g 0 -m

# --- Prepare fieldmap ---
echo "  Preparing fieldmap (dTE=${dTE} ms)..."
FIELDMAP="${FMAP_OUT}fieldmap_rads"
fsl_prepare_fieldmap SIEMENS "${PHASE}" "${MAGN_BRAIN}" "${FIELDMAP}" "${dTE}"

# --- Read EPI params ---
# fugue --dwell needs the EFFECTIVE echo spacing (s). DwellTime in Siemens
# JSONs is the readout sample spacing (~1e-6 s) and must not be used here.
if [ -z "${BOLD_JSON}" ] || [ ! -f "${BOLD_JSON}" ]; then
    echo "ERROR: BOLD_JSON not set or missing: ${BOLD_JSON}"; exit 1
fi

EPI_TE=$(jq -r '.EchoTime // empty'             "${BOLD_JSON}" 2>/dev/null)
DWELL_TIME=$(jq -r '.EffectiveEchoSpacing // empty' "${BOLD_JSON}" 2>/dev/null)
if [ -z "${EPI_TE}" ];     then echo "ERROR: No EchoTime in ${BOLD_JSON}"; exit 1; fi
if [ -z "${DWELL_TIME}" ]; then echo "ERROR: No EffectiveEchoSpacing in ${BOLD_JSON} (needed for fugue --dwell)"; exit 1; fi
echo "  EPI TE=${EPI_TE}s  effective echo spacing=${DWELL_TIME}s"


# --- Unwarp ---
INPUT_EPI="${MC_DIR}${mc_file}_mcf.nii.gz"
OUTPUT_EPI="${MC_DIR}${mc_file}_mcf_unwarp.nii.gz"
if [ ! -f "${INPUT_EPI}" ]; then echo "ERROR: EPI not found: ${INPUT_EPI}"; exit 1; fi

echo "  Unwarping EPI..."
fugue -i "${INPUT_EPI}" \
      --loadfmap="${FIELDMAP}" \
      --dwell="${DWELL_TIME}" \
      --unwarpdir=y- \
      -u "${OUTPUT_EPI}"
echo "  Unwarped EPI saved: ${OUTPUT_EPI}"

# --- Compute mean functional ---
MEAN_OUT="${MC_DIR}meanFunctional_${RUN_NAME}.nii.gz"
fslmaths "${OUTPUT_EPI}" -Tmean "${MEAN_OUT}"
gunzip -c "${MEAN_OUT}"   > "${FIACH_DIR}meanFunctional_${RUN_NAME}.nii"
gunzip -c "${OUTPUT_EPI}" > "${FIACH_DIR}${mc_file}_mcf_unwarp.nii"
echo "  Mean functional saved: ${MEAN_OUT}"
echo "  Unwarping complete for ${SUB} / ${RUN_NAME}"