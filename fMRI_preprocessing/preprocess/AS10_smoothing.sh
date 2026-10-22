#!/bin/bash

if [ -z "${SUB}" ];          then echo "ERROR: SUB not defined.";          exit 1; fi
if [ -z "${smooth_file}" ];  then echo "ERROR: smooth_file not defined.";  exit 1; fi
if [ -z "${FUNC_OUT_DIR}" ]; then echo "ERROR: FUNC_OUT_DIR not defined."; exit 1; fi

if [ -z "${FWHM_MM}" ]; then
    echo "WARNING: FWHM_MM not set, using default of 6mm"
    FWHM_MM=6
fi

sigma=$(echo "scale=6; ${FWHM_MM} / 2.3548" | bc)
echo "  FWHM: ${FWHM_MM}mm  →  sigma: ${sigma}"

# Input may be .nii or .nii.gz
INPUT_NII="${FUNC_OUT_DIR}/${smooth_file}.nii"
INPUT_NIIGZ="${FUNC_OUT_DIR}/${smooth_file}.nii.gz"

if [ -f "${INPUT_NII}" ]; then
    INPUT="${INPUT_NII}"
elif [ -f "${INPUT_NIIGZ}" ]; then
    INPUT="${INPUT_NIIGZ}"
else
    echo "ERROR: Smoothing input not found at either:"
    echo "  ${INPUT_NII}"
    echo "  ${INPUT_NIIGZ}"
    exit 1
fi

echo "  Smoothing: ${INPUT}"
fslmaths "${INPUT}" -s "${sigma}" "${FUNC_OUT_DIR}/${smooth_file}_smooth.nii.gz"