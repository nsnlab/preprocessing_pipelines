#!/bin/bash

if [ -z "${SUB}" ];          then echo "ERROR: SUB not defined.";          exit 1; fi
if [ -z "${RUN_NAME}" ];     then echo "ERROR: RUN_NAME not defined.";     exit 1; fi
if [ -z "${WORK_DIR}" ];     then echo "ERROR: WORK_DIR not defined.";     exit 1; fi
if [ -z "${FUNC_OUT_DIR}" ]; then echo "ERROR: FUNC_OUT_DIR not defined."; exit 1; fi
if [ -z "${T1_DIR}" ];       then echo "ERROR: T1_DIR not defined.";       exit 1; fi

# How many tSNR maps to compute at once. Each job loads a full 4-D run, so keep
# this modest on 7T data. Overridable from the shell: TSNR_JOBS=2 ./wrapper...
TSNR_JOBS="${TSNR_JOBS:-2}"

# --- Input files for tSNR ---
# Work-dir intermediates (native/EPI space)
# Final outputs (T1 and MNI space) are flat in FUNC_OUT_DIR
files=(
    "${WORK_DIR}/NORDIC/NORDIC_Run_BOLD_${SUB}_${RUN_LABEL}.nii"
    "${WORK_DIR}/mc/${mc_file}_mcf.nii.gz"
    "${WORK_DIR}/FIACH/${reg_file}.nii"
    "${FUNC_OUT_DIR}/reg_${reg_file}.nii.gz"
    "${FUNC_OUT_DIR}/norm_${reg_file}.nii.gz"
    "${FUNC_OUT_DIR}/${smooth_file}_smooth.nii.gz"
)

# --- Masks ---
mask_epi="${WORK_DIR}/FIACH/rfBrainMask_${RUN_NAME}.nii"
if [ ! -f "${mask_epi}" ]; then
    echo "ERROR: EPI brain mask not found: ${mask_epi}"
    exit 1
fi
echo "  Using EPI mask: ${mask_epi}"

mask_t1="${T1_DIR}/Masked_UNI.nii"
mask_mni="${FSLDIR}/data/standard/MNI152_T1_1mm_brain_mask.nii.gz"

TSNR_OUT="${FUNC_OUT_DIR}/tsnr"
mkdir -p "${TSNR_OUT}"

# --- One tSNR map (runs as a background job) ---
tsnr_one () {
    local file="$1" base mean std tsnr mask

    base=$(basename "$file" .nii.gz)
    base=$(basename "$base" .nii)

    mean="${TSNR_OUT}/${base}_mean.nii.gz"
    std="${TSNR_OUT}/${base}_std.nii.gz"
    tsnr="${TSNR_OUT}/${base}_tsnr.nii.gz"

    # Select mask based on filename prefix (space-aware)
    # norm_ prefix → MNI space; reg_ prefix → T1 space; else → EPI space
    if [[ "$base" == norm_* ]]; then
        mask="$mask_mni"
    elif [[ "$base" == reg_* ]]; then
        mask="$mask_t1"
    else
        mask="$mask_epi"
    fi

    fslmaths "$file" -mas "$mask" -Tmean "$mean" \
        && fslmaths "$file" -mas "$mask" -Tstd  "$std" \
        && fslmaths "$mean" -div "$std" "$tsnr" \
        && echo "    [${base}] Saved: $tsnr" \
        || echo "    [${base}] WARNING: tSNR failed"
}

echo "Computing tSNR maps for ${SUB} / ${RUN_LABEL} (up to ${TSNR_JOBS} at once):"
pids=()
for file in "${files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "  [SKIP] Not found: $file"
        continue
    fi

    # wait for a free slot before starting the next map
    while (( $(jobs -rp | wc -l) >= TSNR_JOBS )); do sleep 1; done

    echo "  Starting tSNR for: $(basename "$file")"
    tsnr_one "$file" &
    pids+=($!)
done

# wait for every map to finish before the step is marked done
for p in "${pids[@]}"; do wait "$p"; done

echo "tSNR complete for ${SUB} / ${RUN_LABEL}"