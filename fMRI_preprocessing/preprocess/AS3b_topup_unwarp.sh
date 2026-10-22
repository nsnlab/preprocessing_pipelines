#!/bin/bash
# 
# Author: Laith Alexander
# 
# AS3b_topup_unwarp.sh -- reverse-phase-encode (blip-up/blip-down) distortion
# correction, as a drop-in replacement for AS3_unwarping.sh.
#
# WHY THIS EXISTS
#   AS3_unwarping.sh does fugue on a Siemens GRE phase-difference fieldmap
#   (fmap/*_e1.nii.gz + fmap/*_e2_ph.nii.gz), with --unwarpdir=y- hardcoded.
#   The LAPilot has no fmap/ -- it has reverse-PE EPIs in bold_revphase_run-*.
#   Nothing in the existing pipeline can consume those (grep for topup/acqparams/
#   PhaseEncodingDirection across Scripts/ returns nothing), so this step is new.
#
# IT WRITES THE SAME OUTPUT FILENAME AS3 DOES:
#   ${WORK_DIR}/mc/${mc_file}_mcf_unwarp.nii.gz
# plus the mean/uncompressed copies steps 4, 5b and 6 expect. So the rest of the
# pipeline runs unchanged -- swap the call in wrapper_function_preprocess.sh:
#
#     run_or_die 3 "unwarping" bash ${sDIR}AS3b_topup_unwarp.sh
#
# REQUIRED ENVIRONMENT (the wrapper exports all but REVPHASE)
#   SUB RUN_NAME mc_file WORK_DIR BOLD_JSON   and   REVPHASE=/path/to/PA.nii.gz
#
# OPTIONAL
#   PE_DIR_FWD  phase-encode of the task run, as an acqparams triplet.
#               Default "0 -1 0"  (= BIDS j-, i.e. A>>P).
#   PE_DIR_PA   phase-encode of the reverse file. Default "0 1 0" (= j, P>>A).
#   NPA         volumes to average from the PA file (default: all).
#   NFWD        volumes to average from the start of the task run (default 5).
#
# *** POLARITY IS NOT IN THE JSON ***
#   dcm2niix cannot recover PhaseEncodingDirection for the custom dzne_ep3d on
#   XA60, so the defaults above are an ASSUMPTION. Verify empirically first:
#   01_qc_appraise.py reports the y centre-of-mass shift between each run and its
#   PA. If |pa_dy| < 0.5 voxels the pair is NOT opposed and topup will estimate a
#   near-null field and silently do nothing -- exactly the failure mode that bit
#   AS3's dwell-time bug. If the sign of the correction comes out backwards
#   (frontal lobe stretched instead of compressed), swap PE_DIR_FWD and PE_DIR_PA.

set -eo pipefail

for v in SUB RUN_NAME mc_file WORK_DIR; do
    if [ -z "${!v}" ]; then echo "ERROR: ${v} not defined."; exit 1; fi
done

# ---------------------------------------------------------------- find the PA
# The wrapper assigns TARGET_SUB / TARGET_RUN / DIR unconditionally in its config
# block, so passing them as environment variables on the command line does
# NOTHING -- they are overwritten before the loop starts. That rules out driving
# the wrapper one run at a time just to export a per-run REVPHASE.
#
# So resolve it here instead, from the manifest 03_stage_for_pipeline.py writes,
# keyed on SUB and RUN_NAME, which the wrapper DOES export. One wrapper run then
# handles every run with no per-run plumbing. An explicit REVPHASE still wins.
if [ -z "${REVPHASE}" ]; then
    MANIFEST="${REVPHASE_MANIFEST:-${DIR}staged_manifest.tsv}"
    if [ -f "${MANIFEST}" ]; then
        # resolve the column BY HEADER NAME -- the manifest's revphase column is
        # the 10th, not the 8th; a fixed index here once pointed at the TR
        REVPHASE=$(awk -F'\t' -v s="${SUB}" -v r="${RUN_NAME}" \
                   'NR==1 {for (i=1; i<=NF; i++) if ($i=="revphase") c=i; next}
                    $1==s && $2==r && c {print $c}' "${MANIFEST}")
        [ -n "${REVPHASE}" ] && echo "    REVPHASE from ${MANIFEST}"
    fi
fi
if [ -z "${REVPHASE}" ]; then
    echo "ERROR: no reverse-PE file for ${SUB}/${RUN_NAME}."
    echo "       Either export REVPHASE=/path/to/PA.nii.gz, or point"
    echo "       REVPHASE_MANIFEST at the staged_manifest.tsv."
    exit 1
fi
if [ ! -f "${REVPHASE}" ]; then
    echo "ERROR: REVPHASE not found: ${REVPHASE}"; exit 1
fi

MC_DIR="${WORK_DIR}/mc"
FIACH_DIR="${WORK_DIR}/FIACH"
TOPUP_DIR="${WORK_DIR}/topup"
mkdir -p "${TOPUP_DIR}" "${FIACH_DIR}"

INPUT_EPI="${MC_DIR}/${mc_file}_mcf.nii.gz"
OUTPUT_EPI="${MC_DIR}/${mc_file}_mcf_unwarp.nii.gz"
if [ ! -f "${INPUT_EPI}" ]; then
    echo "ERROR: motion-corrected input not found: ${INPUT_EPI}"
    echo "       Run step 2 first."; exit 1
fi

PE_DIR_FWD="${PE_DIR_FWD:-0 -1 0}"
PE_DIR_PA="${PE_DIR_PA:-0 1 0}"
NFWD="${NFWD:-5}"

# ---------------------------------------------------------------- readout time
# topup needs the total readout time (4th column of acqparams.txt). Use
# TotalReadoutTime; fall back to EffectiveEchoSpacing * (PE lines - 1).
# NEVER use DwellTime -- for these data it is ~37x too small, which is the same
# class of error that made AS3's fugue call a silent no-op.
read_trt () {
    local js="$1"
    local trt ees npe
    trt=$(jq -r '.TotalReadoutTime // empty' "${js}" 2>/dev/null)
    if [ -n "${trt}" ]; then echo "${trt}"; return; fi
    ees=$(jq -r '.EffectiveEchoSpacing // empty' "${js}" 2>/dev/null)
    npe=$(jq -r '.AcquisitionMatrixPE // .ReconMatrixPE // empty' "${js}" 2>/dev/null)
    if [ -n "${ees}" ] && [ -n "${npe}" ]; then
        echo "scale=8; ${ees} * (${npe} - 1)" | bc
    fi
}

TRT_FWD=""
[ -n "${BOLD_JSON}" ] && [ -f "${BOLD_JSON}" ] && TRT_FWD=$(read_trt "${BOLD_JSON}")
PA_JSON="${REVPHASE%.nii.gz}.json"
TRT_PA=""
[ -f "${PA_JSON}" ] && TRT_PA=$(read_trt "${PA_JSON}")

if [ -z "${TRT_FWD}" ]; then
    echo "ERROR: no TotalReadoutTime (or EffectiveEchoSpacing x matrix) in the"
    echo "       task-run JSON. Refusing to guess -- the readout time scales the"
    echo "       whole correction. Export BOLD_JSON to a JSON that has it."
    exit 1
fi
[ -z "${TRT_PA}" ] && TRT_PA="${TRT_FWD}"

echo "=== AS3b topup  ${SUB} / ${RUN_NAME}"
echo "    forward : ${INPUT_EPI}"
echo "    reverse : ${REVPHASE}"
echo "    PE fwd  : ${PE_DIR_FWD}   readout ${TRT_FWD} s"
echo "    PE rev  : ${PE_DIR_PA}   readout ${TRT_PA} s"

# ---------------------------------------------------------------- grid check
d_in=$(fslinfo "${INPUT_EPI}"  | awk '/^dim[123]/{printf "%s ", $2}')
d_pa=$(fslinfo "${REVPHASE}"   | awk '/^dim[123]/{printf "%s ", $2}')
if [ "${d_in}" != "${d_pa}" ]; then
    echo "ERROR: grid mismatch -- task run [${d_in}] vs reverse [${d_pa}]."
    echo "       topup requires identical matrices. Do not resample your way out"
    echo "       of this; the PA belongs to a different protocol."
    exit 1
fi

# ---------------------------------------------------------------- b0 pair
echo "--- building the blip pair"
fslroi "${INPUT_EPI}" "${TOPUP_DIR}/fwd_sel" 0 "${NFWD}"
fslmaths "${TOPUP_DIR}/fwd_sel" -Tmean "${TOPUP_DIR}/fwd_mean"
fslmaths "${REVPHASE}"          -Tmean "${TOPUP_DIR}/rev_mean"
fslmerge -t "${TOPUP_DIR}/blip_pair" "${TOPUP_DIR}/fwd_mean" "${TOPUP_DIR}/rev_mean"

printf '%s %s\n%s %s\n' "${PE_DIR_FWD}" "${TRT_FWD}" \
                        "${PE_DIR_PA}"  "${TRT_PA}" > "${TOPUP_DIR}/acqparams.txt"
echo "    acqparams.txt:"
sed 's/^/      /' "${TOPUP_DIR}/acqparams.txt"

# ---------------------------------------------------------------- config
# topup's default b02b0.cnf subsamples by 1,2,4,8 and needs each in-plane
# dimension divisible by the subsampling factor. These matrices (e.g. 180) are
# not divisible by 8, so pick a config that will actually run.
CNF="${FSLDIR}/etc/flirtsch/b02b0.cnf"
dx=$(fslval "${INPUT_EPI}" dim1 | tr -d ' ')
dy=$(fslval "${INPUT_EPI}" dim2 | tr -d ' ')
if [ $((dx % 8)) -ne 0 ] || [ $((dy % 8)) -ne 0 ]; then
    if [ -f "${FSLDIR}/etc/flirtsch/b02b0_1.cnf" ]; then
        CNF="${FSLDIR}/etc/flirtsch/b02b0_1.cnf"
        echo "    matrix ${dx}x${dy} is not divisible by 8 -> using b02b0_1.cnf"
    else
        echo "    WARNING matrix ${dx}x${dy} not divisible by 8 and b02b0_1.cnf is"
        echo "            absent. topup may fail; if it does, zero-pad to a"
        echo "            multiple of 8 or install a newer FSL."
    fi
fi

# ---------------------------------------------------------------- topup
echo "--- topup (this is the slow part)"
topup --imain="${TOPUP_DIR}/blip_pair" \
      --datain="${TOPUP_DIR}/acqparams.txt" \
      --config="${CNF}" \
      --out="${TOPUP_DIR}/topup_results" \
      --fout="${TOPUP_DIR}/topup_field_Hz" \
      --iout="${TOPUP_DIR}/blip_pair_corrected" \
      --verbose

# --- sanity: did topup find a real field, or a null one? ---------------------
FMAX=$(fslstats "${TOPUP_DIR}/topup_field_Hz" -a -P 99 | awk '{print $1}')
echo "    99th-percentile |field| = ${FMAX} Hz"
NULLFIELD=$(echo "${FMAX} < 5" | bc -l)
if [ "${NULLFIELD}" = "1" ]; then
    echo ""
    echo "  *** WARNING: the estimated field is essentially zero. ***"
    echo "      That is what you get when the two inputs have the SAME phase-"
    echo "      encode polarity. Check 01_qc_appraise.py's pa_dy for this run"
    echo "      before using this output -- topup will not error, it will just"
    echo "      hand back your input barely changed."
    echo ""
fi

# ---------------------------------------------------------------- apply
echo "--- applytopup"
applytopup --imain="${INPUT_EPI}" \
           --datain="${TOPUP_DIR}/acqparams.txt" \
           --inindex=1 \
           --topup="${TOPUP_DIR}/topup_results" \
           --method=jac \
           --out="${OUTPUT_EPI}"

if [ ! -f "${OUTPUT_EPI}" ]; then
    echo "ERROR: applytopup produced no output."; exit 1
fi

# ---------------------------------------------------------------- downstream
# Reproduce exactly what AS3 leaves behind, so steps 4, 5b and 6 find their
# inputs where they expect them.
echo "--- mean functional and uncompressed copies"
fslmaths "${OUTPUT_EPI}" -Tmean "${MC_DIR}/meanFunctional_${RUN_NAME}.nii.gz"
gunzip -c "${MC_DIR}/meanFunctional_${RUN_NAME}.nii.gz" \
        > "${FIACH_DIR}/meanFunctional_${RUN_NAME}.nii"
gunzip -c "${OUTPUT_EPI}" > "${FIACH_DIR}/${mc_file}_mcf_unwarp.nii"

echo "=== AS3b done"
echo "    ${OUTPUT_EPI}"
echo "    ${TOPUP_DIR}/topup_field_Hz.nii.gz   <- inspect this before trusting the run"
echo "    ${TOPUP_DIR}/blip_pair_corrected.nii.gz  <- the two volumes should now agree"
