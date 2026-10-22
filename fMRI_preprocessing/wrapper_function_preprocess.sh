#!/bin/bash
# ============================================================
# wrapper_preprocess.sh
# 
# 
# v1.2.0
# Autor: Anna Sadilova
# Last updated: 21.10.2026
# 
# 
# Full preprocessing pipeline wrapper
# ============================================================

# --- Pipeline flags ---
run_NORDIC=true
run_FIACH=true
MNI_space=false
SLICE_TIMING=false
FWHM_MM=3   # smoothing kernel in mm (converted to sigma for fslmaths)
DUMMY_SCANS=3  # number of dummy scans to remove from start of each run (0 = disabled)
PIPELINE_TAG=""   # optional label for variant runs e.g. "smooth6mm" or "noFIACH"
                  # if set: outputs go to func/condition_TAG/ and env saved as pipeline_vars_LABEL_TAG.env
                  # if empty: standard paths, overwrite guard active (won't re-run if smooth output exists)

# --- NORDIC options (AS1_run_NORDIC.m, must sit in preprocess/) ---
NORDIC_PHASE="none"   # "none" = magnitude-only NORDIC
                      # "auto" = use a phase image next to the BIDS bold (part-phase / _ph / _phase)
                      #          if one exists, otherwise magnitude-only



# --- Selective processing ---
# Leave empty ("") to process everything, or set to filter:
#   TARGET_SUB="sub-01"      → only that subject
#   TARGET_CONDITION="cap2"  → only that subfolder (use "" to target direct-in-func files)
#   TARGET_RUN="run1"        → only that run name (as derived after stripping sub-XX_ prefix)
TARGET_SUB=""
TARGET_CONDITION=""
TARGET_RUN=""

# --- Restart / breakpoint control ---
# RESTART_FROM_STEP=0  → resume: use checkpoints, skip already-completed steps
# RESTART_FROM_STEP=N  → clear checkpoints from step N onward, rerun from step N
# RESTART_FROM_STEP=1  → clear all checkpoints, rerun everything from scratch
#
# Step reference:
#   1  Dummy scans + NORDIC      7  Combine regressors
#   2  Motion correction         8  Co-registration + Normalisation
#   3  Unwarping / mean func     9  (marked done with step 8)
#      (topup → GRE fieldmap → none)
#   4  Slice timing             10  Smoothing
#   5  SynthSeg + T1 mask       11  tSNR maps
#   5b Mask → EPI space         12  QC report
#   6  FIACH
RESTART_FROM_STEP=0


# ----- Directory and subjects -----
DIR="/path/to/BIDS/folder/"   # trailing slash required
SUBS=(
    'sub-01'
    'sub-02'
    'sub-03'
)




# ===== Additional settings ====
# --- FIACH options ---
# FIACH_MODE="task" → rclean_ = FIACH output only. Motion + FIACH PCs are NOT
#                     removed; AS7 combines them so they go into the task GLM.
# FIACH_MODE="rs"   → rclean_ = FIACH output with motion + FIACH PCs projected
#                     out voxelwise. Do NOT add these regressors again later.
FIACH_MODE="task"
FIACH_GMM="R"               # "R" (packaged FIACH default, upstream R behaviour) or "ModifiedICH"
FIACH_MASK="user"           # "user" = rfBrainMask_<RUN_NAME>.nii from step 5b, "auto" = FIACH mask
RS_INCLUDE_MOTION=true      # rs mode only (ignored in task mode): true = regress 6 motion + 6 PCs, false = 6 PCs only

# Parallel tSNR processing
TSNR_JOBS=2                 # how many tSNR maps to compute at once (step 11); lower if memory runs out


# ====================================================================================================


# --- Environment --- (need to set up first time running)
export DIR
sDIR="$(dirname $(realpath $0))/preprocess/"

export FSLDIR=$HOME/fsl
source $FSLDIR/etc/fslconf/fsl.sh
export PATH=$FSLDIR/bin:$PATH

export FREESURFER_HOME=${FREESURFER_HOME:-/home/freesurfer/8.1.0}
source "${FREESURFER_HOME}/SetUpFreeSurfer.sh"
export FS_LICENSE=/PROJECTS/scripts/.license


export FIACH_MODE FIACH_GMM FIACH_MASK RS_INCLUDE_MOTION PIPELINE_TAG TSNR_JOBS

export BROWSER=true # to silence zygote error when browser is already open


# === Only need to run the first time ===
chmod +x ${sDIR}AS*.sh

# --- Python virtual environment (Scripts/.venv) ---
# Created once; reused on later launches. Packages are reinstalled only when
# 1_requirements.txt changes (its checksum is stored inside the venv).
SCRIPTS_DIR="$(dirname $(realpath $0))"
VENV_DIR="${SCRIPTS_DIR}/.venv"
REQ_FILE="${SCRIPTS_DIR}/1_requirements.txt"

if [[ ! -f "${VENV_DIR}/bin/activate" ]]; then
    echo "Creating Python virtual environment: ${VENV_DIR}"
    python3 -m venv "${VENV_DIR}" \
        || { echo "ERROR: could not create venv (on Ubuntu: sudo apt install python3-venv)"; exit 1; }
fi
source "${VENV_DIR}/bin/activate"

REQ_HASH=$(sha256sum "${REQ_FILE}" | awk '{print $1}')
if [[ "$(cat "${VENV_DIR}/.requirements.sha256" 2>/dev/null)" != "${REQ_HASH}" ]]; then
    echo "Installing Python dependencies into ${VENV_DIR}..."
    pip install --upgrade pip \
        && pip install -r "${REQ_FILE}" \
        && echo "${REQ_HASH}" > "${VENV_DIR}/.requirements.sha256" \
        || { echo "ERROR: pip install failed"; exit 1; }
else
    echo "Python venv up to date: ${VENV_DIR}"
fi
echo "  Python: $(command -v python3)"




# ====================================================================================================







echo "============================================"
echo "Starting preprocessing pipeline"
echo "  Directory : ${DIR}"
echo "  Subjects  : ${SUBS[*]}"
[[ -n "$PIPELINE_TAG" ]] && echo "  Tag       : ${PIPELINE_TAG}"
[[ "$DUMMY_SCANS" -gt 0 ]] && echo "  Dummies   : removing first ${DUMMY_SCANS} volume(s) per run"
[[ "$run_NORDIC" = true ]] && echo "  NORDIC    : phase=${NORDIC_PHASE}"
[[ "$run_FIACH" = true ]] && echo "  FIACH     : ${FIACH_MODE} mode, ${FIACH_GMM} GMM, ${FIACH_MASK} mask"
[[ -n "$TARGET_SUB" ]] && echo "  Filter SUB: ${TARGET_SUB}"
[[ -n "$TARGET_CONDITION" ]] && echo "  Filter CONDITION: ${TARGET_CONDITION}"
[[ -n "$TARGET_RUN" ]] && echo "  Filter RUN: ${TARGET_RUN}"
[[ "$RESTART_FROM_STEP" -gt 0 ]] && echo "  Restart   : from step ${RESTART_FROM_STEP}"
echo "============================================"

# --- BIDS dataset_description.json for derivatives ---
# Rewritten on every launch so PipelineSettings reflects the flags at the top
# of this script for the most recent run. The per-run record of settings is
# derivatives/stats/<SUB>/pipeline_vars_<LABEL>[_TAG].env (written below).
_generated_on=$(date -u +%Y-%m-%dT%H:%M:%SZ)
for _pipeline in preproc stats work; do
    _root="${DIR}derivatives/${_pipeline}"
    mkdir -p "${_root}"
    cat > "${_root}/dataset_description.json" << JSON
{
  "Name": "${_pipeline}",
  "BIDSVersion": "1.8.0",
  "DatasetType": "derivative",
  "GeneratedBy": [{"Name": "wrapper_preprocess", "Version": "1.2.0"}],
  "PipelineSettings": {
    "GeneratedOn": "${_generated_on}",
    "PIPELINE_TAG": "${PIPELINE_TAG}",
    "DUMMY_SCANS": ${DUMMY_SCANS},
    "run_NORDIC": ${run_NORDIC},
    "NORDIC_PHASE": "${NORDIC_PHASE}",
    "run_FIACH": ${run_FIACH},
    "MNI_space": ${MNI_space},
    "SLICE_TIMING": ${SLICE_TIMING},
    "FWHM_MM": ${FWHM_MM},
    "FIACH_MODE": "${FIACH_MODE}",
    "FIACH_GMM": "${FIACH_GMM}",
    "FIACH_MASK": "${FIACH_MASK}",
    "RS_INCLUDE_MOTION": ${RS_INCLUDE_MOTION}
  }
}
JSON
    echo "Wrote ${_root}/dataset_description.json"
done


# ============================================================
# find_bold_files(func_dir)
#   Returns all *bold.nii.gz sorted:
#     1. Direct files in func/
#     2. Subfolders alphabetically (cap1 -> cap2 ...), runs sorted within each.
#   Empty folders are silently skipped.
#   Phase images (*part-phase*_bold.nii.gz) are not runs and are skipped.
# ============================================================
find_bold_files() {
    local func_dir="${1}"
    local -a results=()

    # Direct files in func/
    while IFS= read -r f; do
        [[ -n "$f" ]] && results+=("$f")
    done < <(find "${func_dir}" -maxdepth 1 -name "*bold.nii.gz" ! -name "*part-phase*" 2>/dev/null | sort)

    # Subfolders sorted (cap1, cap2, ...), runs sorted within each
    while IFS= read -r subdir; do
        [[ -z "$subdir" ]] && continue
        while IFS= read -r f; do
            [[ -n "$f" ]] && results+=("$f")
        done < <(find "${subdir}" -maxdepth 1 -name "*bold.nii.gz" ! -name "*part-phase*" 2>/dev/null | sort)
    done < <(find "${func_dir}" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)

    printf '%s\n' "${results[@]}"
}


# ============================================================
# derive_run_info(func_dir, bold_path)
#   Prints two lines:
#     Line 1 — CONDITION_NAME : subfolder name (e.g. "cap1"), or "" if directly in func/
#     Line 2 — RUN_NAME       : file base stripped of _bold suffix and any leading sub-XX_
#                               prefix, so sub-XX never appears twice in output filenames.
# ============================================================
derive_run_info() {
    local func_dir="${1}"
    local bold_path="${2}"
    local rel="${bold_path#${func_dir}/}"

    local condition_name="" run_file="" run_name=""

    if [[ "$rel" == */* ]]; then
        condition_name="${rel%%/*}"
        run_file="${rel#*/}"
    else
        run_file="$rel"
    fi

    run_name=$(basename "${run_file}" .nii.gz)
    run_name="${run_name%_bold}"
    run_name=$(echo "${run_name}" | sed 's/^sub-[^_]*_//')

    printf '%s\n%s\n' "${condition_name}" "${run_name}"
}


# ============================================================
# Checkpoint helpers
#   Stored in work/ outside derivatives — pipeline scratch only.
#   Keyed by RUN_LABEL so cap1 and cap2 never share checkpoints.
#   Files: work/SUB/checkpoints/RUN_LABEL/stepN.done
# ============================================================
_ckpt_dir() {
    echo "${DIR}derivatives/work/${SUB}/checkpoints/${RUN_LABEL}"
}

step_done() {
    [[ -f "$(_ckpt_dir)/step_${1}.done" ]]
}

mark_done() {
    local cdir
    cdir=$(_ckpt_dir)
    mkdir -p "${cdir}"
    touch "${cdir}/step_${1}.done"
    echo "    v Step ${1} done"
}

clear_from_step() {
    local from="${1}" cdir
    cdir=$(_ckpt_dir)
    mkdir -p "${cdir}"
    [[ "$from" -le 5 ]] && rm -f "${cdir}/step_5b.done"
    for i in $(seq "${from}" 12); do
        rm -f "${cdir}/step_${i}.done"
    done
    echo "  Checkpoints cleared from step ${from} onward (${SUB}/${RUN_LABEL})"
}

run_or_die() {
    local step="${1}" label="${2}"; shift 2
    "$@" && mark_done "${step}" \
         || { echo "  ERROR: Step ${step} (${label}) failed."; \
              echo "         Fix the issue then rerun with RESTART_FROM_STEP=${step}"; \
              exit 1; }
}


# ============================================================
# SUBJECT LOOP
# ============================================================
for SUB in "${SUBS[@]}"; do
    export SUB

    # --- Subject filter ---
    if [[ -n "$TARGET_SUB" && "$SUB" != "$TARGET_SUB" ]]; then
        echo "Skipping ${SUB} (TARGET_SUB=${TARGET_SUB})"
        continue
    fi

    echo ""
    echo "============================================"
    echo "Subject: ${SUB}"
    echo "============================================"

    FUNC_DIR="${DIR}${SUB}/func"
    mapfile -t BOLD_FILES < <(find_bold_files "${FUNC_DIR}")

    if [[ ${#BOLD_FILES[@]} -eq 0 ]]; then
        echo "WARNING: No *bold.nii.gz found under ${FUNC_DIR} — skipping"
        continue
    fi

    echo "Found ${#BOLD_FILES[@]} bold file(s):"
    printf '  %s\n' "${BOLD_FILES[@]}"

    # RUN LOOP
    for BOLD_SRC in "${BOLD_FILES[@]}"; do

        # --- Derive CONDITION_NAME and RUN_NAME separately ---
        mapfile -t _run_info < <(derive_run_info "${FUNC_DIR}" "${BOLD_SRC}")
        CONDITION_NAME="${_run_info[0]}"    # "cap1" or ""
        RUN_NAME="${_run_info[1]}"          # "run1", "task-face", etc. — never contains sub-XX

        # RUN_LABEL: "cap1_run1" when in a subfolder, or just "run1" when direct
        if [[ -n "$CONDITION_NAME" ]]; then
            RUN_LABEL="${CONDITION_NAME}_${RUN_NAME}"
        else
            RUN_LABEL="${RUN_NAME}"
        fi
        # OUT_LABEL: RUN_LABEL plus the variant tag — used for the env file and QC report
        OUT_LABEL="${RUN_LABEL}${PIPELINE_TAG:+_${PIPELINE_TAG}}"
        export BOLD_SRC CONDITION_NAME RUN_NAME RUN_LABEL OUT_LABEL

        # --- Condition filter ---
        if [[ -n "$TARGET_CONDITION" && "$CONDITION_NAME" != "$TARGET_CONDITION" ]]; then
            echo "  Skipping ${RUN_LABEL} (TARGET_CONDITION=${TARGET_CONDITION})"
            continue
        fi

        # --- Run filter ---
        if [[ -n "$TARGET_RUN" && "$RUN_NAME" != "$TARGET_RUN" ]]; then
            echo "  Skipping ${RUN_LABEL} (TARGET_RUN=${TARGET_RUN})"
            continue
        fi

        echo ""
        echo "  ----------------------------------------"
        echo "  ${SUB} / ${OUT_LABEL}"
        [[ -n "$CONDITION_NAME" ]] && echo "  Condition: ${CONDITION_NAME}   Run: ${RUN_NAME}"
        echo "  Source: ${BOLD_SRC}"
        echo "  ----------------------------------------"

        # --- Clear checkpoints from restart point onward ---
        if [[ "$RESTART_FROM_STEP" -gt 0 ]]; then
            clear_from_step "${RESTART_FROM_STEP}"
        fi

        # -------------------------------------------------------
        # Derive output directories
        #
        #   WORK_DIR     — intermediates (mc, FIACH, NORDIC, etc.)
        #                  work/sub-XX/RUN_LABEL/  (shared by standard and tagged runs)
        #
        #   FUNC_OUT_DIR — final func outputs (reg, smooth, tsnr)
        #                  no tag  → func/button/
        #                  tag set → func/button_smooth6mm/   (no condition: func_smooth6mm/)
        #   BASE_FUNC_OUT_DIR — always the untagged equivalent; used to seed tagged runs
        #
        #   T1_DIR       — final T1/anat outputs
        #                  derivatives/preproc/sub-XX/anat/
        #
        #   T1_WORK_DIR  — T1 intermediates (SynthSeg, UNI_T1.nii)
        #                  work/sub-XX/anat/
        # -------------------------------------------------------
        WORK_DIR="${DIR}derivatives/work/${SUB}/${RUN_LABEL}"

        if [[ -n "$CONDITION_NAME" ]]; then
            BASE_FUNC_OUT_DIR="${DIR}derivatives/preproc/${SUB}/func/${CONDITION_NAME}"
        else
            BASE_FUNC_OUT_DIR="${DIR}derivatives/preproc/${SUB}/func"
        fi

        if [[ -n "$PIPELINE_TAG" && -n "$CONDITION_NAME" ]]; then
            FUNC_OUT_DIR="${DIR}derivatives/preproc/${SUB}/func/${CONDITION_NAME}_${PIPELINE_TAG}"
        elif [[ -n "$PIPELINE_TAG" ]]; then
            FUNC_OUT_DIR="${DIR}derivatives/preproc/${SUB}/func_${PIPELINE_TAG}"
        else
            FUNC_OUT_DIR="${BASE_FUNC_OUT_DIR}"
        fi

        T1_DIR="${DIR}derivatives/preproc/${SUB}/anat"
        T1_WORK_DIR="${DIR}derivatives/work/${SUB}/anat"

        export WORK_DIR FUNC_OUT_DIR BASE_FUNC_OUT_DIR T1_DIR T1_WORK_DIR

        # Create directories
        mkdir -p "${WORK_DIR}/mc" \
                 "${WORK_DIR}/FIACH" \
                 "${WORK_DIR}/reg" \
                 "${FUNC_OUT_DIR}" \
                 "${T1_DIR}" \
                 "${T1_WORK_DIR}"
        [ "$run_NORDIC" = true ]   && mkdir -p "${WORK_DIR}/NORDIC"
        [ "$SLICE_TIMING" = true ] && mkdir -p "${WORK_DIR}/slicecor"

        echo "  Work dir    : ${WORK_DIR}"
        echo "  Func out dir: ${FUNC_OUT_DIR}"

        # --- Seed tagged run from base when PIPELINE_TAG is set ---
        # When rerunning from a mid-pipeline step with a tag, earlier outputs
        # (reg, norm) won't exist in the new tagged folder. Copy them over so
        # downstream steps don't fail looking for their inputs.
        # Also stamps checkpoints for all steps before RESTART_FROM_STEP so the
        # pipeline doesn't re-run them — this handles the case where no prior
        # tagged checkpoints exist (fresh tag, first time).
        if [[ -n "$PIPELINE_TAG" ]]; then

            # 1. Copy FUNC_OUT_DIR files from base (reg, norm, tsnr, qc)
            if [[ -d "$BASE_FUNC_OUT_DIR" ]]; then
                echo "  Seeding tagged func dir from base..."
                _seeded=0
                for _f in "${BASE_FUNC_OUT_DIR}"/*.nii.gz "${BASE_FUNC_OUT_DIR}"/*.nii "${BASE_FUNC_OUT_DIR}"/*.html; do
                    [[ -f "$_f" ]] || continue
                    _dest="${FUNC_OUT_DIR}/$(basename "$_f")"
                    if [[ ! -f "$_dest" ]]; then
                        cp "$_f" "$_dest"
                        echo "    Copied: $(basename "$_f")"
                        (( _seeded++ ))
                    fi
                done
                if [[ -d "${BASE_FUNC_OUT_DIR}/tsnr" && ! -d "${FUNC_OUT_DIR}/tsnr" ]]; then
                    cp -r "${BASE_FUNC_OUT_DIR}/tsnr" "${FUNC_OUT_DIR}/tsnr"
                    echo "    Copied: tsnr/"
                    (( _seeded++ ))
                fi
                [[ "$_seeded" -eq 0 ]] && echo "    FUNC_OUT_DIR already seeded"
            else
                echo "  WARNING: BASE_FUNC_OUT_DIR not found — cannot seed: ${BASE_FUNC_OUT_DIR}"
            fi

            # 2. Stamp checkpoints for all steps before RESTART_FROM_STEP
            # so the pipeline skips them even if no prior checkpoints exist.
            if [[ "$RESTART_FROM_STEP" -gt 1 ]]; then
                echo "  Stamping checkpoints for steps 1 to $(( RESTART_FROM_STEP - 1 ))..."
                _cdir=$(_ckpt_dir)
                mkdir -p "${_cdir}"
                # Handle step 5b separately (not a plain integer)
                for _s in $(seq 1 "$(( RESTART_FROM_STEP - 1 ))"); do
                    if [[ ! -f "${_cdir}/step_${_s}.done" ]]; then
                        touch "${_cdir}/step_${_s}.done"
                        echo "    Stamped: step_${_s}"
                    fi
                done
                # Stamp 5b if restarting from step 6 or later
                if [[ "$RESTART_FROM_STEP" -ge 6 && ! -f "${_cdir}/step_5b.done" ]]; then
                    touch "${_cdir}/step_5b.done"
                    echo "    Stamped: step_5b"
                fi
            fi
        fi

        # --- Resolve fmap directory for this run ---
        # Mirrors func folder structure:
        #   func/cap1/run1 → fmap/cap1/   (condition-specific)
        #   func/cap1/run1 → fmap/         (fallback: flat fmap applies to all conditions)
        #   func/run1      → fmap/         (no condition subfolders)
        # If the expected fmap folder doesn't exist → no GRE fieldmap for this run
        FMAP_DIR=""
        FMAP_HAS_CONDITIONS=$(find "${DIR}${SUB}/fmap" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | head -n 1)

        if [[ -n "$CONDITION_NAME" ]]; then
            if [ -d "${DIR}${SUB}/fmap/${CONDITION_NAME}/" ]; then
                FMAP_DIR="${DIR}${SUB}/fmap/${CONDITION_NAME}/"
                echo "  fmap: fmap/${CONDITION_NAME}/"
            elif [ -z "${FMAP_HAS_CONDITIONS}" ] && [ -d "${DIR}${SUB}/fmap/" ]; then
                FMAP_DIR="${DIR}${SUB}/fmap/"
                echo "  fmap: fmap/ (shared, no condition subfolders in fmap)"
            else
                echo "  fmap: fmap/${CONDITION_NAME}/ not found — no GRE fieldmap for ${RUN_LABEL}"
            fi
        else
            if [ -d "${DIR}${SUB}/fmap/" ]; then
                FMAP_DIR="${DIR}${SUB}/fmap/"
                echo "  fmap: fmap/"
            else
                echo "  fmap: none found — no GRE fieldmap"
            fi
        fi
        export FMAP_DIR

        # --- Resolve reverse-PE EPI for this run (primary distortion correction) ---
        # Same lookup AS3b does: staged_manifest.tsv in DIR (or the file given by
        # REVPHASE_MANIFEST from the shell), column "revphase" by header name,
        # matched on SUB (col 1) and RUN_NAME (col 2). Reset every run so one
        # run's file is never reused for the next.
        _manifest="${REVPHASE_MANIFEST:-${DIR}staged_manifest.tsv}"
        REVPHASE=""
        if [[ -f "$_manifest" ]]; then
            REVPHASE=$(awk -F'\t' -v s="${SUB}" -v r="${RUN_NAME}" \
                       'NR==1 {for (i=1; i<=NF; i++) if ($i=="revphase") c=i; next}
                        $1==s && $2==r && c {print $c; exit}' "$_manifest")
        fi
        if [[ -n "$REVPHASE" && ! -f "$REVPHASE" ]]; then
            echo "  WARNING: manifest lists reverse-PE ${REVPHASE} but the file does not exist — ignoring it"
            REVPHASE=""
        fi

        # --- Choose distortion correction: topup → GRE fieldmap → none ---
        if [[ -n "$REVPHASE" ]]; then
            UNWARP_METHOD="topup"
            echo "  Unwarping: topup (reverse-PE: ${REVPHASE})"
        elif [[ -n "$FMAP_DIR" ]]; then
            UNWARP_METHOD="fugue"
            echo "  Unwarping: GRE fieldmap (${FMAP_DIR})"
        else
            UNWARP_METHOD="none"
            echo "  Unwarping: none (no reverse-PE in manifest, no fmap/)"
        fi
        export REVPHASE UNWARP_METHOD
        [[ "$UNWARP_METHOD" = "fugue" ]] && mkdir -p "${WORK_DIR}/fmap"

        # -------------------------------------------------------
        # Pre-compute all pipeline filenames upfront.
        # Deterministic from flags — set once, valid for all steps.
        # -------------------------------------------------------
        if [ "$run_NORDIC" = true ]; then
            if [[ -n "$CONDITION_NAME" ]]; then
                mc_file=NORDIC_Run_BOLD_${SUB}_${CONDITION_NAME}_${RUN_NAME}
            else
                mc_file=NORDIC_Run_BOLD_${SUB}_${RUN_NAME}
            fi
        else
            if [[ -n "$CONDITION_NAME" ]]; then
                mc_file=Run_BOLD_${SUB}_${CONDITION_NAME}_${RUN_NAME}
            else
                mc_file=Run_BOLD_${SUB}_${RUN_NAME}
            fi
        fi

        # AS3 (fugue) and AS3b (topup) both write ${mc_file}_mcf_unwarp
        if [ "$UNWARP_METHOD" != "none" ]; then
            sc_file=${mc_file}_mcf_unwarp
        else
            sc_file=${mc_file}_mcf
        fi

        if [ "$SLICE_TIMING" = true ]; then
            fiach_file=sc${sc_file}
        else
            fiach_file=${sc_file}
        fi

        if [ "$run_FIACH" = true ]; then
            reg_file=rclean_${fiach_file}
        else
            reg_file=${fiach_file}
        fi

        if [ "$MNI_space" = true ]; then
            smooth_file=norm_${reg_file}
        else
            smooth_file=reg_${reg_file}
        fi

        # Motion parameters (mcflirt) — used by FIACH (step 6) and AS7 (step 7)
        MOTION_PAR="${WORK_DIR}/mc/${mc_file}_mcf.par"

        export mc_file sc_file fiach_file reg_file smooth_file MOTION_PAR

        # --- Resolve JSON sidecar early so all steps can use it ---
        BOLD_JSON="${BOLD_SRC%.nii.gz}.json"
        if [[ ! -f "$BOLD_JSON" ]]; then
            echo "  WARNING: No JSON sidecar found at ${BOLD_JSON} — scan params will use defaults"
            BOLD_JSON=""
        else
            echo "  JSON: ${BOLD_JSON}"
        fi
        export BOLD_JSON

        # --- Extract TR from JSON ---
        TR=""
        if [[ -n "$BOLD_JSON" && -f "$BOLD_JSON" ]]; then
            TR=$(jq -r '.RepetitionTime // empty' "${BOLD_JSON}" 2>/dev/null)
        fi
        if [[ -z "$TR" ]]; then
            echo "  WARNING: TR not found in JSON — will be empty in pipeline vars"
        else
            echo "  TR: ${TR}s"
        fi
        export TR DUMMY_SCANS

        # --- Export pipeline variables at start of run ---
        # ENV filename includes tag if set so variant runs don't overwrite each other
        ENV_FILE="${DIR}derivatives/stats/${SUB}/pipeline_vars_${OUT_LABEL}.env"
        mkdir -p "${DIR}derivatives/stats/${SUB}"

        # --- Overwrite guard (no tag only) ---
        # If smooth output already exists and no tag is set, warn and skip rather than clobber.
        # Set RESTART_FROM_STEP or use PIPELINE_TAG to force a rerun.
        if [[ -z "$PIPELINE_TAG" && "$RESTART_FROM_STEP" -eq 0 ]]; then
            _existing_smooth="${FUNC_OUT_DIR}/${smooth_file}_smooth.nii.gz"
            if [[ -f "$_existing_smooth" ]]; then
                echo "  WARNING: Smooth output already exists — skipping to avoid overwrite:"
                echo "           ${_existing_smooth}"
                echo "           Set PIPELINE_TAG or RESTART_FROM_STEP=10 to rerun smoothing."
                continue
            fi
        fi

        {
            echo "SUB=${SUB}"
            echo "DIR=${DIR}"
            echo "CONDITION_NAME=${CONDITION_NAME}"
            echo "RUN_NAME=${RUN_NAME}"
            echo "RUN_LABEL=${RUN_LABEL}"
            echo "PIPELINE_TAG=${PIPELINE_TAG}"
            echo "OUT_LABEL=${OUT_LABEL}"
            echo "BOLD_SRC=${BOLD_SRC}"
            echo "FMAP_DIR=${FMAP_DIR}"
            echo "REVPHASE=${REVPHASE}"
            echo "UNWARP_METHOD=${UNWARP_METHOD}"
            echo "WORK_DIR=${WORK_DIR}"
            echo "FUNC_OUT_DIR=${FUNC_OUT_DIR}"
            echo "T1_DIR=${T1_DIR}"
            echo "T1_WORK_DIR=${T1_WORK_DIR}"
            echo "mc_file=${mc_file}"
            echo "sc_file=${sc_file}"
            echo "fiach_file=${fiach_file}"
            echo "reg_file=${reg_file}"
            echo "smooth_file=${smooth_file}"
            echo "MOTION_PAR=${MOTION_PAR}"
            echo "TR=${TR}"
            echo "DUMMY_SCANS=${DUMMY_SCANS}"
            echo "run_NORDIC=${run_NORDIC}"
            echo "NORDIC_PHASE=${NORDIC_PHASE}"
            echo "MNI_space=${MNI_space}"
            echo "SLICE_TIMING=${SLICE_TIMING}"
            echo "FWHM_MM=${FWHM_MM}"
            echo "run_FIACH=${run_FIACH}"
            echo "FIACH_MODE=${FIACH_MODE}"
            echo "FIACH_GMM=${FIACH_GMM}"
            echo "FIACH_MASK=${FIACH_MASK}"
            echo "RS_INCLUDE_MOTION=${RS_INCLUDE_MOTION}"
        } > "${ENV_FILE}"
        echo "  Vars written: ${ENV_FILE}"

        # -------------------------------------------------------
        # STEP 1: Dummy scan removal + copy to NORDIC folder; apply NORDIC if enabled
        #   DUMMY_SCANS > 0 → first N volumes discarded before NORDIC, so every
        #   later step (motion, FIACH regressors, outputs) uses the trimmed run.
        #   NORDIC_PHASE="none" → magnitude-only.
        #   NORDIC_PHASE="auto" → phase image next to the BIDS bold
        #   (part-phase / _ph / _phase) if one exists, else magnitude-only.
        #   What was actually used is appended to this run's .env file.
        #   temporal_phase=1, phase_filter_width=10 are the fMRI settings
        #   from the AS1_run_NORDIC.m header.
        # -------------------------------------------------------
        if step_done 1; then
            echo "  Step 1: NORDIC/copy [SKIP]"
        else
            NORDIC_DIR="${WORK_DIR}/NORDIC"
            BOLD_COPY="${NORDIC_DIR}/Run_BOLD_${SUB}_${CONDITION_NAME:+${CONDITION_NAME}_}${RUN_NAME}.nii"
            echo "  Step 1: Copying bold to NORDIC folder..."

            if [[ "$DUMMY_SCANS" -gt 0 ]]; then
                echo "  Step 1: Removing ${DUMMY_SCANS} dummy scan(s)..."
                N_VOLS=$(fslnvols "${BOLD_SRC}")
                KEEP_VOLS=$(( N_VOLS - DUMMY_SCANS ))
                if [[ "$KEEP_VOLS" -le 0 ]]; then
                    echo "  ERROR: Step 1 — DUMMY_SCANS (${DUMMY_SCANS}) >= total volumes (${N_VOLS})"
                    exit 1
                fi
                echo "  Step 1: ${N_VOLS} total volumes, keeping ${KEEP_VOLS} (discarding first ${DUMMY_SCANS})"
                TRIMMED="${NORDIC_DIR}/trimmed_${RUN_NAME}.nii.gz"
                fslroi "${BOLD_SRC}" "${TRIMMED}" "${DUMMY_SCANS}" "${KEEP_VOLS}" \
                    || { echo "  ERROR: Step 1 — fslroi failed"; exit 1; }
                gunzip -c "${TRIMMED}" > "${BOLD_COPY}" \
                    && rm -f "${TRIMMED}" \
                    || { echo "  ERROR: Step 1 — failed to copy trimmed bold"; exit 1; }
            else
                gunzip -c "${BOLD_SRC}" > "${BOLD_COPY}" \
                    || { echo "  ERROR: Step 1 — failed to copy bold"; exit 1; }
            fi

            if [ "$run_NORDIC" = true ]; then
                PHASE_SRC=""
                PHASE_COPY=""
                if [ "$NORDIC_PHASE" = "auto" ]; then
                    for cand in "${BOLD_SRC/part-mag/part-phase}" \
                                "${BOLD_SRC%.nii.gz}_ph.nii.gz" \
                                "${BOLD_SRC%_bold.nii.gz}_phase.nii.gz"; do
                        if [ "$cand" != "$BOLD_SRC" ] && [ -f "$cand" ]; then
                            PHASE_SRC="$cand"
                            PHASE_COPY="${BOLD_COPY%.nii}_phase.nii"
                            # Phase gets the same dummy-scan trim as the magnitude so the two still match
                            if [[ "$DUMMY_SCANS" -gt 0 ]]; then
                                PH_TRIMMED="${NORDIC_DIR}/trimmed_phase_${RUN_NAME}.nii.gz"
                                fslroi "$PHASE_SRC" "$PH_TRIMMED" "${DUMMY_SCANS}" "${KEEP_VOLS}" \
                                    && gunzip -c "$PH_TRIMMED" > "$PHASE_COPY" \
                                    && rm -f "$PH_TRIMMED" \
                                    || { echo "  ERROR: Step 1 — failed to trim/copy phase ${PHASE_SRC}"; exit 1; }
                            else
                                gunzip -c "$PHASE_SRC" > "$PHASE_COPY" \
                                    || { echo "  ERROR: Step 1 — failed to copy phase ${PHASE_SRC}"; exit 1; }
                            fi
                            break
                        fi
                    done
                    [[ -z "$PHASE_SRC" ]] && \
                        echo "  WARNING: NORDIC_PHASE=auto but no phase image next to ${BOLD_SRC} — magnitude-only"
                fi

                if [[ -n "$PHASE_SRC" ]]; then
                    echo "  Step 1: NORDIC complex (magnitude + phase: ${PHASE_SRC})"
                    echo "NORDIC_INPUT=complex" >> "${ENV_FILE}"
                else
                    echo "  Step 1: NORDIC magnitude-only"
                    echo "NORDIC_INPUT=magnitude" >> "${ENV_FILE}"
                fi
                echo "NORDIC_PHASE_SRC=${PHASE_SRC}" >> "${ENV_FILE}"

                echo "  Step 1: Running NORDIC (AS1_run_NORDIC)..."
                # Statements go straight to matlab -batch (no temp .m file), so run names
                # with hyphens such as task-face are passed through unchanged.
                # MATLAB calls a function file by its file name, so AS1_run_NORDIC.m is
                # called as AS1_run_NORDIC(...) even though its first line says NIFTI_NORDIC_local.
                # Start a thread pool first: workers share memory. Without it, magnitude-only
                # NORDIC starts the default 'Processes' pool, which copies the full 4-D run to
                # every worker and gets OOM-killed on long runs.
                run_or_die 1 "NORDIC" \
                    matlab -nosplash -batch \
                        "addpath('${sDIR}'); if isempty(gcp('nocreate')), parpool('Threads'); end; ARG = struct('DIROUT', '${NORDIC_DIR}/', 'temporal_phase', 1, 'phase_filter_width', 10); AS1_run_NORDIC('${BOLD_COPY}', '${PHASE_COPY}', '${mc_file}', ARG);"

                # AS1_run_NORDIC returns without writing anything when a run has < 6 volumes
                if [[ ! -f "${NORDIC_DIR}/${mc_file}.nii" ]]; then
                    echo "  ERROR: Step 1 — ${mc_file}.nii not produced by NORDIC. Check the MATLAB output above."
                    rm -f "$(_ckpt_dir)/step_1.done"
                    exit 1
                fi
            else
                echo "  Step 1: NORDIC [disabled by flag] — using raw bold"
                echo "NORDIC_INPUT=none" >> "${ENV_FILE}"
                mark_done 1
            fi
        fi

        # -------------------------------------------------------
        # STEP 2: Motion correction
        # -------------------------------------------------------
        if step_done 2; then
            echo "  Step 2: Motion correction [SKIP]"
        else
            echo "  Step 2: Motion correction..."
            run_or_die 2 "Motion correction" \
                bash ${sDIR}AS2_realign_mcflirt_fsl.sh
            # Verify output was actually created (AS2 can exit 0 on mcflirt failure)
            if [[ ! -f "${WORK_DIR}/mc/${mc_file}_mcf.nii.gz" ]]; then
                echo "  ERROR: Step 2 — mcf output missing despite exit 0. Check AS2 logs."
                echo "         Rerun with RESTART_FROM_STEP=2"
                rm -f "$(_ckpt_dir)/step_2.done"
                exit 1
            fi
        fi

        # -------------------------------------------------------
        # STEP 3: Distortion correction / mean functional
        #   topup  → AS3b_topup_unwarp.sh  (reverse-PE EPI, primary)
        #   fugue  → AS3_unwarping.sh      (GRE fieldmap, fallback)
        #   none   → mean functional only
        #   AS3 and AS3b both write mc/${mc_file}_mcf_unwarp.nii.gz plus the
        #   mean functional and uncompressed copies in FIACH/.
        # -------------------------------------------------------
        if step_done 3; then
            echo "  Step 3: Unwarping/mean func [SKIP]"
        else
            case "$UNWARP_METHOD" in
                topup)
                    echo "  Step 3: Unwarping with topup (reverse-PE: ${REVPHASE})..."
                    run_or_die 3 "Unwarping (topup)" \
                        bash ${sDIR}AS3b_topup_unwarp.sh
                    ;;
                fugue)
                    echo "  Step 3: Unwarping with GRE fieldmap (${FMAP_DIR})..."
                    run_or_die 3 "Unwarping (fugue)" \
                        bash ${sDIR}AS3_unwarping.sh
                    ;;
                none)
                    echo "  Step 3: No reverse-PE or fieldmap — computing mean functional..."
                    run_or_die 3 "Mean functional" bash -c "
                        fslmaths ${WORK_DIR}/mc/${mc_file}_mcf.nii.gz \
                                  -Tmean ${WORK_DIR}/mc/meanFunctional_${RUN_NAME}.nii.gz \
                        && gunzip -c ${WORK_DIR}/mc/meanFunctional_${RUN_NAME}.nii.gz \
                                 > ${WORK_DIR}/FIACH/meanFunctional_${RUN_NAME}.nii \
                        && gunzip -c ${WORK_DIR}/mc/${mc_file}_mcf.nii.gz \
                                 > ${WORK_DIR}/FIACH/${mc_file}_mcf.nii
                    "
                    ;;
            esac
            # Verify the file every later step is named after actually exists
            if [[ ! -f "${WORK_DIR}/mc/${sc_file}.nii.gz" ]]; then
                echo "  ERROR: Step 3 — ${sc_file}.nii.gz not produced (${UNWARP_METHOD})."
                echo "         Rerun with RESTART_FROM_STEP=3"
                rm -f "$(_ckpt_dir)/step_3.done"
                exit 1
            fi
        fi

        # -------------------------------------------------------
        # STEP 4: Slice timing correction
        # -------------------------------------------------------
        if step_done 4; then
            echo "  Step 4: Slice timing [SKIP]"
        else
            if [ "$SLICE_TIMING" = true ]; then
                # AS4 reads mc/${sc_file}.nii, but step 3 writes mc/*.nii.gz and puts the
                # uncompressed copy in FIACH/. Make sure the file AS4 expects exists.
                ST_INPUT="${WORK_DIR}/mc/${sc_file}.nii"
                if [[ ! -f "$ST_INPUT" ]]; then
                    if [[ -f "${WORK_DIR}/mc/${sc_file}.nii.gz" ]]; then
                        gunzip -c "${WORK_DIR}/mc/${sc_file}.nii.gz" > "$ST_INPUT"
                    elif [[ -f "${WORK_DIR}/FIACH/${sc_file}.nii" ]]; then
                        cp "${WORK_DIR}/FIACH/${sc_file}.nii" "$ST_INPUT"
                    else
                        echo "  ERROR: Step 4 — slice-timing input not found. Tried:"
                        echo "         ${WORK_DIR}/mc/${sc_file}.nii(.gz)"
                        echo "         ${WORK_DIR}/FIACH/${sc_file}.nii"
                        exit 1
                    fi
                fi
                echo "  Step 4: Slice timing correction..."
                run_or_die 4 "Slice timing" \
                    bash ${sDIR}AS4_slice_timing_correction.sh
            else
                echo "  Step 4: Slice timing [disabled by flag]"
                mark_done 4
            fi
        fi

        # -------------------------------------------------------
        # STEP 5: Brain mask + segmentation (SynthSeg) — subject-level once
        # -------------------------------------------------------
        T1_CKPT="${T1_WORK_DIR}/.synthseg.done"

        if [[ "$RESTART_FROM_STEP" -gt 0 && "$RESTART_FROM_STEP" -le 5 ]]; then
            rm -f "${T1_CKPT}"
            rm -f "$(_ckpt_dir)/step_5b.done"
        fi

        if [[ -f "$T1_CKPT" ]]; then
            echo "  Step 5: Brain mask/SynthSeg [SKIP — already done for ${SUB}]"
        else
            T1_SRC=$(ls "${DIR}${SUB}/anat/"*T1w.nii.gz 2>/dev/null | head -n 1)
            if [[ ! -f "$T1_SRC" ]]; then
                echo "  ERROR: No T1w.nii.gz found for ${SUB} — cannot continue"
                exit 1
            fi

            SYNTHSEG_DIR="${T1_WORK_DIR}/synthseg"
            mkdir -p "${SYNTHSEG_DIR}"
            export SYNTHSEG_DIR

            echo "  Step 5: Copying T1..."
            gunzip -c "${T1_SRC}" > "${T1_WORK_DIR}/UNI_T1.nii"

            echo "  Step 5: Running SynthSeg + T1 brain mask..."
            run_or_die 5 "SynthSeg + T1 mask" bash -c "
                mri_synthseg \
                    --i '${T1_WORK_DIR}/UNI_T1.nii' \
                    --o '${SYNTHSEG_DIR}/synthseg.nii.gz' \
                    --robust \
                    --vol '${SYNTHSEG_DIR}/volumes.csv' \
                    --qc  '${SYNTHSEG_DIR}/qc.csv' \
                    --threads 8 \
                && python3 ${sDIR}AS5_T1mask.py \
                && touch '${T1_CKPT}'
            "
        fi
        # Always mark run-level step 5 done so pipeline continues
        step_done 5 || mark_done 5

        # -------------------------------------------------------
        # STEP 5b: Register T1 mask to this run's EPI space (run-level)
        # -------------------------------------------------------
        if [[ -f "$(_ckpt_dir)/step_5b.done" ]]; then
            echo "  Step 5b: Mask → EPI space [SKIP]"
        else
            echo "  Step 5b: Registering brain mask to EPI space (${RUN_LABEL})..."
            run_or_die 5b "Mask → EPI space" python3 ${sDIR}AS5b_mask_to_epi.py
        fi

        # -------------------------------------------------------
        # STEP 6: FIACH  (task: no regression | rs: motion + PCs regressed)
        # STEP 7: Combine motion + FIACH regressors
        # -------------------------------------------------------
        if [ "$run_FIACH" = true ]; then
            if step_done 6; then
                echo "  Step 6: FIACH [SKIP]"
            else
                echo "  Step 6: FIACH (${FIACH_MODE} mode)..."

                if [[ -z "$BOLD_JSON" ]]; then
                    echo "  ERROR: No JSON sidecar found — cannot run FIACH"
                    exit 1
                fi

                MEAN_SRC="${WORK_DIR}/mc/meanFunctional_${RUN_NAME}.nii.gz"
                if [[ ! -f "$MEAN_SRC" ]]; then
                    echo "  ERROR: Mean functional not found: ${MEAN_SRC}"
                    echo "         Check that Step 3 completed successfully."
                    exit 1
                fi

                if [[ ! -f "${MOTION_PAR}" ]]; then
                    echo "  ERROR: Motion parameters not found: ${MOTION_PAR}"
                    echo "         Check that Step 2 completed successfully."
                    exit 1
                fi

                run_or_die 6 "FIACH" bash ${sDIR}AS6_call_FIACH.sh

                RCLEAN="${WORK_DIR}/FIACH/rclean_${fiach_file}.nii"
                if [[ ! -f "$RCLEAN" ]]; then
                    echo "  ERROR: Step 6 — rclean_${fiach_file}.nii not produced."
                    echo "         Check MATLAB output above. Rerun with RESTART_FROM_STEP=6"
                    rm -f "$(_ckpt_dir)/step_6.done"
                    exit 1
                fi
            fi

            if step_done 7; then
                echo "  Step 7: Combine regressors [SKIP]"
            else
                echo "  Step 7: Combine motion + FIACH regressors..."
                [[ "$FIACH_MODE" = "rs" ]] && \
                    echo "    NOTE: rs mode — these regressors are already removed from rclean_; do not reuse them in a GLM"
                run_or_die 7 "Combine regressors" \
                    bash ${sDIR}AS7_combine_motion_regs_NORfi.sh
            fi
        else
            echo "  Steps 6-7: FIACH [disabled by flag]"
            step_done 6 || mark_done 6
            step_done 7 || mark_done 7
        fi

        # -------------------------------------------------------
        # STEP 8+9: Co-registration (EPI→T1) and Normalisation (EPI→MNI)
        #   Run as one script — single interpolation per output,
        #   both reg and norm saved. Checkpoints 8 and 9 both marked on success.
        # -------------------------------------------------------
        if step_done 8 && step_done 9; then
            echo "  Steps 8+9: Registration + Normalisation [SKIP]"
        else
            echo "  Steps 8+9: Registration + Normalisation..."
            run_or_die 8 "Registration+Normalisation" \
                python3 ${sDIR}AS8_coreg_and_norm.py
            # Mark step 9 done too since both outputs are produced together
            mark_done 9
        fi

        # -------------------------------------------------------
        # STEP 10: Smoothing
        # -------------------------------------------------------
        if step_done 10; then
            echo "  Step 10: Smoothing [SKIP]"
        else
            echo "  Step 10: Smoothing ${smooth_file} (FWHM=${FWHM_MM}mm)..."
            export FWHM_MM
            run_or_die 10 "Smoothing" \
                bash ${sDIR}AS10_smoothing.sh
        fi

        # -------------------------------------------------------
        # STEP 11: tSNR maps
        # -------------------------------------------------------
        if step_done 11; then
            echo "  Step 11: tSNR [SKIP]"
        else
            echo "  Step 11: tSNR maps..."
            run_or_die 11 "tSNR" \
                bash ${sDIR}AS11_tsnr_maps.sh
        fi

        # -------------------------------------------------------
        # STEP 12: QC report
        # -------------------------------------------------------
        if step_done 12; then
            echo "  Step 12: QC report [SKIP]"
        else
            echo "  Step 12: Generating QC report..."
            run_or_die 12 "QC report" \
                python3 ${sDIR}AS12_preproc_QC.py ${SUB} ${OUT_LABEL}
            xdg-open "${FUNC_OUT_DIR}/${SUB}_${OUT_LABEL}_qc.html" &
        fi

        echo "  Completed: ${SUB}/${OUT_LABEL}  (vars saved to ${ENV_FILE})"

    done  # end run loop

    echo ""
    echo "============================================"
    echo "All runs complete for: ${SUB}"
    echo "============================================"

done  # end subject loop

echo ""
echo "============================================"
echo "Pipeline complete for all subjects"
echo "============================================"