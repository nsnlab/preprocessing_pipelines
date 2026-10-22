# 7T fMRI preprocessing pipeline

**Version:** 1.2.0  
**Author:** Anna Sadilova  
**Last updated:** 22.10.2026

Preprocessing for 7T BOLD fMRI in BIDS layout: NORDIC denoising, motion correction, distortion correction (topup or GRE fieldmap), optional slice timing, SynthSeg brain masking, FIACH noise correction, registration to T1/MNI, smoothing, tSNR maps and an HTML QC report per run.

Everything is driven by one script, `wrapper_function_preprocess.sh`, and configured by the settings at its top.

---

## 1. Requirements

### Software

| Tool | Used for | Notes |
|---|---|---|
| **bash** | wrapper and step scripts | Linux |
| **FSL** | mcflirt, fugue, topup/applytopup, fslmaths, fslroi | expected in `$HOME/fsl` |
| **FreeSurfer 8.1.0** | `mri_synthseg` (step 5); licence needs to be added in the same dir as the wrapper - change wrapper dirs before first use |
| **MATLAB** | NORDIC (step 1), FIACH (step 6) | run headless with `matlab -batch` |
| **SPM12** | NIfTI I/O for FIACH | must be on your saved MATLAB path, or set `SPM_DIR` |
| **Python 3 + venv** | steps 5, 5b, 8, 12 | `python3-venv` must be installed (`sudo apt install python3-venv` on Ubuntu) |
| **jq** | reading JSON sidecars | `sudo apt install jq` |

### MATLAB toolboxes

- **Parallel Computing Toolbox** (NORDIC thread pool)
- **Image Processing Toolbox** (`niftiread` / `niftiwrite` in NORDIC)
- **Signal Processing Toolbox** (`tukeywin` in NORDIC)
- **Statistics and Machine Learning Toolbox**, only if `FIACH_GMM="ModifiedICH"`

Check in MATLAB:
```matlab
license('test','Distrib_Computing_Toolbox')
license('test','Image_Toolbox')
license('test','Signal_Toolbox')
which spm          % must print a path to spm.m
```

### Python packages (virtual environment)

The wrapper manages its own virtual environment in `Scripts/.venv`; no manual install is needed.

- **First launch:** creates the venv and installs `1_requirements.txt`.
- **Later launches:** reuses it, and reinstalls only if `1_requirements.txt` has changed (a checksum is stored in `.venv/.requirements.sha256`).
- **All `python3` steps** (AS5, AS5b, AS8, AS12) run inside the venv.

The QC report needs at least `numpy`, `pandas`, `nilearn` and `python-dotenv`; the rest are listed in `1_requirements.txt`.

To rebuild the venv from scratch (e.g. after a Python upgrade):
```bash
rm -rf /path/to/Scripts/.venv
```
and launch the wrapper again.

---

## 2. Folder layout

### Scripts
```
Scripts/
├── wrapper_function_preprocess.sh
├── 1_requirements.txt
├── .licence (recommended location)
├── README.md
├── .venv/                             # created automatically on first launch
└── preprocess/
    ├── AS1_run_NORDIC.m               # NORDIC (function inside is NIFTI_NORDIC_local)
    ├── AS2_realign_mcflirt_fsl.sh
    ├── AS3_unwarping.sh               # GRE fieldmap (fugue)
    ├── AS3b_topup_unwarp.sh           # reverse-PE (topup)
    ├── AS4_slice_timing_correction.sh
    ├── AS5_T1mask.py
    ├── AS5b_mask_to_epi.py
    ├── AS6_call_FIACH.sh
    ├── AS6_fiach_pipeline.m
    ├── fiach.m                        # FIACH (from FIACH-main/src)
    ├── fiach_bold_contrast.m          # FIACH threshold (from FIACH-main/src)
    ├── AS7_combine_motion_regs_NORfi.sh
    ├── AS8_coreg_and_norm.py
    ├── AS10_smoothing.sh
    ├── AS11_tsnr_maps.sh
    └── AS12_preproc_QC.py
```

If the Scripts folder is under git, add `.venv/` to `.gitignore`.

### Input data (BIDS)
```
DIR/
├── staged_manifest.tsv                # optional: reverse-PE lookup for topup
└── sub-XX/
    ├── anat/sub-XX_*T1w.nii.gz
    ├── func/
    │   ├── sub-XX_<run>_bold.nii.gz + .json          # runs directly in func/, or
    │   └── <condition>/sub-XX_<run>_bold.nii.gz + .json   # one subfolder per condition
    └── fmap/                          # optional GRE fieldmap (or fmap/<condition>/)
        ├── *_e1.nii.gz  *_e1.json
        └── *_e2_ph.nii.gz  *_e2.json
```

The bold JSON sidecar must contain `RepetitionTime` and `EchoTime`. For fugue it also needs `EffectiveEchoSpacing`, and for topup `TotalReadoutTime` or `EffectiveEchoSpacing`.

---

## 3. Configure

Edit the top of `wrapper_function_preprocess.sh`.

### Main settings
| Setting | Meaning |
|---|---|
| `DIR` | BIDS root, **with trailing slash** |
| `SUBS` | subjects to process |
| `run_NORDIC` | `true` / `false` |
| `run_FIACH` | `true` / `false` |
| `MNI_space` | `true` = smooth the MNI-normalised output, `false` = smooth the T1-space output |
| `SLICE_TIMING` | `true` / `false` |
| `FWHM_MM` | smoothing kernel (mm) |
| `DUMMY_SCANS` | volumes to discard from the start of each run (0 = none) |
| `PIPELINE_TAG` | label for a variant run (see section 6) |
| `NORDIC_PHASE` | `"none"` = magnitude-only; `"auto"` = use a phase image next to the bold if one exists |

### Run selection and restarts
| Setting | Meaning |
|---|---|
| `TARGET_SUB` / `TARGET_CONDITION` / `TARGET_RUN` | process only matching runs (`""` = all) |
| `RESTART_FROM_STEP` | `0` = resume from checkpoints; `N` = clear checkpoints from step N and rerun from there |

### FIACH
| Setting | Meaning |
|---|---|
| `FIACH_MODE` | `"task"` = nothing regressed, motion + FIACH PCs go into your task GLM; `"rs"` = motion + PCs regressed out of the data |
| `FIACH_GMM` | `"R"` (default) or `"ModifiedICH"` |
| `FIACH_MASK` | `"user"` = brain mask from step 5b; `"auto"` = FIACH's own mask |
| `RS_INCLUDE_MOTION` | rs mode only: `true` = regress motion + PCs, `false` = PCs only |
| `TSNR_JOBS` | tSNR maps computed in parallel (lower it if memory runs out) |

---

## 4. Run

```bash
cd /path/to/Scripts
bash -n wrapper_function_preprocess.sh          # syntax check, no output = OK
./wrapper_function_preprocess.sh 2>&1 | tee run_$(date +%F_%H%M).log
```

### What the run command does

It runs the pipeline as usual **and** saves everything it prints to a log file.

| Part | Meaning |
|---|---|
| `./wrapper_function_preprocess.sh` | run the pipeline |
| `2>&1` | also capture error messages (they normally go to a separate stream and would be missed) |
| `\|` | pass the output on to the next command |
| `tee ...` | show the output in the terminal **and** write it to a file at the same time |
| `run_$(date +%F_%H%M).log` | the file name, stamped with the current date and time, e.g. `run_2026-10-22_1108.log` |

Why save the log:
- the terminal only keeps a limited amount of scrollback, so early messages from a long run get lost;
- NORDIC's and MATLAB's messages are only in this output;
- each run gets its own file, so earlier logs are never overwritten;
- if something fails, search the log, e.g. `grep -n ERROR run_2026-10-22_1108.log`.

Running plain `./wrapper_function_preprocess.sh` works exactly the same, just without the saved log.

At startup the log shows which Python is used (`Python: .../Scripts/.venv/bin/python3`) and whether the venv was created, updated or reused.

### Override paths for one run
Without editing the script:
```bash
SPM_DIR=/path/to/spm12 ./wrapper_function_preprocess.sh
FREESURFER_HOME=/path/to/freesurfer ./wrapper_function_preprocess.sh
REVPHASE_MANIFEST=/path/to/manifest.tsv ./wrapper_function_preprocess.sh
```

### Test on one run first
```bash
TARGET_SUB="sub-RT"
TARGET_RUN="test_run2"
RESTART_FROM_STEP=1
```

---

## 5. Steps

| Step | What it does | Script |
|---|---|---|
| 1 | copy bold, remove dummy scans, NORDIC | wrapper + `AS1_run_NORDIC.m` |
| 2 | motion correction (mcflirt) | `AS2` |
| 3 | distortion correction: **topup** if a reverse-PE is listed in `staged_manifest.tsv`, else **GRE fieldmap** if `fmap/` exists, else mean functional only | `AS3b` / `AS3` |
| 4 | slice timing (if enabled) | `AS4` |
| 5 | SynthSeg + T1 brain mask (once per subject) | `mri_synthseg` + `AS5` |
| 5b | brain mask → EPI space | `AS5b` |
| 6 | FIACH → `rclean_*.nii` | `AS6_call_FIACH.sh` → `AS6_fiach_pipeline.m` |
| 7 | combine motion + FIACH regressors (12 columns) | `AS7` |
| 8–9 | registration to T1 and MNI | `AS8` |
| 10 | smoothing | `AS10` |
| 11 | tSNR maps (parallel) | `AS11` |
| 12 | HTML QC report | `AS12` |

Each step writes a checkpoint, so a crashed run resumes where it stopped. The error message tells you which `RESTART_FROM_STEP` to use.
If RESTART_FROM_STEP=0, it will continue from last finished step.

---

## 6. Variant runs (`PIPELINE_TAG`)

Set e.g. `PIPELINE_TAG="smooth6mm"` and `RESTART_FROM_STEP=10`:

- Final outputs go to `func/<condition>_smooth6mm/` (no condition: `func_smooth6mm/`).
- The tagged folder is first filled with the untagged outputs, and steps before `RESTART_FROM_STEP` are marked done.
- The env file is saved as `pipeline_vars_<run>_smooth6mm.env`.

Tagged and untagged runs **share the work folder and checkpoints**. Use tags for late-step variants (smoothing, registration). A tag restarting from an early step overwrites the shared intermediates.

Without a tag, a run whose final smoothed file already exists is skipped when `RESTART_FROM_STEP=0` (overwrite guard).

---

## 7. Outputs

```
DIR/derivatives/
├── preproc/sub-XX/
│   ├── anat/                                   # Masked_UNI.nii etc.
│   └── func/<condition>[_TAG]/
│       ├── reg_*.nii.gz  norm_*.nii.gz  *_smooth.nii.gz
│       ├── tsnr/
│       └── sub-XX_<label>_qc.html              # QC report
├── work/sub-XX/<RUN_LABEL>/                     # NORDIC/ mc/ fmap/ topup/ FIACH/ reg/
│   └── FIACH/
│       ├── rclean_*.nii                         # input to registration
│       ├── noisyPCs_*.tsv
│       ├── multi_reg_FIACH_6PCs_sub-XX_<label>.txt   # confounds for the GLM (task mode)
│       ├── FIACH_summary.txt  FIACH_changed_*.nii  FIACH_QC_*.png
│       └── fiach_log_<run>.txt
├── work/sub-XX/checkpoints/<RUN_LABEL>/step_N.done
└── stats/sub-XX/pipeline_vars_<label>[_TAG].env   # every setting and path used for that run
```

`dataset_description.json` in each derivatives folder records the settings of the most recent launch.

### Regenerate a QC report only
Activate the venv first, so the right packages are used:
```bash
source /path/to/Scripts/.venv/bin/activate
DIR=/path/to/BIDS/ python3 /path/to/Scripts/preprocess/AS12_preproc_QC.py sub-XX <label>
deactivate
```

---

## 8. Things to check

- **Dummy scans:** if `DUMMY_SCANS > 0`, shift task onsets by `DUMMY_SCANS × TR`, unless they're already timed from the first kept volume.
- **Distortion correction:** in the QC "EPI over T1" view, frontal and temporal regions should move *towards* the anatomy. If they don't:
  - topup: swap `PE_DIR_FWD` / `PE_DIR_PA`;
  - fugue: check `--unwarpdir` against `PhaseEncodingDirection`.
- **topup:** a warning about a near-zero field means the forward and reverse runs aren't opposed.
- **FIACH task vs rs:** in rs mode the regressors are already removed. Don't add them to a later model.
- **FIACH QC:** a "QC flag count differs" warning in the FIACH log means the changed-voxel overlay shouldn't be trusted.

---

## 9. Troubleshooting

| Message | Cause / fix |
|---|---|
| `could not create venv` | install `python3-venv` (`sudo apt install python3-venv`) |
| `pip install failed` | check the package versions in `1_requirements.txt` and the network; then `rm -rf Scripts/.venv` and relaunch |
| `ModuleNotFoundError` in a Python step | the package is missing from `1_requirements.txt`. Add it there; the next launch installs it |
| `Killed` during NORDIC | out of memory. Already avoided with a thread pool; if it persists, cap workers (`parpool('Threads', 8)` in step 1) |
| `SPM not found` / `NoSPM` | add SPM to your saved MATLAB path (`addpath(...); savepath`) or set `SPM_DIR` |
| `Name is nonexistent or not a directory` (MATLAB) | stale entry in `pathdef.m`; harmless. Remove with `sed -i.bak '/<folder>/d' ~/matlab/pathdef.m` |
| `SetUpFreeSurfer.sh: No such file` | wrong `FREESURFER_HOME` |
| `... not produced` after a step | check that step's output above, then rerun with the `RESTART_FROM_STEP` it suggests |
| zygote / Broken pipe | from opening the QC report in a browser; harmless |