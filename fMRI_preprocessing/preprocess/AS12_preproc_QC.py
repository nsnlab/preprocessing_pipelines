# ------------------------------
# To run interactively:
#   conda activate neuro
#   cd /data/Annie/PROJECTS/scripts
#   python AS12_preproc_QC.py
#
# Called automatically by wrapper:
#   python AS12_preproc_QC.py <SUB> <RUN_LABEL>
# ------------------------------

import sys
import os
import base64
from pathlib import Path
import numpy as np
import pandas as pd
from dotenv import dotenv_values
from nilearn import plotting, image

# ─── CONFIG ───────────────────────────────────────────────────────────────────
# DIR is exported by the wrapper; also readable from env files in interactive mode
DIR = os.environ.get('DIR', '').rstrip('/')
if DIR:
    DIR = DIR + '/'

if len(sys.argv) > 2:
    # Called non-interactively by wrapper
    SUB, RUN_LABEL = sys.argv[1], sys.argv[2]
    selected = [(SUB, RUN_LABEL, None)]
else:
    # Interactive mode — discover env files from stats derivative
    if not DIR:
        raise EnvironmentError("DIR not set. Export DIR or run via wrapper.")

    env_files = sorted(Path(f'{DIR}derivatives/stats').glob('sub-*/pipeline_vars_*.env'))
    if not env_files:
        print('No pipeline_vars_*.env files found.')
        sys.exit(1)

    combos = []
    for f in env_files:
        sub       = f.parts[-2]           # sub-XX
        run_label = f.stem.replace('pipeline_vars_', '')
        combos.append((sub, run_label, f))

    print('Available runs:')
    for i, (sub, run_label, _) in enumerate(combos):
        print(f'  [{i}] {sub}  {run_label}')

    choice = input('\nEnter index (or press Enter to run ALL): ').strip()
    selected = combos if choice == '' else [combos[int(choice)]]

# ─── HTML HELPERS ─────────────────────────────────────────────────────────────
def h1(sections, text):
    sections.append(f'<h1>{text}</h1>')

def h2(sections, text):
    sections.append(f'<h2>{text}</h2>')

def h3(sections, text):
    sections.append(f'<h3>{text}</h3>')

def skip(sections, text):
    sections.append(f'<p style="color:#c0392b;font-family:monospace">[SKIP] {text}</p>')

def preformat(sections, text):
    escaped = text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
    sections.append(f'<pre style="background:#f4f4f4;padding:10px;border-radius:4px">{escaped}</pre>')

def add_view(sections, view, title=''):
    if title:
        h3(sections, title)
    sections.append(view._repr_html_())

def add_png(sections, path, title='', width=700):
    if title:
        h3(sections, title)
    data = base64.b64encode(Path(path).read_bytes()).decode()
    sections.append(f'<img src="data:image/png;base64,{data}" width="{width}" style="display:block;margin:10px 0">')

def add_dataframe(sections, df, title=''):
    if title:
        h3(sections, title)
    sections.append(df.to_html(classes='qc-table', border=0))

# ─── PER-RUN REPORT FUNCTION ──────────────────────────────────────────────────
def run_qc(SUB, RUN_LABEL):
    print(f'\nRunning QC: {SUB} / {RUN_LABEL}')
    sections = []

    env_path = f'{DIR}derivatives/stats/{SUB}/pipeline_vars_{RUN_LABEL}.env'
    config   = dotenv_values(env_path)

    # Read directory variables from env file
    _dir     = config.get('DIR', DIR)
    WORK_DIR     = config.get('WORK_DIR')
    FUNC_OUT_DIR = config.get('FUNC_OUT_DIR')
    T1_DIR       = config.get('T1_DIR')
    T1_WORK_DIR  = config.get('T1_WORK_DIR')

    CONDITION_NAME = config.get('CONDITION_NAME', '')
    RUN_NAME = config.get('RUN_NAME') or (
        RUN_LABEL.split('_', 1)[-1] if '_' in RUN_LABEL else RUN_LABEL
    )

    # Fallbacks if env vars missing (e.g. old env files)
    if not WORK_DIR:
        WORK_DIR = f'{_dir}derivatives/work/{SUB}/{RUN_LABEL}'
        print(f'  WARNING: WORK_DIR not in env, derived as: {WORK_DIR}')
    if not FUNC_OUT_DIR:
        FUNC_OUT_DIR = f'{_dir}derivatives/preproc/{SUB}/func/{CONDITION_NAME}' \
                       if CONDITION_NAME else f'{_dir}derivatives/preproc/{SUB}/func'
        print(f'  WARNING: FUNC_OUT_DIR not in env, derived as: {FUNC_OUT_DIR}')
    if not T1_DIR:
        T1_DIR = f'{_dir}derivatives/preproc/{SUB}/anat'
        print(f'  WARNING: T1_DIR not in env, derived as: {T1_DIR}')
    if not T1_WORK_DIR:
        T1_WORK_DIR = f'{_dir}derivatives/work/{SUB}/anat'
        print(f'  WARNING: T1_WORK_DIR not in env, derived as: {T1_WORK_DIR}')

    mc_file     = config.get('mc_file')
    sc_file     = config.get('sc_file')
    fiach_file  = config.get('fiach_file')
    reg_file    = config.get('reg_file')
    smooth_file = config.get('smooth_file')
    FIACH_MODE  = config.get('FIACH_MODE', 'task')

    # Work-dir paths (intermediates)
    nordic   = f'{WORK_DIR}/NORDIC/{mc_file}.nii'
    mc       = f'{WORK_DIR}/mc/{mc_file}_mcf.nii.gz'
    mc_sc    = f'{WORK_DIR}/mc/{sc_file}.nii.gz'
    fiach    = f'{WORK_DIR}/FIACH/{reg_file}.nii'
    meanfunc = f'{WORK_DIR}/mc/meanFunctional_{RUN_NAME}.nii.gz'
    mask_epi = f'{WORK_DIR}/FIACH/rfBrainMask_{RUN_NAME}.nii'
    tsnr_dir = f'{FUNC_OUT_DIR}/tsnr'

    # FIACH outputs
    fiach_dir    = f'{WORK_DIR}/FIACH'
    fiach_diag   = f'{WORK_DIR}/FIACH_fiach_diagnostics'   # fiach.m writes <input folder>_fiach_diagnostics
    fiach_filt   = f'{fiach_dir}/filt_{fiach_file}.nii'
    changed_pct  = f'{fiach_dir}/FIACH_changed_pct.nii'
    fiach_mask   = f'{fiach_diag}/mask.nii'
    fiach_median = f'{fiach_diag}/median.nii'              # median signal, zero outside mask = masked brain
    fiach_rtsnr  = f'{fiach_diag}/rtsnr.nii'
    fiach_noise  = f'{fiach_diag}/noise_mask.nii'

    # Anat paths
    anat = f'{T1_DIR}/Masked_UNI.nii'
    UNI  = f'{T1_WORK_DIR}/UNI_T1.nii'

    # Final func output paths
    reg    = f'{FUNC_OUT_DIR}/reg_{reg_file}.nii.gz'
    norm   = f'{FUNC_OUT_DIR}/norm_{reg_file}.nii.gz'
    smooth = f'{FUNC_OUT_DIR}/{smooth_file}_smooth.nii.gz'

    # ── Header + path check ───────────────────────────────────────────────────
    h1(sections, f'Preprocessing QC — {SUB} / {RUN_LABEL}')

    lines = [f'RUN_NAME = {RUN_NAME}', f'FIACH_MODE = {FIACH_MODE}', '']
    for label, path in [('anat',       anat),
                        ('nordic',     nordic),
                        ('mc',         mc),
                        ('meanfunc',   meanfunc),
                        ('fiach_filt', fiach_filt),
                        ('fiach',      fiach),
                        ('reg',        reg),
                        ('norm',       norm),
                        ('smooth',     smooth),
                        ('mask_epi',   mask_epi)]:
        mark = '✓' if Path(path).exists() else '✗ MISSING'
        lines.append(f'  [{mark}] {label}: {path}')
    preformat(sections, '\n'.join(lines))


    # ── Realignment ───────────────────────────────────────────────────────────
    h2(sections, 'Realignment and unwarping')

    for label, fname in [('Translation',  f'{mc_file}_mcf_trans.png'),
                         ('Rotation',     f'{mc_file}_mcf_rot.png'),
                         ('Displacement', f'{mc_file}_mcf_disp.png')]:
        p = f'{WORK_DIR}/mc/{fname}'
        if Path(p).exists():
            add_png(sections, p, title=label, width=800)
        else:
            skip(sections, f'{label} plot not found: {p}')

    meanNORfunc    = None
    meanNORmcufunc = None

    if Path(nordic).exists():
        meanNORfunc = image.mean_img(nordic)
    else:
        skip(sections, f'NORDIC not found: {nordic}')

    if Path(mc_sc).exists():
        meanNORmcufunc = image.mean_img(mc_sc)
    else:
        skip(sections, f'mc/unwarp not found: {mc_sc}')

    if meanNORfunc and meanNORmcufunc:
        add_view(sections,
                 plotting.view_img(meanNORfunc, bg_img=False, black_bg=True,
                                   cmap='gray', symmetric_cmap=False,
                                   title='NORDIC'),
                 title='NORDIC')
        add_view(sections,
                 plotting.view_img(meanNORmcufunc, bg_img=False, black_bg=True,
                                   cmap='gray', symmetric_cmap=False,
                                   title='mc/unwarp'),
                 title='mc/unwarp')
        add_view(sections,
                 plotting.view_img(meanNORfunc, bg_img=meanNORmcufunc,
                                   title='NORDIC vs mc/unwarp',
                                   cmap='hot', opacity=0.4, symmetric_cmap=False),
                 title='NORDIC vs mc/unwarp')

    h2(sections, 'EPI overlaid on T1')
    if meanNORfunc and Path(anat).exists():
        add_view(sections,
                 plotting.view_img(meanNORfunc, bg_img=anat, black_bg=True,
                                   title='EPI over T1', opacity=0.2, symmetric_cmap=False),
                 title='EPI over T1')
    else:
        skip(sections, 'meanNORfunc or anat not available')

    if meanNORmcufunc and Path(anat).exists():
        add_view(sections,
                 plotting.view_img(meanNORmcufunc, bg_img=anat, black_bg=True,
                                   title='mc/unwarp over T1', opacity=0.2, symmetric_cmap=False),
                 title='mc/unwarp over T1')
    else:
        skip(sections, 'meanNORmcufunc or anat not available')

    # ── Brain masking ─────────────────────────────────────────────────────────
    h2(sections, 'Brain masking')

    if Path(anat).exists() and Path(UNI).exists():
        add_view(sections,
                 plotting.view_img(anat, bg_img=UNI, black_bg=True,
                                   title='T1 extraction',
                                   opacity=0.3, symmetric_cmap=False),
                 title='T1 extraction (Masked_UNI over UNI_T1)')
    else:
        skip(sections, 'T1 or UNI not found')

    if Path(mask_epi).exists() and Path(meanfunc).exists():
        add_view(sections,
                 plotting.view_img(mask_epi, bg_img=meanfunc, black_bg=True,
                                   title='rfBrainMask over meanFunctional',
                                   cmap='jet', opacity=0.3, symmetric_cmap=False),
                 title='rfBrainMask over meanFunctional')
    else:
        skip(sections, 'mask_epi or meanfunc not found')

    # ── FIACH ─────────────────────────────────────────────────────────────────
    FIACH_MODE     = config.get('FIACH_MODE', 'task')
    BASE_LABEL     = config.get('RUN_LABEL', RUN_LABEL)        # untagged label (AS7 file naming)
    fiach_dir      = f'{WORK_DIR}/FIACH'
    fiach_diag     = f'{WORK_DIR}/FIACH_fiach_diagnostics'      # fiach.m writes <input folder>_fiach_diagnostics
    changed_pct    = f'{fiach_dir}/FIACH_changed_pct.nii'
    fiach_mask     = f'{fiach_diag}/mask.nii'
    fiach_median   = f'{fiach_diag}/median.nii'                 # median signal, zero outside mask = masked brain
    fiach_rtsnr    = f'{fiach_diag}/rtsnr.nii'
    fiach_noise    = f'{fiach_diag}/noise_mask.nii'

    h2(sections, f'FIACH outputs ({FIACH_MODE} mode)')
    preformat(sections, f'FIACH folder: {fiach_dir}')

    if FIACH_MODE == 'rs':
        sections.append('<p>rs mode: motion + FIACH PCs have been regressed out of '
                        '<code>rclean_</code>. Do not include them again in a later model.</p>')
    else:
        sections.append('<p>task mode: nothing regressed. Motion + FIACH PCs '
                        '(confound file below) go into the task GLM.</p>')

    # Summary values written by AS6_fiach_pipeline.m
    summary_path = f'{fiach_dir}/FIACH_summary.txt'
    if Path(summary_path).exists():
        df_summary = pd.read_csv(summary_path, sep='\t', header=None,
                                 names=['Metric', 'Value']).set_index('Metric')
        add_dataframe(sections, df_summary, title='FIACH summary')
    else:
        skip(sections, f'FIACH summary not found: {summary_path}')

    # Changed voxels over the masked brain
    if all(Path(p).exists() for p in [changed_pct, fiach_mask, fiach_median]):
        pct_img  = image.load_img(changed_pct)
        pct_data = pct_img.get_fdata()
        nonzero  = pct_data[pct_data > 0]
        if nonzero.size:
            vmax = float(np.percentile(nonzero, 99))
            overlay_png = f'{fiach_dir}/FIACH_changed_overlay.png'
            display = plotting.plot_stat_map(pct_img, bg_img=fiach_median,
                                             display_mode='z', cut_coords=8,
                                             threshold=1e-6, vmax=vmax, cmap='hot',
                                             symmetric_cbar=False, black_bg=True,
                                             colorbar=True, draw_cross=False,
                                             title='% of frames changed per voxel')
            display.add_contours(fiach_mask, levels=[0.5], colors='cyan', linewidths=0.8)
            display.savefig(overlay_png, dpi=150)
            display.close()
            add_png(sections, overlay_png,
                    title='Changed voxels over masked brain (cyan = FIACH mask)', width=1000)
            add_view(sections,
                     plotting.view_img(pct_img, bg_img=fiach_median, black_bg=True,
                                       threshold=1e-6, vmax=vmax, cmap='hot',
                                       symmetric_cmap=False,
                                       title='% of frames changed per voxel'),
                     title='Changed voxels (interactive)')
        else:
            skip(sections, 'FIACH changed no samples — nothing to overlay')
    else:
        skip(sections, 'Changed-voxel map, FIACH mask or median image not found')

    for title, path, cmap in [('Robust tSNR', fiach_rtsnr, 'jet'),
                              ('Noise mask (PCA voxels)', fiach_noise, 'autumn')]:
        if Path(path).exists() and Path(fiach_median).exists():
            add_view(sections,
                     plotting.view_img(path, bg_img=fiach_median, black_bg=True,
                                       threshold=1e-6, cmap=cmap, symmetric_cmap=False,
                                       opacity=0.6, title=title),
                     title=title)
        else:
            skip(sections, f'{title}: {path}')

    # QC figures written by AS6_fiach_pipeline.m (+ FIACH's own GMM plot)
    for title, p in [('GMM noise segmentation',        f'{fiach_diag}/segmentation.png'),
                     ('Corrections per frame',         f'{fiach_dir}/FIACH_QC_corrections_per_frame.png'),
                     ('Most-corrected voxel 1',        f'{fiach_dir}/FIACH_QC_worstvoxel_1.png'),
                     ('Most-corrected voxel 2',        f'{fiach_dir}/FIACH_QC_worstvoxel_2.png'),
                     ('Most-corrected voxel 3',        f'{fiach_dir}/FIACH_QC_worstvoxel_3.png')]:
        if Path(p).exists():
            add_png(sections, p, title=title, width=700)
        else:
            skip(sections, f'{title}: {p}')

    h3(sections, 'FIACH log')
    log_path = f'{fiach_dir}/fiach_log_{RUN_NAME}.txt'
    if Path(log_path).exists():
        preformat(sections, Path(log_path).read_text())
    else:
        skip(sections, f'FIACH log not found: {log_path}')

    # Confound file from AS7 (named with the untagged RUN_LABEL)
    conf_path = f'{fiach_dir}/multi_reg_FIACH_6PCs_{SUB}_{BASE_LABEL}.txt'
    if Path(conf_path).exists():
        df_conf = pd.DataFrame(np.loadtxt(conf_path))
        n_rows, n_cols = df_conf.shape
        flag = '' if n_cols == 12 else '  ← expected 12 (6 motion + 6 FIACH PCs)'
        add_dataframe(sections, df_conf.head(15),
                      title=f'Confound file — {n_rows} rows × {n_cols} columns{flag} (first 15 rows)')
        preformat(sections, conf_path)
    else:
        skip(sections, f'Confound file not found: {conf_path}')

    # ── Co-registration ───────────────────────────────────────────────────────
    h2(sections, 'Co-registration (EPI to T1)')

    meancoregepi = None
    if Path(reg).exists() and Path(anat).exists():
        meancoregepi = image.mean_img(reg)
        add_view(sections,
                 plotting.view_img(meancoregepi, bg_img=anat, black_bg=True,
                                   title='EPI registered to T1',
                                   opacity=0.3, symmetric_cmap=False),
                 title='EPI registered to T1')
    else:
        skip(sections, 'reg or anat not found')

    # ── Normalisation ─────────────────────────────────────────────────────────
    h2(sections, 'Normalisation')

    if Path(norm).exists():
        meanmnifunc = image.mean_img(norm)
        add_view(sections,
                 plotting.view_img(meanmnifunc, bg_img=anat, black_bg=True,
                                   title='mniFunc over native',
                                   opacity=0.3, symmetric_cmap=False),
                 title='mniFunc over native T1')
        add_view(sections,
                 plotting.view_img(meanmnifunc, black_bg=True,
                                   title='mniFunc over MNI',
                                   opacity=0.3, symmetric_cmap=False),
                 title='mniFunc over MNI template')
    else:
        skip(sections, f'norm not found: {norm}')

    # ── Smoothing ─────────────────────────────────────────────────────────────
    h2(sections, 'Smoothing')

    if Path(smooth).exists():
        meansmooth         = image.mean_img(smooth, copy_header=True)
        meansmooth_clipped = image.math_img("np.clip(img, 0, None)", img=meansmooth)

        if meancoregepi:
            meancoregepi_clipped = image.math_img("np.clip(img, 0, None)", img=meancoregepi)
            add_view(sections,
                     plotting.view_img(meancoregepi_clipped, bg_img=False, black_bg=True,
                                       cmap='Greys_r', symmetric_cmap=False, vmin=0,
                                       title='EPI before smoothing'),
                     title='EPI before smoothing')

        add_view(sections,
                 plotting.view_img(meansmooth_clipped, bg_img=False, black_bg=True,
                                   cmap='Greys_r', symmetric_cmap=False, vmin=0,
                                   title='Smoothed EPI'),
                 title='Smoothed EPI')
    else:
        skip(sections, f'smooth not found: {smooth}')

    # ── tSNR ──────────────────────────────────────────────────────────────────
    h2(sections, 'tSNR maps')

    vmin, vmax = 0, 250

    tsnr_files = {
        'NORDIC (pre-mc)'       : f'{tsnr_dir}/{mc_file}_tsnr.nii.gz',
        'realigned/unwarped'    : f'{tsnr_dir}/{mc_file}_mcf_tsnr.nii.gz',
        'preregistered (FIACH)' : f'{tsnr_dir}/{reg_file}_tsnr.nii.gz',
        'registered (T1)'       : f'{tsnr_dir}/reg_{reg_file}_tsnr.nii.gz',
        'normalised (MNI)'      : f'{tsnr_dir}/norm_{reg_file}_tsnr.nii.gz',
        'smoothed'              : f'{tsnr_dir}/{smooth_file}_smooth_tsnr.nii.gz',
    }

    for title, path in tsnr_files.items():
        if not Path(path).exists():
            skip(sections, f'{title}: {path}')
            continue
        clipped = image.math_img("np.clip(img, 0, None)", img=path)
        add_view(sections,
                 plotting.view_img(clipped, bg_img=False, black_bg=True, cmap='jet',
                                   symmetric_cmap=False, vmin=vmin, vmax=vmax,
                                   colorbar=True, title=title),
                 title=title)

    # ── Write HTML ────────────────────────────────────────────────────────────
    output_path = Path(FUNC_OUT_DIR) / f'{SUB}_{RUN_LABEL}_qc.html'
    output_path.parent.mkdir(parents=True, exist_ok=True)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>QC — {SUB} {RUN_LABEL}</title>
  <style>
    body       {{ font-family: sans-serif; max-width: 1100px; margin: 40px auto; padding: 0 20px; background: #fff; }}
    h1         {{ border-bottom: 2px solid #333; padding-bottom: 8px; }}
    h2         {{ border-bottom: 1px solid #ccc; margin-top: 40px; color: #2c3e50; }}
    h3         {{ color: #555; margin-top: 24px; }}
    .qc-table  {{ border-collapse: collapse; font-size: 12px; }}
    .qc-table td, .qc-table th {{ border: 1px solid #ddd; padding: 4px 8px; }}
    .qc-table th {{ background: #f0f0f0; }}
  </style>
</head>
<body>
{''.join(sections)}
</body>
</html>"""

    output_path.write_text(html, encoding='utf-8')
    print(f'  Saved: {output_path}')


# ─── RUN ──────────────────────────────────────────────────────────────────────
for SUB, RUN_LABEL, _ in selected:
    run_qc(SUB, RUN_LABEL)

print('\nDone.')