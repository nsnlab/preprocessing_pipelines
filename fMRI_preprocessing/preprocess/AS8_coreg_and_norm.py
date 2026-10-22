import ants
import nibabel as nib
import numpy as np
import os
import shutil
import subprocess
from pathlib import Path
from nitransforms.linear import Affine
from nitransforms.nonlinear import DenseFieldTransform
from nitransforms.manip import TransformChain
from nitransforms.resampling import apply as nt_apply
from multiprocessing import Pool
import time

N_CORES = 8

SUB          = os.environ.get('SUB')
RUN_NAME     = os.environ.get('RUN_NAME')
WORK_DIR     = os.environ.get('WORK_DIR')
FUNC_OUT_DIR = os.environ.get('FUNC_OUT_DIR')
T1_DIR       = os.environ.get('T1_DIR')
REG_FILE     = os.environ.get('reg_file')
FSLDIR       = os.environ.get('FSLDIR', os.path.expanduser('~/fsl'))
C3D          = os.environ.get('C3D', '/home/asa25/c3d-1.1.0-Linux-x86_64/bin/c3d_affine_tool')

if not SUB:          raise ValueError("SUB not defined.")
if not RUN_NAME:     raise ValueError("RUN_NAME not defined.")
if not WORK_DIR:     raise ValueError("WORK_DIR not defined.")
if not FUNC_OUT_DIR: raise ValueError("FUNC_OUT_DIR not defined.")
if not T1_DIR:       raise ValueError("T1_DIR not defined.")
if not REG_FILE:     raise ValueError("reg_file not defined.")

# --- Intermediate work paths ---
reg_dir     = f'{WORK_DIR}/reg'
os.makedirs(reg_dir, exist_ok=True)

func_path    = f'{WORK_DIR}/FIACH/{REG_FILE}.nii'
meanfunc_path = f'{WORK_DIR}/mc/meanFunctional_{RUN_NAME}.nii.gz'
init_mat     = f'{reg_dir}/init_{RUN_NAME}.mat'
bbr_mat      = f'{reg_dir}/bbr_{RUN_NAME}.mat'
bbr_mat_itk  = f'{reg_dir}/bbr_{RUN_NAME}.tfm'

# --- T1/anat paths (final outputs from AS5) ---
t1_path   = f'{T1_DIR}/Masked_UNI.nii'
wm_mask   = f'{T1_DIR}/WM_mask.nii.gz'       # saved flat in T1_DIR by AS5
mni_warp  = f'{T1_DIR}/t1_to_mni_warp.nii.gz'
mni_affine = f'{T1_DIR}/t1_to_mni_affine.mat'
mni_path  = f'{FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz'

# --- Final func outputs ---
reg_out  = f'{FUNC_OUT_DIR}/reg_{REG_FILE}.nii.gz'
norm_out = f'{FUNC_OUT_DIR}/norm_{REG_FILE}.nii.gz'

if not os.path.isfile(meanfunc_path):
    raise FileNotFoundError(f'Mean functional not found at:\n  {meanfunc_path}')
print(f'Using mean functional: {meanfunc_path}')


def run(cmd, label):
    print(f'  [{label}] {cmd}')
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.stdout: print(result.stdout)
    if result.stderr: print(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(f'{label} failed (exit {result.returncode})')


def process_volume(args):
    i, vol_data, func_affine, bbr_path, mni_affine_path, mni_warp_path, \
    t1_affine, t1_shape, t1_header, mni_affine_arr, mni_shape, mni_header = args

    bbr      = Affine.from_filename(bbr_path, fmt='itk')
    syn_aff  = Affine.from_filename(mni_affine_path, fmt='itk')
    syn_warp = DenseFieldTransform(mni_warp_path, is_deltas=True)
    chain    = TransformChain([bbr, syn_aff, syn_warp])

    vol_img  = nib.Nifti1Image(vol_data, func_affine)
    t1_ref   = nib.Nifti1Image(np.zeros(t1_shape,  dtype=np.float32), t1_affine, t1_header)
    mni_ref  = nib.Nifti1Image(np.zeros(mni_shape, dtype=np.float32), mni_affine_arr, mni_header)

    warped_t1  = nt_apply(bbr,   vol_img, reference=t1_ref).get_fdata(dtype=np.float32)
    warped_mni = nt_apply(chain, vol_img, reference=mni_ref).get_fdata(dtype=np.float32)
    return i, warped_t1, warped_mni


# Step 1: Check SynthSeg WM mask
if not os.path.isfile(wm_mask):
    raise FileNotFoundError(
        f'WM mask not found: {wm_mask}\n'
        f'Check that Step 5 (SynthSeg) completed successfully.')
print(f'Using SynthSeg WM mask: {wm_mask}')

# Step 2: Initial rigid registration (EPI -> T1, 6 DOF)
print('Initial rigid registration (EPI -> T1)...')
run(f'flirt -in {meanfunc_path} -ref {t1_path} '
    f'-out {reg_dir}/init_meanfunc_{RUN_NAME} -omat {init_mat} -dof 6',
    'FLIRT init')

# Step 3: BBR refinement using SynthSeg WM mask
print('BBR registration (EPI -> T1)...')
run(f'flirt -in {meanfunc_path} -ref {t1_path} '
    f'-out {reg_dir}/bbr_meanfunc_{RUN_NAME} -omat {bbr_mat} '
    f'-wmseg {wm_mask} -cost bbr -init {init_mat} -dof 6',
    'FLIRT BBR')

# Step 4: Convert FSL BBR mat to ITK format
if not os.path.isfile(bbr_mat_itk):
    print('Converting FSL mat to ITK format (c3d)...')
    run(f'{C3D} -ref {t1_path} -src {meanfunc_path} '
        f'{bbr_mat} -fsl2ras -oitk {bbr_mat_itk}',
        'FSL->ITK (c3d)')
    print(f'  ITK mat saved: {bbr_mat_itk}')
else:
    print(f'ITK mat already exists, reusing: {bbr_mat_itk}')

# Step 5: T1 -> MNI nonlinear (ANTs SyN), reused across all runs
img_t1 = ants.image_read(t1_path)
mni    = ants.image_read(mni_path)

if os.path.isfile(mni_warp):
    print(f'T1->MNI warp already exists, reusing: {mni_warp}')
else:
    print('Registering T1 to MNI (ANTs SyN)...')
    reg_t1_mni = ants.registration(
        fixed=mni,
        moving=img_t1,
        type_of_transform='SyN',
        verbose=True
    )
    shutil.copy(reg_t1_mni['fwdtransforms'][0], mni_warp)
    shutil.copy(reg_t1_mni['fwdtransforms'][1], mni_affine)
    print(f'T1->MNI warp saved: {mni_warp}')

# Step 6: Load images and transforms
print('\nLoading images and transforms...')
func_nib    = nib.load(func_path)
t1_img      = nib.load(t1_path)
mni_img     = nib.load(mni_path)
arr_4d      = func_nib.get_fdata(dtype=np.float32)
n_vols      = arr_4d.shape[-1]
func_header = func_nib.header.copy()
print(f'  func shape: {func_nib.shape}')
print(f'  TR: {func_header.get_zooms()[3]:.3f}s')

# Step 7: Resample volume by volume in parallel
args_list = [
    (i, arr_4d[..., i], func_nib.affine,
     bbr_mat_itk, mni_affine, mni_warp,
     t1_img.affine, t1_img.shape[:3], t1_img.header,
     mni_img.affine, mni_img.shape[:3], mni_img.header)
    for i in range(n_vols)
]

result_t1  = np.zeros((*t1_img.shape[:3],  n_vols), dtype=np.float32)
result_mni = np.zeros((*mni_img.shape[:3], n_vols), dtype=np.float32)

print(f'\nTransforming {n_vols} volumes across {N_CORES} cores...')
t_start = time.time()

with Pool(N_CORES) as pool:
    for idx, (i, warped_t1, warped_mni) in enumerate(
            pool.imap_unordered(process_volume, args_list)):
        result_t1[..., i]  = warped_t1
        result_mni[..., i] = warped_mni
        if idx % 20 == 0:
            elapsed = time.time() - t_start
            if idx > 0:
                rate      = elapsed / idx
                remaining = rate * (n_vols - idx)
                print(f'  {idx+1}/{n_vols} done — elapsed: {elapsed/60:.1f}m, ETA: {remaining/60:.1f}m')
            else:
                print(f'  {idx+1}/{n_vols} done')

total = time.time() - t_start
print(f'\nResampling complete — total time: {total/60:.1f}m ({total/n_vols:.1f}s per volume)')

# Step 8: Save outputs
print('\nSaving outputs...')

t1_header = t1_img.header.copy()
t1_header.set_data_dtype(np.float32)
nib.save(nib.Nifti1Image(result_t1, t1_img.affine, t1_header), reg_out)
print(f'  Saved T1-space: {reg_out}')

mni_header = mni_img.header.copy()
mni_header.set_data_dtype(np.float32)
nib.save(nib.Nifti1Image(result_mni, mni_img.affine, mni_header), norm_out)
print(f'  Saved MNI-space: {norm_out}')

print('\nAS8 complete.')