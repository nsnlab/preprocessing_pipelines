import nibabel as nib
import numpy as np
import os
from scipy.ndimage import binary_fill_holes, binary_closing, binary_dilation
from nilearn.image import resample_to_img
from pathlib import Path

SUB          = os.environ.get('SUB')
T1_DIR      = os.environ.get('T1_DIR')       # derivatives/preproc/sub-XX/anat/ — final outputs
T1_WORK_DIR = os.environ.get('T1_WORK_DIR')  # derivatives/work/sub-XX/anat/    — intermediates

if not SUB:          raise ValueError("SUB not defined.")
if not T1_DIR:       raise ValueError("T1_DIR not defined.")
if not T1_WORK_DIR:  raise ValueError("T1_WORK_DIR not defined.")

# --- Inputs (intermediates in work/) ---
t1_path       = f'{T1_WORK_DIR}/UNI_T1.nii'
synthseg_path = f'{T1_WORK_DIR}/synthseg/synthseg.nii.gz'

# --- Outputs (final anat outputs in derivatives/) ---
brain_mask_path = f'{T1_DIR}/brain_mask.nii.gz'
wm_mask_path    = f'{T1_DIR}/WM_mask.nii.gz'
masked_uni_path = f'{T1_DIR}/Masked_UNI.nii'
mask_t1_path    = f'{T1_DIR}/rfBrainMask_T1space.nii'  # shared across all runs

for p in [synthseg_path, t1_path]:
    if not os.path.isfile(p):
        raise FileNotFoundError(f'Required file not found: {p}')

print('Resampling SynthSeg to T1 native space...')
t1_nib      = nib.load(t1_path)
seg_nib     = nib.load(synthseg_path)
seg_t1space = resample_to_img(seg_nib, t1_nib, interpolation='nearest')
seg_data    = seg_t1space.get_fdata()

# Brain mask — binarise segmentation
print('Generating brain mask...')
brain_mask_data = (seg_data >= 1).astype(np.float32)
nib.save(nib.Nifti1Image(brain_mask_data, t1_nib.affine, t1_nib.header), brain_mask_path)
print(f'Brain mask saved: {brain_mask_path}')

# Masked T1
print('Generating masked T1...')
t1_data   = t1_nib.get_fdata().astype(np.float32)
masked_t1 = t1_data * brain_mask_data
nib.save(nib.Nifti1Image(masked_t1, t1_nib.affine, t1_nib.header), masked_uni_path)
print(f'Masked T1 saved: {masked_uni_path}')

# WM mask
wm_mask = ((seg_data == 2) | (seg_data == 41)).astype(np.float32)
nib.save(nib.Nifti1Image(wm_mask, t1_nib.affine, t1_nib.header), wm_mask_path)
print(f'WM mask saved: {wm_mask_path}')

# Morphological expansion — conservative kernel
print('Building expanded brain mask in T1 space...')
mask_data = brain_mask_data.astype(bool)
se_close  = np.ones((3, 3, 3), dtype=bool)
se_dilate = np.ones((1, 1, 1), dtype=bool)
closed    = binary_closing(mask_data, structure=se_close)
dilated   = binary_dilation(closed,   structure=se_dilate)
filled    = binary_fill_holes(dilated).astype(np.float32)
nib.save(nib.Nifti1Image(filled, t1_nib.affine, t1_nib.header), mask_t1_path)
print(f'T1-space brain mask saved: {mask_t1_path}')