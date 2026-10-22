import ants
import nibabel as nib
import numpy as np
import os

SUB      = os.environ.get('SUB')
RUN_NAME = os.environ.get('RUN_NAME')
WORK_DIR = os.environ.get('WORK_DIR')
T1_DIR   = os.environ.get('T1_DIR')

if not SUB:      raise ValueError("SUB not defined.")
if not RUN_NAME: raise ValueError("RUN_NAME not defined.")
if not WORK_DIR: raise ValueError("WORK_DIR not defined.")
if not T1_DIR:   raise ValueError("T1_DIR not defined.")

# --- Inputs ---
t1_path      = f'{T1_DIR}/Masked_UNI.nii'           # final anat output
mask_t1_path = f'{T1_DIR}/rfBrainMask_T1space.nii'  # final anat output
func_path    = f'{WORK_DIR}/mc/meanFunctional_{RUN_NAME}.nii.gz'

# --- Outputs (work — run-specific EPI-space masks) ---
mask_out  = f'{WORK_DIR}/FIACH/rfBrainMask_{RUN_NAME}.nii'
brain_out = f'{WORK_DIR}/FIACH/Masked_UNI_{RUN_NAME}_epispace.nii'

if not os.path.isfile(func_path):
    raise FileNotFoundError(
        f'Mean functional not found: {func_path}\n'
        f'Check that Step 3 completed successfully.')
for p in [t1_path, mask_t1_path]:
    if not os.path.isfile(p):
        raise FileNotFoundError(f'Required file not found: {p}')

print(f'Registering T1 to EPI space for {RUN_NAME} (BOLDRigid)...')
img_t1   = ants.image_read(t1_path)
mask_t1  = ants.image_read(mask_t1_path)
meanFunc = ants.image_read(func_path)

registration = ants.registration(
    fixed=meanFunc, moving=img_t1,
    type_of_transform='BOLDRigid',
    random_seed=100,
    verbose=True
)

mask_epi = ants.apply_transforms(
    fixed=meanFunc, moving=mask_t1,
    transformlist=registration['fwdtransforms'],
    interpolator='nearestNeighbor'
)

meanFunc_nib = nib.load(func_path)
mask_final   = (mask_epi.numpy() > 0.5).astype(np.float32)

# Clip to EPI signal boundary
epi_data   = meanFunc_nib.get_fdata()
epi_thresh = np.percentile(epi_data, 5)
epi_signal = (epi_data > epi_thresh).astype(np.float32)
mask_final = (mask_final * epi_signal).astype(np.float32)

warped_t1    = registration['warpedmovout'].numpy()
brain_masked = (warped_t1 * mask_final).astype(np.float32)

nib.save(nib.Nifti1Image(mask_final,   meanFunc_nib.affine, meanFunc_nib.header), mask_out)
nib.save(nib.Nifti1Image(brain_masked, meanFunc_nib.affine, meanFunc_nib.header), brain_out)

print(f'Unique mask values: {np.unique(mask_final)}')
print(f'Brain voxels: {int(mask_final.sum())}')
print(f'Mask saved:   {mask_out}')
print(f'Brain saved:  {brain_out}')