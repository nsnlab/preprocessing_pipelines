#!/bin/bash
export FSLDIR=$HOME/fsl
source $FSLDIR/etc/fslconf/fsl.sh
export PATH=$FSLDIR/bin:$PATH

if [ -z "${SUB}" ]; then
    echo "ERROR: SUB not defined. Please export SUB from the wrapper."
    exit 1
fi
if [ -z "${WORK_DIR}" ]; then
    echo "ERROR: WORK_DIR not defined. Please export from the wrapper."
    exit 1
fi
if [ -z "${mc_file}" ]; then
    echo "ERROR: mc_file not defined. Please export from the wrapper."
    exit 1
fi

# NORDIC outputs .nii (uncompressed), not .nii.gz
input="${WORK_DIR}/NORDIC/${mc_file}.nii"
output="${WORK_DIR}/mc/${mc_file}_mcf"
refvol="middle"

if [ ! -f "${input}" ]; then
    echo "ERROR: NORDIC input not found: ${input}"
    exit 1
fi

# 1. Motion correction
echo 'Motion correction starting...'
mcflirt \
  -in ${input} \
  -out ${output} \
  -refvol ${refvol} \
  -mats \
  -plots \
  -rmsrel \
  -rmsabs \
  -spline_final \
  -verbose 1

# 2. Generate motion plots
echo 'Generating motion plots...'
fsl_tsplot \
  -i ${output}.par \
  -t 'MCFLIRT estimated rotations (radians)' \
  -u 1 --start=1 --finish=3 \
  -a x,y,z -w 640 -h 144 \
  -o ${output}_rot.png

fsl_tsplot \
  -i ${output}.par \
  -t 'MCFLIRT estimated translations (mm)' \
  -u 1 --start=4 --finish=6 \
  -a x,y,z -w 640 -h 144 \
  -o ${output}_trans.png

fsl_tsplot \
  -i ${output}_abs.rms,${output}_rel.rms \
  -t 'MCFLIRT estimated mean displacement (mm)' \
  -u 1 -w 640 -h 144 \
  -a absolute,relative \
  -o ${output}_disp.png

echo "Motion correction complete. Outputs in mc/ directory:"
echo "  ${output}.nii.gz        - motion corrected BOLD"
echo "  ${output}.par           - motion parameters (6 columns)"
echo "  ${output}_rel.rms       - framewise displacement"
echo "  rot.png, trans.png, disp.png   - QC plots"