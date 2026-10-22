#!/bin/bash

if [ -z "${SUB}" ];      then echo "ERROR: SUB not defined.";      exit 1; fi
if [ -z "${WORK_DIR}" ]; then echo "ERROR: WORK_DIR not defined."; exit 1; fi

N_PCS=6

# Derive motion file from mc_file — always <mc_file>_mcf.par
if [ -n "${MOTION_PAR}" ] && [ -f "${MOTION_PAR}" ]; then
    MOT_FILE="${MOTION_PAR}"
elif [ -n "${mc_file}" ]; then
    MOT_FILE="${WORK_DIR}/mc/${mc_file}_mcf.par"
else
    echo "ERROR: Neither MOTION_PAR nor mc_file defined."
    exit 1
fi

PC_FILE="${WORK_DIR}/FIACH/noisyPCs_${fiach_file}.tsv"
OUT_FILE="${WORK_DIR}/FIACH/multi_reg_FIACH_6PCs_${SUB}_${RUN_LABEL}.txt"

echo "Combining motion + FIACH (${N_PCS}PC)..."
if [ ! -f "${MOT_FILE}" ]; then echo "ERROR: Motion file not found: ${MOT_FILE}"; exit 1; fi
if [ ! -f "${PC_FILE}" ];  then echo "ERROR: PC file not found: ${PC_FILE}";      exit 1; fi

MOT_ROWS=$(wc -l < "${MOT_FILE}")
PC_ROWS=$(wc -l  < "${PC_FILE}")
if [ "${MOT_ROWS}" -ne "${PC_ROWS}" ]; then
    echo "ERROR: Row mismatch (motion=${MOT_ROWS}, PCs=${PC_ROWS})"
    exit 1
fi

echo "  Rows: ${MOT_ROWS} | Motion cols: 6 | PCs: ${N_PCS}"
python3 -c "
import numpy as np
motion = np.loadtxt('${MOT_FILE}')
assert motion.shape[1] == 6, f'Expected 6 motion cols, got {motion.shape[1]}'
with open('${PC_FILE}') as f:
    lines = [l.strip() for l in f if l.strip()]
pcs = np.array([[float(x) for x in l.split('\t') if x.strip()] for l in lines])
assert pcs.shape[1] >= ${N_PCS}, f'Only {pcs.shape[1]} PCs, need ${N_PCS}'
pcs = pcs[:, :${N_PCS}]
assert motion.shape[0] == pcs.shape[0], f'Row mismatch: {motion.shape[0]} vs {pcs.shape[0]}'
print(f'  Motion: {motion.shape}  PCs: {pcs.shape}')
nuis = np.hstack([motion, pcs])
np.savetxt('${OUT_FILE}', nuis, delimiter='\t', fmt='%.6f')
print(f'  Saved: ${OUT_FILE} ({nuis.shape[0]} x {nuis.shape[1]})')
"

echo "Done — nuisance regressors saved for ${SUB}/${RUN_LABEL}"