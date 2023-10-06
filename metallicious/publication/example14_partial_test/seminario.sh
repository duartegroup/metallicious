#!/bin/bash
#$ -cwd
#$ -pe smp 8
#$ -l s_rt=360:00:00
export OMP_NUM_THREADS=1
module load mpi/openmpi3-x86_64
conda activate autode_env
/u/fd/chem1540/miniconda3/envs/autode_env/bin/python3.9 /u/fd/chem1540/github/cgbind2pmd/parametrize_new_sites.py -metal_charge 2 -metal_name Pd -f cage.xyz -keywords PBE0 D3BJ def2-SVP tightOPT freq
