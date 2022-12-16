#!/bin/bash 
#$ -cwd
#$ -l s_rt=300:00:00

#$ -pe smp 4
#$ -l nvidia_p100=1

export ORIG=$PWD
export SCR=$TMPDIR
export LD_LIBRARY_PATH=/usr/local/cuda-10.1/lib64
cp run.tpr $SCR
cd $SCR
export PATH=/usr/local/gromacs-2021.3-GPU-plumed/bin/:$PATH
 module load plumed


gmx mdrun -nt 4 -gpu_id $SGE_GPU -deffnm run

cp -R * $ORIG
rm *.sh.*
