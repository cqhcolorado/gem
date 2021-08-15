#!/bin/bash -l

#SBATCH -q debug
#SBATCH -N 2
#SBATCH -t 00:30:00
#SBATCH -J test
#SBATCH -C knl
#SBATCH -A mp118

cd $SLURM_SUBMIT_DIR
mkdir -p matrix
mkdir -p out
mkdir -p dump
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/project/projectdirs/m499/Software/pspline/DEFAULT/cori_knl/intel/lib

#export OMP_PROC_BIND=true
#export OMP_PLACES=threads
export OMP_NUM_THREADS=16
#export OMP_STACKSIZE=1G
srun -n 32 -c 16 --cpu_bind=cores ./gem > run.out 2> run.err

wait
