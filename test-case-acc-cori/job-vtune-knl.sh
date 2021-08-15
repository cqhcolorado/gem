#!/bin/bash -l

#SBATCH -q debug
#SBATCH -N 2
#SBATCH -t 00:30:00
#SBATCH -J GEMvtune
#SBATCH -C knl
#SBATCH -A mp118
#SBATCH --perf=vtune

cd $SLURM_SUBMIT_DIR
mkdir -p matrix
mkdir -p out
mkdir -p dump
mkdir -p opt
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/project/projectdirs/mp118/jycheng/SOFTWARE/pspline

#export OMP_PROC_BIND=true
#export OMP_PLACES=threads
export OMP_NUM_THREADS=16
#export OMP_STACKSIZE=1G
srun -n 32 -c 16 --cpu_bind=cores vtune -collect hotspots -r opt -finalization-mode=none -trace-mpi --  ./gem > run.out 2> run.err

vtune -archive -r opt
