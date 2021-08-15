#!/bin/bash -l

#SBATCH -A m3502
#SBATCH -C dgx
#SBATCH --reservation=gpu_hackathon_dgx
#SBATCH -N 1
#SBATCH -G 2
#SBATCH -t 0:30:00
#SBATCH -q shared
module list
cd $SLURM_SUBMIT_DIR
mkdir -p matrix
mkdir -p out
mkdir -p dump
c="srun -n 2 --ntasks-per-node=2 -c 1 --cpu-bind=cores --gpus-per-task=1 --gpu-bind=map_gpu:0,1 ./gem_main > run.out 2> run.err"
#c="srun -n 8 --ntasks-per-node=8 -c 1 --cpu-bind=cores --gpus-per-task=1 --gpu-bind=map_gpu:0,1,2,3,4,5,6,7 ./gem_main > run.out 2> run.err"
echo $c
eval $c

wait
