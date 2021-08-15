#!/bin/bash 
### Begin BSUB Options 
#BSUB -P fus123 
#BSUB -J gpu_test
#BSUB -W 00:10 
#BSUB -nnodes 1
#BSUB -q debug 
#BSUB -alloc_flags gpumps
### End BSUB Options and begin shell commands

cd $LS_SUBCWD
module swap xl pgi
module load fftw
module load forge
module load netcdf
module load hdf5

mkdir -p matrix
mkdir -p out
mkdir -p dump

#export PGI_ACC_TIME=1

#setenv OMP_NUM_THREADS 4
setenv OMP_NUM_THREADS 28
jsrun -n 6 -a 1 -c 7 -g 1 -r 6 -EOMP_NUM_THREADS=28 -brs hpcrun -o hpctest.m -e REALTIME -e gpu=nvidia -t ./gem_main >hpcrun.out 2> hpcrun.err &


