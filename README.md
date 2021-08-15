# gem
A. For Cori (8 MPIs, 8 GPUs, 1 MPI/GPU) A.1 Modules: module load cgpu module load nvhpc/21.7 module load openmpi module load cuda A.2 Source files: GEM_OpenACC_Cori GEM_OpenMP_Cori A.3 batch script files: test-case-acc-cori test-case-omp-cori

B. For Summit (6 MPIs, 6GPUs, 1 MPI/GPU) B.1 Modules: module load nvhpc/21.7 module load cuda B.2 Source files: GEM_OpenACC_Summit GEM_OpenMP_Summit B.3 batch script files: test-case-acc-summit test-case-omp-summit

Note that in the gem.in of 'test-case' files, kmx=8 for Cori while kmx=6 for Summit. Additioanlly, the OpenACC codes have the same reresults for CPU and GPU, while the OpenMP codes have different ones.
