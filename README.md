HPC Project: Hybrid Parallel Stencil Method
Overview
This project explores high-performance stencil computations, a core technique in scientific computing for solving partial differential equations (specifically, the 2D Heat Equation). The repository demonstrates a C-based implementation parallelized using a hybrid MPI + OpenMP model. The primary goal is to evaluate the performance and scalability of this approach through rigorous scaling studies on the Leonardo supercomputer.

Project Structure
HPCProject-StencilMethod/
├── include/              # Header files (.h) for the C implementation
├── src/                  # Main C source files (.c)
├── mpi_weak_scaling/     # Scripts and outputs for MPI weak scaling
├── mpi_scaling_results/  # Scripts and outputs for MPI strong scaling
├── omp_scaling_results/  # Scripts and results from OpenMP scaling
├── plot/                 # Scripts and notebooks for analysis and plotting
└── go_scripts/           # SLURM batch scripts for job submission
└── README.md             # This file
How to Compile and Run
Dependencies
An MPI implementation (e.g., OpenMPI)

A C compiler with OpenMP support (e.g., GCC)

Compilation
The hybrid MPI+OpenMP code is compiled using the mpicc wrapper, which correctly links the necessary MPI and OpenMP libraries. From the root directory, run:

Bash

# Load necessary modules on the HPC system first
module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0

# Compile the code
mpicc -fopenmp -o stencil_hybrid src/stencil_template_parallel.c -Iinclude -lm
Running
Jobs are submitted to the SLURM scheduler using batch scripts. An example script to run a hybrid job on one full node would contain:

Bash

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=7   # 7 MPI tasks per node
#SBATCH --cpus-per-task=16  # 16 OpenMP threads per task (7*16=112 cores)
#SBATCH --partition=your-partition
#SBATCH --account=your-account
#SBATCH --time=00:10:00

# Set environment variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# Run the executable
mpirun -np ${SLURM_NTASKS} ./stencil_hybrid -x 1000 -y 1000 -n 1000 -p 1
Submit the job with: sbatch your_script_name.sh

Performance Analysis and Results
OpenMP Scaling
An initial study was performed on a single node with one MPI task to determine the optimal number of OpenMP threads.

Finding: The analysis showed that 16 threads per MPI task provided the best balance of speedup and efficiency. This configuration was used for all subsequent MPI scaling tests.

Note: The initial implementation, which placed the OpenMP parallel directive inside the main loop, demonstrated poor scaling due to high fork-join overhead. The final version uses a "Single Parallel Region" pattern to achieve excellent OpenMP performance.

MPI Strong Scaling
The strong scaling study was performed with a fixed total problem size, using 7 MPI tasks per node and 16 threads per task to fully utilize the 112 cores on each node.

Result: The code demonstrated good speedup up to 4 nodes.

Analysis: As the node count increased to 8 and beyond, a growing gap appeared between the measured speedup and the ideal linear speedup. The instrumentation showed that this was caused by the MPI communication overhead, which began to represent a larger fraction of the total runtime as the problem was divided among more processes.

MPI Weak Scaling
The weak scaling study was performed by keeping the problem size per process constant while increasing the number of nodes.

Result: The total computation time remained nearly constant as expected.

Analysis: A slight increase in total runtime was observed as the node count grew. This is attributed to the increase in the total surface-to-volume ratio of the decomposed domain, which leads to a larger total amount of data being communicated across the network.

Conclusion
This project successfully demonstrates the implementation and analysis of a hybrid MPI+OpenMP stencil code.

The hybrid model effectively utilizes the resources of a modern HPC node by mapping MPI tasks to sockets and OpenMP threads to the cores within those sockets.

Strong scaling is effective but is ultimately limited by communication overhead, as described by Amdahl's Law.

Weak scaling proves the application's capability to solve increasingly larger problems by adding more resources, with performance limited by the total communication cost.

The results highlight the critical trade-off between parallel computation and inter-process communication that governs the performance of distributed-memory applications.
