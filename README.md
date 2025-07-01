---
# HPC Project: The Stencil Method

## Overview

This project explores high-performance stencil computations—a core technique in scientific computing for solving partial differential equations (specifically, the 2D Heat Equation). The repository demonstrates the implementation and analysis of a hybrid MPI+OpenMP stencil code.

---

## Project Structure

```
HPCProject-StencilMethod/
├── include/               # Header files (.h) for the C implementation
├── src/                   # Main C source files (.c)
├── mpi_weak_scaling/      # Scripts and outputs for MPI weak scaling
├── mpi_scaling_results/   # Scripts and outputs for MPI strong scaling
├── omp_scaling_results/   # Scripts and results from OpenMP scaling
├── plot/                  # Scripts and notebooks for analysis and plotting
├── go_scripts/            # SLURM batch scripts for job submission
└── README.md              # This file
```

---

## How to Compile and Run

### Dependencies

- An MPI implementation (e.g., OpenMPI)
- A C compiler with OpenMP support (e.g., GCC)

### Compilation

Compile the hybrid MPI+OpenMP code using the `mpicc` wrapper to link MPI and OpenMP libraries:

```bash
# Load necessary modules on the HPC system
module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0

# Compile the code
mpicc -fopenmp -o stencil_hybrid src/stencil_template_parallel.c -Iinclude -lm
```

### Running

Jobs are submitted to the SLURM scheduler using batch scripts. Example script for a hybrid job on one full node:

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=7     # 7 MPI tasks per node
#SBATCH --cpus-per-task=16      # 16 OpenMP threads per task (7*16=112 cores)
#SBATCH --partition=your-partition
#SBATCH --account=your-account
#SBATCH --time=00:10:00

# Set environment variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# Run the executable
mpirun -np ${SLURM_NTASKS} ./stencil_hybrid -x 1000 -y 1000 -n 1000 -p 1
```

Submit the job with:

```bash
sbatch your_script_name.sh
```

---

## Performance Analysis and Results

### OpenMP Scaling

- Initial studies on a single node (one MPI task) determined the optimal number of OpenMP threads.
- **Finding:** 16 threads per MPI task provided the best balance of speedup and efficiency.

### MPI Strong Scaling

- Conducted with a fixed total problem size, using 7 MPI tasks per node and 16 threads per task (112 cores/node).
- **Result:** Good speedup up to 4 nodes.
- **Analysis:** Beyond 4 nodes, speedup deviates from ideal due to MPI communication overhead.

### MPI Weak Scaling

- Performed by keeping the problem size per process constant while increasing the number of nodes.
- **Result:** Total computation time remained nearly constant.
- **Analysis:** Slight runtime increases with more nodes, attributed to increased communication cost.

---

## Conclusion

- The hybrid model efficiently utilizes modern HPC resources, mapping MPI tasks to sockets and OpenMP threads to cores.
- **Strong scaling** is effective but ultimately limited by the serial portion of the code(Amdahl's Law).
- **Weak scaling** shows the application's capacity to solve larger problems by adding resources, with performance limited by communication cost.
- The results highlight the trade-off between parallel computation and inter-process communication in distributed-memory applications.

---
