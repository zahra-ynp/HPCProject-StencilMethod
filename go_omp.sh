#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH --partition=dcgp_usr_prod
#SBATCH --account=uTS25_Tornator_0
#SBATCH --time=01:00:00
#SBATCH --job-name=omp_scaling_112_threads
#SBATCH --output=omp_scaling_results/log_%j.out
#SBATCH --error=omp_scaling_results/err_%j.err

# === Setup ===
module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0


EXEC="stencil_hybrid"
ARGS="-x 1000 -y 1000 -n 1000"

THREAD_COUNTS=(1 2 4 8 16 32 56 84 112)

echo "Starting OpenMP scaling tests..."

# === Loop over thread counts ===
for t in "${THREAD_COUNTS[@]}"; do
    export OMP_NUM_THREADS=$t
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close

    echo "Running with $t OpenMP threads..."
    OUT_FILE="omp_scaling_results/out_omp_threads_${t}.txt
    mpirun -np 1 ./${EXEC} ${ARGS} > ${OUT_FILE}
done
echo "All runs completed."
