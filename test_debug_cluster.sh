#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=7
#SBATCH --mem=0
#SBATCH --partition=dcgp_usr_prod
#SBATCH --account=uTS25_Tornator_0
#SBATCH --time=00:15:00
#SBATCH --job-name=stencil_debug_test
#SBATCH --output=debug_test_%j.out
#SBATCH --error=debug_test_%j.err

# === Setup ===
module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0

echo "=== MPI Stencil Debug Test on Leonardo Cluster ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $SLURM_NODELIST"
echo "Number of MPI tasks: $SLURM_NTASKS"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Total threads: $((SLURM_NTASKS * SLURM_CPUS_PER_TASK))"

# Set OpenMP environment variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

# Compile the code first
echo "Compiling stencil code..."
mpicc -fopenmp -O3 -o stencil_hybrid src/stencil_template_parallel.c -Iinclude -lm
if [ $? -ne 0 ]; then
    echo "ERROR: Compilation failed!"
    exit 1
fi
echo "Compilation successful!"

# Test 1: Very small problem (100x100, 10 iterations)
echo ""
echo "Test 1: Very small problem (100x100, 10 iterations)"
srun --nodes=1 --ntasks=4 --cpus-per-task=7 ./stencil_hybrid -x 100 -y 100 -n 10 -v 1
TEST1_EXIT=$?

if [ $TEST1_EXIT -eq 0 ]; then
    echo "✓ Test 1 PASSED"
else
    echo "✗ Test 1 FAILED (exit code: $TEST1_EXIT)"
    echo "Check the output files for details"
    exit 1
fi

# Test 2: Small problem (500x500, 50 iterations)
echo ""
echo "Test 2: Small problem (500x500, 50 iterations)"
srun --nodes=1 --ntasks=4 --cpus-per-task=7 ./stencil_hybrid -x 500 -y 500 -n 50 -v 1
TEST2_EXIT=$?

if [ $TEST2_EXIT -eq 0 ]; then
    echo "✓ Test 2 PASSED"
else
    echo "✗ Test 2 FAILED (exit code: $TEST2_EXIT)"
    echo "Check the output files for details"
    exit 1
fi

# Test 3: Medium problem (1000x1000, 100 iterations)
echo ""
echo "Test 3: Medium problem (1000x1000, 100 iterations)"
srun --nodes=1 --ntasks=4 --cpus-per-task=7 ./stencil_hybrid -x 1000 -y 1000 -n 100 -v 1
TEST3_EXIT=$?

if [ $TEST3_EXIT -eq 0 ]; then
    echo "✓ Test 3 PASSED"
else
    echo "✗ Test 3 FAILED (exit code: $TEST3_EXIT)"
    echo "Check the output files for details"
    exit 1
fi

# Test 4: Your original problem size (1000x1000, 1000 iterations)
echo ""
echo "Test 4: Original problem size (1000x1000, 1000 iterations)"
srun --nodes=1 --ntasks=4 --cpus-per-task=7 ./stencil_hybrid -x 1000 -y 1000 -n 1000
TEST4_EXIT=$?

if [ $TEST4_EXIT -eq 0 ]; then
    echo "✓ Test 4 PASSED"
else
    echo "✗ Test 4 FAILED (exit code: $TEST4_EXIT)"
    echo "Check the output files for details"
    exit 1
fi

echo ""
echo "=== All tests completed successfully! ==="
echo "Your code is ready for strong scaling experiments."
echo "You can now run: sbatch go_mpi_scaling" 