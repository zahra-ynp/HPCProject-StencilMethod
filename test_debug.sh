#!/bin/bash

# Debug test script for MPI stencil code
# This script helps identify issues with the MPI implementation

echo "=== MPI Stencil Debug Test ==="
echo "Testing with small problem size and 4 MPI processes..."

# Set OpenMP environment variables
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

# Test with very small problem size first
echo "Test 1: Very small problem (100x100, 10 iterations)"
mpirun -np 4 ./stencil_hybrid -x 100 -y 100 -n 10 -v 1 2>&1 | tee debug_test1.log

if [ $? -eq 0 ]; then
    echo "✓ Test 1 PASSED"
else
    echo "✗ Test 1 FAILED - check debug_test1.log"
    exit 1
fi

echo ""
echo "Test 2: Small problem (500x500, 50 iterations)"
mpirun -np 4 ./stencil_hybrid -x 500 -y 500 -n 50 -v 1 2>&1 | tee debug_test2.log

if [ $? -eq 0 ]; then
    echo "✓ Test 2 PASSED"
else
    echo "✗ Test 2 FAILED - check debug_test2.log"
    exit 1
fi

echo ""
echo "Test 3: Medium problem (1000x1000, 100 iterations)"
mpirun -np 4 ./stencil_hybrid -x 1000 -y 1000 -n 100 -v 1 2>&1 | tee debug_test3.log

if [ $? -eq 0 ]; then
    echo "✓ Test 3 PASSED"
else
    echo "✗ Test 3 FAILED - check debug_test3.log"
    exit 1
fi

echo ""
echo "All tests completed successfully!"
echo "You can now run your scaling tests with confidence." 