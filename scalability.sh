#!/bin/bash

# Test scalability of the Laplace solver

# Ensure the script exits if any command fails
set -e

# Compile the program with the SCALABILITY_TEST environment variable set
make clean

make SCALABILITY_TEST=0

# Define the number of processors and grid sizes for testing
processor_counts=(1 2 4 8)

# Output file for results
output_file="output/scalabilityResults.txt"

if [ -f "$output_file" ] ; then
    rm "$output_file"
fi
echo "*******************************" >> $output_file
echo "Scalability Test" >> $output_file
echo "*******************************" >> $output_file
echo "" >> $output_file

# Run the tests and capture the output
for np in "${processor_counts[@]}"; do  
        echo "Running test with ${np} processors"
        echo "-------------------------------------" >> $output_file
        mpirun -np ${np} ./main >> $output_file
        echo "" >> $output_file
done

# Summarize results
echo "Scalability test completed. Results are stored in $output_file."