#!/bin/bash

# Define test configurations
particle_counts=(100000 1000000)
thread_counts=(2 4 6 8 12 24 48)
steps=10

# Output file
output_file="optimization_results.txt"

# Ensure the binary is built with the latest optimizations
echo "Building optimized binary..."
make clean && make -j8

# Create or clear the output file
echo "# N-body simulation benchmark results" > $output_file
echo "# Particle counts: ${particle_counts[*]}" >> $output_file
echo "# Thread counts: ${thread_counts[*]}" >> $output_file
echo "# Steps: $steps" >> $output_file
echo "# Format: particles threads run_number time(seconds) particles_per_second" >> $output_file
echo "" >> $output_file

# Run benchmarks
for parts in "${particle_counts[@]}"
do
	for threads in "${thread_counts[@]}"
	do
		echo "Running benchmark with $parts particles using $threads threads..."
		
		# Run each test configuration multiple times for statistical validity
		for run in {1..5}
		do
			echo "  Run $run of 5"
			
			# Capture the output with time information
			output=$(OMP_NUM_THREADS=$threads ./kdtree-sim $steps $parts 2>&1)
			
			# Extract the time and performance metrics
			time_seconds=$(echo "$output" | grep "Simulation completed in" | awk '{print $4}')
			particles_per_second=$(echo "$output" | grep "Performance" | awk '{print $2}')
			
			# Log the results
			echo "$parts $threads $run $time_seconds $particles_per_second" >> $output_file
		done
	done
done

# Generate summary statistics
echo -e "\nGenerating summary statistics..."
echo -e "\n# Summary Statistics (average time per configuration)" >> $output_file

for parts in "${particle_counts[@]}"
do
	for threads in "${thread_counts[@]}"
	do
		# Calculate average time and performance
		avg_time=$(grep "^$parts $threads" $output_file | awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }')
		avg_perf=$(grep "^$parts $threads" $output_file | awk '{ sum += $5; n++ } END { if (n > 0) print sum / n; }')
		
		echo "$parts particles, $threads threads: $avg_time seconds, $avg_perf particles/sec" >> $output_file
	done
done

echo "Benchmark completed. Results saved to $output_file"
