#define _POSIX_C_SOURCE 199309L  // Required for clock_gettime and CLOCK_MONOTONIC

#include "kdtree.h"
#include "particle.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <string.h>
#include <math.h>  // Add for fabs function

// Measure elapsed time using high-resolution clock
double get_time() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
}

// Print energy details with nice formatting
void print_energy_info(const char *label, const SystemEnergy *energy) {
  printf("=== %s ===\n", label);
  printf("  Kinetic Energy:   %.10e\n", energy->kinetic);
  printf("  Potential Energy: %.10e\n", energy->potential);
  printf("  Total Energy:     %.10e\n", energy->total);
}

int main(int argc, char *argv[]) {
  if (argc < 3) {
    printf("Usage: %s <steps> <particles> [threads]\n", argv[0]);
    printf("  steps     : Number of simulation steps to run\n");
    printf("  particles : Number of particles to simulate\n");
    printf("  threads   : Number of OpenMP threads (optional, default: max available)\n");
    return 1;
  }
  
  int steps = atoi(argv[1]);
  int n = atoi(argv[2]);
  
  // Set thread count if specified, otherwise use all available cores
  if (argc >= 4) {
    int threads = atoi(argv[3]);
    if (threads > 0) {
      omp_set_num_threads(threads);
    }
  }
  
  printf("Running N-body simulation with %d particles for %d steps using %d threads.\n", 
         n, steps, omp_get_max_threads());

  double dt = 1e-3; // Time step

  // Initialize and measure time
  double start_time = get_time();
  
  // Generate particles in circular orbits
  Particle_array_t particles = circular_orbits(n);
  
  // Calculate initial energy to verify accuracy
  SystemEnergy initial_energy = calculate_system_energy(&particles);
  print_energy_info("Initial Energy", &initial_energy);
  
  // Run simulation
  simple_sim(&particles, dt, steps);
  
  // Calculate final energy and verify conservation
  SystemEnergy final_energy = calculate_system_energy(&particles);
  print_energy_info("Final Energy", &final_energy);
  
  // Calculate energy difference (should be close to zero for a good simulation)
  double energy_diff = fabs(final_energy.total - initial_energy.total);
  double energy_rel_error = energy_diff / fabs(initial_energy.total);
  
  printf("\n=== Energy Conservation ===\n");
  printf("  Absolute difference: %.10e\n", energy_diff);
  printf("  Relative error:      %.10e (%.6f%%)\n", 
         energy_rel_error, energy_rel_error * 100.0);
  
  // Assess accuracy
  if (energy_rel_error < 1e-3) {
    printf("  Status: GOOD - Energy is well conserved\n");
  } else if (energy_rel_error < 1e-2) {
    printf("  Status: ACCEPTABLE - Minor energy drift\n");
  } else {
    printf("  Status: POOR - Significant energy drift, simulation accuracy concerns\n");
  }
  
  // Calculate and print elapsed time
  double elapsed = get_time() - start_time;
  printf("\nSimulation completed in %.6f seconds\n", elapsed);
  printf("Performance: %.2f particle-steps per second\n", 
         (double)n * steps / elapsed);

  // Clean up
  FREE_ARRAY(particles);

  return 0;
}
