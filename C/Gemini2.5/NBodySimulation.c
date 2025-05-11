#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> // For srand and timing (optional)

// --- Constants ---
#define G              1.0     // Gravitational constant (simplified units)
#define DT             0.001   // Time step (adjust based on system dynamics)
#define NUM_STEPS      100    // Number of simulation steps
#define N_ORBITING     10000 // Number of orbiting bodies <<< WARNING: VERY LARGE for sequential!
#define N_TOTAL        (N_ORBITING + 1) // Total number of bodies (central + orbiting)
#define EPSILON_SQ     1e-6    // Softening factor squared to prevent singularities
#define M_PI           3.1415927 // Gemini 2.5 didn't create this.

// Central Body Properties
#define CENTRAL_MASS   1000.0

// Orbiting Body Properties
#define ORBITING_MASS  0.01     // Mass of smaller bodies (should be << CENTRAL_MASS)
#define MIN_RADIUS     5.0      // Minimum orbital radius
#define MAX_RADIUS     50.0     // Maximum orbital radius

// --- Data Structure ---
typedef struct {
    double *posX, *posY, *posZ;
    double *velX, *velY, *velZ;
    double *mass;
    int N; // Number of bodies in the system
} System;

// --- Function Prototypes ---
void allocate_system(System *sys, int N);
void free_system(System *sys);
void initialize_system(System *sys, double central_mass, double orbiting_mass,
                       double min_radius, double max_radius, int n_orbiting);
// double calculate_energy(const System *sys);
void calculate_forces(const System *sys, double *forceX, double *forceY, double *forceZ);
void simulate_step(System *sys);

// --- Main Function ---
int main() {
    printf("N-Body Simulation\n");
    printf("-------------------\n");
    printf("Bodies: %d (%d orbiting + 1 central)\n", N_TOTAL, N_ORBITING);
    printf("Steps:  %d\n", NUM_STEPS);
    printf("DT:     %.4f\n", DT);
    printf("G:      %.2f\n", G);
    printf("WARNING: N=%d is very large for a sequential O(N^2) simulation.\n", N_TOTAL);
    printf("         Execution will likely take an extremely long time.\n\n");

    System sys;
    allocate_system(&sys, N_TOTAL);

    // Initialize the system
    printf("Initializing system...\n");
    initialize_system(&sys, CENTRAL_MASS, ORBITING_MASS, MIN_RADIUS, MAX_RADIUS, N_ORBITING);
    printf("Initialization complete.\n");

    // Calculate initial energy
    // printf("Calculating initial energy...\n");
    // double initial_energy = calculate_energy(&sys);
    // printf("Initial System Energy: %.6e\n", initial_energy);

    // --- Simulation Loop ---
    printf("Starting simulation...\n");
    time_t start_time = time(NULL); // Optional: track time

    for (int step = 0; step < NUM_STEPS; ++step) {
        simulate_step(&sys);

        // Optional: Print progress
        if ((step + 1) % (NUM_STEPS / 10) == 0 || step == NUM_STEPS - 1) {
            time_t current_time = time(NULL);
            double elapsed_sec = difftime(current_time, start_time);
            printf("Step %d / %d completed. Elapsed time: %.2f s\n",
                   step + 1, NUM_STEPS, elapsed_sec);
        }
    }

    time_t end_time = time(NULL);
    double total_time = difftime(end_time, start_time);
    printf("Simulation finished.\n");
    printf("Total simulation time: %.2f seconds\n", total_time);

    // Calculate final energy
    // printf("Calculating final energy...\n");
    // double final_energy = calculate_energy(&sys);
    // printf("Final System Energy:   %.6e\n", final_energy);

    // // Verify energy conservation
    // double energy_diff = final_energy - initial_energy;
    // double energy_rel_diff = (initial_energy != 0.0) ? (energy_diff / initial_energy) : 0.0;
    // printf("Energy Difference (Final - Initial): %.6e\n", energy_diff);
    // printf("Relative Energy Difference: %.6e (%.4f%%)\n", energy_rel_diff, energy_rel_diff * 100.0);

    // Clean up
    free_system(&sys);

    return 0;
}

// --- Function Implementations ---

// Allocate memory for system data
void allocate_system(System *sys, int N) {
    sys->N = N;
    sys->posX = (double*)malloc(N * sizeof(double));
    sys->posY = (double*)malloc(N * sizeof(double));
    sys->posZ = (double*)malloc(N * sizeof(double));
    sys->velX = (double*)malloc(N * sizeof(double));
    sys->velY = (double*)malloc(N * sizeof(double));
    sys->velZ = (double*)malloc(N * sizeof(double));
    sys->mass = (double*)malloc(N * sizeof(double));

    // Check for allocation errors
    if (!sys->posX || !sys->posY || !sys->posZ || !sys->velX || !sys->velY || !sys->velZ || !sys->mass) {
        fprintf(stderr, "Error: Memory allocation failed!\n");
        // Free any already allocated memory before exiting
        free(sys->posX); free(sys->posY); free(sys->posZ);
        free(sys->velX); free(sys->velY); free(sys->velZ);
        free(sys->mass);
        exit(EXIT_FAILURE);
    }
}

// Free memory used by the system
void free_system(System *sys) {
    free(sys->posX);
    free(sys->posY);
    free(sys->posZ);
    free(sys->velX);
    free(sys->velY);
    free(sys->velZ);
    free(sys->mass);
    sys->N = 0; // Reset N
}

// Initialize the system with a central body and orbiting bodies
void initialize_system(System *sys, double central_mass, double orbiting_mass,
                       double min_radius, double max_radius, int n_orbiting) {
    if (sys->N != n_orbiting + 1) {
        fprintf(stderr, "Error: System size mismatch during initialization!\n");
        exit(EXIT_FAILURE);
    }

    srand(time(NULL)); // Seed random number generator

    // --- Central Body (index 0) ---
    sys->mass[0] = central_mass;
    sys->posX[0] = 0.0;
    sys->posY[0] = 0.0;
    sys->posZ[0] = 0.0;
    sys->velX[0] = 0.0;
    sys->velY[0] = 0.0;
    sys->velZ[0] = 0.0;

    // --- Orbiting Bodies (indices 1 to N-1) ---
    for (int i = 1; i < sys->N; ++i) {
        sys->mass[i] = orbiting_mass;

        // Generate random radius and angle for position (in XY plane for simplicity)
        double r = min_radius + (max_radius - min_radius) * ((double)rand() / RAND_MAX);
        double angle = 2.0 * M_PI * ((double)rand() / RAND_MAX); // 0 to 2*pi

        sys->posX[i] = r * cos(angle);
        sys->posY[i] = r * sin(angle);
        sys->posZ[i] = 0.0; // Keep orbits initially in the XY plane

        // Calculate velocity for a circular orbit around the central mass
        // v = sqrt(G * M_central / r)
        // Velocity vector is perpendicular to position vector (-y, x)
        if (r > 0) { // Avoid division by zero if r happens to be 0
            double orbital_speed = sqrt(G * sys->mass[0] / r);
            sys->velX[i] = -orbital_speed * sin(angle); // -v * (y/r)
            sys->velY[i] =  orbital_speed * cos(angle); //  v * (x/r)
            sys->velZ[i] =  0.0;
        } else {
            sys->velX[i] = 0.0;
            sys->velY[i] = 0.0;
            sys->velZ[i] = 0.0;
        }
         // Optional: Add slight random perturbation to Z position/velocity
         // sys->posZ[i] = ( (double)rand() / RAND_MAX - 0.5 ) * 0.1 * r; // Small z offset
         // sys->velZ[i] = ( (double)rand() / RAND_MAX - 0.5 ) * 0.1 * orbital_speed; // Small z velocity
    }
}

// Calculate the total energy (Kinetic + Potential) of the system
// double calculate_energy(const System *sys) {
//     double kinetic_energy = 0.0;
//     double potential_energy = 0.0;

//     // Kinetic Energy: Sum(0.5 * m_i * v_i^2)
//     for (int i = 0; i < sys->N; ++i) {
//         double v_sq = sys->velX[i] * sys->velX[i] +
//                       sys->velY[i] * sys->velY[i] +
//                       sys->velZ[i] * sys->velZ[i];
//         kinetic_energy += 0.5 * sys->mass[i] * v_sq;
//     }

//     // Potential Energy: Sum(-G * m_i * m_j / |r_i - r_j|) for i < j
//     for (int i = 0; i < sys->N; ++i) {
//         for (int j = i + 1; j < sys->N; ++j) {
//             double dx = sys->posX[j] - sys->posX[i];
//             double dy = sys->posY[j] - sys->posY[i];
//             double dz = sys->posZ[j] - sys->posZ[i];

//             double r_sq = dx * dx + dy * dy + dz * dz;
//             double r = sqrt(r_sq + EPSILON_SQ); // Add softening

//             if (r > sqrt(EPSILON_SQ)) { // Avoid potential division by zero if r is extremely close to 0
//                  potential_energy -= G * sys->mass[i] * sys->mass[j] / r;
//             }
//             // If r is effectively zero due to softening, potential is large but finite.
//             // A different handling might be needed for specific physical scenarios.
//         }
//     }

//     return kinetic_energy + potential_energy;
// }

// Calculate gravitational forces between all pairs of bodies
void calculate_forces(const System *sys, double *forceX, double *forceY, double *forceZ) {
    // Initialize forces to zero
    for (int i = 0; i < sys->N; ++i) {
        forceX[i] = 0.0;
        forceY[i] = 0.0;
        forceZ[i] = 0.0;
    }

    // Calculate pairwise forces (Newton's Law of Gravitation)
    // F_ij = G * m_i * m_j / |r_ij|^2 * (r_j - r_i) / |r_ij|
    for (int i = 0; i < sys->N; ++i) {
        for (int j = i + 1; j < sys->N; ++j) { // Loop j from i+1 to avoid double counting and self-interaction
            double dx = sys->posX[j] - sys->posX[i];
            double dy = sys->posY[j] - sys->posY[i];
            double dz = sys->posZ[j] - sys->posZ[i];

            double r_sq = dx * dx + dy * dy + dz * dz;
            double dist = sqrt(r_sq + EPSILON_SQ); // Add softening

            // Force magnitude: G * m1 * m2 / (r^2 + eps^2)
            // Avoid division by zero if dist is extremely small
            if (dist > 1e-10) { // Use a small threshold instead of exact zero check
                 double force_mag_over_dist = G * sys->mass[i] * sys->mass[j] / (r_sq + EPSILON_SQ) / dist;

                 // Calculate force components F_x = F * (dx/dist), etc.
                 double Fx = force_mag_over_dist * dx;
                 double Fy = force_mag_over_dist * dy;
                 double Fz = force_mag_over_dist * dz;

                 // Add force to particle i (Force exerted by j on i)
                 forceX[i] += Fx;
                 forceY[i] += Fy;
                 forceZ[i] += Fz;

                 // Add opposite force to particle j (Newton's 3rd Law)
                 forceX[j] -= Fx;
                 forceY[j] -= Fy;
                 forceZ[j] -= Fz;
            }
        }
    }
}

// Perform one step of the simulation using Kick-Step method
void simulate_step(System *sys) {
    // Allocate temporary arrays for forces
    // (Could be allocated once outside the loop if preferred, but N might change)
    double *forceX = (double*)malloc(sys->N * sizeof(double));
    double *forceY = (double*)malloc(sys->N * sizeof(double));
    double *forceZ = (double*)malloc(sys->N * sizeof(double));
    if (!forceX || !forceY || !forceZ) {
        fprintf(stderr, "Error: Force array allocation failed in simulate_step!\n");
        // Consider more robust error handling or freeing system memory
        exit(EXIT_FAILURE);
    }

    // 1. Calculate forces based on current positions
    calculate_forces(sys, forceX, forceY, forceZ);

    // 2. Update velocities ("Kick") based on forces
    // a = F/m  =>  v_new = v_old + a * dt = v_old + (F/m) * dt
    for (int i = 0; i < sys->N; ++i) {
        if (sys->mass[i] > 0) { // Avoid division by zero for massless particles (if any)
            sys->velX[i] += (forceX[i] / sys->mass[i]) * DT;
            sys->velY[i] += (forceY[i] / sys->mass[i]) * DT;
            sys->velZ[i] += (forceZ[i] / sys->mass[i]) * DT;
        }
    }

    // 3. Update positions ("Drift") based on *new* velocities
    // x_new = x_old + v_new * dt
    for (int i = 0; i < sys->N; ++i) {
        sys->posX[i] += sys->velX[i] * DT;
        sys->posY[i] += sys->velY[i] * DT;
        sys->posZ[i] += sys->velZ[i] * DT;
    }

    // Free temporary force arrays
    free(forceX);
    free(forceY);
    free(forceZ);
}