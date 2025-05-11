#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Define the structure for a body
typedef struct {
    double mass;
    double position[3];
    double velocity[3];
} Body;

// Define the structure for the system
typedef struct {
    int num_bodies;
    Body *bodies;
} System;

// Gravitational constant
const double G = 6.67430e-11;

// Function to calculate the distance squared between two bodies
double distance_squared(const Body *b1, const Body *b2) {
    double dx = b1->position[0] - b2->position[0];
    double dy = b1->position[1] - b2->position[1];
    double dz = b1->position[2] - b2->position[2];
    return dx * dx + dy * dy + dz * dz;
}

// Function to calculate the force on body i due to body j
void calculate_force(const Body *body_i, const Body *body_j, double force[3]) {
    double r2 = distance_squared(body_i, body_j);
    if (r2 > 1e-9) { // Avoid division by zero for very close bodies
        double r = sqrt(r2);
        double magnitude = (G * body_i->mass * body_j->mass) / r2;
        force[0] = magnitude * (body_j->position[0] - body_i->position[0]) / r;
        force[1] = magnitude * (body_j->position[1] - body_i->position[1]) / r;
        force[2] = magnitude * (body_j->position[2] - body_i->position[2]) / r;
    } else {
        force[0] = force[1] = force[2] = 0.0;
    }
}

// Function to calculate the total energy of the system
// double calculate_total_energy(const System *system) {
//     double kinetic_energy = 0.0;
//     double potential_energy = 0.0;

//     for (int i = 0; i < system->num_bodies; i++) {
//         // Kinetic energy
//         kinetic_energy += 0.5 * system->bodies[i].mass * (system->bodies[i].velocity[0] * system->bodies[i].velocity[0] +
//                                                           system->bodies[i].velocity[1] * system->bodies[i].velocity[1] +
//                                                           system->bodies[i].velocity[2] * system->bodies[i].velocity[2]);

//         // Potential energy
//         for (int j = i + 1; j < system->num_bodies; j++) {
//             potential_energy -= (G * system->bodies[i].mass * system->bodies[j].mass) / sqrt(distance_squared(&system->bodies[i], &system->bodies[j]));
//         }
//     }

//     return kinetic_energy + potential_energy;
// }

// Function to initialize a system with a central body and orbiting smaller bodies
System* initialize_circular_orbits(int num_orbiting_bodies, double central_mass, double orbit_radius, double orbiting_mass) {
    System *system = (System*)malloc(sizeof(System));
    if (system == NULL) {
        fprintf(stderr, "Memory allocation failed for system.\n");
        exit(EXIT_FAILURE);
    }
    system->num_bodies = num_orbiting_bodies + 1;
    system->bodies = (Body*)malloc(sizeof(Body) * system->num_bodies);
    if (system->bodies == NULL) {
        fprintf(stderr, "Memory allocation failed for bodies.\n");
        free(system);
        exit(EXIT_FAILURE);
    }

    // Initialize the central body
    system->bodies[0].mass = central_mass;
    system->bodies[0].position[0] = 0.0;
    system->bodies[0].position[1] = 0.0;
    system->bodies[0].position[2] = 0.0;
    system->bodies[0].velocity[0] = 0.0;
    system->bodies[0].velocity[1] = 0.0;
    system->bodies[0].velocity[2] = 0.0;

    // Initialize the orbiting bodies
    for (int i = 1; i < system->num_bodies; i++) {
        system->bodies[i].mass = orbiting_mass;
        system->bodies[i].position[0] = orbit_radius * cos(2 * M_PI * (i - 1) / num_orbiting_bodies);
        system->bodies[i].position[1] = orbit_radius * sin(2 * M_PI * (i - 1) / num_orbiting_bodies);
        system->bodies[i].position[2] = 0.0; // Orbiting in the xy-plane

        // Calculate the orbital velocity for a circular orbit
        double orbital_speed = sqrt(G * central_mass / orbit_radius);
        system->bodies[i].velocity[0] = -orbital_speed * sin(2 * M_PI * (i - 1) / num_orbiting_bodies);
        system->bodies[i].velocity[1] = orbital_speed * cos(2 * M_PI * (i - 1) / num_orbiting_bodies);
        system->bodies[i].velocity[2] = 0.0;
    }

    return system;
}

// Function to perform one kick-step
void kick_step(System *system, double dt) {
    int n = system->num_bodies;
    double (*forces)[3] = malloc(sizeof(double[n][3]));
    if (forces == NULL) {
        fprintf(stderr, "Memory allocation failed for forces.\n");
        exit(EXIT_FAILURE);
    }

    // Calculate forces
    for (int i = 0; i < n; i++) {
        forces[i][0] = forces[i][1] = forces[i][2] = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                double f[3];
                calculate_force(&system->bodies[i], &system->bodies[j], f);
                forces[i][0] += f[0];
                forces[i][1] += f[1];
                forces[i][2] += f[2];
            }
        }
    }

    // Update velocities (kick)
    for (int i = 0; i < n; i++) {
        system->bodies[i].velocity[0] += (forces[i][0] / system->bodies[i].mass) * dt * 0.5;
        system->bodies[i].velocity[1] += (forces[i][1] / system->bodies[i].mass) * dt * 0.5;
        system->bodies[i].velocity[2] += (forces[i][2] / system->bodies[i].mass) * dt * 0.5;
    }

    free(forces);
}

// Function to perform one drift-step
void drift_step(System *system, double dt) {
    for (int i = 0; i < system->num_bodies; i++) {
        system->bodies[i].position[0] += system->bodies[i].velocity[0] * dt;
        system->bodies[i].position[1] += system->bodies[i].velocity[1] * dt;
        system->bodies[i].position[2] += system->bodies[i].velocity[2] * dt;
    }
}

// Function to perform one full first-order kick-step
void first_order_kick_step(System *system, double dt) {
    kick_step(system, dt);
    drift_step(system, dt);
    kick_step(system, dt);
}

// Function to free the memory allocated for the system
void free_system(System *system) {
    if (system != NULL) {
        free(system->bodies);
        free(system);
    }
}

int main() {
    int num_orbiting_bodies = 10000;
    double central_mass = 1.989e30; // Mass of the Sun (kg)
    double orbit_radius = 1.496e11; // 1 AU (m)
    double orbiting_mass = 5.972e24 / num_orbiting_bodies; // Mass of the Earth (kg)
    int num_steps = 100;
    double time_step = 3600 * 24 * 7; // 1 week in seconds

    // Initialize the system
    System *system = initialize_circular_orbits(num_orbiting_bodies, central_mass, orbit_radius, orbiting_mass);
    if (system == NULL) {
        return 1;
    }

    printf("Initial number of bodies: %d\n", system->num_bodies);

    // Calculate initial energy
    // double initial_energy = calculate_total_energy(system);
    // printf("Initial total energy: %e J\n", initial_energy);

    // Run the simulation
    printf("Running simulation for %d steps...\n", num_steps);
    time_t start_time = time(NULL);
    for (int step = 0; step < num_steps; step++) {
        first_order_kick_step(system, time_step);
        if ((step + 1) % 100 == 0) {
            printf("Step %d completed.\n", step + 1);
        }
    }
    time_t end_time = time(NULL);
    double elapsed_time = difftime(end_time, start_time);
    printf("Simulation finished in %.2f seconds.\n", elapsed_time);

    // Calculate final energy
    // double final_energy = calculate_total_energy(system);
    // printf("Final total energy: %e J\n", final_energy);

    // // Calculate the energy difference
    // double energy_difference = fabs(final_energy - initial_energy);
    // double relative_energy_difference = energy_difference / fabs(initial_energy);
    // printf("Absolute energy difference: %e J\n", energy_difference);
    // printf("Relative energy difference: %e\n", relative_energy_difference);

    // Free the memory
    free_system(system);

    return 0;
}