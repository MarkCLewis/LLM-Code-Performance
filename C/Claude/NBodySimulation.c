#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define G 6.67430e-11  // Gravitational constant
#define SOLAR_MASS 1.989e30  // Mass of the Sun in kg
#define AU 1.496e11  // Astronomical Unit in meters
#define YEAR 3.154e7  // Year in seconds

typedef struct {
    double x, y, z;       // Position
    double vx, vy, vz;    // Velocity
    double mass;          // Mass
} Body;

// Vector operations
double vector_distance(double x1, double y1, double z1, double x2, double y2, double z2) {
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

// Calculate total energy of the system (kinetic + potential)
double calculate_energy(Body *bodies, int n) {
    double energy = 0.0;
    
    // Calculate kinetic and potential energy
    for (int i = 0; i < n; i++) {
        // Kinetic energy
        double v_squared = bodies[i].vx * bodies[i].vx + 
                           bodies[i].vy * bodies[i].vy + 
                           bodies[i].vz * bodies[i].vz;
        energy += 0.5 * bodies[i].mass * v_squared;
        
        // Potential energy with all other bodies
        for (int j = i + 1; j < n; j++) {
            double distance = vector_distance(
                bodies[i].x, bodies[i].y, bodies[i].z,
                bodies[j].x, bodies[j].y, bodies[j].z
            );
            energy -= G * bodies[i].mass * bodies[j].mass / distance;
        }
    }
    
    return energy;
}

// Update velocities of all bodies
void update_velocities(Body *bodies, int n, double dt) {
    for (int i = 0; i < n; i++) {
        double fx = 0.0, fy = 0.0, fz = 0.0;
        
        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            
            double dx = bodies[j].x - bodies[i].x;
            double dy = bodies[j].y - bodies[i].y;
            double dz = bodies[j].z - bodies[i].z;
            
            double distance = sqrt(dx*dx + dy*dy + dz*dz);
            double distance_cubed = distance * distance * distance;
            
            // Skip if bodies are too close (to avoid numerical instability)
            if (distance < 1e-10) continue;
            
            // Calculate gravitational force
            double force = G * bodies[i].mass * bodies[j].mass / distance_cubed;
            
            fx += force * dx;
            fy += force * dy;
            fz += force * dz;
        }
        
        // Update velocities (F = ma => a = F/m)
        bodies[i].vx += fx / bodies[i].mass * dt;
        bodies[i].vy += fy / bodies[i].mass * dt;
        bodies[i].vz += fz / bodies[i].mass * dt;
    }
}

// Update positions of all bodies
void update_positions(Body *bodies, int n, double dt) {
    for (int i = 0; i < n; i++) {
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
        bodies[i].z += bodies[i].vz * dt;
    }
}

// Kick-step integration (first-order symplectic)
void simulate_step(Body *bodies, int n, double dt) {
    update_velocities(bodies, n, dt);  // Kick
    update_positions(bodies, n, dt);   // Drift
}

// Initialize system with central body and orbiting bodies
void initialize_system(Body *bodies, int n_bodies, double central_mass, double orbit_radius_min, double orbit_radius_max) {
    if (n_bodies < 1) return;
    
    // Set central body at origin
    bodies[0].x = 0.0;
    bodies[0].y = 0.0;
    bodies[0].z = 0.0;
    bodies[0].vx = 0.0;
    bodies[0].vy = 0.0;
    bodies[0].vz = 0.0;
    bodies[0].mass = central_mass;
    
    srand(time(NULL));
    
    // Create orbiting bodies
    for (int i = 1; i < n_bodies; i++) {
        // Random angle for position
        double phi = ((double)rand() / RAND_MAX) * 2.0 * M_PI;
        double theta = ((double)rand() / RAND_MAX) * M_PI - M_PI/2.0;
        
        // Random distance within range
        double radius = orbit_radius_min + ((double)rand() / RAND_MAX) * (orbit_radius_max - orbit_radius_min);
        
        // Position on sphere
        bodies[i].x = radius * cos(theta) * cos(phi);
        bodies[i].y = radius * cos(theta) * sin(phi);
        bodies[i].z = radius * sin(theta);
        
        // Velocity for circular orbit
        double orbit_speed = sqrt(G * central_mass / radius);
        
        // Cross product with up vector to get perpendicular velocity
        bodies[i].vx = -orbit_speed * sin(phi);
        bodies[i].vy = orbit_speed * cos(phi);
        bodies[i].vz = 0.0;
        
        // Small mass compared to central body
        bodies[i].mass = central_mass * 1e-6;
    }
}

void print_system_stats(Body *bodies, int n) {
    printf("System statistics:\n");
    printf("Central body: position (%.2e, %.2e, %.2e), mass: %.2e\n", 
           bodies[0].x, bodies[0].y, bodies[0].z, bodies[0].mass);
    
    // Calculate system boundaries
    double min_x = bodies[0].x, max_x = bodies[0].x;
    double min_y = bodies[0].y, max_y = bodies[0].y;
    double min_z = bodies[0].z, max_z = bodies[0].z;
    
    for (int i = 1; i < n; i++) {
        if (bodies[i].x < min_x) min_x = bodies[i].x;
        if (bodies[i].x > max_x) max_x = bodies[i].x;
        if (bodies[i].y < min_y) min_y = bodies[i].y;
        if (bodies[i].y > max_y) max_y = bodies[i].y;
        if (bodies[i].z < min_z) min_z = bodies[i].z;
        if (bodies[i].z > max_z) max_z = bodies[i].z;
    }
    
    printf("System boundaries:\n");
    printf("  X: [%.2e, %.2e]\n", min_x, max_x);
    printf("  Y: [%.2e, %.2e]\n", min_y, max_y);
    printf("  Z: [%.2e, %.2e]\n", min_z, max_z);
}

int main() {
    // Problem parameters
    const int N = 10000 + 1;  // 1 central + 1 million orbiting bodies
    const int STEPS = 100;
    const double dt = 0.01 * YEAR;  // Timestep of 0.01 years
    
    printf("Starting N-body simulation with %d bodies for %d steps\n", N, STEPS);
    printf("Allocating memory...\n");
    
    // Allocate memory for bodies
    Body *bodies = (Body *)malloc(N * sizeof(Body));
    if (!bodies) {
        printf("Memory allocation failed!\n");
        return 1;
    }
    
    printf("Initializing system...\n");
    
    // Initialize system with Sun-like central body and Earth-like orbits
    initialize_system(bodies, N, SOLAR_MASS, 0.8 * AU, 5.0 * AU);
    
    print_system_stats(bodies, N);
    
    // Calculate initial energy
    printf("Calculating initial energy...\n");
    double initial_energy = calculate_energy(bodies, N);
    printf("Initial energy: %.10e\n", initial_energy);
    
    printf("Starting simulation...\n");
    clock_t start = clock();
    
    // Run simulation
    for (int step = 0; step < STEPS; step++) {
        if (step % 100 == 0) {
            printf("Step %d of %d (%.1f%%)\n", step, STEPS, 100.0 * step / STEPS);
        }
        simulate_step(bodies, N, dt);
    }
    
    clock_t end = clock();
    double cpu_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Simulation completed in %.2f seconds\n", cpu_time);
    
    // Calculate final energy
    printf("Calculating final energy...\n");
    double final_energy = calculate_energy(bodies, N);
    printf("Final energy: %.10e\n", final_energy);
    
    // Calculate energy error
    double energy_error = (final_energy - initial_energy) / initial_energy;
    printf("Relative energy error: %.10e (%.10f%%)\n", energy_error, energy_error * 100.0);
    
    print_system_stats(bodies, N);
    
    // Free memory
    free(bodies);
    
    return 0;
}