#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <string.h>

#define G 6.67430e-11  // Gravitational constant
#define SOLAR_MASS 1.989e30  // Mass of the Sun in kg
#define AU 1.496e11  // Astronomical Unit in meters
#define YEAR 3.154e7  // Year in seconds
#define THETA 0.3     // Barnes-Hut opening angle threshold

typedef struct {
    double x, y, z;       // Position
    double vx, vy, vz;    // Velocity
    double mass;          // Mass
} Body;

// Node structure for kD-tree / octree
typedef struct OctNode {
    double x_min, x_max;  // Boundaries
    double y_min, y_max;
    double z_min, z_max;
    
    double center_x, center_y, center_z;  // Center of mass
    double total_mass;                    // Total mass
    
    int body_index;       // Index of body if this is a leaf (-1 otherwise)
    int is_leaf;          // Flag for leaf nodes
    int num_bodies;       // Number of bodies in this node
    
    struct OctNode* children[8];  // Octants
} OctNode;

// Vector operations
double vector_distance(double x1, double y1, double z1, double x2, double y2, double z2) {
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

// Calculate width of node
double node_width(OctNode* node) {
    return node->x_max - node->x_min;
}

// Create a new node
OctNode* create_node(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max) {
    OctNode* node = (OctNode*)malloc(sizeof(OctNode));
    if (!node) {
        fprintf(stderr, "Memory allocation failed for octree node!\n");
        exit(1);
    }
    
    node->x_min = x_min;
    node->x_max = x_max;
    node->y_min = y_min;
    node->y_max = y_max;
    node->z_min = z_min;
    node->z_max = z_max;
    
    node->center_x = 0.0;
    node->center_y = 0.0;
    node->center_z = 0.0;
    node->total_mass = 0.0;
    
    node->body_index = -1;
    node->is_leaf = 1;
    node->num_bodies = 0;
    
    for (int i = 0; i < 8; i++) {
        node->children[i] = NULL;
    }
    
    return node;
}

// Determine which octant a body belongs to
int get_octant(OctNode* node, Body* body) {
    int octant = 0;
    double mid_x = (node->x_min + node->x_max) / 2.0;
    double mid_y = (node->y_min + node->y_max) / 2.0;
    double mid_z = (node->z_min + node->z_max) / 2.0;
    
    if (body->x >= mid_x) octant |= 1;
    if (body->y >= mid_y) octant |= 2;
    if (body->z >= mid_z) octant |= 4;
    
    return octant;
}

// Create child nodes for each octant
void subdivide_node(OctNode* node) {
    double mid_x = (node->x_min + node->x_max) / 2.0;
    double mid_y = (node->y_min + node->y_max) / 2.0;
    double mid_z = (node->z_min + node->z_max) / 2.0;
    
    node->children[0] = create_node(node->x_min, mid_x, node->y_min, mid_y, node->z_min, mid_z);
    node->children[1] = create_node(mid_x, node->x_max, node->y_min, mid_y, node->z_min, mid_z);
    node->children[2] = create_node(node->x_min, mid_x, mid_y, node->y_max, node->z_min, mid_z);
    node->children[3] = create_node(mid_x, node->x_max, mid_y, node->y_max, node->z_min, mid_z);
    node->children[4] = create_node(node->x_min, mid_x, node->y_min, mid_y, mid_z, node->z_max);
    node->children[5] = create_node(mid_x, node->x_max, node->y_min, mid_y, mid_z, node->z_max);
    node->children[6] = create_node(node->x_min, mid_x, mid_y, node->y_max, mid_z, node->z_max);
    node->children[7] = create_node(mid_x, node->x_max, mid_y, node->y_max, mid_z, node->z_max);
    
    node->is_leaf = 0;
}

// Insert a body into the octree
void insert_body(OctNode* node, Body* bodies, int body_index) {
    // Update center of mass and total mass
    double new_mass = node->total_mass + bodies[body_index].mass;
    
    // Weighted average for center of mass
    if (new_mass > 0) {
        node->center_x = (node->center_x * node->total_mass + bodies[body_index].x * bodies[body_index].mass) / new_mass;
        node->center_y = (node->center_y * node->total_mass + bodies[body_index].y * bodies[body_index].mass) / new_mass;
        node->center_z = (node->center_z * node->total_mass + bodies[body_index].z * bodies[body_index].mass) / new_mass;
        node->total_mass = new_mass;
    }
    
    node->num_bodies++;
    
    // If node is empty, store body index
    if (node->num_bodies == 1) {
        node->body_index = body_index;
        return;
    }
    
    // If this was a leaf with one body, need to move that body before adding another
    if (node->is_leaf && node->num_bodies == 2) {
        int old_body_index = node->body_index;
        subdivide_node(node);
        
        // Move existing body to the appropriate child
        int octant = get_octant(node, &bodies[old_body_index]);
        insert_body(node->children[octant], bodies, old_body_index);
        
        node->body_index = -1;  // No longer a leaf containing a single body
    }
    
    // Insert new body into appropriate child if not a leaf
    if (!node->is_leaf) {
        int octant = get_octant(node, &bodies[body_index]);
        insert_body(node->children[octant], bodies, body_index);
    }
}

// Free the octree memory
void free_octree(OctNode* node) {
    if (!node) return;
    
    for (int i = 0; i < 8; i++) {
        if (node->children[i]) {
            free_octree(node->children[i]);
        }
    }
    
    free(node);
}

// Calculate boundaries for the octree
void calculate_boundaries(Body* bodies, int n, double* x_min, double* x_max, 
                         double* y_min, double* y_max, double* z_min, double* z_max) {
    *x_min = *y_min = *z_min = 1e20;
    *x_max = *y_max = *z_max = -1e20;
    
    for (int i = 0; i < n; i++) {
        if (bodies[i].x < *x_min) *x_min = bodies[i].x;
        if (bodies[i].x > *x_max) *x_max = bodies[i].x;
        if (bodies[i].y < *y_min) *y_min = bodies[i].y;
        if (bodies[i].y > *y_max) *y_max = bodies[i].y;
        if (bodies[i].z < *z_min) *z_min = bodies[i].z;
        if (bodies[i].z > *z_max) *z_max = bodies[i].z;
    }
    
    // Add some padding
    double padding = 0.01 * (*x_max - *x_min);
    *x_min -= padding;
    *x_max += padding;
    *y_min -= padding;
    *y_max += padding;
    *z_min -= padding;
    *z_max += padding;
    
    // Make sure boundaries are cubic for octree
    double max_range = fmax(*x_max - *x_min, fmax(*y_max - *y_min, *z_max - *z_min));
    double x_center = (*x_min + *x_max) / 2.0;
    double y_center = (*y_min + *y_max) / 2.0;
    double z_center = (*z_min + *z_max) / 2.0;
    
    *x_min = x_center - max_range / 2.0;
    *x_max = x_center + max_range / 2.0;
    *y_min = y_center - max_range / 2.0;
    *y_max = y_center + max_range / 2.0;
    *z_min = z_center - max_range / 2.0;
    *z_max = z_center + max_range / 2.0;
}

// Build the octree for a set of bodies
OctNode* build_octree(Body* bodies, int n) {
    double x_min, x_max, y_min, y_max, z_min, z_max;
    calculate_boundaries(bodies, n, &x_min, &x_max, &y_min, &y_max, &z_min, &z_max);
    
    OctNode* root = create_node(x_min, x_max, y_min, y_max, z_min, z_max);
    
    for (int i = 0; i < n; i++) {
        insert_body(root, bodies, i);
    }
    
    return root;
}

// Compute force on body using Barnes-Hut approximation
void compute_force(OctNode* node, Body* bodies, int body_index, double* fx, double* fy, double* fz) {
    // Skip if node is empty or contains the body itself
    if (node->num_bodies == 0 || (node->is_leaf && node->body_index == body_index)) {
        return;
    }
    
    double dx = node->center_x - bodies[body_index].x;
    double dy = node->center_y - bodies[body_index].y;
    double dz = node->center_z - bodies[body_index].z;
    double distance = sqrt(dx*dx + dy*dy + dz*dz);
    
    // If node is a leaf with a single body, calculate direct force
    if (node->is_leaf && node->body_index != -1) {
        // Skip if bodies are too close (to avoid numerical instability)
        if (distance < 1e-10) return;
        
        double distance_cubed = distance * distance * distance;
        double force = G * bodies[body_index].mass * bodies[node->body_index].mass / distance_cubed;
        
        *fx += force * dx;
        *fy += force * dy;
        *fz += force * dz;
        return;
    }
    
    // Check if node is far enough for approximation
    double s = node_width(node);
    if (s / distance < THETA) {
        // Use approximation - treat node as a single body at center of mass
        // Skip if too close
        if (distance < 1e-10) return;
        
        double distance_cubed = distance * distance * distance;
        double force = G * bodies[body_index].mass * node->total_mass / distance_cubed;
        
        *fx += force * dx;
        *fy += force * dy;
        *fz += force * dz;
    } else {
        // Node is too close, recurse into children
        for (int i = 0; i < 8; i++) {
            if (node->children[i]) {
                compute_force(node->children[i], bodies, body_index, fx, fy, fz);
            }
        }
    }
}

// Update velocities of all bodies using the Barnes-Hut approximation
void update_velocities_barnes_hut(Body* bodies, int n, double dt, OctNode* root) {
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; i++) {
        double fx = 0.0, fy = 0.0, fz = 0.0;
        
        // Compute force using Barnes-Hut approximation
        compute_force(root, bodies, i, &fx, &fy, &fz);
        
        // Update velocities (F = ma => a = F/m)
        bodies[i].vx += fx / bodies[i].mass * dt;
        bodies[i].vy += fy / bodies[i].mass * dt;
        bodies[i].vz += fz / bodies[i].mass * dt;
    }
}

// Update positions of all bodies
void update_positions(Body* bodies, int n, double dt) {
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
        bodies[i].z += bodies[i].vz * dt;
    }
}

// Calculate total energy of the system (kinetic + potential)
// double calculate_energy(Body* bodies, int n) {
//     double total_energy = 0.0;
    
//     #pragma omp parallel
//     {
//         double energy = 0.0;
        
//         // Calculate kinetic and potential energy
//         #pragma omp for schedule(dynamic)
//         for (int i = 0; i < n; i++) {
//             // Kinetic energy
//             double v_squared = bodies[i].vx * bodies[i].vx + 
//                                bodies[i].vy * bodies[i].vy + 
//                                bodies[i].vz * bodies[i].vz;
//             energy += 0.5 * bodies[i].mass * v_squared;
            
//             // Potential energy with all other bodies
//             for (int j = i + 1; j < n; j++) {
//                 double distance = vector_distance(
//                     bodies[i].x, bodies[i].y, bodies[i].z,
//                     bodies[j].x, bodies[j].y, bodies[j].z
//                 );
//                 energy -= G * bodies[i].mass * bodies[j].mass / distance;
//             }
//         }
        
//         // Sum up the energy from all threads
//         #pragma omp critical
//         {
//             total_energy += energy;
//         }
//     }
    
//     return total_energy;
// }

// Kick-step integration with Barnes-Hut approximation
void simulate_step_barnes_hut(Body* bodies, int n, double dt) {
    // Build octree
    OctNode* root = build_octree(bodies, n);
    
    // Update velocities using Barnes-Hut approximation
    update_velocities_barnes_hut(bodies, n, dt, root);
    
    // Update positions
    update_positions(bodies, n, dt);
    
    // Free octree
    free_octree(root);
}

// Initialize system with central body and orbiting bodies
void initialize_system(Body* bodies, int n_bodies, double central_mass, double orbit_radius_min, double orbit_radius_max) {
    if (n_bodies < 1) return;
    
    // Set central body at origin
    bodies[0].x = 0.0;
    bodies[0].y = 0.0;
    bodies[0].z = 0.0;
    bodies[0].vx = 0.0;
    bodies[0].vy = 0.0;
    bodies[0].vz = 0.0;
    bodies[0].mass = central_mass;
    
    // Use thread-safe random number generation
    unsigned int seed = time(NULL);
    
    // Create orbiting bodies in parallel
    #pragma omp parallel
    {
        // Each thread needs its own random seed
        unsigned int local_seed = seed + omp_get_thread_num();
        
        #pragma omp for
        for (int i = 1; i < n_bodies; i++) {
            // Thread-safe random number generation
            double phi = ((double)rand_r(&local_seed) / RAND_MAX) * 2.0 * M_PI;
            double theta = ((double)rand_r(&local_seed) / RAND_MAX) * M_PI - M_PI/2.0;
            
            // Random distance within range
            double radius = orbit_radius_min + ((double)rand_r(&local_seed) / RAND_MAX) * (orbit_radius_max - orbit_radius_min);
            
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
}

void print_system_stats(Body* bodies, int n) {
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
    const int N = 100000 + 1;  // 1 central + 1 million orbiting bodies
    const int STEPS = 10;
    const double dt = 0.01 * YEAR;  // Timestep of 0.01 years
    
    // OpenMP setup
    int num_threads = omp_get_max_threads();
    printf("Starting N-body simulation with %d bodies for %d steps using %d threads\n", 
           N, STEPS, num_threads);
    printf("Using Barnes-Hut approximation with theta = %.2f\n", THETA);
    printf("Allocating memory...\n");
    
    // Allocate memory for bodies
    Body* bodies = (Body*)malloc(N * sizeof(Body));
    if (!bodies) {
        printf("Memory allocation failed!\n");
        return 1;
    }
    
    printf("Initializing system...\n");
    
    // Initialize system with Sun-like central body and Earth-like orbits
    initialize_system(bodies, N, SOLAR_MASS, 0.8 * AU, 5.0 * AU);
    
    print_system_stats(bodies, N);
    
    // Calculate initial energy
    // printf("Calculating initial energy...\n");
    // double initial_energy = calculate_energy(bodies, N);
    // printf("Initial energy: %.10e\n", initial_energy);
    
    printf("Starting simulation...\n");
    double start_time = omp_get_wtime();
    
    // Run simulation
    for (int step = 0; step < STEPS; step++) {
        if (step % 100 == 0) {
            printf("Step %d of %d (%.1f%%)\n", step, STEPS, 100.0 * step / STEPS);
        }
        simulate_step_barnes_hut(bodies, N, dt);
    }
    
    double end_time = omp_get_wtime();
    double elapsed_time = end_time - start_time;
    printf("Simulation completed in %.2f seconds\n", elapsed_time);
    
    // Calculate final energy
    // printf("Calculating final energy...\n");
    // double final_energy = calculate_energy(bodies, N);
    // printf("Final energy: %.10e\n", final_energy);
    
    // // Calculate energy error
    // double energy_error = (final_energy - initial_energy) / initial_energy;
    // printf("Relative energy error: %.10e (%.10f%%)\n", energy_error, energy_error * 100.0);
    
    print_system_stats(bodies, N);
    
    // Free memory
    free(bodies);
    
    return 0;
}