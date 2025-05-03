//Okay, let's implement some of the more straightforward optimizations that can be done without a complete architectural overhaul. We'll focus on:

//Optimized kD-tree Construction (using an index array to avoid body swaps).
//Experimenting with OpenMP scheduling in the force calculation.
//Here's the modified code:

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

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

// Define the structure for a kD-tree node
typedef struct KDNode {
    int body_index;     // Index of the body in this leaf node (-1 if internal)
    double center_of_mass[3];
    double total_mass;
    double min[3];      // Minimum bounds of the node
    double max[3];      // Maximum bounds of the node
    struct KDNode *left;
    struct KDNode *right;
} KDNode;

// Gravitational constant
const double G = 6.67430e-11;
const double THETA = 0.3; // Barnes-Hut opening angle

// Function to calculate the distance squared between two points
double distance_squared_points(const double p1[3], const double p2[3]) {
    double dx = p1[0] - p2[0];
    double dy = p1[1] - p2[1];
    double dz = p1[2] - p2[2];
    return dx * dx + dy * dy + dz * dz;
}

// Function to calculate the distance between two points
double distance_points(const double p1[3], const double p2[3]) {
    return sqrt(distance_squared_points(p1, p2));
}

// Function to calculate the force on body i due to a mass element (node or body)
void calculate_force_element(const Body *body_i, const double mass_element_cm[3], double mass_element_mass, double force[3]) {
    double r2 = distance_squared_points(body_i->position, mass_element_cm);
    if (r2 > 1e-9) { // Avoid division by zero
        double r = sqrt(r2);
        double magnitude = (G * body_i->mass * mass_element_mass) / r2;
        force[0] += magnitude * (mass_element_cm[0] - body_i->position[0]) / r;
        force[1] += magnitude * (mass_element_cm[1] - body_i->position[1]) / r;
        force[2] += magnitude * (mass_element_cm[2] - body_i->position[2]) / r;
    }
}

// Function to calculate the total energy of the system (parallelized)
double calculate_total_energy(const System *system) {
    double kinetic_energy = 0.0;
    double potential_energy = 0.0;
    int n = system->num_bodies;

    #pragma omp parallel for reduction(+:kinetic_energy)
    for (int i = 0; i < n; i++) {
        kinetic_energy += 0.5 * system->bodies[i].mass * (system->bodies[i].velocity[0] * system->bodies[i].velocity[0] +
                                                          system->bodies[i].velocity[1] * system->bodies[i].velocity[1] +
                                                          system->bodies[i].velocity[2] * system->bodies[i].velocity[2]);
    }

    #pragma omp parallel for reduction(+:potential_energy)
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            potential_energy -= (G * system->bodies[i].mass * system->bodies[j].mass) / sqrt(distance_squared_points(system->bodies[i].position, system->bodies[j].position));
        }
    }

    return kinetic_energy + potential_energy;
}

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

// Function to calculate the bounding box of a set of bodies (using indices)
void calculate_bounds_indexed(const System *system, const int *indices, int start, int end, double min[3], double max[3]) {
    if (start > end) return;
    for (int d = 0; d < 3; d++) {
        min[d] = system->bodies[indices[start]].position[d];
        max[d] = system->bodies[indices[start]].position[d];
    }
    for (int i = start + 1; i <= end; i++) {
        for (int d = 0; d < 3; d++) {
            if (system->bodies[indices[i]].position[d] < min[d]) min[d] = system->bodies[indices[i]].position[d];
            if (system->bodies[indices[i]].position[d] > max[d]) max[d] = system->bodies[indices[i]].position[d];
        }
    }
}

// Function to create a new kD-tree node
KDNode* create_kd_node(int body_index) {
    KDNode *node = (KDNode*)malloc(sizeof(KDNode));
    if (node == NULL) {
        fprintf(stderr, "Memory allocation failed for kD-tree node.\n");
        exit(EXIT_FAILURE);
    }
    node->body_index = body_index;
    node->total_mass = 0.0;
    node->center_of_mass[0] = node->center_of_mass[1] = node->center_of_mass[2] = 0.0;
    node->left = NULL;
    node->right = NULL;
    return node;
}

// Recursive function to build the kD-tree (using indices)
KDNode* build_kd_tree_recursive(System *system, int *indices, int start_index, int end_index, double min[3], double max[3], int depth) {
    if (start_index > end_index) {
        return NULL;
    }

    KDNode *node = create_kd_node(-1); // Internal node initially
    for (int d = 0; d < 3; d++) {
        node->min[d] = min[d];
        node->max[d] = max[d];
    }

    if (start_index == end_index) {
        node->body_index = indices[start_index];
        node->total_mass = system->bodies[indices[start_index]].mass;
        for (int d = 0; d < 3; d++) {
            node->center_of_mass[d] = system->bodies[indices[start_index]].position[d];
        }
        return node;
    }

    // Calculate center of mass and total mass for the current node
    node->total_mass = 0.0;
    for (int i = start_index; i <= end_index; i++) {
        int index = indices[i];
        node->total_mass += system->bodies[index].mass;
        for (int d = 0; d < 3; d++) {
            node->center_of_mass[d] += system->bodies[index].mass * system->bodies[index].position[d];
        }
    }
    for (int d = 0; d < 3; d++) {
        node->center_of_mass[d] /= node->total_mass;
    }

    // Determine the splitting dimension
    int split_dim = depth % 3;
    double split_val = (min[split_dim] + max[split_dim]) / 2.0;

    // Partition the indices based on the splitting dimension
    int left_end = start_index - 1;
    int right_start = end_index + 1;
    int i = start_index;
    while (i < right_start) {
        if (system->bodies[indices[i]].position[split_dim] <= split_val) {
            left_end++;
            // Swap indices[i] and indices[left_end]
            int temp = indices[i];
            indices[i] = indices[left_end];
            indices[left_end] = temp;
            i++;
        } else {
            right_start--;
            // Swap indices[i] and indices[right_start]
            int temp = indices[i];
            indices[i] = indices[right_start];
            indices[right_start] = temp;
        }
    }

    // Create new bounds for the left and right children
    double left_max[3], right_min[3];
    for (int d = 0; d < 3; d++) {
        left_max[d] = max[d];
        right_min[d] = min[d];
    }
    left_max[split_dim] = split_val;
    right_min[split_dim] = split_val;

    // Recursively build the left and right subtrees
    node->left = build_kd_tree_recursive(system, indices, start_index, left_end, min, left_max, depth + 1);
    node->right = build_kd_tree_recursive(system, indices, left_end + 1, end_index, right_min, max, depth + 1);

    return node;
}

// Function to build the kD-tree (using indices)
KDNode* build_kd_tree(System *system) {
    int *indices = malloc(sizeof(int) * system->num_bodies);
    if (indices == NULL) {
        fprintf(stderr, "Memory allocation failed for indices.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < system->num_bodies; i++) {
        indices[i] = i;
    }
    double min_bounds[3], max_bounds[3];
    calculate_bounds_indexed(system, indices, 0, system->num_bodies - 1, min_bounds, max_bounds);
    KDNode* root = build_kd_tree_recursive(system, indices, 0, system->num_bodies - 1, min_bounds, max_bounds, 0);
    free(indices);
    return root;
}

// Recursive function to calculate the force on a body using the kD-tree
void calculate_force_kd_tree_recursive(const Body *body_i, const KDNode *node, const System *system, double force[3]) {
    if (node == NULL) {
        return;
    }

    if (node->body_index != -1) { // Leaf node (contains a single body)
        if (node->body_index != (body_i - &system->bodies[0])) { // Don't interact with itself
            calculate_force_element(body_i, node->center_of_mass, node->total_mass, force);
        }
    } else { // Internal node
        double s = node->max[0] - node->min[0]; // Size of the node (assuming cubic-like)
        double d = distance_points(body_i->position, node->center_of_mass);

        if (s / d < THETA || d == 0.0) { // Far enough or at the same location (shouldn't happen)
            calculate_force_element(body_i, node->center_of_mass, node->total_mass, force);
        } else { // Need to go deeper
            calculate_force_kd_tree_recursive(body_i, node->left, system, force);
            calculate_force_kd_tree_recursive(body_i, node->right, system, force);
        }
    }
}

// Function to calculate the force on a body using the kD-tree
void calculate_force_kd_tree(const Body *body_i, const KDNode *root, const System *system, double force[3]) {
    force[0] = force[1] = force[2] = 0.0;
    calculate_force_kd_tree_recursive(body_i, root, system, force);
}

// Function to free the memory allocated for the kD-tree
void free_kd_tree(KDNode *node) {
    if (node == NULL) {
        return;
    }
    free_kd_tree(node->left);
    free_kd_tree(node->right);
    free(node);
}

// Function to perform one kick-step using the kD-tree for force calculation (parallelized)
void kick_step_kd_tree(System *system, double dt) {
    int n = system->num_bodies;
    KDNode *root = build_kd_tree(system);
    double (*forces)[3] = malloc(sizeof(double[n][3]));
    if (forces == NULL) {
        fprintf(stderr, "Memory allocation failed for forces.\n");
        free_kd_tree(root);
        exit(EXIT_FAILURE);
    }

    // Experiment with different OpenMP scheduling
    #pragma omp parallel for schedule(dynamic, 1000)
    for (int i = 0; i < n; i++) {
        forces[i][0] = forces[i][1] = forces[i][2] = 0.0;
        calculate_force_kd_tree(&system->bodies[i], root, system, forces[i]);
    }

    // Update velocities (kick)
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        system->bodies[i].velocity[0] += (forces[i][0] / system->bodies[i].mass) * dt * 0.5;
        system->bodies[i].velocity[1] += (forces[i][1] / system->bodies[i].mass) * dt * 0.5;
        system->bodies[i].velocity[2] += (forces[i][2] / system->bodies[i].mass) * dt * 0.5;
    }

    free(forces);
    free_kd_tree(root);
}

// Function to perform one drift-step (parallelized)
void drift_step(System *system, double dt) {
    int n = system->num_bodies;
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        system->bodies[i].position[0] += system->bodies