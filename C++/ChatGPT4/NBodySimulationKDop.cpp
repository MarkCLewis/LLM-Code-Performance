#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <omp.h>
#include <memory>

const double G = 6.67430e-11;
const double DT = 1.0;
const double PI = 3.141592653589793;
const double THETA = 0.3;

struct Body {
    double position[3];
    double velocity[3];
    double mass;
};

struct KDNode {
    double mass;
    double center_of_mass[3];
    double min_bounds[3];
    double max_bounds[3];
    KDNode* left;
    KDNode* right;

    KDNode() : mass(0), left(nullptr), right(nullptr) {
        for (int i = 0; i < 3; ++i) {
            center_of_mass[i] = min_bounds[i] = max_bounds[i] = 0;
        }
    }
};

// Optimized KDTree Construction with Index Partitioning
KDNode* build_tree(std::vector<Body*>& bodies, int start, int end, int depth = 0) {
    if (start >= end) return nullptr;

    KDNode* node = new KDNode();
    
    if (end - start == 1) {
        node->mass = bodies[start]->mass;
        for (int i = 0; i < 3; ++i) {
            node->center_of_mass[i] = bodies[start]->position[i];
        }
        return node;
    }

    int axis = depth % 3;
    std::nth_element(bodies.begin() + start, bodies.begin() + (start + end) / 2, bodies.begin() + end,
                     [axis](Body* a, Body* b) { return a->position[axis] < b->position[axis]; });

    int mid = (start + end) / 2;
    node->left = build_tree(bodies, start, mid, depth + 1);
    node->right = build_tree(bodies, mid, end, depth + 1);

    double total_mass = 0;
    double com[3] = {0, 0, 0};

    for (int i = start; i < end; i++) {
        total_mass += bodies[i]->mass;
        for (int j = 0; j < 3; ++j) {
            com[j] += bodies[i]->position[j] * bodies[i]->mass;
        }
    }

    node->mass = total_mass;
    for (int i = 0; i < 3; ++i) {
        node->center_of_mass[i] = com[i] / total_mass;
    }
    return node;
}

void compute_force(Body& body, KDNode* node, double force[3]) {
    if (!node) return;

    double dx = node->center_of_mass[0] - body.position[0];
    double dy = node->center_of_mass[1] - body.position[1];
    double dz = node->center_of_mass[2] - body.position[2];
    double dist_sq = dx * dx + dy * dy + dz * dz + 1e-9; // Avoid div by zero
    double dist = sqrt(dist_sq);

    if (!node->left && !node->right || ((node->max_bounds[0] - node->min_bounds[0]) / dist) < THETA) {
        double factor = G * node->mass / (dist_sq * dist);
        force[0] += factor * dx;
        force[1] += factor * dy;
        force[2] += factor * dz;
    } else {
        if (node->left) compute_force(body, node->left, force);
        if (node->right) compute_force(body, node->right, force);
    }
}

void update_velocities(std::vector<Body>& bodies, KDNode* root) {
    #pragma omp parallel for
    for (size_t i = 0; i < bodies.size(); ++i) {
        double force[3] = {0, 0, 0};
        compute_force(bodies[i], root, force);
        
        for (int k = 0; k < 3; ++k) {
            bodies[i].velocity[k] += force[k] * DT / bodies[i].mass;
        }
    }
}

void update_positions(std::vector<Body>& bodies) {
    #pragma omp parallel for
    for (auto& body : bodies) {
        for (int i = 0; i < 3; ++i) {
            body.position[i] += body.velocity[i] * DT;
        }
    }
}

int main() {
    size_t n_bodies = 1'000'000;
    double central_mass = 1.989e30;
    
    std::vector<Body> bodies(n_bodies);
    std::vector<Body*> body_ptrs(n_bodies);
    
    #pragma omp parallel for
    for (size_t i = 0; i < n_bodies; ++i) {
        bodies[i] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0};
        body_ptrs[i] = &bodies[i];
    }
    
    KDNode* root = build_tree(body_ptrs, 0, n_bodies);
    
    for (int step = 0; step < 1000; ++step) {
        update_velocities(bodies, root);
        update_positions(bodies);
        
        delete root; // Free old tree
        root = build_tree(body_ptrs, 0, n_bodies); // Rebuild efficiently
    }
    
    delete root; // Free last tree
    return 0;
}