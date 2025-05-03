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
    std::vector<Body*> bodies;
    std::unique_ptr<KDNode> left;
    std::unique_ptr<KDNode> right;

    KDNode() : mass(0) {
        for (int i = 0; i < 3; ++i) {
            center_of_mass[i] = 0;
            min_bounds[i] = 0;
            max_bounds[i] = 0;
        }
    }
};

std::unique_ptr<KDNode> build_tree(std::vector<Body*>& bodies, int depth = 0) {
    if (bodies.empty()) return nullptr;
    
    auto node = std::make_unique<KDNode>();
    node->bodies = bodies;
    
    if (bodies.size() == 1) {
        node->mass = bodies[0]->mass;
        for (int i = 0; i < 3; ++i) {
            node->center_of_mass[i] = bodies[0]->position[i];
        }
        return node;
    }
    
    int axis = depth % 3;
    std::sort(bodies.begin(), bodies.end(), [axis](Body* a, Body* b) {
        return a->position[axis] < b->position[axis];
    });
    
    size_t mid = bodies.size() / 2;
    std::vector<Body*> left_bodies(bodies.begin(), bodies.begin() + mid);
    std::vector<Body*> right_bodies(bodies.begin() + mid, bodies.end());
    
    node->left = build_tree(left_bodies, depth + 1);
    node->right = build_tree(right_bodies, depth + 1);
    
    double total_mass = 0;
    double com[3] = {0, 0, 0};
    for (Body* body : bodies) {
        total_mass += body->mass;
        for (int i = 0; i < 3; ++i) {
            com[i] += body->position[i] * body->mass;
        }
    }
    node->mass = total_mass;
    for (int i = 0; i < 3; ++i) {
        node->center_of_mass[i] = com[i] / total_mass;
    }
    return node;
}

void compute_force(Body& body, KDNode* node, double force[3]) {
    if (!node || node->bodies.empty()) return;
    
    double dx = node->center_of_mass[0] - body.position[0];
    double dy = node->center_of_mass[1] - body.position[1];
    double dz = node->center_of_mass[2] - body.position[2];
    double dist_sq = dx * dx + dy * dy + dz * dz;
    double dist = sqrt(dist_sq);
    
    if (node->bodies.size() == 1 || (dist > 0 && ((node->max_bounds[0] - node->min_bounds[0]) / dist) < THETA)) {
        double factor = G * node->mass / (dist_sq * dist);
        force[0] += factor * dx;
        force[1] += factor * dy;
        force[2] += factor * dz;
    } else {
        compute_force(body, node->left.get(), force);
        compute_force(body, node->right.get(), force);
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
    for (size_t i = 0; i < n_bodies; ++i) {
        bodies[i] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0};
    }
    
    std::vector<Body*> body_ptrs;
    for (auto& body : bodies) {
        body_ptrs.push_back(&body);
    }
    
    auto root = build_tree(body_ptrs);
    
    for (int step = 0; step < 1000; ++step) {
        update_velocities(bodies, root.get());
        update_positions(bodies);
        root = build_tree(body_ptrs);
    }
    
    return 0;
}
