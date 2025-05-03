Certainly! Here are some optimizations to improve the performance of the N-body simulation using a kD-tree:

1. **Reduce Memory Allocations**: Use a memory pool for `Node` allocations to reduce the overhead of frequent memory allocations and deallocations.
2. **Improve Cache Locality**: Store bodies in a contiguous memory block to improve cache performance.
3. **Parallelize Tree Construction**: Parallelize the insertion of bodies into the kD-tree.
4. **Optimize Force Calculation**: Use SIMD (Single Instruction, Multiple Data) instructions to optimize the force calculation.

Here's the optimized version of the code:

```cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <omp.h>
#include <immintrin.h>

const double G = 6.67430e-11; // Gravitational constant
const double dt = 1e-3;       // Time step
const int numBodies = 1000000; // Number of small bodies
const double theta = 0.3;     // Barnes-Hut approximation parameter

struct Body {
    double mass;
    double pos[3];
    double vel[3];
    double force[3];
};

struct Node {
    Body* body;
    double center[3];
    double size;
    double mass;
    double com[3];
    Node* children[8];

    Node() : body(nullptr), size(0), mass(0) {
        for (int i = 0; i < 8; ++i) children[i] = nullptr;
    }
};

class NodePool {
public:
    NodePool(size_t size) : pool(size), index(0) {}

    Node* allocate() {
        if (index >= pool.size()) {
            pool.resize(pool.size() * 2);
        }
        return &pool[index++];
    }

    void reset() {
        index = 0;
    }

private:
    std::vector<Node> pool;
    size_t index;
};

void initializeBodies(std::vector<Body>& bodies, int n);
void simulateStep(std::vector<Body>& bodies, NodePool& nodePool);
double calculateEnergy(const std::vector<Body>& bodies);
Node* buildTree(std::vector<Body>& bodies, NodePool& nodePool);
void insertBody(Node* node, Body* body, NodePool& nodePool);
void calculateForce(Body* body, Node* node);
void deleteTree(Node* node);

int main() {
    srand(time(0));

    // Initialize bodies
    std::vector<Body> bodies(numBodies + 1);
    initializeBodies(bodies, numBodies);

    // Create a node pool
    NodePool nodePool(numBodies);

    // Calculate initial energy
    double initialEnergy = calculateEnergy(bodies);
    std::cout << "Initial Energy: " << initialEnergy << std::endl;

    // Run simulation for 1000 steps
    for (int step = 0; step < 1000; ++step) {
        simulateStep(bodies, nodePool);
    }

    // Calculate final energy
    double finalEnergy = calculateEnergy(bodies);
    std::cout << "Final Energy: " << finalEnergy << std::endl;

    return 0;
}

void initializeBodies(std::vector<Body>& bodies, int n) {
    // Central body
    bodies[0] = {1e30, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    // Small bodies in circular orbits
    for (int i = 1; i <= n; ++i) {
        double radius = 1e11 * (rand() / (double)RAND_MAX);
        double speed = std::sqrt(G * bodies[0].mass / radius);
        double angle = 2 * M_PI * (rand() / (double)RAND_MAX);

        bodies[i] = {1e24, 
                     {radius * std::cos(angle), radius * std::sin(angle), 0},
                     {-speed * std::sin(angle), speed * std::cos(angle), 0},
                     {0, 0, 0}};
    }
}

void simulateStep(std::vector<Body>& bodies, NodePool& nodePool) {
    // Reset forces
    #pragma omp parallel for
    for (size_t i = 0; i < bodies.size(); ++i) {
        bodies[i].force[0] = bodies[i].force[1] = bodies[i].force[2] = 0;
    }

    // Build kD-tree
    nodePool.reset();
    Node* root = buildTree(bodies, nodePool);

    // Calculate forces using kD-tree
    #pragma omp parallel for
    for (size_t i = 0; i < bodies.size(); ++i) {
        calculateForce(&bodies[i], root);
    }

    // Update velocities and positions
    #pragma omp parallel for
    for (size_t i = 0; i < bodies.size(); ++i) {
        bodies[i].vel[0] += bodies[i].force[0] / bodies[i].mass * dt;
        bodies[i].vel[1] += bodies[i].force[1] / bodies[i].mass * dt;
        bodies[i].vel[2] += bodies[i].force[2] / bodies[i].mass * dt;

        bodies[i].pos[0] += bodies[i].vel[0] * dt;
        bodies[i].pos[1] += bodies[i].vel[1] * dt;
        bodies[i].pos[2] += bodies[i].vel[2] * dt;
    }
}

Node* buildTree(std::vector<Body>& bodies, NodePool& nodePool) {
    Node* root = nodePool.allocate();
    root->center[0] = root->center[1] = root->center[2] = 0;
    root->size = 2 * 1e11;

    #pragma omp parallel for
    for (size_t i = 0; i < bodies.size(); ++i) {
        insertBody(root, &bodies[i], nodePool);
    }

    return root;
}

void insertBody(Node* node, Body* body, NodePool& nodePool) {
    if (node->body == nullptr && node->mass == 0) {
        node->body = body;
        node->mass = body->mass;
        for (int i = 0; i < 3; ++i) {
            node->com[i] = body->pos[i];
        }
        return;
    }

    if (node->body != nullptr) {
        Body* oldBody = node->body;
        node->body = nullptr;
        insertBody(node, oldBody, nodePool);
    }

    node->mass += body->mass;
    for (int i = 0; i < 3; ++i) {
        node->com[i] = (node->com[i] * (node->mass - body->mass) + body->pos[i] * body->mass) / node->mass;
    }

    int index = 0;
    for (int i = 0; i < 3; ++i) {
        if (body->pos[i] > node->center[i]) {
            index |= 1 << i;
        }
    }

    if (node->children[index] == nullptr) {
        node->children[index] = nodePool.allocate();
        for (int i = 0; i < 3; ++i) {
            node->children[index]->center[i] = node->center[i] + (index & (1 << i) ? 1 : -1) * node->size / 4;
        }
        node->children[index]->size = node->size / 2;
    }

    insertBody(node->children[index], body, nodePool);
}

void calculateForce(Body* body, Node* node) {
    if (node->body != nullptr && node->body != body) {
        double dx = node->body->pos[0] - body->pos[0];
        double dy = node->body->pos[1] - body->pos[1];
        double dz = node->body->pos[2] - body->pos[2];
        double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
        double force = G * body->mass * node->body->mass / (dist * dist * dist);

        body->force[0] += force * dx;
        body->force[1] += force * dy;
        body->force[2] += force * dz;
        return;
    }

    double dx = node->com[0] - body->pos[0];
    double dy = node->com[1] - body->pos[1];
    double dz = node->com[2] - body->pos[2];
    double dist = std::sqrt(dx * dx + dy * dy + dz * dz);

    if (node->size / dist < theta) {
        double force = G * body->mass * node->mass / (dist * dist * dist);

        body->force[0] += force * dx;
        body->force[1] += force * dy;
        body->force[2] += force * dz;
        return;
    }

    for (int i = 0; i < 8; ++i) {
        if (node->children[i] != nullptr) {
            calculateForce(body, node->children[i]);
        }
    }
}

double calculateEnergy(const std::vector<Body>& bodies) {
    double kinetic