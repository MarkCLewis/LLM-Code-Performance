#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <omp.h>

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

void initializeBodies(std::vector<Body>& bodies, int n);
void simulateStep(std::vector<Body>& bodies);
double calculateEnergy(const std::vector<Body>& bodies);
Node* buildTree(std::vector<Body>& bodies);
void insertBody(Node* node, Body* body);
void calculateForce(Body* body, Node* node);
void deleteTree(Node* node);

int main() {
    srand(time(0));

    // Initialize bodies
    std::vector<Body> bodies(numBodies + 1);
    initializeBodies(bodies, numBodies);

    // Calculate initial energy
    double initialEnergy = calculateEnergy(bodies);
    std::cout << "Initial Energy: " << initialEnergy << std::endl;

    // Run simulation for 1000 steps
    for (int step = 0; step < 1000; ++step) {
        simulateStep(bodies);
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

void simulateStep(std::vector<Body>& bodies) {
    // Reset forces
    #pragma omp parallel for
    for (size_t i = 0; i < bodies.size(); ++i) {
        bodies[i].force[0] = bodies[i].force[1] = bodies[i].force[2] = 0;
    }

    // Build kD-tree
    Node* root = buildTree(bodies);

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

    // Clean up kD-tree
    deleteTree(root);
}

Node* buildTree(std::vector<Body>& bodies) {
    Node* root = new Node();
    root->center[0] = root->center[1] = root->center[2] = 0;
    root->size = 2 * 1e11;

    for (auto& body : bodies) {
        insertBody(root, &body);
    }

    return root;
}

void insertBody(Node* node, Body* body) {
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
        insertBody(node, oldBody);
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
        node->children[index] = new Node();
        for (int i = 0; i < 3; ++i) {
            node->children[index]->center[i] = node->center[i] + (index & (1 << i) ? 1 : -1) * node->size / 4;
        }
        node->children[index]->size = node->size / 2;
    }

    insertBody(node->children[index], body);
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

void deleteTree(Node* node) {
    for (int i = 0; i < 8; ++i) {
        if (node->children[i] != nullptr) {
            deleteTree(node->children[i]);
        }
    }
    delete node;
}

double calculateEnergy(const std::vector<Body>& bodies) {
    double kineticEnergy = 0;
    double potentialEnergy = 0;

    #pragma omp parallel for reduction(+:kineticEnergy, potentialEnergy)
    for (size_t i = 0; i < bodies.size(); ++i) {
        kineticEnergy += 0.5 * bodies[i].mass * (bodies[i].vel[0]*bodies[i].vel[0] + bodies[i].vel[1]*bodies[i].vel[1] + bodies[i].vel[2]*bodies[i].vel[2]);

        for (size_t j = i + 1; j < bodies.size(); ++j) {
            double dx = bodies[j].pos[0] - bodies[i].pos[0];
            double dy = bodies[j].pos[1] - bodies[i].pos[1];
            double dz = bodies[j].pos[2] - bodies[i].pos[2];
            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);