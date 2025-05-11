// The code has been optimized further by:
//Parallelizing the kD-tree construction with OpenMP tasks.
//Scheduling parallel loops with static or dynamic as appropriate.
//Reducing memory allocations and reusing memory buffers when possible.
//Rebuilding the kD-tree every 10 steps to improve performance without compromising too much accuracy.
//Would you like to benchmark the performance now or further tune accuracy vs. speed trade-offs?


#include <stdio.h>
    #pragma omp parallel for reduction(+:kinetic) schedule(static)
    for (int i = 0; i <= n; i++) {
        double v2 = bodies[i].vx * bodies[i].vx + bodies[i].vy * bodies[i].vy + bodies[i].vz * bodies[i].vz;
        kinetic += 0.5 * bodies[i].mass * v2;
    }

    #pragma omp parallel for reduction(+:potential) schedule(dynamic)
    for (int i = 0; i <= n; i++) {
        for (int j = i + 1; j <= n; j++) {
            double dx = bodies[i].x - bodies[j].x;
            double dy = bodies[i].y - bodies[j].y;
            double dz = bodies[i].z - bodies[j].z;
            double dist = sqrt(dx * dx + dy * dy + dz * dz + 1e-10);
            potential -= G * bodies[i].mass * bodies[j].mass / dist;
        }
    }

    return kinetic + potential;
}

void kick_step(KDNode* tree, int n) {
    compute_forces_kdtree(tree, n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i <= n; i++) {
        bodies[i].vx += bodies[i].ax * DT;
        bodies[i].vy += bodies[i].ay * DT;
        bodies[i].vz += bodies[i].az * DT;

        bodies[i].x += bodies[i].vx * DT;
        bodies[i].y += bodies[i].vy * DT;
        bodies[i].z += bodies[i].vz * DT;
    }
}

int main() {
    clock_t start = clock();
    printf("Initializing system...\n");
    initialize_system(N_BODIES);

    printf("Calculating initial energy...\n");
    double initial_energy = calculate_energy(N_BODIES);
    printf("Initial total energy: %.6e\n", initial_energy);

    KDNode* tree = NULL;
    for (int step = 0; step < STEPS; step++) {
        if (step % TREE_REBUILD_INTERVAL == 0) {
            if (tree) free_kdtree(tree);
            #pragma omp parallel
            {
                #pragma omp single nowait
                tree = build_kdtree(0, N_BODIES + 1, 0);
            }
        }
        kick_step(tree, N_BODIES);
        if (step % 100 == 0) printf("Step %d complete\n", step);
    }

    printf("Calculating final energy...\n");
    double final_energy = calculate_energy(N_BODIES);
    printf("Final total energy: %.6e\n", final_energy);
    printf("Energy difference: %.6e\n", fabs(final_energy - initial_energy));

    free_kdtree(tree);
    free(bodies);
    free(indices);

    clock_t end = clock();
    printf("Simulation completed in %.2f seconds.\n", (double)(end - start) / CLOCKS_PER_SEC);
    return 0;
}