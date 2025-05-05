// Below are only the modified sections from the previous code (to keep it concise):
// Replace computeForces() with:
public static void computeForces() {
    for (Body b : bodies) {
        b.ax = b.ay = b.az = 0;
    }

    int n = bodies.length;

    java.util.Arrays.stream(bodies).parallel().forEach(bi -> {
        for (Body bj : bodies) {
            if (bi == bj) continue;

            double dx = bj.x - bi.x;
            double dy = bj.y - bi.y;
            double dz = bj.z - bi.z;

            double distSqr = dx * dx + dy * dy + dz * dz + 1e-10;
            double dist = Math.sqrt(distSqr);
            double force = G * bj.mass / (distSqr * dist);

            bi.ax += dx * force;
            bi.ay += dy * force;
            bi.az += dz * force;
        }
    });
}
// This parallel version is not symmetric, so it's about twice as expensive as the original i < j version, but parallel and simpler to implement. Further optimization may require fine-grained task splitting.

// Replace computeEnergy() with:
public static double computeEnergy() {
    double kinetic = java.util.Arrays.stream(bodies).parallel().mapToDouble(b -> {
        double v2 = b.vx * b.vx + b.vy * b.vy + b.vz * b.vz;
        return 0.5 * b.mass * v2;
    }).sum();

    int n = bodies.length;
    double potential = java.util.stream.IntStream.range(0, n).parallel().mapToDouble(i -> {
        Body bi = bodies[i];
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            Body bj = bodies[j];
            double dx = bj.x - bi.x;
            double dy = bj.y - bi.y;
            double dz = bj.z - bi.z;
            double dist = Math.sqrt(dx * dx + dy * dy + dz * dz + 1e-10);
            sum -= G * bi.mass * bj.mass / dist;
        }
        return sum;
    }).sum();

    return kinetic + potential;
}

// You must run this with adequate heap space, like:
// java -Xmx16g NBodySimulation
// Java Parallel Streams scale well up to a few dozen cores. Beyond that, you'll want finer control (e.g., with ForkJoinPool or ExecutorService).
// This is still O(N²). For real performance at scale, consider a Barnes-Hut tree or Fast Multipole Method.
