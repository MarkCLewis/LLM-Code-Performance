import java.util.Random;

public class NBodySimulationPar {

    static final int N_BODIES = 10_000;
    static final int STEPS = 100;
    static final double G = 6.67430e-11;
    static final double DT = 1.0;

    static class Body {
        double x, y, z;
        double vx, vy, vz;
        double ax, ay, az;
        double mass;

        public Body(double x, double y, double z, double vx, double vy, double vz, double mass) {
            this.x = x; this.y = y; this.z = z;
            this.vx = vx; this.vy = vy; this.vz = vz;
            this.mass = mass;
        }
    }

    static Body[] bodies;

    public static void initializeSystem(int n) {
        bodies = new Body[n + 1];
        Random rand = new Random(42);

        // Central body
        double centralMass = 1e20;
        bodies[0] = new Body(0, 0, 0, 0, 0, 0, centralMass);

        double radius = 1e7;

        for (int i = 1; i <= n; i++) {
            double angle = 2 * Math.PI * i / n;
            double r = radius * (1.0 + 0.1 * rand.nextDouble());

            double x = r * Math.cos(angle);
            double y = r * Math.sin(angle);
            double z = 0;

            double v = Math.sqrt(G * centralMass / r);
            double vx = -v * Math.sin(angle);
            double vy = v * Math.cos(angle);
            double vz = 0;

            bodies[i] = new Body(x, y, z, vx, vy, vz, 1.0);
        }
    }

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

    public static void kickStep() {
        for (Body b : bodies) {
            b.vx += b.ax * DT;
            b.vy += b.ay * DT;
            b.vz += b.az * DT;

            b.x += b.vx * DT;
            b.y += b.vy * DT;
            b.z += b.vz * DT;
        }
    }

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
    
    public static void main(String[] args) {
        long start = System.currentTimeMillis();
        System.out.println("Initializing...");
        initializeSystem(N_BODIES);

        System.out.println("Computing initial energy...");
        double initialEnergy = computeEnergy();
        System.out.printf("Initial total energy: %.6e\n", initialEnergy);

        for (int step = 0; step < STEPS; step++) {
            computeForces();
            kickStep();

            if (step % 100 == 0) {
                System.out.println("Step " + step + " complete");
            }
        }

        System.out.println("Computing final energy...");
        double finalEnergy = computeEnergy();
        System.out.printf("Final total energy: %.6e\n", finalEnergy);
        System.out.printf("Energy difference: %.6e\n", Math.abs(finalEnergy - initialEnergy));

        long end = System.currentTimeMillis();
        System.out.printf("Simulation completed in %.2f seconds.\n", (end - start) / 1000.0);
    }
}


// Below are only the modified sections from the previous code (to keep it concise):
// Replace computeForces() with:


// Replace computeEnergy() with:


// You must run this with adequate heap space, like:
// java -Xmx16g NBodySimulation
// Java Parallel Streams scale well up to a few dozen cores. Beyond that, you'll want finer control (e.g., with ForkJoinPool or ExecutorService).
// This is still O(N²). For real performance at scale, consider a Barnes-Hut tree or Fast Multipole Method.
