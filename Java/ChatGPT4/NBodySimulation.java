import java.util.Random;

public class NBodySimulation {

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
        for (int i = 0; i < n; i++) {
            Body bi = bodies[i];
            for (int j = i + 1; j < n; j++) {
                Body bj = bodies[j];

                double dx = bj.x - bi.x;
                double dy = bj.y - bi.y;
                double dz = bj.z - bi.z;

                double distSqr = dx * dx + dy * dy + dz * dz + 1e-10;
                double dist = Math.sqrt(distSqr);
                double force = G / (distSqr * dist);

                double fx = force * dx;
                double fy = force * dy;
                double fz = force * dz;

                bi.ax += fx * bj.mass;
                bi.ay += fy * bj.mass;
                bi.az += fz * bj.mass;

                bj.ax -= fx * bi.mass;
                bj.ay -= fy * bi.mass;
                bj.az -= fz * bi.mass;
            }
        }
    }

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

    // public static double computeEnergy() {
    //     double kinetic = 0.0;
    //     double potential = 0.0;

    //     for (Body b : bodies) {
    //         double v2 = b.vx * b.vx + b.vy * b.vy + b.vz * b.vz;
    //         kinetic += 0.5 * b.mass * v2;
    //     }

    //     int n = bodies.length;
    //     for (int i = 0; i < n; i++) {
    //         Body bi = bodies[i];
    //         for (int j = i + 1; j < n; j++) {
    //             Body bj = bodies[j];
    //             double dx = bj.x - bi.x;
    //             double dy = bj.y - bi.y;
    //             double dz = bj.z - bi.z;
    //             double dist = Math.sqrt(dx * dx + dy * dy + dz * dz + 1e-10);
    //             potential -= G * bi.mass * bj.mass / dist;
    //         }
    //     }

    //     return kinetic + potential;
    // }

    public static void main(String[] args) {
        long start = System.currentTimeMillis();
        System.out.println("Initializing...");
        initializeSystem(N_BODIES);

        // System.out.println("Computing initial energy...");
        // double initialEnergy = computeEnergy();
        // System.out.printf("Initial total energy: %.6e\n", initialEnergy);

        for (int step = 0; step < STEPS; step++) {
            computeForces();
            kickStep();

            if (step % 100 == 0) {
                System.out.println("Step " + step + " complete");
            }
        }

        // System.out.println("Computing final energy...");
        // double finalEnergy = computeEnergy();
        // System.out.printf("Final total energy: %.6e\n", finalEnergy);
        // System.out.printf("Energy difference: %.6e\n", Math.abs(finalEnergy - initialEnergy));

        long end = System.currentTimeMillis();
        System.out.printf("Simulation completed in %.2f seconds.\n", (end - start) / 1000.0);
    }
}
