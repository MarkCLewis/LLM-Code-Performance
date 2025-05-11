import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

public class NBodySimulationPar {
    private static final double G = 6.67430e-11; // Gravitational constant
    private static final double DT = 3600*24*365*1e-3; // Time step
    private static final int NUM_BODIES = 10000; // Number of bodies

    static class Body {
        double x, y, z;
        double vx, vy, vz;
        double mass;

        Body(double x, double y, double z, double vx, double vy, double vz, double mass) {
            this.x = x;
            this.y = y;
            this.z = z;
            this.vx = vx;
            this.vy = vy;
            this.vz = vz;
            this.mass = mass;
        }
    }

    private static void initializeBodies(List<Body> bodies) {
        bodies.add(new Body(0, 0, 0, 0, 0, 0, 1e30)); // Central body mass

        IntStream.range(1, NUM_BODIES).parallel().forEach(i -> {
            double angle = 2 * Math.PI * i / (NUM_BODIES - 1);
            bodies.add(new Body(
                Math.cos(angle) * 1e11,
                Math.sin(angle) * 1e11,
                0,
                -Math.sin(angle) * Math.sqrt(G * bodies.get(0).mass / 1e11),
                Math.cos(angle) * Math.sqrt(G * bodies.get(0).mass / 1e11),
                0,
                1e24 / NUM_BODIES // Small body mass
            ));
        });
    }

    // private static double calculateEnergy(List<Body> bodies) {
    //     return bodies.parallelStream().mapToDouble(body -> {
    //         double kinetic = 0.5 * body.mass * (body.vx * body.vx + body.vy * body.vy + body.vz * body.vz);
    //         double potential = bodies.stream().filter(other -> body != other).mapToDouble(other -> {
    //             double dx = body.x - other.x;
    //             double dy = body.y - other.y;
    //             double dz = body.z - other.z;
    //             double distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
    //             return -G * body.mass * other.mass / distance;
    //         }).sum();
    //         return kinetic + 0.5 * potential;
    //     }).sum();
    // }

    private static void kickStep(List<Body> bodies) {
        bodies.parallelStream().forEach(body -> {
            double ax = 0.0, ay = 0.0, az = 0.0;
            for (Body other : bodies) {
                if (body != other) {
                    double dx = other.x - body.x;
                    double dy = other.y - body.y;
                    double dz = other.z - body.z;
                    double distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
                    double force = G * other.mass / (distance * distance * distance);
                    ax += force * dx;
                    ay += force * dy;
                    az += force * dz;
                }
            }
            body.vx += ax * DT;
            body.vy += ay * DT;
            body.vz += az * DT;
        });

        bodies.parallelStream().forEach(body -> {
            body.x += body.vx * DT;
            body.y += body.vy * DT;
            body.z += body.vz * DT;
        });
    }

    public static void main(String[] args) {
        List<Body> bodies = new ArrayList<>(NUM_BODIES);
        initializeBodies(bodies);

        // double initialEnergy = calculateEnergy(bodies);
        // System.out.printf("Initial energy: %e%n", initialEnergy);

        for (int step = 0; step < 100; step++) {
            kickStep(bodies);
        }

        // double finalEnergy = calculateEnergy(bodies);
        // System.out.printf("Final energy: %e%n", finalEnergy);
    }
}
