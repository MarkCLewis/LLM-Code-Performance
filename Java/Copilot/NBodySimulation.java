import java.util.ArrayList;
import java.util.List;

public class NBodySimulation {
    private static final double G = 6.67430e-11; // Gravitational constant
    private static final double DT = 1e-3; // Time step
    private static final int NUM_BODIES = 1000000; // Number of bodies

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

        for (int i = 1; i < NUM_BODIES; i++) {
            double angle = 2 * Math.PI * i / (NUM_BODIES - 1);
            bodies.add(new Body(
                Math.cos(angle) * 1e11,
                Math.sin(angle) * 1e11,
                0,
                -Math.sin(angle) * Math.sqrt(G * bodies.get(0).mass / 1e11),
                Math.cos(angle) * Math.sqrt(G * bodies.get(0).mass / 1e11),
                0,
                1e24 // Small body mass
            ));
        }
    }

    private static double calculateEnergy(List<Body> bodies) {
        double energy = 0.0;
        for (Body body : bodies) {
            double kinetic = 0.5 * body.mass * (body.vx * body.vx + body.vy * body.vy + body.vz * body.vz);
            double potential = 0.0;
            for (Body other : bodies) {
                if (body != other) {
                    double dx = body.x - other.x;
                    double dy = body.y - other.y;
                    double dz = body.z - other.z;
                    double distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
                    potential -= G * body.mass * other.mass / distance;
                }
            }
            energy += kinetic + 0.5 * potential;
        }
        return energy;
    }

    private static void kickStep(List<Body> bodies) {
        for (Body body : bodies) {
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
        }

        for (Body body : bodies) {
            body.x += body.vx * DT;
            body.y += body.vy * DT;
            body.z += body.vz * DT;
        }
    }

    public static void main(String[] args) {
        List<Body> bodies = new ArrayList<>(NUM_BODIES);
        initializeBodies(bodies);

        double initialEnergy = calculateEnergy(bodies);
        System.out.printf("Initial energy: %e%n", initialEnergy);

        for (int step = 0; step < 1000; step++) {
            kickStep(bodies);
        }

        double finalEnergy = calculateEnergy(bodies);
        System.out.printf("Final energy: %e%n", finalEnergy);
    }
}
