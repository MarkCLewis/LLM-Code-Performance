import java.util.Random;

public class NBodySimulation {

    private static final double GRAVITATIONAL_CONSTANT = 6.674e-11;
    private static final double TIME_STEP = 0.001; // Adjust as needed
    private static final int NUM_STEPS = 1000;

    private Body[] bodies;

    public NBodySimulation(int numBodies) {
        bodies = new Body[numBodies + 1]; // +1 for the central body
        initializeSystem();
    }

    private void initializeSystem() {
        Random random = new Random(0); // Seed for reproducibility

        // Initialize central body (Sun-like)
        bodies[0] = new Body(1.989e30, 0, 0, 0, 0, 0, 0); // Mass, x, y, z, vx, vy, vz

        // Initialize smaller bodies (planet-like) in circular orbits
        for (int i = 1; i < bodies.length; i++) {
            double mass = 1e25 + random.nextDouble() * 1e26; // Varying mass
            double radius = 1e11 + random.nextDouble() * 1e12; // Varying orbital radius
            double speed = Math.sqrt(GRAVITATIONAL_CONSTANT * bodies[0].mass / radius);

            double angle = random.nextDouble() * 2 * Math.PI;
            double x = radius * Math.cos(angle);
            double y = radius * Math.sin(angle);
            double z = 0; // Orbit in x-y plane for simplicity

            double vx = -speed * Math.sin(angle);
            double vy = speed * Math.cos(angle);
            double vz = 0;

            bodies[i] = new Body(mass, x, y, z, vx, vy, vz);
        }
    }


    public void simulate() {
        double initialEnergy = calculateTotalEnergy();
        System.out.println("Initial Energy: " + initialEnergy);

        for (int step = 0; step < NUM_STEPS; step++) {
            // Kick step (update velocities)
            for (int i = 0; i < bodies.length; i++) {
                double fx = 0, fy = 0, fz = 0;
                for (int j = 0; j < bodies.length; j++) {
                    if (i == j) continue; // Don't interact with self

                    double dx = bodies[j].x - bodies[i].x;
                    double dy = bodies[j].y - bodies[i].y;
                    double dz = bodies[j].z - bodies[i].z;
                    double distSq = dx * dx + dy * dy + dz * dz;
                    double forceMagnitude = GRAVITATIONAL_CONSTANT * bodies[i].mass * bodies[j].mass / distSq;
                    
                    fx += forceMagnitude * dx / Math.sqrt(distSq);
                    fy += forceMagnitude * dy / Math.sqrt(distSq);
                    fz += forceMagnitude * dz / Math.sqrt(distSq);
                }

                bodies[i].vx += fx / bodies[i].mass * TIME_STEP;
                bodies[i].vy += fy / bodies[i].mass * TIME_STEP;
                bodies[i].vz += fz / bodies[i].mass * TIME_STEP;
            }

            // Drift step (update positions)
            for (int i = 0; i < bodies.length; i++) {
                bodies[i].x += bodies[i].vx * TIME_STEP;
                bodies[i].y += bodies[i].vy * TIME_STEP;
                bodies[i].z += bodies[i].vz * TIME_STEP;
            }
        }
        double finalEnergy = calculateTotalEnergy();
        System.out.println("Final Energy: " + finalEnergy);
        System.out.println("Energy Change: " + (finalEnergy - initialEnergy));

    }

    private double calculateTotalEnergy() {
        double totalEnergy = 0;

        // Kinetic Energy
        for (Body body : bodies) {
            totalEnergy += 0.5 * body.mass * (body.vx * body.vx + body.vy * body.vy + body.vz * body.vz);
        }

        // Potential Energy
        for (int i = 0; i < bodies.length; i++) {
            for (int j = i + 1; j < bodies.length; j++) { // Avoid double-counting
                double dx = bodies[j].x - bodies[i].x;
                double dy = bodies[j].y - bodies[i].y;
                double dz = bodies[j].z - bodies[i].z;
                double dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
                totalEnergy -= GRAVITATIONAL_CONSTANT * bodies[i].mass * bodies[j].mass / dist;
            }
        }

        return totalEnergy;
    }


    private static class Body {
        double mass, x, y, z, vx, vy, vz;

        public Body(double mass, double x, double y, double z, double vx, double vy, double vz) {
            this.mass = mass;
            this.x = x;
            this.y = y;
            this.z = z;
            this.vx = vx;
            this.vy = vy;
            this.vz = vz;
        }
    }

    public static void main(String[] args) {
        int numSmallBodies = 100000;  // One million small bodies
        NBodySimulation simulation = new NBodySimulation(numSmallBodies);
        simulation.simulate();

    }
}
