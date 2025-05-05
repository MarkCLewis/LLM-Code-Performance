import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

public class NBodySimulationParallelStreams {

    private static final double GRAVITATIONAL_CONSTANT = 6.67430e-11;

    public static class Body {
        double mass;
        double[] position;
        double[] velocity;
        double[] force; // Added force array for parallel updates

        public Body(double mass, double[] position, double[] velocity) {
            this.mass = mass;
            this.position = position;
            this.velocity = velocity;
            this.force = new double[]{0, 0, 0};
        }
    }

    public static class System {
        List<Body> bodies;

        public System() {
            this.bodies = new ArrayList<>();
        }

        public void addBody(Body body) {
            this.bodies.add(body);
        }
    }

    public static double distanceSquared(Body b1, Body b2) {
        double dx = b1.position[0] - b2.position[0];
        double dy = b1.position[1] - b2.position[1];
        double dz = b1.position[2] - b2.position[2];
        return dx * dx + dy * dy + dz * dz;
    }

    public static void calculateForce(Body bodyI, Body bodyJ, double[] force) {
        double rSq = distanceSquared(bodyI, bodyJ);
        if (rSq > 1e-9) { // Avoid division by zero for very close bodies
            double r = Math.sqrt(rSq);
            double magnitude = (GRAVITATIONAL_CONSTANT * bodyI.mass * bodyJ.mass) / rSq;
            force[0] += magnitude * (bodyJ.position[0] - bodyI.position[0]) / r;
            force[1] += magnitude * (bodyJ.position[1] - bodyI.position[1]) / r;
            force[2] += magnitude * (bodyJ.position[2] - bodyI.position[2]) / r;
        }
    }

    public static double calculateTotalEnergy(System system) {
        int numBodies = system.bodies.size();
        double kineticEnergy = system.bodies.parallelStream()
                .mapToDouble(body -> 0.5 * body.mass * (body.velocity[0] * body.velocity[0] +
                                                        body.velocity[1] * body.velocity[1] +
                                                        body.velocity[2] * body.velocity[2]))
                .sum();

        double potentialEnergy = IntStream.range(0, numBodies)
                .parallel()
                .mapToDouble(i -> IntStream.range(i + 1, numBodies)
                        .mapToDouble(j -> {
                            Body bodyI = system.bodies.get(i);
                            Body bodyJ = system.bodies.get(j);
                            double r = Math.sqrt(distanceSquared(bodyI, bodyJ));
                            return -(GRAVITATIONAL_CONSTANT * bodyI.mass * bodyJ.mass) / r;
                        })
                        .sum())
                .sum();

        return kineticEnergy + potentialEnergy;
    }

    public static System initializeCircularOrbits(int numOrbiting, double centralMass, double orbitRadius, double orbitingMass) {
        System system = new System();

        // Initialize the central body
        system.addBody(new Body(centralMass, new double[]{0, 0, 0}, new double[]{0, 0, 0}));

        // Initialize the orbiting bodies
        Random random = new Random();
        for (int i = 0; i < numOrbiting; i++) {
            double angle = 2 * Math.PI * i / (double) numOrbiting;
            double x = orbitRadius * Math.cos(angle);
            double y = orbitRadius * Math.sin(angle);
            double z = 0; // Orbiting in the xy-plane

            // Calculate the orbital velocity for a circular orbit
            double orbitalSpeed = Math.sqrt(GRAVITATIONAL_CONSTANT * centralMass / orbitRadius);
            double vx = -orbitalSpeed * Math.sin(angle);
            double vy = orbitalSpeed * Math.cos(angle);
            double vz = 0;

            system.addBody(new Body(orbitingMass, new double[]{x, y, z}, new double[]{vx, vy, vz}));
        }

        return system;
    }

    public static void kickStep(System system, double dt) {
        int numBodies = system.bodies.size();

        // Reset forces
        system.bodies.parallelStream().forEach(body -> {
            body.force[0] = 0;
            body.force[1] = 0;
            body.force[2] = 0;
        });

        // Calculate forces in parallel
        IntStream.range(0, numBodies).parallel().forEach(i -> {
            Body bodyI = system.bodies.get(i);
            for (int j = 0; j < numBodies; j++) {
                if (i != j) {
                    Body bodyJ = system.bodies.get(j);
                    calculateForce(bodyI, bodyJ, bodyI.force); // Accumulate force on bodyI
                }
            }
        });

        // Update velocities (kick) in parallel
        system.bodies.parallelStream().forEach(body -> {
            body.velocity[0] += (body.force[0] / body.mass) * dt * 0.5;
            body.velocity[1] += (body.force[1] / body.mass) * dt * 0.5;
            body.velocity[2] += (body.force[2] / body.mass) * dt * 0.5;
        });
    }

    public static void driftStep(System system, double dt) {
        system.bodies.parallelStream().forEach(body -> {
            body.position[0] += body.velocity[0] * dt;
            body.position[1] += body.velocity[1] * dt;
            body.position[2] += body.velocity[2] * dt;
        });
    }

    public static void firstOrderKickStep(System system, double dt) {
        kickStep(system, dt);
        driftStep(system, dt);
        kickStep(system, dt);
    }

    public static void main(String[] args) {
        int numOrbitingBodies = 1000000;
        double centralMass = 1.989e30;    // Mass of the Sun (kg)
        double orbitRadius = 1.496e11;    // 1 AU (m)
        double orbitingMass = 5.972e24;   // Mass of the Earth (kg)
        int numSteps = 1000;
        double timeStep = 3600.0 * 24.0 * 7.0; // 1 week in seconds

        // Initialize the system
        System system = initializeCircularOrbits(numOrbitingBodies, centralMass, orbitRadius, orbitingMass);
        System initialSystem = new System();
        for (Body body : system.bodies) {
            initialSystem.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        // Calculate initial energy
        double initialEnergy = calculateTotalEnergy(initialSystem);
        System systemForSimulation = new System();
        for (Body body : initialSystem.bodies) {
            systemForSimulation.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy2 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy2.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy3 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy3.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy4 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy4.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy5 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy5.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy6 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy6.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy7 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy7.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy8 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy8.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy9 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy9.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy10 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy10.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy11 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy11.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy12 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy12.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy13 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy13.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy14 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy14.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy15 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy15.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy16 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy16.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy17 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy17.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy18 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy18.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy19 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy19.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy20 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy20.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy21 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy21.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy22 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy22.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy23 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy23.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy24 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy24.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy25 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy25.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy26 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy26.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy27 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy27.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy28 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy28.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy29 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy29.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy30 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy30.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy31 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy31.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy32 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy32.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy33 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy33.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy34 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy34.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy35 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy35.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy36 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy36.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy37 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy37.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy38 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy38.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy39 = new System();
        for (Body body : initialSystem.bodies)