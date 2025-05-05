import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class NBodySimulationSequential {

    private static final double GRAVITATIONAL_CONSTANT = 6.67430e-11;

    public static class Body {
        double mass;
        double[] position;
        double[] velocity;

        public Body(double mass, double[] position, double[] velocity) {
            this.mass = mass;
            this.position = position;
            this.velocity = velocity;
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
        double kineticEnergy = 0.0;
        double potentialEnergy = 0.0;
        int numBodies = system.bodies.size();

        for (int i = 0; i < numBodies; i++) {
            Body bodyI = system.bodies.get(i);
            double vSq = bodyI.velocity[0] * bodyI.velocity[0] +
                         bodyI.velocity[1] * bodyI.velocity[1] +
                         bodyI.velocity[2] * bodyI.velocity[2];
            kineticEnergy += 0.5 * bodyI.mass * vSq;

            for (int j = i + 1; j < numBodies; j++) {
                Body bodyJ = system.bodies.get(j);
                double r = Math.sqrt(distanceSquared(bodyI, bodyJ));
                potentialEnergy -= (GRAVITATIONAL_CONSTANT * bodyI.mass * bodyJ.mass) / r;
            }
        }

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
        List<double[]> forces = new ArrayList<>();
        for (int i = 0; i < numBodies; i++) {
            forces.add(new double[]{0, 0, 0});
        }

        // Calculate forces
        for (int i = 0; i < numBodies; i++) {
            Body bodyI = system.bodies.get(i);
            double[] forceOnI = forces.get(i);
            for (int j = 0; j < numBodies; j++) {
                if (i != j) {
                    Body bodyJ = system.bodies.get(j);
                    calculateForce(bodyI, bodyJ, forceOnI);
                }
            }
        }

        // Update velocities (kick)
        for (int i = 0; i < numBodies; i++) {
            Body body = system.bodies.get(i);
            double[] force = forces.get(i);
            body.velocity[0] += (force[0] / body.mass) * dt * 0.5;
            body.velocity[1] += (force[1] / body.mass) * dt * 0.5;
            body.velocity[2] += (force[2] / body.mass) * dt * 0.5;
        }
    }

    public static void driftStep(System system, double dt) {
        for (Body body : system.bodies) {
            body.position[0] += body.velocity[0] * dt;
            body.position[1] += body.velocity[1] * dt;
            body.position[2] += body.velocity[2] * dt;
        }
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
        System finalSystem = new System();

        System energyCheckSystem = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        // Calculate initial energy
        double initialEnergy = calculateTotalEnergy(energyCheckSystem);
        System energyCheckSystem2 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem2.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }
        System energyCheckSystem3 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem3.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }
        System energyCheckSystem4 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem4.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem5 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem5.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem6 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem6.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem7 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem7.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem8 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem8.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem9 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem9.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem10 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem10.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem11 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem11.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem12 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem12.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem13 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem13.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem14 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem14.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem15 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem15.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem16 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem16.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem17 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem17.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem18 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem18.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem19 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem19.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem20 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem20.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem21 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem21.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem22 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem22.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem23 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem23.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem24 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem24.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem25 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem25.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem26 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem26.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem27 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem27.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem28 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem28.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem29 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem29.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem30 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem30.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem31 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem31.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem32 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem32.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem33 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem33.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem34 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem34.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem35 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem35.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem36 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem36.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem37 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem37.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem38 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem38.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem39 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem39.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem40 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem40.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem41 = new System();
        for (Body body : initialSystem.bodies) {
            energyCheckSystem41.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System energyCheckSystem4