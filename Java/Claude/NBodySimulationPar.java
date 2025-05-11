import java.util.Random;
import java.util.Arrays;
import java.util.stream.IntStream;

public class NBodySimulationPar {
    // Gravitational constant (using an arbitrary value for simulation purposes)
    private static final double G = 6.67430e-11;
    
    // Simulation bodies
    private Body[] bodies;
    
    // Simulation parameters
    private double timeStep;
    private int numSteps;
    
    // Class to represent a body in the simulation
    static class Body {
        // Position components
        double x, y, z;
        // Velocity components
        double vx, vy, vz;
        // Mass
        double mass;
        
        public Body(double x, double y, double z, double vx, double vy, double vz, double mass) {
            this.x = x;
            this.y = y;
            this.z = z;
            this.vx = vx;
            this.vy = vy;
            this.vz = vz;
            this.mass = mass;
        }
        
        // Copy constructor
        public Body(Body other) {
            this.x = other.x;
            this.y = other.y;
            this.z = other.z;
            this.vx = other.vx;
            this.vy = other.vy;
            this.vz = other.vz;
            this.mass = other.mass;
        }
    }
    
    public NBodySimulationPar(int numBodies, double timeStep, int numSteps) {
        this.bodies = new Body[numBodies];
        this.timeStep = timeStep;
        this.numSteps = numSteps;
    }
    
    /**
     * Initialize a central body with specified number of smaller bodies on circular orbits
     * @param centralMass Mass of the central body
     * @param bodyMass Mass of each orbiting body
     * @param minRadius Minimum orbit radius
     * @param maxRadius Maximum orbit radius
     */
    public void initializeSystem(double centralMass, double bodyMass, double minRadius, double maxRadius) {
        // Central body at origin with no initial velocity
        bodies[0] = new Body(0, 0, 0, 0, 0, 0, centralMass);
        
        // Create a thread-safe random number generator with fixed seed for reproducibility
        Random seedGenerator = new Random(42);
        long seed = seedGenerator.nextLong();
        
        // Use parallel streams to initialize the orbiting bodies
        IntStream.range(1, bodies.length)
            .parallel()
            .forEach(i -> {
                // Each thread needs its own random generator to avoid contention
                // Use a different seed for each thread based on the index
                Random random = new Random(seed + i);
                
                // Random radius between minRadius and maxRadius
                double radius = minRadius + random.nextDouble() * (maxRadius - minRadius);
                
                // Random position on a sphere of the given radius
                double theta = random.nextDouble() * 2 * Math.PI;
                double phi = Math.acos(2 * random.nextDouble() - 1);
                
                double x = radius * Math.sin(phi) * Math.cos(theta);
                double y = radius * Math.sin(phi) * Math.sin(theta);
                double z = radius * Math.cos(phi);
                
                // Calculate velocity for a circular orbit
                // v = sqrt(G * M / r)
                double velocity = Math.sqrt(G * centralMass / radius);
                
                // Direction perpendicular to radius vector (for circular orbit)
                // We need a vector perpendicular to (x, y, z)
                // One approach is to cross product with a reference vector
                double[] perpDirection = crossProduct(x, y, z, 0, 0, 1);
                
                // If the result is close to zero, use a different reference vector
                if (magnitude(perpDirection) < 1e-10) {
                    perpDirection = crossProduct(x, y, z, 0, 1, 0);
                }
                
                // Normalize the perpendicular direction
                double perpMagnitude = magnitude(perpDirection);
                perpDirection[0] /= perpMagnitude;
                perpDirection[1] /= perpMagnitude;
                perpDirection[2] /= perpMagnitude;
                
                // Set velocity components
                double vx = velocity * perpDirection[0];
                double vy = velocity * perpDirection[1];
                double vz = velocity * perpDirection[2];
                
                bodies[i] = new Body(x, y, z, vx, vy, vz, bodyMass);
            });
    }
    
    // Helper function to compute cross product
    private double[] crossProduct(double x1, double y1, double z1, double x2, double y2, double z2) {
        double[] result = new double[3];
        result[0] = y1 * z2 - z1 * y2;
        result[1] = z1 * x2 - x1 * z2;
        result[2] = x1 * y2 - y1 * x2;
        return result;
    }
    
    // Helper function to compute magnitude of a vector
    private double magnitude(double[] vector) {
        return Math.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
    }
    
    // Calculate total energy (kinetic + potential) of the system
    public double calculateTotalEnergy() {
        // Calculate kinetic energy for each body using parallel streams
        double kineticEnergy = Arrays.stream(bodies)
            .parallel()
            .mapToDouble(body -> {
                double speedSquared = body.vx * body.vx + body.vy * body.vy + body.vz * body.vz;
                return 0.5 * body.mass * speedSquared;
            })
            .sum();
        
        // Calculate potential energy using parallel streams
        double potentialEnergy = IntStream.range(0, bodies.length)
            .parallel()
            .mapToDouble(i -> {
                Body bi = bodies[i];
                return IntStream.range(i + 1, bodies.length)
                    .mapToDouble(j -> {
                        Body bj = bodies[j];
                        
                        double dx = bi.x - bj.x;
                        double dy = bi.y - bj.y;
                        double dz = bi.z - bj.z;
                        
                        double distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
                        
                        // Avoid division by zero
                        if (distance > 0) {
                            return -G * bi.mass * bj.mass / distance;
                        }
                        return 0.0;
                    })
                    .sum();
            })
            .sum();
            
        return kineticEnergy + potentialEnergy;
    }
    
    // Run the simulation using a first-order kick-step method (leapfrog)
    public void runSimulation() {
        // Display initial energy
        double initialEnergy = calculateTotalEnergy();
        System.out.println("Initial total energy: " + initialEnergy);
        
        // Run simulation steps
        for (int step = 0; step < numSteps; step++) {
            if (step % 100 == 0) {
                System.out.println("Step " + step + " of " + numSteps);
            }
            
            updateVelocities();
            updatePositions();
        }
        
        // Display final energy
        double finalEnergy = calculateTotalEnergy();
        System.out.println("Final total energy: " + finalEnergy);
        System.out.println("Energy change: " + (finalEnergy - initialEnergy));
        System.out.println("Relative energy change: " + ((finalEnergy - initialEnergy) / initialEnergy));
    }
    
    // Update velocities using forces (kick step)
    private void updateVelocities() {
        // Use parallel streams for multithreading
        java.util.stream.IntStream.range(0, bodies.length)
            .parallel()
            .forEach(i -> {
                Body bi = bodies[i];
                double fx = 0, fy = 0, fz = 0;
                
                // Calculate net force on body i from all other bodies
                for (int j = 0; j < bodies.length; j++) {
                    if (i == j) continue;
                    
                    Body bj = bodies[j];
                    
                    double dx = bj.x - bi.x;
                    double dy = bj.y - bi.y;
                    double dz = bj.z - bi.z;
                    
                    double distanceSquared = dx * dx + dy * dy + dz * dz;
                    double distance = Math.sqrt(distanceSquared);
                    
                    // Add small softening parameter to avoid numerical instability
                    double softening = 1e-10;
                    distanceSquared += softening;
                    
                    // Calculate gravitational force magnitude (F = G * m1 * m2 / r^2)
                    double forceMagnitude = G * bi.mass * bj.mass / distanceSquared;
                    
                    // Force components (F * direction)
                    fx += forceMagnitude * dx / distance;
                    fy += forceMagnitude * dy / distance;
                    fz += forceMagnitude * dz / distance;
                }
                
                // Update velocity using F = ma => a = F/m => Δv = a * Δt
                bi.vx += fx / bi.mass * timeStep;
                bi.vy += fy / bi.mass * timeStep;
                bi.vz += fz / bi.mass * timeStep;
            });
    }
    
    // Update positions using velocities (drift step)
    private void updatePositions() {
        // Use parallel streams for updating positions
        java.util.Arrays.stream(bodies)
            .parallel()
            .forEach(body -> {
                body.x += body.vx * timeStep;
                body.y += body.vy * timeStep;
                body.z += body.vz * timeStep;
            });
    }
    
    // Main method to run the simulation
    public static void main(String[] args) {
        // Number of bodies (1 central + orbiting bodies)
        final int NUM_BODIES = 10_001; // 1 million small bodies + 1 central body
        
        // Time step (smaller is more accurate but slower)
        final double TIME_STEP = 0.001 * 3600*24*365;
        
        // Number of simulation steps
        final int NUM_STEPS = 100;
        
        // Create simulation
        NBodySimulationPar simulation = new NBodySimulationPar(NUM_BODIES, TIME_STEP, NUM_STEPS);
        
        // Central mass (e.g., like a star)
        double centralMass = 1.0e30;
        
        // Orbiting body mass (e.g., like planets)
        double bodyMass = 1.0e20;
        
        // Orbit radius range
        double minRadius = 1.0e8;
        double maxRadius = 5.0e8;
        
        System.out.println("Initializing " + (NUM_BODIES - 1) + " bodies around central mass...");
        long startTime = System.currentTimeMillis();
        
        // Initialize the system
        simulation.initializeSystem(centralMass, bodyMass, minRadius, maxRadius);
        
        long initTime = System.currentTimeMillis() - startTime;
        System.out.println("Initialization completed in " + initTime + " ms");
        
        // Run the simulation
        System.out.println("Running simulation for " + NUM_STEPS + " steps...");
        startTime = System.currentTimeMillis();
        
        simulation.runSimulation();
        
        long simTime = System.currentTimeMillis() - startTime;
        System.out.println("Simulation completed in " + simTime + " ms");
    }
}