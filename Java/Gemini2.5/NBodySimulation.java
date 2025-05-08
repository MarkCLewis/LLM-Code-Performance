import java.util.Random;
import java.util.Locale;
import java.util.concurrent.TimeUnit;

/**
 * Sequential 3D N-Body Simulation using a first-order Kick-Step (Euler-Cromer) method.
 */
public class NBodySimulation {

    // Constants
    static final double G = 1.0; // Gravitational constant (simulation units)
    // Softening factor squared to avoid division by zero and mitigate large forces during close encounters
    // A small value is necessary for numerical stability, especially with discrete time steps.
    static final double SOFTENING_SQUARED = 1e-6;

    // System state arrays (using separate arrays for potential cache/memory efficiency with large N)
    static double[] masses;
    static double[] posX, posY, posZ; // Position components
    static double[] velX, velY, velZ; // Velocity components
    // Temporary storage for forces calculated in each step (avoids recalculation)
    static double[] forceX, forceY, forceZ;
    static int N; // Total number of bodies

    /**
     * Initializes the simulation system with a central body and orbiting smaller bodies.
     *
     * @param centralMass Mass of the central body.
     * @param numSmallBodies Number of smaller orbiting bodies.
     * @param smallMass Mass of each small body.
     * @param minDist Minimum orbital radius for small bodies.
     * @param maxDist Maximum orbital radius for small bodies.
     * @param seed Seed for the random number generator for reproducibility.
     */
    public static void initializeSystem(double centralMass, int numSmallBodies, double smallMass,
                                        double minDist, double maxDist, long seed) {
        if (minDist <= 0 || maxDist <= minDist) {
             throw new IllegalArgumentException("Radii must be positive and maxDist > minDist.");
        }

        N = numSmallBodies + 1; // Total number of bodies
        System.out.printf("Allocating memory for %,d bodies...\n", N);
        try {
            masses = new double[N];
            posX = new double[N]; posY = new double[N]; posZ = new double[N];
            velX = new double[N]; velY = new double[N]; velZ = new double[N];
            // Allocate force arrays once
            forceX = new double[N]; forceY = new double[N]; forceZ = new double[N];
        } catch (OutOfMemoryError e) {
             System.err.println("**************************************************************");
             System.err.println("ERROR: Failed to allocate memory for " + N + " bodies.");
             System.err.println("Java Heap Space insufficient. You may need to increase the heap size");
             System.err.println("using the -Xmx flag (e.g., java -Xmx8G NBodySimulation)");
             System.err.println("**************************************************************");
             throw e; // Re-throw the error after logging
        }
        System.out.println("Memory allocation successful.");

        Random rand = new Random(seed);

        // 1. Initialize Central Body (at index 0)
        masses[0] = centralMass;
        posX[0] = 0.0; posY[0] = 0.0; posZ[0] = 0.0;
        velX[0] = 0.0; velY[0] = 0.0; velZ[0] = 0.0;

        // 2. Initialize Small Bodies (indices 1 to N-1)
        for (int i = 1; i < N; i++) {
            masses[i] = smallMass;

            // --- Position ---
            // Choose a random radius
            double r = minDist + (maxDist - minDist) * rand.nextDouble();
            // Choose a random direction in 3D space (uniform distribution on a sphere)
            double phi = 2.0 * Math.PI * rand.nextDouble();       // Azimuthal angle (0 to 2*pi)
            double costheta = 2.0 * rand.nextDouble() - 1.0;      // Cosine of polar angle (-1 to 1)
            double sintheta = Math.sqrt(1.0 - costheta * costheta); // Sine of polar angle

            // Convert spherical to Cartesian coordinates for position
            double px = r * sintheta * Math.cos(phi);
            double py = r * sintheta * Math.sin(phi);
            double pz = r * costheta;
            posX[i] = px;
            posY[i] = py;
            posZ[i] = pz;

            // --- Velocity for Circular Orbit ---
            // Speed required for a circular orbit at radius r around the central mass M
            double orbitalSpeed = Math.sqrt(G * centralMass / r);

            // Velocity vector must be perpendicular to the position vector (relative to central body at origin)
            // We need a direction vector `v_dir` such that dot(pos_vec, v_dir) = 0.
            // A robust way is to take the cross product of the position vector `p = (px, py, pz)`
            // with *any* vector not parallel to it. A random vector works well usually.
            double[] posVec = {px, py, pz};
            double[] randomVec = {rand.nextDouble() - 0.5, rand.nextDouble() - 0.5, rand.nextDouble() - 0.5};

            double crossX = posVec[1] * randomVec[2] - posVec[2] * randomVec[1];
            double crossY = posVec[2] * randomVec[0] - posVec[0] * randomVec[2];
            double crossZ = posVec[0] * randomVec[1] - posVec[1] * randomVec[0];

            double crossMag = Math.sqrt(crossX * crossX + crossY * crossY + crossZ * crossZ);

            // Handle the (rare) case where the random vector is parallel to the position vector
            if (crossMag < 1e-10) {
                // If parallel, try a standard axis like Z, unless position is along Z, then use X.
                if (Math.abs(posVec[0]) < 1e-9 && Math.abs(posVec[1]) < 1e-9) { // Along Z-axis
                    randomVec[0] = 1; randomVec[1] = 0; randomVec[2] = 0; // Use X-axis
                } else {
                    randomVec[0] = 0; randomVec[1] = 0; randomVec[2] = 1; // Use Z-axis
                }
                // Recalculate cross product
                crossX = posVec[1] * randomVec[2] - posVec[2] * randomVec[1];
                crossY = posVec[2] * randomVec[0] - posVec[0] * randomVec[2];
                crossZ = posVec[0] * randomVec[1] - posVec[1] * randomVec[0];
                crossMag = Math.sqrt(crossX * crossX + crossY * crossY + crossZ * crossZ);
            }

            // Normalize the cross product to get the velocity direction
            double vx_dir = crossX / crossMag;
            double vy_dir = crossY / crossMag;
            double vz_dir = crossZ / crossMag;

            // Final velocity vector
            velX[i] = orbitalSpeed * vx_dir;
            velY[i] = orbitalSpeed * vy_dir;
            velZ[i] = orbitalSpeed * vz_dir;
        }
         System.out.printf("Finished initializing %,d bodies.\n", N);
    }

    /**
     * Calculates the total energy (Kinetic + Potential) of the system.
     * This is an O(N^2) calculation.
     * @return The total energy of the system.
     */
    public static double calculateTotalEnergy() {
        double kineticEnergy = 0.0;
        double potentialEnergy = 0.0;

        // --- Kinetic Energy ---
        // KE = sum(0.5 * m_i * v_i^2)
        for (int i = 0; i < N; i++) {
            double vSq = velX[i] * velX[i] + velY[i] * velY[i] + velZ[i] * velZ[i];
            kineticEnergy += 0.5 * masses[i] * vSq;
        }

        // --- Potential Energy ---
        // PE = sum_{i < j} (-G * m_i * m_j / |r_i - r_j|)
        // This requires iterating through all unique pairs (O(N^2))
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) { // Start j from i+1 to avoid double counting and self-interaction
                double dx = posX[i] - posX[j];
                double dy = posY[i] - posY[j];
                double dz = posZ[i] - posZ[j];

                double distSq = dx * dx + dy * dy + dz * dz;
                // Use softening factor in distance calculation for potential energy as well for consistency
                double dist = Math.sqrt(distSq + SOFTENING_SQUARED);

                potentialEnergy -= G * masses[i] * masses[j] / dist;
            }
        }
        return kineticEnergy + potentialEnergy;
    }

    /**
     * Performs a single simulation time step using the first-order Kick-Step method.
     * 1. Calculate forces based on current positions.
     * 2. Update velocities based on forces (Kick).
     * 3. Update positions based on *new* velocities (Step).
     * This involves an O(N^2) force calculation.
     *
     * @param dt The time step duration.
     */
    public static void simulationStep(double dt) {
        // 1. Calculate Forces (O(N^2))
        // Reset forces from previous step
        for (int i = 0; i < N; i++) {
            forceX[i] = 0.0; forceY[i] = 0.0; forceZ[i] = 0.0;
        }

        // Iterate through all unique pairs of bodies (i, j) where i < j
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                // Calculate distance vector from i to j
                double dx = posX[j] - posX[i];
                double dy = posY[j] - posY[i];
                double dz = posZ[j] - posZ[i];

                // Calculate squared distance, including the softening factor
                double distSq = dx * dx + dy * dy + dz * dz + SOFTENING_SQUARED;

                // Calculate 1 / distance^3 = 1 / (distSq)^(3/2)
                // distSq = (r^2 + eps^2)
                // (distSq)^(3/2) = sqrt(distSq * distSq * distSq)
                double distSixth = distSq * distSq * distSq; // (r^2 + eps^2)^3
                double invDistCube = 1.0 / Math.sqrt(distSixth); // 1 / (r^2 + eps^2)^(3/2)

                // Calculate scalar force magnitude component F_ij = G * m_i * m_j / (r^2 + eps^2)
                // Force vector F_vec = F_ij * (r_vec / r) = G * m_i * m_j * r_vec / r^3
                // Including softening: F_vec = G * m_i * m_j * r_vec / (r^2 + eps^2)^(3/2)
                double forceScalar = G * masses[i] * masses[j] * invDistCube;

                // Calculate force components
                double fx = forceScalar * dx;
                double fy = forceScalar * dy;
                double fz = forceScalar * dz;

                // Add force to body i (force exerted by j on i)
                forceX[i] += fx; forceY[i] += fy; forceZ[i] += fz;

                // Subtract force from body j (Newton's 3rd Law: F_ji = -F_ij)
                forceX[j] -= fx; forceY[j] -= fy; forceZ[j] -= fz;
            }
        }

        // 2. Update Velocities (Kick) & 3. Update Positions (Step)
        for (int i = 0; i < N; i++) {
            // Calculate acceleration: a = F / m
            double accX = forceX[i] / masses[i];
            double accY = forceY[i] / masses[i];
            double accZ = forceZ[i] / masses[i];

            // Kick: Update velocity: v(t + dt) = v(t) + a(t) * dt
            velX[i] += accX * dt;
            velY[i] += accY * dt;
            velZ[i] += accZ * dt;

            // Step: Update position: r(t + dt) = r(t) + v(t + dt) * dt
            // Note: Using the *updated* velocity for the position step (Euler-Cromer)
            posX[i] += velX[i] * dt;
            posY[i] += velY[i] * dt;
            posZ[i] += velZ[i] * dt;
        }
    }

    /**
     * Helper method to format time duration.
     */
    private static String formatDuration(long millis) {
        long hours = TimeUnit.MILLISECONDS.toHours(millis);
        long minutes = TimeUnit.MILLISECONDS.toMinutes(millis) % 60;
        long seconds = TimeUnit.MILLISECONDS.toSeconds(millis) % 60;
        long milliseconds = millis % 1000;
        if (hours > 0) {
             return String.format("%dh %02dm %02ds", hours, minutes, seconds);
        } else if (minutes > 0) {
            return String.format("%dm %02ds %03dms", minutes, seconds, milliseconds);
        } else {
             return String.format("%ds %03dms", seconds, milliseconds);
        }
    }


    // Main simulation driver method
    public static void main(String[] args) {
        // Set Locale for consistent number formatting (e.g., using '.' for decimals)
        Locale.setDefault(Locale.US);
        System.out.println("Starting Sequential 3D N-Body Simulation in Java...");

        // --- Simulation Parameters ---
        int numSmallBodies = 10000; // ONE MILLION orbiting bodies as requested
        int numSteps = 100;            // Number of simulation steps as requested
        double centralMass = 1.0e6;      // Mass of the central body (arbitrary simulation units)
        double smallMass = 1.0;          // Mass of the small bodies
        double minDist = 5.0;           // Minimum orbital radius
        double maxDist = 50.0;          // Maximum orbital radius
        double dt = 0.005;              // Time step duration. Smaller is more accurate but slower.
                                        // Needs to be small enough to resolve closest orbits/interactions.
        long seed = 123456789L;         // Seed for random initialization

        System.out.println("-------------------- Parameters --------------------");
        System.out.printf("Number of small bodies: %,d\n", numSmallBodies);
        System.out.printf("Total bodies (N):       %,d\n", numSmallBodies + 1);
        System.out.printf("Number of steps:        %,d\n", numSteps);
        System.out.printf("Time step (dt):         %.5f\n", dt);
        System.out.printf("Central Mass:           %.2e\n", centralMass);
        System.out.printf("Small Body Mass:        %.2e\n", smallMass);
        System.out.printf("Orbital Radii Range:    [%.1f, %.1f]\n", minDist, maxDist);
        System.out.printf("Softening^2:            %.2e\n", SOFTENING_SQUARED);
        System.out.printf("Random Seed:            %d\n", seed);
        System.out.println("--------------------------------------------------");

        // --- *********************************************************** ---
        // ---                        MAJOR WARNING                        ---
        // --- *********************************************************** ---
        long estimatedOpsPerStep = (long) (numSmallBodies + 1) * (numSmallBodies) / 2; // Approx N^2/2 interactions
        long totalEstimatedOps = estimatedOpsPerStep * numSteps;
        System.out.println("\n************************** WARNING **************************");
        System.out.printf ("Simulating N = %,d bodies sequentially requires O(N^2) operations per step.\n", (numSmallBodies+1));
        System.out.printf ("Estimated pairwise interactions per step: ~%,d\n", estimatedOpsPerStep);
        System.out.printf ("Estimated total pairwise interactions:    ~%,d\n", totalEstimatedOps);
        System.out.println("This computation will be EXTREMELY SLOW.");
        System.out.println("Expect runtimes potentially spanning MANY HOURS or even DAYS on typical hardware.");
        System.out.println("Consider significantly reducing 'numSmallBodies' (e.g., to 1000) for testing.");
        System.out.println("Running with 1 Million bodies is often impractical without parallelization");
        System.out.println("and/or more advanced algorithms (like Barnes-Hut O(N log N)).");
        System.out.println("Ensure you have sufficient RAM allocated (use -Xmx JVM flag if needed).");
        System.out.println("***********************************************************\n");
        // --- End Warning ---


        // 1. Initialize the System
        System.out.println("Initializing system...");
        long initStartTime = System.currentTimeMillis();
        try {
             initializeSystem(centralMass, numSmallBodies, smallMass, minDist, maxDist, seed);
        } catch (OutOfMemoryError e) {
            System.err.println("Exiting due to OutOfMemoryError during initialization.");
            return; // Stop execution if memory allocation failed
        }
        long initEndTime = System.currentTimeMillis();
        System.out.printf("System initialization complete. Time: %s\n", formatDuration(initEndTime - initStartTime));


        // 2. Calculate Initial Energy
        System.out.println("Calculating initial energy...");
        long energyStartTime = System.currentTimeMillis();
        double initialEnergy = calculateTotalEnergy();
        long energyEndTime = System.currentTimeMillis();
        System.out.printf("Initial Total Energy: %.8e (Calculation time: %s)\n",
                          initialEnergy, formatDuration(energyEndTime - energyStartTime));


        // 3. Run Simulation Loop
        System.out.printf("Starting simulation loop for %,d steps...\n", numSteps);
        long simStartTime = System.currentTimeMillis();
        long lastReportTime = simStartTime;
        // Determine a reasonable reporting interval based on N
        int reportInterval = (N > 10000) ? 1 : (N > 1000 ? 10 : 100); // Report more often for large N
        if(numSteps < 100) reportInterval = 1;

        for (int step = 0; step < numSteps; step++) {
            simulationStep(dt);

            // --- Progress Reporting ---
            if ((step + 1) % reportInterval == 0 || step == numSteps - 1) {
                long currentTime = System.currentTimeMillis();
                double timeElapsedTotal = (currentTime - simStartTime) / 1000.0;
                // Calculate average time for the *last interval* of steps
                double intervalDurationMillis = (currentTime - lastReportTime);
                double avgStepTimeMillis = intervalDurationMillis / reportInterval;
                // Estimate remaining time
                int stepsRemaining = numSteps - (step + 1);
                double estimatedRemainingMillis = stepsRemaining * avgStepTimeMillis;

                System.out.printf("  Step %d/%d completed. Total time: %.2f s. Avg step time (last %d): %.3f ms. Est. remaining: %s\n",
                                  step + 1, numSteps, timeElapsedTotal,
                                  reportInterval, avgStepTimeMillis, formatDuration((long)estimatedRemainingMillis) );
                lastReportTime = currentTime;

                // Adjust report interval dynamically if steps are very slow
                if (avgStepTimeMillis > 30000 && reportInterval>1) { // If steps take > 30s, report every step
                    reportInterval = 1;
                    System.out.println("  (Step time > 30s, adjusting report interval to 1)");
                } else if (avgStepTimeMillis > 5000 && reportInterval > 10) { // If steps take > 5s, report every 10 steps
                     reportInterval = 10;
                     System.out.println("  (Step time > 5s, adjusting report interval to 10)");
                }
            }
        }
        long simEndTime = System.currentTimeMillis();
        System.out.printf("Simulation loop complete. Total simulation time: %s\n", formatDuration(simEndTime - simStartTime));


        // 4. Calculate Final Energy
        System.out.println("Calculating final energy...");
        energyStartTime = System.currentTimeMillis();
        double finalEnergy = calculateTotalEnergy();
        energyEndTime = System.currentTimeMillis();
        System.out.printf("Final Total Energy:   %.8e (Calculation time: %s)\n",
                          finalEnergy, formatDuration(energyEndTime - energyStartTime));


        // 5. Verify Energy Conservation
        System.out.println("-------------------- Energy Check --------------------");
        System.out.printf("Initial Energy:        %.8e\n", initialEnergy);
        System.out.printf("Final Energy:          %.8e\n", finalEnergy);
        // Avoid division by zero if initial energy is zero (unlikely but possible)
        if (Math.abs(initialEnergy) < 1e-15) {
             System.out.printf("Energy Change (Absolute): %.4e\n", finalEnergy - initialEnergy);
        } else {
            double relativeEnergyChange = Math.abs((finalEnergy - initialEnergy) / initialEnergy);
            System.out.printf("Relative Energy Change: %.4e (%.4f%%)\n", relativeEnergyChange, relativeEnergyChange * 100.0);
             // Assess conservation based on relative change
             if (relativeEnergyChange < 1e-3) { // Threshold depends on dt, steps, softening
                  System.out.println("Energy conservation appears reasonable for this method and parameters.");
             } else if (relativeEnergyChange < 1e-1){
                 System.out.println("Noticeable energy drift detected. Consider reducing dt or using a higher-order integrator.");
             } else {
                  System.out.println("WARNING: Significant energy drift detected! Accuracy is questionable.");
                  System.out.println("         This might be due to large dt, insufficient softening, or limitations of the 1st-order method.");
             }
        }
        System.out.println("----------------------------------------------------");
        System.out.println("Simulation Finished.");
    }
}