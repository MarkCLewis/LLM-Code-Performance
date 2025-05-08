import java.util.Random;
import java.util.Locale;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream; // Import for streams
import java.util.Arrays; // For parallel stream operations

/**
 * Multi-threaded 3D N-Body Simulation using Java Parallel Streams and
 * a first-order Kick-Step (Euler-Cromer) method.
 */
public class NBodySimulationParallel {

    // Constants (same as before)
    static final double G = 1.0;
    static final double SOFTENING_SQUARED = 1e-6;

    // System state arrays (same as before)
    static double[] masses;
    static double[] posX, posY, posZ;
    static double[] velX, velY, velZ;
    static double[] forceX, forceY, forceZ; // Final forces for the step
    // Temporary arrays for parallel force calculation to avoid race conditions
    static double[] tempForceX, tempForceY, tempForceZ;
    static int N; // Total number of bodies

    /**
     * Initializes the simulation system (remains sequential for determinism).
     * Allocates main arrays and temporary force arrays.
     */
    public static void initializeSystem(double centralMass, int numSmallBodies, double smallMass,
                                        double minDist, double maxDist, long seed) {
        if (minDist <= 0 || maxDist <= minDist) {
            throw new IllegalArgumentException("Radii must be positive and maxDist > minDist.");
        }

        N = numSmallBodies + 1;
        System.out.printf("Allocating memory for %,d bodies...\n", N);
        try {
            masses = new double[N];
            posX = new double[N]; posY = new double[N]; posZ = new double[N];
            velX = new double[N]; velY = new double[N]; velZ = new double[N];
            forceX = new double[N]; forceY = new double[N]; forceZ = new double[N];
            // Allocate temporary force arrays needed for parallel calculation
            tempForceX = new double[N]; tempForceY = new double[N]; tempForceZ = new double[N];
        } catch (OutOfMemoryError e) {
            System.err.println("**************************************************************");
            System.err.println("ERROR: Failed to allocate memory for " + N + " bodies.");
            System.err.println("Java Heap Space insufficient. Increase heap size (-Xmx flag).");
            System.err.println("**************************************************************");
            throw e;
        }
        System.out.println("Memory allocation successful.");

        Random rand = new Random(seed);

        // Central Body (index 0)
        masses[0] = centralMass;
        posX[0] = 0.0; posY[0] = 0.0; posZ[0] = 0.0;
        velX[0] = 0.0; velY[0] = 0.0; velZ[0] = 0.0;

        // Small Bodies (indices 1 to N-1) - Kept sequential for deterministic init
        for (int i = 1; i < N; i++) {
            masses[i] = smallMass;
            double r = minDist + (maxDist - minDist) * rand.nextDouble();
            double phi = 2.0 * Math.PI * rand.nextDouble();
            double costheta = 2.0 * rand.nextDouble() - 1.0;
            double sintheta = Math.sqrt(1.0 - costheta * costheta);
            double px = r * sintheta * Math.cos(phi);
            double py = r * sintheta * Math.sin(phi);
            double pz = r * costheta;
            posX[i] = px; posY[i] = py; posZ[i] = pz;

            double orbitalSpeed = Math.sqrt(G * centralMass / r);
            double[] posVec = {px, py, pz};
            double[] randomVec = {rand.nextDouble() - 0.5, rand.nextDouble() - 0.5, rand.nextDouble() - 0.5};
            double crossX = posVec[1] * randomVec[2] - posVec[2] * randomVec[1];
            double crossY = posVec[2] * randomVec[0] - posVec[0] * randomVec[2];
            double crossZ = posVec[0] * randomVec[1] - posVec[1] * randomVec[0];
            double crossMag = Math.sqrt(crossX * crossX + crossY * crossY + crossZ * crossZ);
            if (crossMag < 1e-10) {
                if (Math.abs(posVec[0]) < 1e-9 && Math.abs(posVec[1]) < 1e-9) {
                    randomVec[0] = 1; randomVec[1] = 0; randomVec[2] = 0;
                } else {
                    randomVec[0] = 0; randomVec[1] = 0; randomVec[2] = 1;
                }
                crossX = posVec[1] * randomVec[2] - posVec[2] * randomVec[1];
                crossY = posVec[2] * randomVec[0] - posVec[0] * randomVec[2];
                crossZ = posVec[0] * randomVec[1] - posVec[1] * randomVec[0];
                crossMag = Math.sqrt(crossX * crossX + crossY * crossY + crossZ * crossZ);
            }
            velX[i] = orbitalSpeed * crossX / crossMag;
            velY[i] = orbitalSpeed * crossY / crossMag;
            velZ[i] = orbitalSpeed * crossZ / crossMag;
        }
        System.out.printf("Finished initializing %,d bodies.\n", N);
    }

    /**
     * Calculates the total energy (Kinetic + Potential) using parallel streams.
     * @return The total energy of the system.
     */
    public static double calculateTotalEnergyParallel() {
        // --- Kinetic Energy (Parallel Summation) ---
        double kineticEnergy = IntStream.range(0, N).parallel()
            .mapToDouble(i -> {
                double vSq = velX[i] * velX[i] + velY[i] * velY[i] + velZ[i] * velZ[i];
                return 0.5 * masses[i] * vSq;
            })
            .sum();

        // --- Potential Energy (Parallel Summation over pairs i < j) ---
        // Parallelize the outer loop (i), calculate sum for j > i sequentially within the thread
        double potentialEnergy = IntStream.range(0, N).parallel()
            .mapToDouble(i -> {
                double pe_i = 0.0; // Potential energy contribution related to body i
                // Inner loop iterates j > i
                for (int j = i + 1; j < N; j++) {
                    double dx = posX[i] - posX[j];
                    double dy = posY[i] - posY[j];
                    double dz = posZ[i] - posZ[j];
                    double distSq = dx * dx + dy * dy + dz * dz;
                    double dist = Math.sqrt(distSq + SOFTENING_SQUARED);
                    pe_i -= G * masses[i] * masses[j] / dist;
                }
                return pe_i;
            })
            .sum(); // Sum up contributions from all i

        return kineticEnergy + potentialEnergy;
    }

    /**
     * Performs a single simulation time step using parallel streams.
     * 1. Calculate forces in parallel (O(N^2) work).
     * 2. Update velocities and positions in parallel (O(N) work).
     * @param dt The time step duration.
     */
    public static void simulationStepParallel(double dt) {

        // 1. Calculate Forces (Parallel O(N^2))
        // Parallelize the outer loop (over body i). Each thread calculates the total force on body i.
        // Forces are stored in temporary arrays to prevent race conditions when applying Newton's 3rd law implicitly.
        IntStream.range(0, N).parallel().forEach(i -> {
            double fx_i = 0.0; // Accumulated force x on body i for this step
            double fy_i = 0.0; // Accumulated force y on body i for this step
            double fz_i = 0.0; // Accumulated force z on body i for this step

            // Inner loop: iterate through all other bodies j to calculate force exerted *by* j *on* i
            for (int j = 0; j < N; j++) {
                if (i == j) continue; // Skip self-interaction

                // Calculate distance vector from i to j
                double dx = posX[j] - posX[i];
                double dy = posY[j] - posY[i];
                double dz = posZ[j] - posZ[i];

                double distSq = dx * dx + dy * dy + dz * dz + SOFTENING_SQUARED;
                double invDistCube = 1.0 / Math.sqrt(distSq * distSq * distSq);
                double forceScalar = G * masses[j] * masses[i] * invDistCube; // Note: masses[j]*masses[i]

                // Add force component exerted by j on i
                fx_i += forceScalar * dx;
                fy_i += forceScalar * dy;
                fz_i += forceScalar * dz;
            }
            // Store the calculated total force for body i in temporary arrays
            tempForceX[i] = fx_i;
            tempForceY[i] = fy_i;
            tempForceZ[i] = fz_i;
        });

        // Copy forces from temporary arrays to the main force arrays - acts as a synchronization point
        System.arraycopy(tempForceX, 0, forceX, 0, N);
        System.arraycopy(tempForceY, 0, forceY, 0, N);
        System.arraycopy(tempForceZ, 0, forceZ, 0, N);


        // 2. Update Velocities (Kick) & 3. Update Positions (Step) (Parallel O(N))
        // This loop is independent for each body i, so it's perfectly parallelizable.
        IntStream.range(0, N).parallel().forEach(i -> {
            // Calculate acceleration: a = F / m
            double accX = forceX[i] / masses[i];
            double accY = forceY[i] / masses[i];
            double accZ = forceZ[i] / masses[i];

            // Kick: Update velocity: v(t + dt) = v(t) + a(t) * dt
            velX[i] += accX * dt;
            velY[i] += accY * dt;
            velZ[i] += accZ * dt;

            // Step: Update position: r(t + dt) = r(t) + v(t + dt) * dt
            posX[i] += velX[i] * dt;
            posY[i] += velY[i] * dt;
            posZ[i] += velZ[i] * dt;
        });
    }

    /** Helper method to format time duration. */
    private static String formatDuration(long millis) {
        long hours = TimeUnit.MILLISECONDS.toHours(millis);
        long minutes = TimeUnit.MILLISECONDS.toMinutes(millis) % 60;
        long seconds = TimeUnit.MILLISECONDS.toSeconds(millis) % 60;
        long milliseconds = millis % 1000;
        if (hours > 0) return String.format("%dh %02dm %02ds", hours, minutes, seconds);
        if (minutes > 0) return String.format("%dm %02ds %03dms", minutes, seconds, milliseconds);
        return String.format("%ds %03dms", seconds, milliseconds);
    }


    // Main simulation driver method
    public static void main(String[] args) {
        Locale.setDefault(Locale.US);
        System.out.println("Starting Parallel (Java Streams) 3D N-Body Simulation...");
        // Display the number of processors available for parallel streams
        int processors = Runtime.getRuntime().availableProcessors();
        System.out.println("Number of available processors: " + processors);
        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", Integer.toString(processors));
        System.out.println("Using common ForkJoinPool with parallelism: " +
            java.util.concurrent.ForkJoinPool.commonPool().getParallelism());


        // --- Simulation Parameters --- (Same as before)
        int numSmallBodies = 10000; // ONE MILLION orbiting bodies as requested
        int numSteps = 100;            // Number of simulation steps as requested
        double centralMass = 1.0e6;
        double smallMass = 1.0;
        double minDist = 5.0;
        double maxDist = 50.0;
        double dt = 0.005;
        long seed = 123456789L;

        System.out.println("-------------------- Parameters --------------------");
        System.out.printf("Number of small bodies: %,d\n", numSmallBodies);
        System.out.printf("Total bodies (N):       %,d\n", numSmallBodies + 1);
        System.out.printf("Number of steps:        %,d\n", numSteps);
        System.out.printf("Time step (dt):         %.5f\n", dt);
        // ... (print other parameters)
        System.out.println("--------------------------------------------------");

        // --- WARNING --- (Still relevant, parallel just makes it faster, not instantaneous)
        System.out.println("\n************************** WARNING **************************");
        System.out.printf ("Simulating N = %,d bodies requires O(N^2) work per step.\n", (numSmallBodies+1));
        System.out.println("Parallel streams will utilize multiple CPU cores, but the computation remains intensive.");
        System.out.println("Expect significant runtime, though potentially much faster than sequential.");
        System.out.println("Ensure sufficient RAM allocated (use -Xmx JVM flag).");
        System.out.println("***********************************************************\n");


        // 1. Initialize the System (Sequential)
        System.out.println("Initializing system (sequentially)...");
        long initStartTime = System.currentTimeMillis();
        try {
            initializeSystem(centralMass, numSmallBodies, smallMass, minDist, maxDist, seed);
        } catch (OutOfMemoryError e) {
            System.err.println("Exiting due to OutOfMemoryError during initialization.");
            return;
        }
        long initEndTime = System.currentTimeMillis();
        System.out.printf("System initialization complete. Time: %s\n", formatDuration(initEndTime - initStartTime));


        // 2. Calculate Initial Energy (Parallel)
        System.out.println("Calculating initial energy (parallel)...");
        long energyStartTime = System.currentTimeMillis();
        double initialEnergy = calculateTotalEnergyParallel(); // Use parallel version
        long energyEndTime = System.currentTimeMillis();
        System.out.printf("Initial Total Energy: %.8e (Calculation time: %s)\n",
                          initialEnergy, formatDuration(energyEndTime - energyStartTime));


        // 3. Run Simulation Loop (Parallel Steps)
        System.out.printf("Starting simulation loop for %,d steps (using parallel steps)...\n", numSteps);
        long simStartTime = System.currentTimeMillis();
        long lastReportTime = simStartTime;
        int reportInterval = (N > 10000) ? 1 : (N > 1000 ? 10 : 100);
        if(numSteps < 100) reportInterval = 1;

        for (int step = 0; step < numSteps; step++) {
            simulationStepParallel(dt); // Use parallel version

            // --- Progress Reporting --- (Same logic as before)
            if ((step + 1) % reportInterval == 0 || step == numSteps - 1) {
                long currentTime = System.currentTimeMillis();
                double timeElapsedTotal = (currentTime - simStartTime) / 1000.0;
                double intervalDurationMillis = (currentTime - lastReportTime);
                double avgStepTimeMillis = intervalDurationMillis / reportInterval;
                int stepsRemaining = numSteps - (step + 1);
                double estimatedRemainingMillis = stepsRemaining * avgStepTimeMillis;

                System.out.printf("  Step %d/%d completed. Total time: %.2f s. Avg step time (last %d): %.3f ms. Est. remaining: %s\n",
                                  step + 1, numSteps, timeElapsedTotal,
                                  reportInterval, avgStepTimeMillis, formatDuration((long)estimatedRemainingMillis) );
                lastReportTime = currentTime;
                 // Adjust report interval dynamically
                 if (avgStepTimeMillis > 30000 && reportInterval > 1) { reportInterval = 1; System.out.println("  (Step time > 30s, adjusting report interval to 1)");}
                 else if (avgStepTimeMillis > 5000 && reportInterval > 10) { reportInterval = 10; System.out.println("  (Step time > 5s, adjusting report interval to 10)");}
            }
        }
        long simEndTime = System.currentTimeMillis();
        System.out.printf("Simulation loop complete. Total simulation time: %s\n", formatDuration(simEndTime - simStartTime));


        // 4. Calculate Final Energy (Parallel)
        System.out.println("Calculating final energy (parallel)...");
        energyStartTime = System.currentTimeMillis();
        double finalEnergy = calculateTotalEnergyParallel(); // Use parallel version
        energyEndTime = System.currentTimeMillis();
        System.out.printf("Final Total Energy:   %.8e (Calculation time: %s)\n",
                          finalEnergy, formatDuration(energyEndTime - energyStartTime));


        // 5. Verify Energy Conservation (Same logic as before)
        System.out.println("-------------------- Energy Check --------------------");
        System.out.printf("Initial Energy:        %.8e\n", initialEnergy);
        System.out.printf("Final Energy:          %.8e\n", finalEnergy);
        if (Math.abs(initialEnergy) < 1e-15) {
            System.out.printf("Energy Change (Absolute): %.4e\n", finalEnergy - initialEnergy);
        } else {
            double relativeEnergyChange = Math.abs((finalEnergy - initialEnergy) / initialEnergy);
            System.out.printf("Relative Energy Change: %.4e (%.4f%%)\n", relativeEnergyChange, relativeEnergyChange * 100.0);
            if (relativeEnergyChange < 1e-3) System.out.println("Energy conservation appears reasonable.");
            else if (relativeEnergyChange < 1e-1) System.out.println("Noticeable energy drift detected.");
            else System.out.println("WARNING: Significant energy drift detected!");
        }
        System.out.println("----------------------------------------------------");
        System.out.println("Parallel Simulation Finished.");
    }
}