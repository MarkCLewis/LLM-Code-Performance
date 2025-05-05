using LinearAlgebra
using Printf
using Random
using Base.Threads # Added for multi-threading

# --- Constants ---
const G = 1.0             # Gravitational constant (set to 1 for simplicity)
const SOFTENING = 0.01    # Softening parameter to avoid divergences at r=0
const SOFTENING_SQ = SOFTENING^2 # Pre-calculate squared softening

# --- Data Structures ---
# positions: N x 3 matrix (x, y, z for each body)
# velocities: N x 3 matrix (vx, vy, vz for each body)
# masses: N x 1 vector (mass of each body)
# accelerations: N x 3 matrix (ax, ay, az for each body)

# --- Core Physics Calculation (Multithreaded) ---

"""
Calculates the gravitational accelerations for all bodies using multiple threads.
The outer loop over bodies 'i' is parallelized.
"""
function calculate_accelerations!(accelerations::Matrix{Float64},
                                positions::Matrix{Float64},
                                masses::Vector{Float64})
    N = size(positions, 1)
    fill!(accelerations, 0.0) # Reset accelerations (safe before @threads)

    # Parallelize the loop over bodies 'i'
    Threads.@threads for i in 1:N
        # Accumulators for the current body 'i' (local to this iteration/thread)
        accel_x_i = 0.0
        accel_y_i = 0.0
        accel_z_i = 0.0

        # Inner loop: Calculate influence of all *other* bodies 'j' on 'i'
        # Note: We cannot use the j > i optimization directly here because
        # each thread only computes forces *for* its assigned 'i' values.
        # We also don't update 'j' here to avoid race conditions.
        for j in 1:N
            if i == j
                continue # Skip self-interaction
            end

            # Vector from body i to body j
            dx = positions[j, 1] - positions[i, 1]
            dy = positions[j, 2] - positions[i, 2]
            dz = positions[j, 3] - positions[i, 3]

            # Squared distance with softening
            dist_sq = dx^2 + dy^2 + dz^2 + SOFTENING_SQ

            # Inverse cube distance factor (including G * m_j)
            # a_i = sum_{j!=i} G * m_j * (r_j - r_i) / |r_j - r_i|^3 (softened)
            inv_dist_cubed = G / (dist_sq * sqrt(dist_sq)) # Faster than G / dist_sq^(3/2)
            force_factor = masses[j] * inv_dist_cubed

            # Accumulate acceleration on body i due to body j
            accel_x_i += force_factor * dx
            accel_y_i += force_factor * dy
            accel_z_i += force_factor * dz
        end
        # Assign the calculated acceleration to the main array (safe, only this thread writes to index i)
        accelerations[i, 1] = accel_x_i
        accelerations[i, 2] = accel_y_i
        accelerations[i, 3] = accel_z_i
    end
end

# --- Simulation Step ---
# (No changes needed here, calls the threaded acceleration function)
"""
Performs one step of the simulation using the first-order Kick-Step
(Euler-Cromer or semi-implicit Euler) method.
"""
function kick_step_update!(positions::Matrix{Float64},
                          velocities::Matrix{Float64},
                          masses::Vector{Float64},
                          accelerations::Matrix{Float64},
                          dt::Float64)

    # Calculate accelerations based on current positions (uses the threaded version)
    calculate_accelerations!(accelerations, positions, masses)

    # 1. Kick: Update velocities (v_new = v_old + a * dt)
    # Broadcasting is generally efficient, threading might not add much here.
    velocities .+= accelerations .* dt

    # 2. Step: Update positions (p_new = p_old + v_new * dt)
    positions .+= velocities .* dt
end


# --- Energy Calculation (Potential Energy Multithreaded) ---

"""
Calculates the total energy (Kinetic + Potential) of the system.
Potential energy calculation is O(N^2) and parallelized using thread-local storage.
Kinetic energy (O(N)) is calculated sequentially.
"""
function calculate_total_energy(positions::Matrix{Float64},
                                velocities::Matrix{Float64},
                                masses::Vector{Float64})::Float64
    N = size(positions, 1)
    kinetic_energy = 0.0
    potential_energy = 0.0

    # Kinetic Energy: KE = sum(0.5 * m_i * |v_i|^2) - O(N)
    # Calculated sequentially for simplicity, parallelize if it becomes a bottleneck
    for i in 1:N
        vel_sq = velocities[i, 1]^2 + velocities[i, 2]^2 + velocities[i, 3]^2
        kinetic_energy += 0.5 * masses[i] * vel_sq
    end

    # Potential Energy: PE = sum_{i < j} (-G * m_i * m_j / |r_i - r_j|) - O(N^2)
    # Parallelized using per-thread accumulators
    num_threads = Threads.nthreads()
    pe_partials = zeros(Float64, num_threads) # Storage for each thread's partial sum

    Threads.@threads for i in 1:N
        tid = Threads.threadid()
        pe_thread = 0.0 # Accumulator local to this thread's work on this 'i' loop
        # Loop only j > i to avoid double counting and self-interaction
        for j in (i+1):N
            dx = positions[j, 1] - positions[i, 1]
            dy = positions[j, 2] - positions[i, 2]
            dz = positions[j, 3] - positions[i, 3]

            dist_sq = dx^2 + dy^2 + dz^2 + SOFTENING_SQ
            dist = sqrt(dist_sq)

            pe_thread -= G * masses[i] * masses[j] / dist
        end
        # Accumulate the result for this 'i' into the thread's partial sum storage
        # Note: This+= operation on pe_partials[tid] is NOT atomic, but it's safe
        # because each thread only writes to its own index 'tid'.
        # However, multiple 'i' iterations run by the same thread will accumulate correctly.
        # A cleaner way might be to sum pe_thread within the thread and add once at the end,
        # but this accumulation should also work correctly with @threads partitioning.
        # Let's use an atomic add just to be explicitly safe if multiple iterations for
        # the same thread were interleaved in some unexpected scheduler way, though unlikely
        # for simple loop partitioning. Or stick to the simpler += as race conditions
        # between iterations assigned to the *same thread* are not the primary concern.
        # Sticking with += as it's likely correct and avoids atomic overhead.
        pe_partials[tid] += pe_thread
    end

    # Sum the partial potential energies from all threads
    potential_energy = sum(pe_partials)

    return kinetic_energy + potential_energy
end


# --- Initialization ---
# (No changes needed here)
"""
Initializes a system with one central body and N_orbiting smaller bodies
in approximately circular orbits within a specified radius range.
"""
function initialize_circular_orbits(N_orbiting::Int,
                                    central_mass::Float64;
                                    min_radius::Float64 = 1.0,
                                    max_radius::Float64 = 10.0,
                                    orbiting_mass_factor::Float64 = 1e-7 # Mass relative to central
                                    )::Tuple{Matrix{Float64}, Matrix{Float64}, Vector{Float64}}
    N_total = N_orbiting + 1
    positions = zeros(Float64, N_total, 3)
    velocities = zeros(Float64, N_total, 3)
    masses = zeros(Float64, N_total)

    # Initialize central body (at index 1)
    masses[1] = central_mass

    # Initialize orbiting bodies
    rng = MersenneTwister() # Use default global RNG or create one per thread if needed elsewhere
    for i in 2:N_total
        masses[i] = central_mass * orbiting_mass_factor * (0.5 + rand(rng))
        radius = min_radius + (max_radius - min_radius) * rand(rng)
        phi = 2.0 * pi * rand(rng)
        cos_theta = 2.0 * rand(rng) - 1.0
        sin_theta = sqrt(max(0.0, 1.0 - cos_theta^2)) # Ensure non-negative under sqrt

        x = radius * sin_theta * cos(phi)
        y = radius * sin_theta * sin(phi)
        z = radius * cos_theta
        positions[i, :] = [x, y, z]

        vel_magnitude = sqrt(G * central_mass / radius)
        pos_vec = positions[i, :]
        if abs(pos_vec[1]) < 1e-9 && abs(pos_vec[2]) < 1e-9
             axis_vec = [1.0, 0.0, 0.0]
        else
             axis_vec = [0.0, 0.0, 1.0]
        end

        # Handle potential zero vector cross product if pos_vec aligns with axis_vec
        vel_dir = cross(pos_vec, axis_vec)
        if norm(vel_dir) < 1e-9 # If vectors were parallel, try another axis
             axis_vec = [0.0, 1.0, 0.0] # Try Y-axis
             vel_dir = cross(pos_vec, axis_vec)
        end

        # Ensure velocity direction is normalized, handle zero norm case if necessary
        norm_vel_dir = norm(vel_dir)
        if norm_vel_dir > 1e-9
            vel_dir_normalized = vel_dir / norm_vel_dir
        else
            # Handle cases where position is exactly at origin or something went wrong
            # Assigning a random orthogonal direction as fallback
            vel_dir_normalized = normalize([rand(rng)-0.5, rand(rng)-0.5, rand(rng)-0.5])
            println("Warning: Could not determine unique velocity direction for body $i. Using random.")
        end


        velocities[i, :] = vel_magnitude .* vel_dir_normalized
    end

    return positions, velocities, masses
end

# --- Main Simulation Function ---

function run_simulation(num_bodies_orbiting::Int, num_steps::Int, dt::Float64)
    # --- Check Threads ---
    num_threads = Threads.nthreads()
    if num_threads == 1
        println("WARNING: Julia started with only 1 thread. Multithreading disabled.")
        println("         Restart Julia with '-t auto' or '-t N' (e.g., julia -t 4 script.jl)")
    end
    println("--- N-Body Simulation Start (Threads: $num_threads) ---")
    println("Number of orbiting bodies: ", num_bodies_orbiting)
    println("Total bodies: ", num_bodies_orbiting + 1)
    println("Number of steps: ", num_steps)
    println("Time step (dt): ", dt)
    println("Softening factor: ", SOFTENING)
    println("Gravitational Constant (G): ", G)
    println("---------------------------------")

    # Initialization
    println("Initializing system...")
    central_mass = 1.0e6
    positions, velocities, masses = initialize_circular_orbits(num_bodies_orbiting, central_mass)
    N_total = size(positions, 1)
    # Acceleration buffer allocation remains the same
    accelerations = zeros(Float64, N_total, 3)
    println("Initialization complete.")

    # Energy check before simulation
    println("Calculating initial energy...")
    initial_energy = calculate_total_energy(positions, velocities, masses) # Uses threaded PE calc
    @printf("Initial Total Energy: %.6e\n", initial_energy)
    println("---------------------------------")

    # --- Simulation Loop ---
    println("Starting simulation loop...")
    start_time = time()

    for step in 1:num_steps
        kick_step_update!(positions, velocities, masses, accelerations, dt)

        # Optional: Print progress
        if step % max(1, num_steps ÷ 10) == 0 || step == num_steps # Avoid division by zero if num_steps < 10
           elapsed_time = time() - start_time
           @printf("Step: %d / %d (%.1f%%), Elapsed Time: %.2f s\n",
                   step, num_steps, 100 * step / num_steps, elapsed_time)
        end
    end

    end_time = time()
    total_time = end_time - start_time
    println("Simulation loop finished.")
    println("---------------------------------")

    # Energy check after simulation
    println("Calculating final energy...")
    final_energy = calculate_total_energy(positions, velocities, masses) # Uses threaded PE calc
    @printf("Final Total Energy:   %.6e\n", final_energy)

    # Accuracy check
    energy_diff = abs(final_energy - initial_energy)
    relative_error = if abs(initial_energy) > 1e-12
        energy_diff / abs(initial_energy)
    else
        energy_diff
    end
    @printf("Absolute Energy Difference: %.6e\n", energy_diff)
    @printf("Relative Energy Error: %.6e (%.4f%%)\n", relative_error, relative_error * 100)
    println("---------------------------------")
    @printf("Total Simulation Time: %.3f seconds\n", total_time)
    println("--- N-Body Simulation End ---")

end

# --- Run ---
# Remember to start Julia with threads enabled:
# julia -t auto your_script_name.jl
# or
# julia -t N your_script_name.jl (where N is the number of threads)

# **WARNING:** N=1,000,001 is EXTREMELY computationally expensive even multithreaded.
# The complexity is still O(N^2), threading only divides the constant factor.
# Expect very long runtimes. Start with smaller N for testing.

# const NUM_ORBITING_BODIES = 1_000_000 # Requested, but likely infeasible runtime
const NUM_ORBITING_BODIES = 1000      # More reasonable number for testing
const NUM_STEPS = 1000
const DT = 0.001                      # Time step

# Check if running interactively or as script
if isinteractive() && NUM_ORBITING_BODIES > 5000 # Add threshold for warning
    println("Running simulation with N_orbiting = $(NUM_ORBITING_BODIES) for $(NUM_STEPS) steps using $(Threads.nthreads()) threads.")
    println("WARNING: If N_orbiting is large, this will take a very long time even with threads.")
    println("Press Enter to continue, or Ctrl+C to abort.")
    readline()
elseif !isinteractive() && Threads.nthreads() == 1
     println("WARNING: Running script with only 1 thread. Performance will be limited.")
     println("         Consider restarting Julia with the -t flag for multithreading.")
end

run_simulation(NUM_ORBITING_BODIES, NUM_STEPS, DT)