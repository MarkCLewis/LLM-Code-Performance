using LinearAlgebra
using Printf
using Random

# --- Constants ---
const G = 1.0             # Gravitational constant (set to 1 for simplicity)
const SOFTENING = 0.01    # Softening parameter to avoid divergences at r=0
const SOFTENING_SQ = SOFTENING^2 # Pre-calculate squared softening

# --- Data Structures ---
# We use a Structure of Arrays (SoA) approach for potentially better cache performance
# positions: N x 3 matrix (x, y, z for each body)
# velocities: N x 3 matrix (vx, vy, vz for each body)
# masses: N x 1 vector (mass of each body)
# accelerations: N x 3 matrix (ax, ay, az for each body)

# --- Core Physics Calculation ---

"""
Calculates the gravitational accelerations for all bodies.
This is the most computationally expensive part (O(N^2)).
"""
function calculate_accelerations!(accelerations::Matrix{Float64},
                                positions::Matrix{Float64},
                                masses::Vector{Float64})
    N = size(positions, 1)
    fill!(accelerations, 0.0) # Reset accelerations

    for i in 1:N
        for j in (i+1):N # Loop j > i to calculate interactions only once
            # Vector from body i to body j
            dx = positions[j, 1] - positions[i, 1]
            dy = positions[j, 2] - positions[i, 2]
            dz = positions[j, 3] - positions[i, 3]

            # Squared distance with softening
            dist_sq = dx^2 + dy^2 + dz^2 + SOFTENING_SQ

            # Inverse cube distance factor (including G)
            # F_ij = G * m_i * m_j * (r_j - r_i) / |r_j - r_i|^3
            # a_i = F_ij / m_i = G * m_j * (r_j - r_i) / |r_j - r_i|^3
            # a_j = F_ji / m_j = -G * m_i * (r_j - r_i) / |r_j - r_i|^3
            inv_dist_cubed = G / (dist_sq * sqrt(dist_sq)) # Faster than G / dist_sq^(3/2)

            accel_factor_i = masses[j] * inv_dist_cubed
            accel_factor_j = masses[i] * inv_dist_cubed

            # Accumulate acceleration on body i due to body j
            accelerations[i, 1] += accel_factor_i * dx
            accelerations[i, 2] += accel_factor_i * dy
            accelerations[i, 3] += accel_factor_i * dz

            # Accumulate acceleration on body j due to body i (Newton's 3rd Law)
            accelerations[j, 1] -= accel_factor_j * dx
            accelerations[j, 2] -= accel_factor_j * dy
            accelerations[j, 3] -= accel_factor_j * dz
        end
    end
end

# --- Simulation Step ---

"""
Performs one step of the simulation using the first-order Kick-Step
(also known as Euler-Cromer or semi-implicit Euler) method.
1. Kick: Update velocities based on accelerations at time t.
2. Step: Update positions based on *new* velocities at time t+dt.
"""
function kick_step_update!(positions::Matrix{Float64},
                          velocities::Matrix{Float64},
                          masses::Vector{Float64},
                          accelerations::Matrix{Float64},
                          dt::Float64)

    # Calculate accelerations based on current positions
    calculate_accelerations!(accelerations, positions, masses)

    # 1. Kick: Update velocities (v_new = v_old + a * dt)
    # Use broadcasting for efficiency
    velocities .+= accelerations .* dt

    # 2. Step: Update positions (p_new = p_old + v_new * dt)
    # Use broadcasting for efficiency
    positions .+= velocities .* dt
end


# --- Energy Calculation ---

"""
Calculates the total energy (Kinetic + Potential) of the system.
Potential energy calculation is O(N^2).
"""
function calculate_total_energy(positions::Matrix{Float64},
                                velocities::Matrix{Float64},
                                masses::Vector{Float64})::Float64
    N = size(positions, 1)
    kinetic_energy = 0.0
    potential_energy = 0.0

    # Kinetic Energy: KE = sum(0.5 * m_i * |v_i|^2)
    for i in 1:N
        vel_sq = velocities[i, 1]^2 + velocities[i, 2]^2 + velocities[i, 3]^2
        kinetic_energy += 0.5 * masses[i] * vel_sq
    end

    # Potential Energy: PE = sum_{i < j} (-G * m_i * m_j / |r_i - r_j|)
    # Includes softening in distance calculation for consistency
    for i in 1:N
        for j in (i+1):N
            dx = positions[j, 1] - positions[i, 1]
            dy = positions[j, 2] - positions[i, 2]
            dz = positions[j, 3] - positions[i, 3]

            dist_sq = dx^2 + dy^2 + dz^2 + SOFTENING_SQ
            dist = sqrt(dist_sq)

            potential_energy -= G * masses[i] * masses[j] / dist
        end
    end

    return kinetic_energy + potential_energy
end


# --- Initialization ---

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
    # positions[1, :] is already [0, 0, 0]
    # velocities[1, :] is already [0, 0, 0]

    # Initialize orbiting bodies
    rng = MersenneTwister() # Random number generator
    for i in 2:N_total
        # Assign mass
        masses[i] = central_mass * orbiting_mass_factor * (0.5 + rand(rng)) # Add some variation

        # Generate random position in 3D space within the radius range
        radius = min_radius + (max_radius - min_radius) * rand(rng) # Uniform radius distribution
        
        # Use spherical coordinates for uniform distribution on sphere surface, then scale by radius
        phi = 2.0 * pi * rand(rng)       # Azimuthal angle (0 to 2pi)
        cos_theta = 2.0 * rand(rng) - 1.0 # Cosine of polar angle (-1 to 1) -> uniform on sphere
        sin_theta = sqrt(1.0 - cos_theta^2)

        x = radius * sin_theta * cos(phi)
        y = radius * sin_theta * sin(phi)
        z = radius * cos_theta
        positions[i, :] = [x, y, z]

        # Calculate velocity for a circular orbit (assuming central mass dominates)
        # v = sqrt(G * M_central / r)
        vel_magnitude = sqrt(G * central_mass / radius)

        # Velocity vector must be perpendicular to the position vector.
        # A simple way to get *a* perpendicular vector is to cross the position
        # vector with a non-parallel vector (e.g., Z-axis [0,0,1]), unless the position
        # itself is along the Z-axis.
        pos_vec = positions[i, :]
        if abs(pos_vec[1]) < 1e-9 && abs(pos_vec[2]) < 1e-9 # If position is along Z-axis
             axis_vec = [1.0, 0.0, 0.0] # Use X-axis instead
        else
             axis_vec = [0.0, 0.0, 1.0] # Use Z-axis
        end

        vel_dir = cross(pos_vec, axis_vec) # Direction is perpendicular to position and axis
        vel_dir_normalized = normalize(vel_dir) # Make it a unit vector

        velocities[i, :] = vel_magnitude .* vel_dir_normalized
    end

    return positions, velocities, masses
end

# --- Main Simulation Function ---

function run_simulation(num_bodies_orbiting::Int, num_steps::Int, dt::Float64)
    println("--- N-Body Simulation Start ---")
    println("Number of orbiting bodies: ", num_bodies_orbiting)
    println("Total bodies: ", num_bodies_orbiting + 1)
    println("Number of steps: ", num_steps)
    println("Time step (dt): ", dt)
    println("Softening factor: ", SOFTENING)
    println("Gravitational Constant (G): ", G)
    println("---------------------------------")

    # Initialization
    println("Initializing system...")
    central_mass = 1.0e6 # Example central mass (much larger than orbiting)
    positions, velocities, masses = initialize_circular_orbits(num_bodies_orbiting, central_mass)
    N_total = size(positions, 1)
    accelerations = zeros(Float64, N_total, 3) # Allocate acceleration buffer
    println("Initialization complete.")

    # Energy check before simulation
    println("Calculating initial energy...")
    initial_energy = calculate_total_energy(positions, velocities, masses)
    @printf("Initial Total Energy: %.6e\n", initial_energy)
    println("---------------------------------")

    # --- Simulation Loop ---
    println("Starting simulation loop...")
    start_time = time()

    for step in 1:num_steps
        kick_step_update!(positions, velocities, masses, accelerations, dt)

        # Optional: Print progress
        if step % (num_steps ÷ 10) == 0 || step == num_steps
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
    final_energy = calculate_total_energy(positions, velocities, masses)
    @printf("Final Total Energy:   %.6e\n", final_energy)

    # Accuracy check
    energy_diff = abs(final_energy - initial_energy)
    relative_error = if abs(initial_energy) > 1e-12 # Avoid division by zero if initial E is tiny
        energy_diff / abs(initial_energy)
    else
        energy_diff # Absolute difference if initial energy is near zero
    end
    @printf("Absolute Energy Difference: %.6e\n", energy_diff)
    @printf("Relative Energy Error: %.6e (%.4f%%)\n", relative_error, relative_error * 100)
    println("---------------------------------")
    @printf("Total Simulation Time: %.3f seconds\n", total_time)
    println("--- N-Body Simulation End ---")

end

# --- Run ---
# **WARNING:** N=1,000,001 is EXTREMELY computationally expensive for O(N^2)
# A single step requires ~ N^2 / 2 force calculations.
# For N=1e6, one step is ~ 0.5 * (1e6)^2 = 5e11 calculations.
# 1000 steps require ~ 5e14 calculations.
# On a modern CPU (~10 GFLOPS), this could take conceptually:
# 5e14 Flops / (10e9 Flops/sec) = 50,000 seconds ~= 14 hours PER STEP (very rough estimate)
# The actual time will vary hugely based on CPU, memory bandwidth, Julia overhead etc.
# Running with 1 million bodies sequentially is likely infeasible.
# Starting with a smaller number (e.g., 1000) is recommended for testing.

# const NUM_ORBITING_BODIES = 1_000_000 # As requested, but likely infeasible
const NUM_ORBITING_BODIES = 1000      # More reasonable number for testing
const NUM_STEPS = 1000
const DT = 0.001                      # Time step

# Check if running interactively or as script to avoid running huge sim automatically
if isinteractive()
    println("Running simulation with N_orbiting = $(NUM_ORBITING_BODIES) for $(NUM_STEPS) steps.")
    println("WARNING: If N_orbiting = 1,000,000 this will take an extremely long time.")
    println("Press Enter to continue, or Ctrl+C to abort.")
    readline()
end

run_simulation(NUM_ORBITING_BODIES, NUM_STEPS, DT)