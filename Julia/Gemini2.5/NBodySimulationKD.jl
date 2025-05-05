using LinearAlgebra
using Printf
using Random
using Base.Threads
using Statistics # For median

# --- Constants ---
const G = 1.0
const SOFTENING = 0.01
const SOFTENING_SQ = SOFTENING^2
const THETA = 0.3       # Opening angle for Barnes-Hut approximation
const THETA_SQ = THETA^2
const MAX_LEAF_SIZE = 8 # Maximum particles in a leaf node

# --- Data Structures for k-D Tree ---

struct AxisAlignedBoundingBox
    min_coord::Vector{Float64}
    max_coord::Vector{Float64}
end

# Helper function to get the center of a box
center(box::AxisAlignedBoundingBox) = 0.5 .* (box.min_coord .+ box.max_coord)

# Helper function to get the size (max dimension) of a box
box_size(box::AxisAlignedBoundingBox) = maximum(box.max_coord .- box.min_coord)

# Helper function to compute bounding box for a set of particles
function compute_bounds(indices::AbstractVector{Int}, positions::Matrix{Float64})::AxisAlignedBoundingBox
    if isempty(indices)
        # Return an empty/invalid box? Or handle upstream? Let's return zero box.
        return AxisAlignedBoundingBox(zeros(3), zeros(3))
    end
    # Use view for efficiency
    particle_positions = view(positions, indices, :)
    min_c = vec(minimum(particle_positions, dims=1))
    max_c = vec(maximum(particle_positions, dims=1))
    # Add small epsilon to avoid zero-size boxes if all particles coincide
    epsilon = 1e-9
    max_c .+= ifelse.(min_c .== max_c, epsilon, 0.0)
    return AxisAlignedBoundingBox(min_c, max_c)
end

# Helper function to merge two bounding boxes
function merge_bounds(b1::AxisAlignedBoundingBox, b2::AxisAlignedBoundingBox)::AxisAlignedBoundingBox
    min_c = min.(b1.min_coord, b2.min_coord)
    max_c = max.(b1.max_coord, b2.max_coord)
    return AxisAlignedBoundingBox(min_c, max_c)
end

mutable struct KDTreeNode
    bounds::AxisAlignedBoundingBox
    center_of_mass::Vector{Float64}
    total_mass::Float64
    particle_indices::Union{Vector{Int}, Nothing} # Indices if leaf, nothing otherwise
    children::Union{Vector{KDTreeNode}, Nothing}  # Child nodes if internal, nothing otherwise
    # Debug/Info (Optional)
    axis::Int # Splitting axis (0, 1, or 2) -> 1,2,3 in Julia
    split_value::Float64
    is_leaf::Bool

    # Constructor for internal node
    KDTreeNode(bounds, com, mass, children, axis, split_val) =
        new(bounds, com, mass, nothing, children, axis, split_val, false)

    # Constructor for leaf node
    KDTreeNode(bounds, com, mass, indices) =
        new(bounds, com, mass, indices, nothing, -1, NaN, true)
end


# --- k-D Tree Construction ---

"""
Builds the k-D tree recursively.
"""
function build_kdtree(indices::Vector{Int}, positions::Matrix{Float64}, masses::Vector{Float64}; depth::Int = 0)::KDTreeNode
    num_particles = length(indices)
    current_bounds = compute_bounds(indices, positions)

    # Base case: If few particles, create a leaf node
    if num_particles <= MAX_LEAF_SIZE
        com = zeros(3)
        total_mass = 0.0
        if num_particles > 0
            # Use views for particle data within the leaf
            leaf_positions = view(positions, indices, :)
            leaf_masses = view(masses, indices)
            total_mass = sum(leaf_masses)
            if total_mass > 1e-12 # Avoid division by zero if masses are tiny/zero
                 # Weighted average for center of mass: sum(m_i * r_i) / sum(m_i)
                 com = vec(sum(leaf_positions .* leaf_masses, dims=1) ./ total_mass)
            else # Assign geometric center if total mass is negligible
                 com = center(current_bounds)
            end
        end
        # Ensure indices are copied if needed downstream
        return KDTreeNode(current_bounds, com, total_mass, copy(indices))
    end

    # Recursive step: Choose axis and split
    axis = (depth % 3) + 1 # Cycle through axes 1, 2, 3 (x, y, z)

    # Find median coordinate along the axis
    # Using quickselect (via partialsort!) is faster O(N) than full sort O(N log N)
    coords = view(positions, indices, axis)
    median_idx = div(num_particles, 2) + 1 # Index for median value
    # partialsort! rearranges 'indices' such that the median element's original index
    # is placed correctly if indices were sorted by positions[:, axis].
    # It partitions 'indices' based on the values in 'coords'.
    median_split_value = partialsort!(copy(coords), median_idx) # Value used for split
    
    # Partition indices based on the median value
    left_indices = filter(idx -> positions[idx, axis] < median_split_value, indices)
    right_indices = filter(idx -> positions[idx, axis] >= median_split_value, indices)

    # Handle cases where median split doesn't partition well (e.g., many identical coords)
     if isempty(left_indices) || isempty(right_indices)
         # Fallback: simple split in the middle, might not be balanced
         mid = div(num_particles, 2)
         sorted_indices = sort(indices, by=idx -> positions[idx, axis])
         left_indices = sorted_indices[1:mid]
         right_indices = sorted_indices[mid+1:end]
         # Recalculate median_split_value if needed for node info, though less critical now
         if !isempty(right_indices)
            median_split_value = positions[right_indices[1], axis]
         end
     end

    # Build children recursively
    left_child = build_kdtree(left_indices, positions, masses, depth=depth+1)
    right_child = build_kdtree(right_indices, positions, masses, depth=depth+1)

    # Create internal node
    children = [left_child, right_child]
    merged_bounds = merge_bounds(left_child.bounds, right_child.bounds)
    total_mass = left_child.total_mass + right_child.total_mass
    com = zeros(3)
    if total_mass > 1e-12
         com = (left_child.center_of_mass .* left_child.total_mass .+
                right_child.center_of_mass .* right_child.total_mass) ./ total_mass
    else
         com = center(merged_bounds)
    end

    return KDTreeNode(merged_bounds, com, total_mass, children, axis, median_split_value)
end


# --- Force Calculation using k-D Tree ---

"""
Recursive helper function to traverse the tree and calculate acceleration
for a target particle. Accumulates acceleration directly into `accel_accum`.
"""
function force_traverse_recursive!(accel_accum::Vector{Float64}, # Pass accumulator by reference
                                   target_idx::Int,
                                   target_pos::SubArray{Float64}, # Use view for target position
                                   node::KDTreeNode,
                                   positions::Matrix{Float64},
                                   masses::Vector{Float64})

    # Vector from target particle to node's center of mass
    dx = node.center_of_mass[1] - target_pos[1]
    dy = node.center_of_mass[2] - target_pos[2]
    dz = node.center_of_mass[3] - target_pos[3]

    dist_sq = dx^2 + dy^2 + dz^2

    # Avoid self-interaction if target is the *only* particle in a node
    # (More robust check needed if leaf nodes can contain the target particle)
    # This check is implicitly handled by leaf node iteration later.
    # Need check to avoid zero distance if target is exactly at node COM
    if dist_sq < 1e-18 # Effectively zero distance
        return # Skip interaction
    end

    # Barnes-Hut Criterion: s^2 / d^2 < theta^2
    node_size = box_size(node.bounds)
    if (node_size^2 / dist_sq < THETA_SQ) || node.is_leaf
        # If node is far enough OR it's a leaf node:

        if node.is_leaf
            # Iterate through particles in the leaf
            for leaf_idx in node.particle_indices
                if leaf_idx == target_idx
                    continue # Skip self-interaction
                end
                # Direct calculation for particles in the leaf
                dx_leaf = positions[leaf_idx, 1] - target_pos[1]
                dy_leaf = positions[leaf_idx, 2] - target_pos[2]
                dz_leaf = positions[leaf_idx, 3] - target_pos[3]
                dist_sq_leaf = dx_leaf^2 + dy_leaf^2 + dz_leaf^2 + SOFTENING_SQ
                inv_dist_cubed_leaf = G / (dist_sq_leaf * sqrt(dist_sq_leaf))
                force_factor_leaf = masses[leaf_idx] * inv_dist_cubed_leaf

                accel_accum[1] += force_factor_leaf * dx_leaf
                accel_accum[2] += force_factor_leaf * dy_leaf
                accel_accum[3] += force_factor_leaf * dz_leaf
            end
        else
            # Treat internal node as a single pseudo-particle (if far enough)
            dist_sq_soft = dist_sq + SOFTENING_SQ
            inv_dist_cubed = G / (dist_sq_soft * sqrt(dist_sq_soft))
            force_factor = node.total_mass * inv_dist_cubed

            accel_accum[1] += force_factor * dx
            accel_accum[2] += force_factor * dy
            accel_accum[3] += force_factor * dz
        end
    else
        # Node is too close, traverse children (if internal node)
        if !node.is_leaf
            for child_node in node.children
                # Only traverse children that might contain interactions (optional optimization: check bounds intersection)
                # Basic version: always traverse children if node is too close
                 force_traverse_recursive!(accel_accum, target_idx, target_pos, child_node, positions, masses)
            end
        # If it's a leaf but somehow failed the s/d < theta check (shouldn't happen with || node.is_leaf condition),
        # the leaf particle interactions would be handled here too, but the logic above covers it.
        end
    end
end


"""
Calculates accelerations for all bodies using the k-D tree (Barnes-Hut).
Multithreaded loop over particles, each traversing the tree.
"""
function calculate_accelerations_kdtree!(accelerations::Matrix{Float64},
                                         positions::Matrix{Float64},
                                         masses::Vector{Float64})
    N = size(positions, 1)
    particle_indices = collect(1:N)

    # --- Build Tree (Sequential) ---
    # println("Building k-D tree...") # Optional debug message
    # build_time = @elapsed root_node = build_kdtree(particle_indices, positions, masses)
    root_node = build_kdtree(particle_indices, positions, masses)
    # println("Tree build time: $build_time seconds") # Optional debug message

    # --- Calculate Forces (Parallel using Tree) ---
    fill!(accelerations, 0.0) # Reset accelerations

    Threads.@threads for i in 1:N
        # Each thread calculates acceleration for its assigned particle 'i'
        # Need thread-local accumulator or direct write to accelerations[i,:]
        # Direct write is safe as only one thread writes to accelerations[i,:]
        accel_i = zeros(3) # Thread-local accumulator vector
        target_pos_view = view(positions, i, :) # Use view for efficiency
        force_traverse_recursive!(accel_i, i, target_pos_view, root_node, positions, masses)
        accelerations[i, :] = accel_i
    end
end


# --- Simulation Step ---
"""
Performs one step using Kick-Step, calling the k-D tree acceleration function.
"""
function kick_step_update!(positions::Matrix{Float64},
                          velocities::Matrix{Float64},
                          masses::Vector{Float64},
                          accelerations::Matrix{Float64},
                          dt::Float64)

    # Calculate accelerations using k-D Tree approximation
    calculate_accelerations_kdtree!(accelerations, positions, masses)

    # 1. Kick: Update velocities (v_new = v_old + a * dt)
    velocities .+= accelerations .* dt

    # 2. Step: Update positions (p_new = p_old + v_new * dt)
    positions .+= velocities .* dt
end


# --- Energy Calculation (Exact O(N^2) - FOR REFERENCE ONLY with Tree Code) ---
# NOTE: This calculates the EXACT potential energy. The tree code uses approximate forces,
# so the energy calculated by THIS function will NOT be conserved by the simulation dynamics.
# It serves as a baseline reference.

function calculate_total_energy(positions::Matrix{Float64},
                                velocities::Matrix{Float64},
                                masses::Vector{Float64})::Float64
    N = size(positions, 1)
    kinetic_energy = 0.0

    # Kinetic Energy (O(N)) - Sequential
    for i in 1:N
        vel_sq = velocities[i, 1]^2 + velocities[i, 2]^2 + velocities[i, 3]^2
        kinetic_energy += 0.5 * masses[i] * vel_sq
    end

    # Exact Potential Energy (O(N^2)) - Parallelized
    num_threads = Threads.nthreads()
    pe_partials = zeros(Float64, num_threads)

    Threads.@threads for i in 1:N
        tid = Threads.threadid()
        pe_thread = 0.0
        for j in (i+1):N
            dx = positions[j, 1] - positions[i, 1]
            dy = positions[j, 2] - positions[i, 2]
            dz = positions[j, 3] - positions[i, 3]
            dist_sq = dx^2 + dy^2 + dz^2 + SOFTENING_SQ
            dist = sqrt(dist_sq)
            pe_thread -= G * masses[i] * masses[j] / dist
        end
        pe_partials[tid] += pe_thread
    end
    potential_energy = sum(pe_partials)

    return kinetic_energy + potential_energy
end


# --- Initialization --- (No changes needed here)
function initialize_circular_orbits(N_orbiting::Int,
                                    central_mass::Float64;
                                    min_radius::Float64 = 1.0,
                                    max_radius::Float64 = 10.0,
                                    orbiting_mass_factor::Float64 = 1e-7
                                    )::Tuple{Matrix{Float64}, Matrix{Float64}, Vector{Float64}}
    N_total = N_orbiting + 1
    positions = zeros(Float64, N_total, 3)
    velocities = zeros(Float64, N_total, 3)
    masses = zeros(Float64, N_total)
    masses[1] = central_mass
    rng = MersenneTwister()
    for i in 2:N_total
        masses[i] = central_mass * orbiting_mass_factor * (0.5 + rand(rng))
        radius = min_radius + (max_radius - min_radius) * rand(rng)
        phi = 2.0 * pi * rand(rng)
        cos_theta = 2.0 * rand(rng) - 1.0
        sin_theta = sqrt(max(0.0, 1.0 - cos_theta^2))
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
        vel_dir = cross(pos_vec, axis_vec)
        if norm(vel_dir) < 1e-9
             axis_vec = [0.0, 1.0, 0.0]
             vel_dir = cross(pos_vec, axis_vec)
        end
        norm_vel_dir = norm(vel_dir)
        if norm_vel_dir > 1e-9
            vel_dir_normalized = vel_dir / norm_vel_dir
        else
            vel_dir_normalized = normalize([rand(rng)-0.5, rand(rng)-0.5, rand(rng)-0.5])
            # println("Warning: Could not determine unique velocity direction for body $i. Using random.")
        end
        velocities[i, :] = vel_magnitude .* vel_dir_normalized
    end
    return positions, velocities, masses
end

# --- Main Simulation Function ---
function run_simulation(num_bodies_orbiting::Int, num_steps::Int, dt::Float64)
    num_threads = Threads.nthreads()
    if num_threads == 1
        println("WARNING: Julia started with only 1 thread. Multithreading will be limited.")
        println("         Restart Julia with '-t auto' or '-t N' (e.g., julia -t 4 script.jl)")
    end
    println("--- N-Body Simulation Start (k-D Tree, Threads: $num_threads) ---")
    println("Algorithm: Barnes-Hut approximation with k-D Tree")
    println("Opening Angle (theta): ", THETA)
    println("Number of orbiting bodies: ", num_bodies_orbiting)
    println("Total bodies: ", num_bodies_orbiting + 1)
    println("Number of steps: ", num_steps)
    println("Time step (dt): ", dt)
    println("---------------------------------")

    println("Initializing system...")
    central_mass = 1.0e6
    positions, velocities, masses = initialize_circular_orbits(num_bodies_orbiting, central_mass)
    N_total = size(positions, 1)
    accelerations = zeros(Float64, N_total, 3)
    println("Initialization complete.")

    println("Calculating initial exact energy (Reference)...")
    initial_energy = calculate_total_energy(positions, velocities, masses) # Exact O(N^2) calc
    @printf("Initial Exact Energy: %.6e\n", initial_energy)
    println("---------------------------------")

    println("Starting simulation loop...")
    start_time = time()

    for step in 1:num_steps
        # Store positions before step for energy calc if needed, or calculate E less often
        # pos_before = copy(positions) # Example if needed

        kick_step_update!(positions, velocities, masses, accelerations, dt) # Uses k-D Tree

        if step % max(1, num_steps ÷ 10) == 0 || step == num_steps
           elapsed_time = time() - start_time
           # Optional: Calculate exact energy periodically to see drift
           # current_exact_energy = calculate_total_energy(positions, velocities, masses)
           # @printf("Step: %d / %d (%.1f%%), Elapsed: %.2f s, Exact Energy: %.6e\n",
           #         step, num_steps, 100 * step / num_steps, elapsed_time, current_exact_energy)
           @printf("Step: %d / %d (%.1f%%), Elapsed Time: %.2f s\n",
                   step, num_steps, 100 * step / num_steps, elapsed_time)
        end
    end

    end_time = time()
    total_time = end_time - start_time
    println("Simulation loop finished.")
    println("---------------------------------")

    println("Calculating final exact energy (Reference)...")
    final_energy = calculate_total_energy(positions, velocities, masses) # Exact O(N^2) calc
    @printf("Final Exact Energy:   %.6e\n", final_energy)

    energy_diff = abs(final_energy - initial_energy)
    relative_error = if abs(initial_energy) > 1e-12
        energy_diff / abs(initial_energy)
    else
        energy_diff
    end
    println("** NOTE: Energy calculated using exact O(N^2) method. **")
    println("** Simulation used approximate O(N log N) tree forces. **")
    println("** Therefore, this energy difference reflects the **")
    println("** approximation error, not just numerical integration error. **")
    @printf("Absolute Exact Energy Difference: %.6e\n", energy_diff)
    @printf("Relative Exact Energy Error: %.6e (%.4f%%)\n", relative_error, relative_error * 100)
    println("---------------------------------")
    @printf("Total Simulation Time: %.3f seconds\n", total_time)
    println("--- N-Body Simulation End ---")
end

# --- Run ---
# Remember to start Julia with threads: julia -t auto nbody_kdtree.jl

# Now N=1,000,000 might be feasible (though still demanding) due to O(N log N)
const NUM_ORBITING_BODIES = 100_000 # Reduced from 1M for faster testing
# const NUM_ORBITING_BODIES = 1_000_000 # Uncomment for full scale (will take time!)
const NUM_STEPS = 100 # Reduced steps for faster testing
# const NUM_STEPS = 1000 # Uncomment for full run
const DT = 0.001

if isinteractive() && NUM_ORBITING_BODIES > 50000
    println("Running k-D Tree simulation with N_orbiting = $(NUM_ORBITING_BODIES) for $(NUM_STEPS) steps using $(Threads.nthreads()) threads.")
    println("Theta = $THETA. This may still take significant time and memory.")
    println("Press Enter to continue, or Ctrl+C to abort.")
    readline()
elseif !isinteractive() && Threads.nthreads() == 1
     println("WARNING: Running script with only 1 thread. Performance will be limited.")
     println("         Consider restarting Julia with the -t flag for multithreading.")
end

run_simulation(NUM_ORBITING_BODIES, NUM_STEPS, DT)