# Import necessary packages
using LinearAlgebra
using Printf
using Random
using Base.Threads
using Statistics
using StaticArrays # Added for SVector/MVector

# --- Constants ---
const G = 1.0
const SOFTENING = 0.01
const SOFTENING_SQ = SOFTENING^2
const THETA = 0.3       # Opening angle for Barnes-Hut approximation
const THETA_SQ = THETA^2
# --- Tunable Parameters ---
const MAX_LEAF_SIZE = 16 # Increased leaf size << TUNABLE PARAMETER >>

# --- Data Types using StaticArrays ---
const Vec3 = SVector{3, Float64} # Use immutable StaticVector for positions/CoM
const MVec3 = MVector{3, Float64} # Use mutable StaticVector for accumulators

struct AxisAlignedBoundingBox
    min_coord::Vec3
    max_coord::Vec3
end

center(box::AxisAlignedBoundingBox) = 0.5 * (box.min_coord + box.max_coord)
box_size(box::AxisAlignedBoundingBox) = maximum(box.max_coord - box.min_coord)

# Recompute bounds using StaticArrays more efficiently
function compute_bounds(indices::AbstractVector{Int}, positions::Matrix{Float64})::AxisAlignedBoundingBox
    if isempty(indices)
        return AxisAlignedBoundingBox(Vec3(Inf, Inf, Inf), Vec3(-Inf, -Inf, -Inf)) # Empty box
    end
    # Initialize with the first particle's position
    # Note: Accessing rows and converting to SVector
    first_pos = Vec3(@view positions[indices[1], :])
    min_c = MVec3(first_pos) # Mutable copy
    max_c = MVec3(first_pos) # Mutable copy

    @inbounds for i in 2:length(indices)
        idx = indices[i]
        pos_i = Vec3(@view positions[idx, :])
        min_c .= min.(min_c, pos_i)
        max_c .= max.(max_c, pos_i)
    end

    # Add small epsilon for robustness against coincident particles
    epsilon = 1e-9
    max_c .+= ifelse.(min_c .== max_c, epsilon, 0.0)

    return AxisAlignedBoundingBox(Vec3(min_c), Vec3(max_c)) # Convert back to SVector
end

function merge_bounds(b1::AxisAlignedBoundingBox, b2::AxisAlignedBoundingBox)::AxisAlignedBoundingBox
    min_c = min.(b1.min_coord, b2.min_coord)
    max_c = max.(b1.max_coord, b2.max_coord)
    return AxisAlignedBoundingBox(min_c, max_c)
end

mutable struct KDTreeNode
    bounds::AxisAlignedBoundingBox
    center_of_mass::Vec3
    total_mass::Float64
    particle_indices::Union{UnitRange{Int}, Nothing} # Store range within sorted indices if leaf
    children::Union{Vector{KDTreeNode}, Nothing}
    # Debug/Info (Optional)
    # axis::Int # Splitting axis (1, 2, or 3)
    # split_value::Float64
    is_leaf::Bool

    # Simplified constructors
    KDTreeNode(bounds, com, mass, children) = new(bounds, com, mass, nothing, children, false) # Internal
    KDTreeNode(bounds, com, mass, indices_range) = new(bounds, com, mass, indices_range, nothing, true) # Leaf
end


# --- k-D Tree Construction (Optimized Allocations) ---

"""
Builds the k-D tree recursively using in-place partitioning (implicitly via sorting views).
Operates on a shared `sorted_indices` array, passing views down.
"""
function build_kdtree_recursive!(
    node_idx_range::UnitRange{Int}, # Range within sorted_indices array this node represents
    sorted_indices::Vector{Int},    # Shared array of particle indices, sorted along axes
    positions::Matrix{Float64},
    masses::Vector{Float64};
    depth::Int = 0
)::KDTreeNode

    num_particles = length(node_idx_range)
    current_indices_view = view(sorted_indices, node_idx_range)
    current_bounds = compute_bounds(current_indices_view, positions)

    # Base case: Leaf node
    if num_particles <= MAX_LEAF_SIZE
        com = Vec3(0.0, 0.0, 0.0)
        total_mass = 0.0
        if num_particles > 0
            # Calculate CoM and total mass for the leaf
            m_com_accum = MVec3(0.0, 0.0, 0.0)
            @inbounds for i in node_idx_range
                idx = sorted_indices[i]
                m = masses[idx]
                total_mass += m
                pos_idx = Vec3(@view positions[idx, :])
                m_com_accum .+= m .* pos_idx
            end
            if total_mass > 1e-12
                com = Vec3(m_com_accum ./ total_mass)
            else
                com = center(current_bounds)
            end
        end
        # Leaf node stores the *range* within the sorted_indices array
        return KDTreeNode(current_bounds, com, total_mass, node_idx_range)
    end

    # Recursive step: Choose axis and partition
    axis = (depth % 3) + 1
    
    # Sort the current slice of 'sorted_indices' *in-place* along the chosen axis
    # This is the main partitioning step - avoids filter allocations
    sort!(current_indices_view, by = idx -> @inbounds positions[idx, axis])

    # Find median index within the current range
    median_local_idx = div(num_particles, 2) # Use integer division for split point
    median_global_idx = node_idx_range.start + median_local_idx -1 # Adjust to global index in sorted_indices

    # Define index ranges for children (views into the now-sorted range)
    left_range = node_idx_range.start : median_global_idx
    right_range = (median_global_idx + 1) : node_idx_range.stop

    # Build children recursively (passing views/ranges)
    left_child = build_kdtree_recursive!(left_range, sorted_indices, positions, masses, depth=depth+1)
    right_child = build_kdtree_recursive!(right_range, sorted_indices, positions, masses, depth=depth+1)

    # Create internal node
    children = [left_child, right_child]
    merged_bounds = merge_bounds(left_child.bounds, right_child.bounds)
    total_mass = left_child.total_mass + right_child.total_mass
    com = Vec3(0.0, 0.0, 0.0)
    if total_mass > 1e-12
         # CoM weighted average using StaticArrays
         com = (left_child.center_of_mass * left_child.total_mass +
                right_child.center_of_mass * right_child.total_mass) / total_mass
    else
         com = center(merged_bounds)
    end

    return KDTreeNode(merged_bounds, com, total_mass, children)
end

# Wrapper to initialize tree building
function build_kdtree(N::Int, positions::Matrix{Float64}, masses::Vector{Float64})::KDTreeNode
    sorted_indices = collect(1:N) # Initial indices, will be sorted in place recursively
    return build_kdtree_recursive!(1:N, sorted_indices, positions, masses)
end


# --- Force Calculation using k-D Tree (Optimized) ---

"""
Recursive helper using StaticArrays and @inbounds. Modifies MVec3 accumulator.
"""
@inline function force_traverse_recursive!(
    accel_accum::MVec3, # Mutable accumulator
    target_idx::Int,
    target_pos::Vec3,   # Immutable position
    node::KDTreeNode,
    positions::Matrix{Float64},
    masses::Vector{Float64},
    sorted_indices::Vector{Int} # Needed to access leaf particles
)

    # Vector from target particle to node's center of mass (StaticArray subtraction)
    dr = node.center_of_mass - target_pos
    dist_sq = sum(abs2, dr) # Faster than dx^2+dy^2+dz^2? Equivalent. dot(dr, dr) also good.

    if dist_sq < 1e-18 # Avoid division by zero / huge forces if exactly at CoM
        return
    end

    node_size = box_size(node.bounds)

    # Criterion: s^2 / d^2 < theta^2
    if node.is_leaf || (node_size^2 < THETA_SQ * dist_sq)
        if node.is_leaf
            # Iterate through particles in the leaf using the stored range
            @inbounds for i in node.particle_indices
                leaf_idx = sorted_indices[i]
                if leaf_idx == target_idx
                    continue # Skip self-interaction
                end
                # Direct calculation
                leaf_pos = Vec3(@view positions[leaf_idx, :])
                dr_leaf = leaf_pos - target_pos
                dist_sq_leaf = sum(abs2, dr_leaf) + SOFTENING_SQ
                
                # Check distance again (optional, softening might handle it)
                if dist_sq_leaf < 1e-18 continue end

                inv_dist = 1.0 / sqrt(dist_sq_leaf) # Calculate 1/r
                inv_dist_cubed = inv_dist * inv_dist * inv_dist # 1/r^3
                
                # a_i += G * m_j * (r_j - r_i) / |r_j - r_i|^3 (softened)
                force_factor = G * masses[leaf_idx] * inv_dist_cubed
                accel_accum .+= force_factor .* dr_leaf # Use broadcasting '.'
            end
        else # Treat internal node as pseudo-particle (if far enough)
            dist_sq_soft = dist_sq + SOFTENING_SQ
            inv_dist = 1.0 / sqrt(dist_sq_soft)
            inv_dist_cubed = inv_dist * inv_dist * inv_dist
            force_factor = G * node.total_mass * inv_dist_cubed
            accel_accum .+= force_factor .* dr
        end
    else # Node is too close (and not a leaf), traverse children
         # No need for !node.is_leaf check due to logic order
        @inbounds begin # Assuming children is always length 2 for binary tree
             force_traverse_recursive!(accel_accum, target_idx, target_pos, node.children[1], positions, masses, sorted_indices)
             force_traverse_recursive!(accel_accum, target_idx, target_pos, node.children[2], positions, masses, sorted_indices)
        end
    end
end


"""
Calculates accelerations using k-D tree, optimized with pre-allocation.
"""
function calculate_accelerations_kdtree!(
    accelerations::Matrix{Float64},
    positions::Matrix{Float64},
    masses::Vector{Float64},
    accel_partials::Vector{MVec3} # Pre-allocated thread-local storage
)
    N = size(positions, 1)

    # --- Build Tree (Optimized) ---
    sorted_indices = collect(1:N) # Prepare index array for sorting during build
    root_node = build_kdtree_recursive!(1:N, sorted_indices, positions, masses)

    # --- Calculate Forces (Parallel using Tree & Pre-allocated Accumulators) ---
    Threads.@threads for i in 1:N
        tid = Threads.threadid()
        # Reset pre-allocated accumulator for this thread's particle
        accel_accum = accel_partials[tid]
        accel_accum .= 0.0 # Reset using broadcast

        @inbounds target_pos = Vec3(@view positions[i, :]) # Get target pos as SVector

        # Call traversal, writing result into accel_accum
        force_traverse_recursive!(accel_accum, i, target_pos, root_node, positions, masses, sorted_indices)

        # Write result back to the main accelerations matrix
        @inbounds accelerations[i, :] = accel_accum
    end
end


# --- Simulation Step (Unchanged Logic) ---
function kick_step_update!(positions::Matrix{Float64},
                          velocities::Matrix{Float64},
                          masses::Vector{Float64},
                          accelerations::Matrix{Float64},
                          dt::Float64,
                          accel_partials::Vector{MVec3}) # Pass storage

    # Calculate accelerations using k-D Tree approximation
    calculate_accelerations_kdtree!(accelerations, positions, masses, accel_partials)

    # 1. Kick & 2. Step (Vectorized updates)
    # Use @.. macro for fusion, potentially faster for simple broadcasts
    @.. velocities += accelerations * dt
    @.. positions += velocities * dt
end

# --- Energy Calculation (Exact O(N^2) - Optimized with StaticArrays) ---
function calculate_total_energy(positions::Matrix{Float64},
                                velocities::Matrix{Float64},
                                masses::Vector{Float64})::Float64
    N = size(positions, 1)
    kinetic_energy = 0.0

    # Kinetic Energy (O(N)) - Sequential, using StaticArrays internally might help slightly
    @inbounds for i in 1:N
        # vel_sq = sum(abs2, @view velocities[i, :]) # Slightly slower due to view?
        vel_sq = velocities[i, 1]^2 + velocities[i, 2]^2 + velocities[i, 3]^2
        kinetic_energy += 0.5 * masses[i] * vel_sq
    end

    # Exact Potential Energy (O(N^2)) - Parallelized, maybe StaticArrays help slightly
    num_threads = Threads.nthreads()
    pe_partials = zeros(Float64, num_threads)

    Threads.@threads for i in 1:N
        tid = Threads.threadid()
        pe_thread = 0.0
        # Convert target position to SVector once per outer loop
        @inbounds target_pos_i = Vec3(@view positions[i, :])
        @inbounds mi = masses[i]

        @inbounds for j in (i+1):N
            target_pos_j = Vec3(@view positions[j, :])
            dr = target_pos_j - target_pos_i
            dist_sq = sum(abs2, dr) + SOFTENING_SQ

            # Avoid sqrt if possible, or ensure non-zero before sqrt
            if dist_sq > 1e-18
                 dist = sqrt(dist_sq)
                 pe_thread -= G * mi * masses[j] / dist
            end
        end
        # This += is safe as each thread writes to its own accumulator index
        pe_partials[tid] += pe_thread
    end
    potential_energy = sum(pe_partials)

    return kinetic_energy + potential_energy
end


function initialize_circular_orbits(N_orbiting::Int,
                                central_mass::Float64;
                                min_radius::Float64 = 1.0,
                                max_radius::Float64 = 10.0,
                                orbiting_mass_factor::Float64 = 1e-7
                                )::Tuple{Matrix{Float64}, Matrix{Float64}, Vector{Float64}}
    N_total = N_orbiting + 1
    # Keep main storage as Matrix/Vector for compatibility, convert to SVector inside loops
    positions = zeros(Float64, N_total, 3)
    velocities = zeros(Float64, N_total, 3)
    masses = zeros(Float64, N_total)
    masses[1] = central_mass
    rng = MersenneTwister() # Or use default RNG

    for i in 2:N_total
        masses[i] = central_mass * orbiting_mass_factor * (0.5 + rand(rng))
        radius = min_radius + (max_radius - min_radius) * rand(rng)
        phi = 2.0 * pi * rand(rng)
        cos_theta = 2.0 * rand(rng) - 1.0
        sin_theta = sqrt(max(0.0, 1.0 - cos_theta^2))

        # Calculate position components
        x = radius * sin_theta * cos(phi)
        y = radius * sin_theta * sin(phi)
        z = radius * cos_theta
        pos_vec = Vec3(x, y, z) # Create StaticVector
        positions[i, :] = pos_vec # Assign back to Matrix row

        # Calculate velocity
        vel_magnitude = sqrt(G * central_mass / max(radius, 1e-9)) # Avoid division by zero

        # Determine velocity direction (perpendicular)
        axis_vec = (abs(x) < 1e-9 && abs(y) < 1e-9) ? Vec3(1.0, 0.0, 0.0) : Vec3(0.0, 0.0, 1.0)
        vel_dir_cross = cross(pos_vec, axis_vec)

        if sum(abs2, vel_dir_cross) < 1e-18 # If parallel, try another axis
            axis_vec = Vec3(0.0, 1.0, 0.0)
            vel_dir_cross = cross(pos_vec, axis_vec)
        end

        # Normalize velocity direction
        norm_vel_dir = norm(vel_dir_cross)
        vel_dir_normalized = if norm_vel_dir > 1e-9
            vel_dir_cross / norm_vel_dir
        else # Fallback for zero position vector etc.
            normalize(Vec3(rand(rng)-0.5, rand(rng)-0.5, rand(rng)-0.5))
        end

        vel_vec = vel_magnitude * vel_dir_normalized
        velocities[i, :] = vel_vec # Assign back to Matrix row
    end
    return positions, velocities, masses
end

# --- Main Simulation Function (Adjusted for pre-allocation) ---
function run_simulation(num_bodies_orbiting::Int, num_steps::Int, dt::Float64)
    num_threads = Threads.nthreads()
    # ... (Startup messages as before) ...
    println("--- N-Body Simulation Start (k-D Tree Optimized, Threads: $num_threads) ---")
    println("Algorithm: Barnes-Hut (theta = $THETA), MAX_LEAF_SIZE = $MAX_LEAF_SIZE")
    # ... (Parameter printing) ...

    println("Initializing system...")
    central_mass = 1.0e6
    positions, velocities, masses = initialize_circular_orbits(num_bodies_orbiting, central_mass)
    N_total = size(positions, 1)
    accelerations = zeros(Float64, N_total, 3)

    # --- Pre-allocate thread-local storage for acceleration accumulators ---
    accel_partials = [zeros(MVec3) for _ in 1:num_threads]
    println("Initialization complete.")

    println("Calculating initial exact energy (Reference)...")
    initial_energy = calculate_total_energy(positions, velocities, masses)
    @printf("Initial Exact Energy: %.6e\n", initial_energy)
    println("---------------------------------")

    println("Starting simulation loop...")
    start_time = time()

    for step in 1:num_steps
        # Pass pre-allocated storage to the update function
        kick_step_update!(positions, velocities, masses, accelerations, dt, accel_partials)

        # ... (Progress printing as before) ...
    end

    end_time = time()
    # ... (Final energy calculation and output as before) ...
    println("** NOTE: Energy calculation uses exact O(N^2) method for reference. **")
    println("** Simulation dynamics used O(N log N) tree approximation. **")
    # ... (Rest of the final output) ...

end

# --- Run ---
# ... (Run parameters and execution logic as before) ...
# Remember: julia -t auto script.jl

const NUM_ORBITING_BODIES = 100_000 # Keep reduced for testing
const NUM_STEPS = 10             # Keep reduced for testing
const DT = 0.001

# ... (Interactive check / run call as before) ...
if isinteractive() && NUM_ORBITING_BODIES > 50000
    println("Running k-D Tree simulation with N_orbiting = $(NUM_ORBITING_BODIES) for $(NUM_STEPS) steps using $(Threads.nthreads()) threads.")
    println("Theta = $THETA, MAX_LEAF_SIZE = $MAX_LEAF_SIZE. This may still take significant time and memory.")
    println("Press Enter to continue, or Ctrl+C to abort.")
    readline()
elseif !isinteractive() && Threads.nthreads() == 1
     println("WARNING: Running script with only 1 thread. Performance will be limited.")
     println("         Consider restarting Julia with the -t flag for multithreading.")
end

run_simulation(NUM_ORBITING_BODIES, NUM_STEPS, DT)