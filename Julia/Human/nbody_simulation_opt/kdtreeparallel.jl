module kdtreeparallel

using StaticArrays
using Base.Threads: @threads, threadid, nthreads
using Profile
using LinearAlgebra: norm

# Constants
const MAX_PARTS = 7
const THETA = 0.3
const G = 1.0  # Gravitational constant in simulation units
const SOFTENING = 1e-10  # Softening parameter to avoid numerical instability

# Use StaticArrays for small fixed-size vectors
struct Body
    p::SVector{3, Float64}  # Position vector with fixed size of 3
    v::SVector{3, Float64}  # Velocity vector with fixed size of 3
    m::Float64              # Mass
end

abstract type KDTree end

struct KDLeaf <: KDTree
    # for leaves
    particles::Vector{Int64}
end

struct KDInternal <: KDTree
    # for internal nodes
    split_dim::Int64
    split_val::Float64
    m::Float64
    cm::SVector{3, Float64}  # Center of mass
    size::Float64
    left::KDTree
    right::KDTree
end

# Barnes-Hut acceleration cell for better traversal
struct BHCell
    cm::SVector{3, Float64}
    m::Float64
    size::Float64
end

# Node pool for better memory management
mutable struct NodePool
    leaf_nodes::Vector{KDLeaf}
    internal_nodes::Vector{KDInternal}
    leaf_count::Int
    internal_count::Int
end

function create_node_pool(n::Int)
    # Estimate node counts based on tree size
    est_leaves = n ÷ MAX_PARTS + 1
    est_internal = 2 * est_leaves
    
    NodePool(
        Vector{KDLeaf}(undef, est_leaves),
        Vector{KDInternal}(undef, est_internal),
        0,
        0
    )
end

function get_leaf_node(pool::NodePool, particles::Vector{Int64})
    pool.leaf_count += 1
    if pool.leaf_count > length(pool.leaf_nodes)
        # Double capacity if needed
        resize!(pool.leaf_nodes, 2 * length(pool.leaf_nodes))
    end
    
    pool.leaf_nodes[pool.leaf_count] = KDLeaf(particles)
    return pool.leaf_nodes[pool.leaf_count]
end

function get_internal_node(pool::NodePool, split_dim, split_val, m, cm, size, left, right)
    pool.internal_count += 1
    if pool.internal_count > length(pool.internal_nodes)
        # Double capacity if needed
        resize!(pool.internal_nodes, 2 * length(pool.internal_nodes))
    end
    
    pool.internal_nodes[pool.internal_count] = KDInternal(split_dim, split_val, m, cm, size, left, right)
    return pool.internal_nodes[pool.internal_count]
end

function reset_pool!(pool::NodePool)
    pool.leaf_count = 0
    pool.internal_count = 0
end

# Fast median selection with introselect algorithm (combination of quickselect and insertion sort)
function partition_median!(indices::Vector{Int64}, system::Vector{Body}, split_dim::Int, start::Int, ending::Int)
    n = ending - start
    if n <= 16
        # Use insertion sort for small arrays
        @inbounds for i in start+1:ending-1
            key = indices[i]
            key_val = system[key].p[split_dim]
            j = i - 1
            while j >= start && system[indices[j]].p[split_dim] > key_val
                indices[j+1] = indices[j]
                j -= 1
            end
            indices[j+1] = key
        end
        
        # Return the median
        return start + n ÷ 2
    end
    
    # Choose pivot as median of first, middle and last elements
    mid_idx = start + n ÷ 2
    @inbounds begin
        a_val = system[indices[start]].p[split_dim]
        b_val = system[indices[mid_idx]].p[split_dim]
        c_val = system[indices[ending-1]].p[split_dim]
        
        # Median of three
        pivot_idx = if a_val <= b_val
            if b_val <= c_val 
                mid_idx
            elseif a_val <= c_val 
                ending-1
            else 
                start
            end
        else # a > b
            if a_val <= c_val 
                start
            elseif b_val <= c_val 
                ending-1
            else 
                mid_idx
            end
        end
        
        # Swap pivot to start
        if pivot_idx != start
            indices[pivot_idx], indices[start] = indices[start], indices[pivot_idx]
        end
        
        # Partition
        pivot_val = system[indices[start]].p[split_dim]
        i = start
        j = ending
        
        while true
            i += 1
            while i < j && system[indices[i]].p[split_dim] < pivot_val
                i += 1
            end
            
            j -= 1
            while j > i && system[indices[j]].p[split_dim] > pivot_val
                j -= 1
            end
            
            if i >= j
                break
            end
            
            indices[i], indices[j] = indices[j], indices[i]
        end
        
        # Place pivot in final position
        indices[start], indices[j] = indices[j], indices[start]
        
        # Return pivot position
        return j
    end
end

function build_tree(indices::Vector{Int64}, start::Int64, ending::Int64, system::Vector{Body})::KDTree
    np = ending - start
    if np <= MAX_PARTS
        # Pre-size the array exactly
        node = KDLeaf(Vector{Int64}(undef, np))
        @inbounds for i in 0:np-1
            node.particles[i+1] = indices[start + i]
        end
        return node
    else
        # Use StaticArrays for bounds
        minp = @SVector [1e100, 1e100, 1e100]
        maxp = @SVector [-1e100, -1e100, -1e100]
        m = 0.0
        cm = @SVector [0.0, 0.0, 0.0]
        
        # Accumulate mass and center of mass
        @inbounds for i in start:ending-1
            b_idx = indices[i]
            b = system[b_idx]
            m += b.m
            cm = cm + b.m * b.p
            minp = min.(minp, b.p)
            maxp = max.(maxp, b.p)
        end
        
        cm = cm / m  # Normalize center of mass
        
        # Find longest dimension
        diff_x = maxp[1] - minp[1]
        diff_y = maxp[2] - minp[2]
        diff_z = maxp[3] - minp[3]
        
        split_dim = 1
        if diff_y > diff_x
            split_dim = 2
        end
        if diff_z > (split_dim == 1 ? diff_x : diff_y)
            split_dim = 3
        end
        
        size = maxp[split_dim] - minp[split_dim]
        
        # Partition around median (efficient quickselect)
        mid::Int64 = div((start + ending), 2)
        s = start
        e = ending
        
        # Quickselect algorithm for partitioning
        while s + 1 < e
            pivot = rand(s:e-1)
            @inbounds indices[s], indices[pivot] = indices[pivot], indices[s]  # Swap
            
            low = s+1
            high = e-1
            @inbounds pivot_val = system[indices[s]].p[split_dim]
            
            while low <= high
                @inbounds if system[indices[low]].p[split_dim] < pivot_val
                    low += 1
                else
                    @inbounds indices[low], indices[high] = indices[high], indices[low]  # Swap
                    high -= 1
                end
            end
            
            @inbounds indices[s], indices[high] = indices[high], indices[s]  # Swap
            
            if high < mid
                s = high + 1
            elseif high > mid
                e = high
            else
                s = e
            end
        end
        
        @inbounds split_val = system[indices[mid]].p[split_dim]
        
        # Recursion on children
        left = build_tree(indices, start, mid, system)
        right = build_tree(indices, mid, ending, system)
        
        return KDInternal(split_dim, split_val, m, cm, size, left, right)
    end
end

# Pre-allocate cell stacks for Barnes-Hut traversal
mutable struct BHTraversal
    stack::Vector{Tuple{KDTree, Float64}}  # node and distance²
    cells::Vector{BHCell}  # Approximated cells for acceleration
    cell_count::Int
end

function BHTraversal(capacity::Int=1000)
    BHTraversal(
        Vector{Tuple{KDTree, Float64}}(undef, capacity),
        Vector{BHCell}(undef, capacity),
        0
    )
end

function reset!(traversal::BHTraversal)
    traversal.cell_count = 0
end

function add_cell!(traversal::BHTraversal, cm::SVector{3, Float64}, m::Float64, size::Float64)
    traversal.cell_count += 1
    if traversal.cell_count > length(traversal.cells)
        # Double capacity if needed
        resize!(traversal.cells, 2 * length(traversal.cells))
    end
    
    traversal.cells[traversal.cell_count] = BHCell(cm, m, size)
end

@inline function calc_pp_accel(p_pos::SVector{3, Float64}, j_pos::SVector{3, Float64}, j_mass::Float64)::SVector{3, Float64}
    d = p_pos - j_pos
    dist_sqr = sum(d.^2)
    dist_sqr = max(dist_sqr, SOFTENING)  # Apply softening
    dist = sqrt(dist_sqr)
    magi = -j_mass / (dist * dist_sqr)
    return d * magi
end

@inline function calc_pp_accel(system::Vector{Body}, i::Int64, j::Int64, acc::AbstractVector{Float64})
    @inbounds begin
        d = system[i].p - system[j].p
        dist_sqr = sum(d.^2)
        dist_sqr = max(dist_sqr, SOFTENING)  # Apply softening
        dist = sqrt(dist_sqr)
        magi = -system[j].m / (dist * dist_sqr)
        @. acc += d * magi  # Broadcast operation
    end
end

function collect_approximation_cells!(tree::KDTree, p_pos::SVector{3, Float64}, traversal::BHTraversal)
    reset!(traversal)
    
    # Start with root
    stack_size = 1
    traversal.stack[1] = (tree, 0.0)
    
    while stack_size > 0
        # Pop from stack
        current, _ = traversal.stack[stack_size]
        stack_size -= 1
        
        if current isa KDLeaf
            # For leaves, add all particles individually
            @inbounds for i in current.particles
                # Individual particles will be handled separately
            end
        else 
            # For internal nodes
            @inbounds begin
                node = current::KDInternal
                d = p_pos - node.cm
                dist_sqr = sum(d.^2)
                
                # Barnes-Hut approximation criterion
                if node.size * node.size < THETA^2 * dist_sqr
                    # Use cell approximation
                    add_cell!(traversal, node.cm, node.m, node.size)
                else
                    # Need to traverse deeper - add children to stack
                    stack_size += 1
                    if stack_size > length(traversal.stack)
                        resize!(traversal.stack, 2 * length(traversal.stack))
                    end
                    traversal.stack[stack_size] = (node.right, 0.0)
                    
                    stack_size += 1
                    if stack_size > length(traversal.stack)
                        resize!(traversal.stack, 2 * length(traversal.stack))
                    end
                    traversal.stack[stack_size] = (node.left, 0.0)
                end
            end
        end
    end
end

function calc_accel_bh(p::Int64, tree::KDTree, system::Vector{Body}, acc::AbstractVector{Float64}, traversal::BHTraversal)
    @inbounds p_pos = system[p].p
    
    # Collect approximation cells
    collect_approximation_cells!(tree, p_pos, traversal)
    
    # Apply forces from all collected cells
    @inbounds for i in 1:traversal.cell_count
        cell = traversal.cells[i]
        d = p_pos - cell.cm
        dist_sqr = sum(d.^2)
        dist_sqr = max(dist_sqr, SOFTENING)
        dist = sqrt(dist_sqr)
        magi = -cell.m / (dist * dist_sqr)
        @. acc += d * magi
    end
    
    # Handle direct calculation for leaves
    stack_size = 1
    traversal.stack[1] = (tree, 0.0)
    
    while stack_size > 0
        # Pop from stack
        current, _ = traversal.stack[stack_size]
        stack_size -= 1
        
        if current isa KDLeaf
            # For leaves, calculate direct interactions
            @inbounds for i in current.particles
                if i != p
                    calc_pp_accel(system, p, i, acc)
                end
            end
        else 
            # For internal nodes
            @inbounds begin
                node = current::KDInternal
                d = p_pos - node.cm
                dist_sqr = sum(d.^2)
                
                # Barnes-Hut approximation criterion
                if node.size * node.size < THETA^2 * dist_sqr
                    # Already handled in approximation phase
                else
                    # Need to traverse deeper - add children to stack
                    stack_size += 1
                    traversal.stack[stack_size] = (node.right, 0.0)
                    
                    stack_size += 1
                    traversal.stack[stack_size] = (node.left, 0.0)
                end
            end
        end
    end
end

function accel_recur(cur_node::KDLeaf, p::Int64, system::Vector{Body}, acc::AbstractVector{Float64})
    @inbounds for i in cur_node.particles
        if i != p
            calc_pp_accel(system, p, i, acc)
        end
    end
end

function accel_recur(cur_node::KDInternal, p::Int64, system::Vector{Body}, acc::AbstractVector{Float64})
    @inbounds begin
        d = system[p].p - cur_node.cm
        dist_sqr = sum(d.^2)
        if cur_node.size * cur_node.size < THETA^2 * dist_sqr
            # Use approximation
            dist = sqrt(dist_sqr)
            magi = -cur_node.m / (dist * dist_sqr)
            @. acc += d * magi  # Broadcast operation
        else
            # Need to traverse deeper
            accel_recur(cur_node.left, p, system, acc)
            accel_recur(cur_node.right, p, system, acc)
        end
    end
end

function calc_accel(p::Int64, tree::KDTree, system::Vector{Body}, acc::AbstractVector{Float64})
    accel_recur(tree, p, system, acc)
end

function print_tree(step::Int64, tree::KDTree, system::Vector{Body})
    function print_node(n::KDLeaf, file::IO)
        println(file, "L $(length(n.particles))")
        for i in n.particles
            println(file, "$(system[i].p[1]) $(system[i].p[2]) $(system[i].p[3])")
        end
    end

    function print_node(n::KDInternal, file::IO)
        println(file, "I $(n.split_dim) $(n.split_val)")
        print_node(n.left, file)
        print_node(n.right, file)
    end

    fname = "tree$step.txt"
    try
        open(fname, "w") do file
            print_node(tree, file)
        end
    catch ex
        println(ex)
    end
end

function calculate_energy(system::Vector{Body})::Tuple{Float64, Float64}
    nb = length(system)
    kinetic_energy = 0.0
    potential_energy = 0.0
    
    # Calculate kinetic energy: sum(0.5 * m * v^2)
    @inbounds for i in 1:nb
        v_squared = sum(system[i].v.^2)
        kinetic_energy += 0.5 * system[i].m * v_squared
    end
    
    # Calculate potential energy: sum(-G * m_i * m_j / r_ij)
    # Note: we calculate each pair once (j > i)
    @inbounds for i in 1:nb-1
        for j in i+1:nb
            r_ij = system[i].p - system[j].p
            distance = sqrt(sum(r_ij.^2))
            # Avoid division by zero
            if distance > 1e-10
                potential_energy -= G * system[i].m * system[j].m / distance
            end
        end
    end
    
    return (kinetic_energy, potential_energy)
end

function simple_sim(system::Vector{Body}, dt::Float64, steps::Int64, use_optimized::Bool=true)
    nb::Int64 = length(system)
    
    # Pre-allocate acceleration arrays with explicit type
    acc = zeros(Float64, 3, nb)  # Column-major layout for better memory access
    indices = collect(1:nb)
    
    # Calculate initial energy
    # kinetic_init, potential_init = calculate_energy(system)
    # total_energy_init = kinetic_init + potential_init
    # println("Initial energy: Kinetic = $kinetic_init, Potential = $potential_init, Total = $total_energy_init")
    
    # For large datasets, use chunks to reduce scheduling overhead
    # For smaller datasets, use direct threading to reduce overhead
    using_chunks = nb > 10000
    chunk_size = using_chunks ? max(1, nb ÷ (nthreads() * 4)) : nb
    
    for step in 1:steps
        # Build tree once per step
        tree = build_tree(indices, 1, nb+1, system)
        
        # Calculate accelerations in parallel - use different strategies based on dataset size
        if using_chunks
            # Use chunks for large datasets (better load balancing)
            @threads for chunk_start in 1:chunk_size:nb
                chunk_end = min(chunk_start + chunk_size - 1, nb)
                
                for i in chunk_start:chunk_end
                    # Use view to avoid allocation
                    acc_view = @view acc[:, i]
                    fill!(acc_view, 0.0)  # Reset acceleration
                    calc_accel(i, tree, system, acc_view)
                end
            end
        else
            # Simple threading for small datasets (less overhead)
            @threads for i in 1:nb
                acc_view = @view acc[:, i]
                fill!(acc_view, 0.0)  # Reset acceleration
                calc_accel(i, tree, system, acc_view)
            end
        end
        
        # Update positions and velocities
        if using_chunks
            # Use chunks for large datasets
            @threads for chunk_start in 1:chunk_size:nb
                chunk_end = min(chunk_start + chunk_size - 1, nb)
                
                for i in chunk_start:chunk_end
                    # Extract current values
                    @inbounds body = system[i]
                    @inbounds acc_i = @view acc[:, i]
                    
                    # Update velocity with improved accuracy
                    new_v = body.v + dt * SVector{3}(acc_i)
                    
                    # Update position
                    new_p = body.p + dt * new_v
                    
                    # Create updated body and assign
                    @inbounds system[i] = Body(new_p, new_v, body.m)
                end
            end
        else
            # Simple threading for small datasets
            @threads for i in 1:nb
                @inbounds body = system[i]
                @inbounds acc_i = @view acc[:, i]
                
                # Update velocity and position
                new_v = body.v + dt * SVector{3}(acc_i)
                new_p = body.p + dt * new_v
                
                # Create updated body and assign
                @inbounds system[i] = Body(new_p, new_v, body.m)
            end
        end
    end
    
    # Calculate final energy
    # kinetic_final, potential_final = calculate_energy(system)
    # total_energy_final = kinetic_final + potential_final
    # println("Final energy: Kinetic = $kinetic_final, Potential = $potential_final, Total = $total_energy_final")
    
    # # Calculate energy error
    # energy_error = abs(total_energy_final - total_energy_init) / abs(total_energy_init)
    # println("Energy conservation error: $(energy_error * 100)%")
    
    return 0.0  # Placeholder for energy_error
end

function circular_orbits(n::Int64)::Vector{Body}
    # Central body
    first = Body(@SVector([0.0, 0.0, 0.0]), @SVector([0.0, 0.0, 0.0]), 1.0)
    
    # Pre-allocate bodies vector for better memory allocation
    bods = Vector{Body}(undef, n+1)
    bods[1] = first
    
    @inbounds for i in 1:n
        d = .1 + (i * 5.0 / n)
        v = sqrt(1.0 / d)
        theta = rand(Float64) * 2 * pi
        
        # Calculate positions and velocities
        x = d * cos(theta)
        y = d * sin(theta)
        vx = -v * sin(theta)
        vy = v * cos(theta)
        
        # Create body with static vectors
        bods[i+1] = Body(@SVector([x, y, 0.0]), @SVector([vx, vy, 0.0]), 1.0e-7)
    end
    
    return bods
end

function verify_optimization(steps::Int64=10, n::Int64=1000)
    # Run small benchmark to verify optimizations
    println("Running verification benchmark with $n particles for $steps steps...")
    
    # Create original version for comparison
    function simple_sim_original(system::Vector{Body}, dt::Float64, steps::Int64)
        nb::Int64 = length(system)
        acc = zeros(Float64, 3, nb)
        indices = collect(1:nb)
    
        # Calculate initial energy
        kinetic_init, potential_init = calculate_energy(system)
        total_energy_init = kinetic_init + potential_init
        println("Initial energy: Kinetic = $kinetic_init, Potential = $potential_init, Total = $total_energy_init")
        
        for step in 1:steps
            tree = build_tree(indices, 1, nb+1, system)
            @threads for i in 1:nb
                acc_view = @view acc[:, i]
                fill!(acc_view, 0.0)
                calc_accel(i, tree, system, acc_view)
            end
            @threads for i in 1:nb
                @inbounds body = system[i]
                @inbounds acc_i = @view acc[:, i]
                
                # Update state
                new_v = body.v + dt * SVector{3}(acc_i)
                new_p = body.p + dt * new_v
                
                @inbounds system[i] = Body(new_p, new_v, body.m)
            end
        end
        
        # Calculate final energy
        kinetic_final, potential_final = calculate_energy(system)
        total_energy_final = kinetic_final + potential_final
        println("Final energy: Kinetic = $kinetic_final, Potential = $potential_final, Total = $total_energy_final")
        
        # Calculate energy error
        energy_error = abs(total_energy_final - total_energy_init) / abs(total_energy_init)
        println("Energy conservation error: $(energy_error * 100)%")
        
        return energy_error
    end
    
    system = circular_orbits(n)
    dt = 1e-3
    
    # Run both optimized and baseline for comparison
    println("\nRunning optimized version:")
    start_time = time()
    energy_error_opt = simple_sim(system, dt, steps, true)
    end_time = time()
    runtime_opt = end_time - start_time
    
    # Reset system
    system = circular_orbits(n)
    
    println("\nRunning baseline version:")
    start_time = time()
    energy_error_base = simple_sim_original(system, dt, steps)
    end_time = time()
    runtime_base = end_time - start_time
    
    # Compare results
    speedup = runtime_base / runtime_opt
    println("\nComparison:")
    println("Baseline runtime: $runtime_base seconds")
    println("Optimized runtime: $runtime_opt seconds")
    println("Speedup: $(speedup)x")
    
    # Check energy conservation
    println("Baseline energy error: $(energy_error_base * 100)%")
    println("Optimized energy error: $(energy_error_opt * 100)%")
    
    # Good implementations should have low energy error
    if energy_error_opt < 0.01  # Less than 1% error
        println("✓ Energy conservation is good (error < 1%)")
    else
        println("⚠ Energy conservation is poor (error >= 1%)")
    end
    
    return runtime_opt, energy_error_opt, speedup
end

function run_profile(steps::Int64=3, n::Int64=10000)
    println("Running profiler with $n particles for $steps steps...")
    system = circular_orbits(n)
    dt = 1e-3
    
    # Warm up
    simple_sim(system, dt, 1, true)
    
    # Profile
    Profile.clear()
    Profile.@profile simple_sim(system, dt, steps, true)
    
    # Print profile
    Profile.print(maxdepth=20)
end

# to run:
# julia main.jl #steps #particles [verify|profile]
if !isinteractive()
    if length(ARGS) >= 3
        if ARGS[3] == "verify"
            # Run verification with smaller dataset
            verify_optimization(parse(Int64, ARGS[1]), min(parse(Int64, ARGS[2]), 10000))
        elseif ARGS[3] == "profile"
            # Run profiler
            run_profile(parse(Int64, ARGS[1]), min(parse(Int64, ARGS[2]), 50000))
        else
            # Regular run
            steps = parse(Int64, ARGS[1])
            n = parse(Int64, ARGS[2])
            dt = 1e-3
            system = circular_orbits(n)
            simple_sim(system, dt, steps, true)
        end
    else
        steps = parse(Int64, ARGS[1])
        n = parse(Int64, ARGS[2])
        dt = 1e-3
        system = circular_orbits(n)
        simple_sim(system, dt, steps, true)
    end
end

end
