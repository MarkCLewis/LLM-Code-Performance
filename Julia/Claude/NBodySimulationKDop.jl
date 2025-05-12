module NBodySimulation

using LinearAlgebra
using Random
using Printf
using Base.Threads
using StaticArrays
using SIMD  # Add SIMD for explicit vectorization
using Profile # For profiling

# Constants
const G = 6.67430e-11  # Gravitational constant
const EPSILON = 1e-10   # Softening parameter to avoid division by zero

"""
    Body

Represents a body in 3D space with position, velocity, and mass.
Using StaticArrays for better performance.
"""
struct Body
    position::SVector{3, Float64}  # 3D position vector
    velocity::SVector{3, Float64}  # 3D velocity vector
    mass::Float64                  # Mass of the body
end

# Constructor for Body with regular arrays
Body(position::Vector{Float64}, velocity::Vector{Float64}, mass::Float64) = 
    Body(SVector{3, Float64}(position), SVector{3, Float64}(velocity), mass)

"""
    Cell

Represents a cell in the Barnes-Hut octree.
"""
mutable struct Cell
    # Geometric properties
    center::SVector{3, Float64}       # Center of the cell
    size::Float64                     # Side length of the cell
    
    # Physical properties
    total_mass::Float64               # Total mass of all bodies in this cell
    com::SVector{3, Float64}          # Center of mass of all bodies in this cell
    
    # Tree structure
    children::Vector{Union{Cell, Body, Nothing}}  # 8 children cells (octants) or body or nothing
    num_bodies::Int                   # Number of bodies in this cell and all subcells
    
    # Constructors
    function Cell(center::SVector{3, Float64}, size::Float64)
        new(
            center, 
            size, 
            0.0, 
            SVector{3, Float64}(0.0, 0.0, 0.0), 
            [nothing for _ in 1:8], 
            0
        )
    end
    
    Cell(center::Vector{Float64}, size::Float64) = Cell(SVector{3, Float64}(center), size)
end

"""
    System

Contains all bodies in the simulation and the octree for force calculations.
"""
struct System
    bodies::Vector{Body}
    octree::Ref{Union{Cell, Nothing}}  # Use Ref for mutable reference
    
    # Pre-allocated memory for parallelization
    accelerations::Vector{SVector{3, Float64}}
    energy_chunks::Vector{Float64}
    
    # Constructor
    function System(bodies::Vector{Body})
        n = length(bodies)
        accelerations = Vector{SVector{3, Float64}}(undef, n)
        energy_chunks = Vector{Float64}(undef, nthreads())
        return new(bodies, Ref{Union{Cell, Nothing}}(nothing), accelerations, energy_chunks)
    end
end

"""
    update_octree!(system)

Update the octree for a system.
"""
function update_octree!(system::System)
    system.octree[] = build_octree(system.bodies)
end

# Bit operations for octants are faster than conditionals
@inline function determine_octant(position::SVector{3, Float64}, center::SVector{3, Float64})
    octant = 0
    octant |= position[1] >= center[1] ? 1 : 0
    octant |= position[2] >= center[2] ? 2 : 0
    octant |= position[3] >= center[3] ? 4 : 0
    return octant + 1  # 1-based indexing
end

"""
    get_octant_center(cell_center, cell_size, octant)

Get the center of a specific octant of a cell.
"""
@inline function get_octant_center(cell_center::SVector{3, Float64}, cell_size::Float64, octant::Int)
    # Pre-calculated offset values
    quarter_size = cell_size * 0.25
    offsets = (
        SVector{3, Float64}(-quarter_size, -quarter_size, -quarter_size), # octant 1
        SVector{3, Float64}(quarter_size, -quarter_size, -quarter_size),  # octant 2
        SVector{3, Float64}(-quarter_size, quarter_size, -quarter_size),  # octant 3
        SVector{3, Float64}(quarter_size, quarter_size, -quarter_size),   # octant 4
        SVector{3, Float64}(-quarter_size, -quarter_size, quarter_size),  # octant 5
        SVector{3, Float64}(quarter_size, -quarter_size, quarter_size),   # octant 6
        SVector{3, Float64}(-quarter_size, quarter_size, quarter_size),   # octant 7
        SVector{3, Float64}(quarter_size, quarter_size, quarter_size)     # octant 8
    )
    return cell_center + offsets[octant]
end

"""
    insert_body!(cell, body)

Insert a body into the octree.
"""
function insert_body!(cell::Cell, body::Body)
    # If this cell is empty, just put the body here
    if cell.num_bodies == 0
        # Find the octant
        octant = determine_octant(body.position, cell.center)
        cell.children[octant] = body
        cell.total_mass = body.mass
        cell.com = body.position
        cell.num_bodies = 1
        return
    end
    
    # If this cell has just one body
    if cell.num_bodies == 1
        # Find the single body
        for i in 1:8
            if isa(cell.children[i], Body)
                existing_body = cell.children[i]
                # Remove the body from the cell
                cell.children[i] = nothing
                
                # Create subcells and insert both bodies
                new_size = cell.size / 2
                
                # Insert the existing body
                octant = determine_octant(existing_body.position, cell.center)
                if cell.children[octant] === nothing || !isa(cell.children[octant], Cell)
                    cell.children[octant] = Cell(get_octant_center(cell.center, cell.size, octant), new_size)
                end
                insert_body!(cell.children[octant]::Cell, existing_body)
                
                # Insert the new body
                octant = determine_octant(body.position, cell.center)
                if cell.children[octant] === nothing || !isa(cell.children[octant], Cell)
                    cell.children[octant] = Cell(get_octant_center(cell.center, cell.size, octant), new_size)
                end
                insert_body!(cell.children[octant]::Cell, body)
                
                break
            end
        end
    else
        # If this cell already has multiple bodies, just insert the new body in the appropriate octant
        octant = determine_octant(body.position, cell.center)
        
        if cell.children[octant] === nothing
            # Create a new subcell
            cell.children[octant] = Cell(get_octant_center(cell.center, cell.size, octant), cell.size / 2)
        elseif isa(cell.children[octant], Body)
            # Convert this child to a cell
            existing_body = cell.children[octant]
            cell.children[octant] = Cell(get_octant_center(cell.center, cell.size, octant), cell.size / 2)
            insert_body!(cell.children[octant]::Cell, existing_body)
        end
        
        # Insert the body into the subcell
        if isa(cell.children[octant], Cell)
            insert_body!(cell.children[octant]::Cell, body)
        end
    end
    
    # Update the properties of this cell
    cell.num_bodies += 1
    
    # Recalculate center of mass
    cell.com = (cell.com * cell.total_mass + body.position * body.mass) / (cell.total_mass + body.mass)
    cell.total_mass += body.mass
end

"""
    build_octree(bodies)

Build an octree from a list of bodies.
"""
function build_octree(bodies::Vector{Body})
    # Find the bounding box using SIMD for faster min/max operations
    min_coords = SVector{3, Float64}(Inf, Inf, Inf)
    max_coords = SVector{3, Float64}(-Inf, -Inf, -Inf)
    
    @inbounds for body in bodies
        min_coords = min.(min_coords, body.position)
        max_coords = max.(max_coords, body.position)
    end
    
    # Add a small padding to ensure all bodies are inside
    min_coords = min_coords .- EPSILON
    max_coords = max_coords .+ EPSILON
    
    # Calculate the center and size of the root cell
    center = (min_coords + max_coords) / 2
    size = maximum(max_coords - min_coords) * 1.01  # 1% larger to be safe
    
    # Create the root cell
    root = Cell(center, size)
    
    # Insert all bodies
    @inbounds for body in bodies
        insert_body!(root, body)
    end
    
    return root
end

"""
    compute_acceleration_bh(octree, position, mass, theta, ε)

Compute the acceleration of a body at a given position due to gravitational forces
using the Barnes-Hut approximation with parameter theta.
"""
@inline function compute_acceleration_bh(
    cell::Union{Cell, Body, Nothing}, 
    position::SVector{3, Float64}, 
    mass::Float64, 
    theta::Float64
)
    # Base case: empty cell
    if cell === nothing
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    
    # If this is a body, calculate direct force
    if isa(cell, Body)
        # Skip self-interaction
        r_vec = cell.position - position
        r_squared = sum(r_vec.^2)
        
        if r_squared < EPSILON^2
            return SVector{3, Float64}(0.0, 0.0, 0.0)
        end
        
        # Softened distance squared
        r_squared += EPSILON^2  
        r = sqrt(r_squared)
        
        # Gravitational acceleration: a = G * m / r^2 * (r_vec / r)
        return G * cell.mass / r_squared * (r_vec / r)
    end
    
    # At this point, cell must be a Cell
    # Calculate the distance from the body to the cell's center of mass
    r_vec = cell.com - position
    r_squared = sum(r_vec.^2)
    r = sqrt(r_squared)
    
    # If the cell is far enough away, use the approximation
    if cell.size / r < theta || cell.num_bodies == 1
        # Softened distance squared
        r_squared += EPSILON^2
        r = sqrt(r_squared)
        
        # Gravitational acceleration: a = G * M / r^2 * (r_vec / r)
        return G * cell.total_mass / r_squared * (r_vec / r)
    else
        # Otherwise, recursively compute forces from children
        acceleration = SVector{3, Float64}(0.0, 0.0, 0.0)
        
        @inbounds for child in cell.children
            acceleration += compute_acceleration_bh(child, position, mass, theta)
        end
        
        return acceleration
    end
end

"""
    compute_acceleration_bh(system, body_index, theta)

Compute the acceleration of a body due to gravitational forces from all other bodies
using the Barnes-Hut approximation with parameter theta.
"""
@inline function compute_acceleration_bh(system::System, body_index::Int, theta::Float64)
    body = system.bodies[body_index]
    
    # Compute acceleration using Barnes-Hut approximation
    return compute_acceleration_bh(system.octree[], body.position, body.mass, theta)
end

"""
    kick_step!(system, dt, theta)

Update velocities of all bodies using computed accelerations (kick step).
Uses multithreading and Barnes-Hut algorithm with parameter theta.
"""
function kick_step!(system::System, dt::Float64, theta::Float64)
    n = length(system.bodies)
    
    # Update the octree only once per kick step
    update_octree!(system)
    
    # Compute accelerations in parallel using Barnes-Hut
    @threads for i in 1:n
        system.accelerations[i] = compute_acceleration_bh(system, i, theta)
    end
    
    # Update velocities - this is a good candidate for SIMD operations
    @inbounds for i in 1:n
        system.bodies[i] = Body(
            system.bodies[i].position,
            system.bodies[i].velocity + system.accelerations[i] * dt,
            system.bodies[i].mass
        )
    end
end

"""
    drift_step!(system, dt)

Update positions of all bodies using their velocities (drift step).
Uses multithreading for updating positions.
"""
function drift_step!(system::System, dt::Float64)
    n = length(system.bodies)
    
    # Using in-place updates
    @threads for i in 1:n
        @inbounds body = system.bodies[i]
        @inbounds system.bodies[i] = Body(
            body.position + body.velocity * dt,
            body.velocity,
            body.mass
        )
    end
end

"""
    calculate_energy(system)

Calculate the total energy (kinetic + potential) of the system.
Uses multithreading and avoids unnecessary calculations.
"""
# function calculate_energy(system::System)
#     n = length(system.bodies)
    
#     # Reset energy chunks
#     fill!(system.energy_chunks, 0.0)
    
#     # Calculate kinetic and potential energy in one pass for better cache efficiency
#     @threads for i in 1:n
#         thread_id = threadid()
        
#         # Calculate kinetic energy for this body
#         @inbounds body_i = system.bodies[i]
#         v_squared = sum(body_i.velocity.^2)
#         system.energy_chunks[thread_id] += 0.5 * body_i.mass * v_squared
        
#         # Calculate potential energy with all remaining bodies
#         @inbounds for j in (i+1):n
#             body_j = system.bodies[j]
            
#             r_vec = body_j.position - body_i.position
#             r = sqrt(sum(r_vec.^2) + EPSILON^2)  # Softened distance
            
#             system.energy_chunks[thread_id] -= G * body_i.mass * body_j.mass / r
#         end
#     end
    
#     # Sum up all energy chunks
#     return sum(system.energy_chunks)
# end

"""
    run_simulation!(system, dt, num_steps, theta)

Run the N-body simulation for a specified number of steps using
Barnes-Hut approximation with parameter theta.
Reports thread usage and performance metrics.
"""
function run_simulation!(system::System, dt::Float64, num_steps::Int, theta::Float64)
    println("Running simulation with $(nthreads()) threads")
    println("Using Barnes-Hut approximation with theta = $theta")
    
    # initial_energy = calculate_energy(system)
    # println("Initial total energy: ", initial_energy)
    
    # Pre-compute half time step and store frequently used values
    half_dt = dt/2
    report_interval = max(1, num_steps ÷ 10)
    
    # Start timing
    start_time = time()
    
    # Run the simulation for num_steps
    for step in 1:num_steps
        if step % report_interval == 0
            elapsed_so_far = time() - start_time
            steps_per_second = step / elapsed_so_far
            estimated_total = num_steps / steps_per_second
            
            println("Completed step $step of $num_steps ($(round(100.0 * step / num_steps, digits=1))%)")
            println("Current performance: $(round(steps_per_second, digits=2)) steps/second")
            println("Estimated total runtime: $(round(estimated_total / 60, digits=1)) minutes")
            
            # Optional: force garbage collection during status reports to reduce memory pressure
            GC.gc()
        end
        
        # Kick-drift-kick leapfrog integrator
        kick_step!(system, half_dt, theta)  # Half kick
        drift_step!(system, dt)             # Full drift
        kick_step!(system, half_dt, theta)  # Half kick
    end
    
    # End timing
    end_time = time()
    elapsed = end_time - start_time
    
    # final_energy = calculate_energy(system)
    # println("Final total energy: ", final_energy)
    # println("Energy change: ", final_energy - initial_energy)
    # println("Relative energy change: ", (final_energy - initial_energy) / initial_energy)
    println("bodie[1] %e %e %e", system.bodies[1].position[1], system.bodies[1].position[2], system.bodies[1].position[3])
    println("Simulation completed in $(elapsed) seconds")
    println("Average performance: $(round(num_steps / elapsed, digits=2)) steps/second")
    
    return system
end

"""
    initialize_system(central_mass, central_position, central_velocity, 
                     num_orbiting, orbiting_mass, 
                     min_radius, max_radius)

Initialize a system with a central body and `num_orbiting` smaller bodies
in circular orbits at randomly distributed radii. Uses batch processing for
better performance with large numbers of bodies.
"""
function initialize_system(central_mass, central_position, central_velocity, 
                         num_orbiting, orbiting_mass, 
                         min_radius, max_radius)
    # Convert inputs to SVectors for performance
    central_pos_sv = SVector{3, Float64}(central_position)
    central_vel_sv = SVector{3, Float64}(central_velocity)
    
    # Create central body
    central_body = Body(central_pos_sv, central_vel_sv, central_mass)
    
    # Initialize array of bodies with pre-allocation
    bodies = Vector{Body}(undef, num_orbiting + 1)
    bodies[1] = central_body
    
    # Create orbiting bodies with random positions on circular orbits
    rng = MersenneTwister(1234)  # For reproducibility
    
    # Batch size for orbit generation (improves memory access patterns)
    batch_size = min(10000, num_orbiting)
    
    # Pre-allocate arrays for batch processing
    radii = Vector{Float64}(undef, batch_size)
    thetas = Vector{Float64}(undef, batch_size)
    phis = Vector{Float64}(undef, batch_size)
    v_orbits = Vector{Float64}(undef, batch_size)
    
    # Process in batches
    for batch_start in 1:batch_size:num_orbiting
        batch_end = min(batch_start + batch_size - 1, num_orbiting)
        batch_count = batch_end - batch_start + 1
        
        # Generate batch of random values
        for i in 1:batch_count
            radii[i] = min_radius + (max_radius - min_radius) * rand(rng)
            thetas[i] = 2π * rand(rng)  # Azimuthal angle
            phis[i] = acos(2 * rand(rng) - 1)  # Polar angle
            v_orbits[i] = sqrt(G * central_mass / radii[i])
        end
        
        # Process batch
        for i in 1:batch_count
            idx = batch_start + i - 1
            r = radii[i]
            θ = thetas[i]
            φ = phis[i]
            
            # Convert to Cartesian coordinates
            x = r * sin(φ) * cos(θ)
            y = r * sin(φ) * sin(θ)
            z = r * cos(φ)
            
            position = central_pos_sv + SVector{3, Float64}(x, y, z)
            
            # Velocity vector perpendicular to position vector (for circular orbit)
            # We need a vector perpendicular to the radial vector
            rel_pos = SVector{3, Float64}(x, y, z)
            
            if abs(z) < 0.9 * r  # Avoid numerical issues near poles
                perp_vector = cross(SVector{3, Float64}(0.0, 0.0, 1.0), rel_pos)
            else
                perp_vector = cross(SVector{3, Float64}(1.0, 0.0, 0.0), rel_pos)
            end
            
            # Normalize and scale to orbital velocity
            velocity_direction = perp_vector / norm(perp_vector)
            velocity = central_vel_sv + velocity_direction * v_orbits[i]
            
            # Create and add the orbiting body
            bodies[idx + 1] = Body(position, velocity, orbiting_mass)
        end
    end
    
    return System(bodies)
end

"""
    main(num_orbiting)

Create and run an N-body simulation with a central body and a specified number of orbiting bodies.
Uses Barnes-Hut algorithm with theta=0.3.
Shows thread information and performance metrics.
"""
function main(num_orbiting=100000)
    # Barnes-Hut theta parameter
    theta = 0.3
    
    # Display thread information
    println("Julia is running with $(nthreads()) threads")
    println("To increase threads, set the JULIA_NUM_THREADS environment variable")
    println("or start Julia with the --threads=auto or --threads=N flag")
    println("Using Barnes-Hut approximation with theta = $theta")
    
    # Parameters
    central_mass = 1.0e30      # Mass of central body (approximately solar mass)
    central_position = [0.0, 0.0, 0.0]
    central_velocity = [0.0, 0.0, 0.0]
    orbiting_mass = 1.0e20     # Mass of orbiting bodies
    min_radius = 1.0e10        # Minimum orbital radius
    max_radius = 5.0e10        # Maximum orbital radius
    dt = 3600.0                # Time step (seconds)
    num_steps = 10           # Number of simulation steps
    
    println("Initializing system with $num_orbiting orbiting bodies...")
    init_start = time()
    
    # Initialize the system
    system = initialize_system(
        central_mass, central_position, central_velocity,
        num_orbiting, orbiting_mass, min_radius, max_radius
    )
    
    init_time = time() - init_start
    println("Initialization completed in $(round(init_time, digits=2)) seconds")
    println("Starting simulation with $(length(system.bodies)) bodies for $num_steps steps...")
    
    # Run the simulation
    run_simulation!(system, dt, num_steps, theta)
    
    println("Simulation complete!")
    return system
end

"""
    run_with_profile(num_orbiting=10000, num_steps=10)

Run a smaller simulation with profiling enabled to identify bottlenecks.
"""
function run_with_profile(num_orbiting=100000, num_steps=10)
    # Barnes-Hut theta parameter
    theta = 0.3
    
    # Parameters
    central_mass = 1.0e30      # Mass of central body (approximately solar mass)
    central_position = [0.0, 0.0, 0.0]
    central_velocity = [0.0, 0.0, 0.0]
    orbiting_mass = 1.0e20     # Mass of orbiting bodies
    min_radius = 1.0e10        # Minimum orbital radius
    max_radius = 5.0e10        # Maximum orbital radius
    dt = 3600.0                # Time step (seconds)
    
    # Initialize the system
    system = initialize_system(
        central_mass, central_position, central_velocity,
        num_orbiting, orbiting_mass, min_radius, max_radius
    )
    
    # Run profile
    Profile.clear()
    Profile.init(n = 10^7, delay = 0.01)
    @profile run_simulation!(system, dt, num_steps, theta)
    Profile.print(format=:flat, sortedby=:count)
    
    return system
end

"""
    reduced_main(num_orbiting=1000)

Run the simulation with a smaller number of particles for testing.
"""
function reduced_main(num_orbiting=10)
    main(num_orbiting)
end

end  # module NBodySimulation

# Execute the simulation with reduced particle count for testing
# NBodySimulation.reduced_main()

# Uncomment the next line to run the full simulation with 1 million bodies
NBodySimulation.main()

# To profile and identify bottlenecks, uncomment this line:
# NBodySimulation.run_with_profile()