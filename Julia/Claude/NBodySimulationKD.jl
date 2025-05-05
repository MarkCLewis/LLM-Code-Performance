"""
    System

Contains all bodies in the simulation and the octree for force calculations.
"""
struct System
    bodies::Vector{Body}
    octree::Union{Cell, Nothing}
    
    # Constructor without octree (will be built later)
    System(bodies::Vector{Body}) = new(bodies, nothing)
end

"""
    update_octree!(system)

Update the octree for a system.
"""
function update_octree!(system::System)
    system.octree = build_octree(system.bodies)
endmodule NBodySimulation

using LinearAlgebra
using Random
using Printf
using Base.Threads
using StaticArrays

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
    
    function Cell(center::Vector{Float64}, size::Float64)
        new(
            SVector{3, Float64}(center), 
            size, 
            0.0, 
            SVector{3, Float64}(0.0, 0.0, 0.0), 
            [nothing for _ in 1:8], 
            0
        )
    end
end

"""
    determine_octant(position, center)

Determine which octant (0-7) a position falls into relative to the center.
Octants are numbered as follows:
    0: (-, -, -)   1: (+, -, -)
    2: (-, +, -)   3: (+, +, -)
    4: (-, -, +)   5: (+, -, +)
    6: (-, +, +)   7: (+, +, +)
"""
function determine_octant(position::SVector{3, Float64}, center::SVector{3, Float64})
    # Determine octant based on position relative to center
    octant = 0
    if position[1] >= center[1]
        octant |= 1
    end
    if position[2] >= center[2]
        octant |= 2
    end
    if position[3] >= center[3]
        octant |= 4
    end
    return octant + 1  # 1-based indexing
end

"""
    get_octant_center(cell_center, cell_size, octant)

Get the center of a specific octant of a cell.
"""
function get_octant_center(cell_center::SVector{3, Float64}, cell_size::Float64, octant::Int)
    # Calculate the offset for each dimension
    offset_x = ((octant - 1) & 1) == 0 ? -0.25 : 0.25
    offset_y = ((octant - 1) & 2) == 0 ? -0.25 : 0.25
    offset_z = ((octant - 1) & 4) == 0 ? -0.25 : 0.25
    
    # Calculate the center of the octant
    return cell_center + cell_size * SVector{3, Float64}(offset_x, offset_y, offset_z)
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
    # Find the bounding box
    min_coords = SVector{3, Float64}(Inf, Inf, Inf)
    max_coords = SVector{3, Float64}(-Inf, -Inf, -Inf)
    
    for body in bodies
        min_coords = min.(min_coords, body.position)
        max_coords = max.(max_coords, body.position)
    end
    
    # Add a small padding to ensure all bodies are inside
    min_coords = min_coords .- 1e-10
    max_coords = max_coords .+ 1e-10
    
    # Calculate the center and size of the root cell
    center = (min_coords + max_coords) / 2
    size = maximum(max_coords - min_coords) * 1.01  # 1% larger to be safe
    
    # Create the root cell
    root = Cell(center, size)
    
    # Insert all bodies
    for body in bodies
        insert_body!(root, body)
    end
    
    return root
end

"""
    initialize_system(central_mass, central_position, central_velocity, 
                     num_orbiting, orbiting_mass, 
                     min_radius, max_radius)

Initialize a system with a central body and `num_orbiting` smaller bodies
in circular orbits at randomly distributed radii.
"""
function initialize_system(central_mass, central_position, central_velocity, 
                         num_orbiting, orbiting_mass, 
                         min_radius, max_radius)
    # Convert inputs to SVectors for performance
    central_pos_sv = SVector{3, Float64}(central_position)
    central_vel_sv = SVector{3, Float64}(central_velocity)
    
    # Create central body
    central_body = Body(central_pos_sv, central_vel_sv, central_mass)
    
    # Initialize system with central body
    bodies = [central_body]
    
    # Create orbiting bodies with random positions on circular orbits
    rng = MersenneTwister(1234)  # For reproducibility
    
    for i in 1:num_orbiting
        # Random radius between min_radius and max_radius
        r = min_radius + (max_radius - min_radius) * rand(rng)
        
        # Random point on a sphere
        θ = 2π * rand(rng)  # Azimuthal angle
        φ = acos(2 * rand(rng) - 1)  # Polar angle
        
        # Convert to Cartesian coordinates
        x = r * sin(φ) * cos(θ)
        y = r * sin(φ) * sin(θ)
        z = r * cos(φ)
        
        position = central_pos_sv + SVector{3, Float64}(x, y, z)
        
        # Calculate velocity for a circular orbit
        # v = sqrt(G * M / r)
        G = 6.67430e-11  # Gravitational constant
        v_orbit = sqrt(G * central_mass / r)
        
        # Velocity vector perpendicular to position vector (for circular orbit)
        # We need a vector perpendicular to the radial vector
        if abs(z) < 0.9 * r  # Avoid numerical issues near poles
            perp_vector = cross(SVector{3, Float64}(0.0, 0.0, 1.0), SVector{3, Float64}(x, y, z))
        else
            perp_vector = cross(SVector{3, Float64}(1.0, 0.0, 0.0), SVector{3, Float64}(x, y, z))
        end
        
        # Normalize and scale to orbital velocity
        velocity_direction = perp_vector / norm(perp_vector)
        velocity = central_vel_sv + velocity_direction * v_orbit
        
        # Create and add the orbiting body
        push!(bodies, Body(position, velocity, orbiting_mass))
    end
    
    return System(bodies)
end

"""
    compute_acceleration_bh(octree, position, mass, theta, G, ε)

Compute the acceleration of a body at a given position due to gravitational forces
using the Barnes-Hut approximation with parameter theta.
"""
function compute_acceleration_bh(cell::Union{Cell, Body, Nothing}, position::SVector{3, Float64}, 
                               mass::Float64, theta::Float64, G::Float64, ε::Float64)
    # Base case: empty cell
    if cell === nothing
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    
    # If this is a body, calculate direct force
    if isa(cell, Body)
        # Skip self-interaction
        if norm(position - cell.position) < 1e-10
            return SVector{3, Float64}(0.0, 0.0, 0.0)
        end
        
        # Calculate gravitational force
        r_vec = cell.position - position
        r_squared = sum(r_vec.^2) + ε^2  # Softened distance squared
        r = sqrt(r_squared)
        
        # Gravitational acceleration: a = G * m / r^2 * (r_vec / r)
        return G * cell.mass / r_squared * (r_vec / r)
    end
    
    # At this point, cell must be a Cell
    # Calculate the distance from the body to the cell's center of mass
    r_vec = cell.com - position
    r = norm(r_vec)
    
    # If the cell is far enough away, use the approximation
    if cell.size / r < theta || cell.num_bodies == 1
        # Softened distance squared
        r_squared = r^2 + ε^2
        
        # Gravitational acceleration: a = G * M / r^2 * (r_vec / r)
        return G * cell.total_mass / r_squared * (r_vec / r)
    else
        # Otherwise, recursively compute forces from children
        acceleration = SVector{3, Float64}(0.0, 0.0, 0.0)
        
        for child in cell.children
            acceleration += compute_acceleration_bh(child, position, mass, theta, G, ε)
        end
        
        return acceleration
    end
end

"""
    compute_acceleration_bh(system, body_index, theta)

Compute the acceleration of a body due to gravitational forces from all other bodies
using the Barnes-Hut approximation with parameter theta.
"""
function compute_acceleration_bh(system::System, body_index::Int, theta::Float64)
    G = 6.67430e-11  # Gravitational constant
    ε = 1e-10        # Softening parameter to avoid division by zero
    
    body = system.bodies[body_index]
    
    # Make sure the octree is up to date
    if system.octree === nothing
        update_octree!(system)
    end
    
    # Compute acceleration using Barnes-Hut approximation
    return compute_acceleration_bh(system.octree, body.position, body.mass, theta, G, ε)
end

"""
    kick_step!(system, dt, theta)

Update velocities of all bodies using computed accelerations (kick step).
Uses multithreading and Barnes-Hut algorithm with parameter theta.
"""
function kick_step!(system::System, dt::Float64, theta::Float64)
    n = length(system.bodies)
    accelerations = Vector{SVector{3, Float64}}(undef, n)
    
    # Update the octree
    update_octree!(system)
    
    # Compute accelerations in parallel using Barnes-Hut
    @threads for i in 1:n
        accelerations[i] = compute_acceleration_bh(system, i, theta)
    end
    
    # Update velocities
    for i in 1:n
        system.bodies[i] = Body(
            system.bodies[i].position,
            system.bodies[i].velocity + accelerations[i] * dt,
            system.bodies[i].mass
        )
    end
end

"""
    drift_step!(system, dt)

Update positions of all bodies using their velocities (drift step).
Uses multithreading for updating positions.
"""
function drift_step!(system, dt)
    n = length(system.bodies)
    new_bodies = Vector{Body}(undef, n)
    
    @threads for i in 1:n
        body = system.bodies[i]
        new_bodies[i] = Body(
            body.position + body.velocity * dt,
            body.velocity,
            body.mass
        )
    end
    
    # Copy the new bodies back to the system
    for i in 1:n
        system.bodies[i] = new_bodies[i]
    end
end

"""
    calculate_energy(system)

Calculate the total energy (kinetic + potential) of the system.
Uses multithreading for kinetic energy calculation and atomic operations
for thread-safe summation of potential energy.
"""
function calculate_energy(system)
    G = 6.67430e-11  # Gravitational constant
    ε = 1e-10        # Softening parameter
    
    n = length(system.bodies)
    
    # Calculate kinetic energy in parallel: KE = 0.5 * m * v^2
    kinetic_energies = zeros(Float64, nthreads())
    
    @threads for i in 1:n
        thread_id = threadid()
        body = system.bodies[i]
        v_squared = sum(body.velocity.^2)
        kinetic_energies[thread_id] += 0.5 * body.mass * v_squared
    end
    
    kinetic_energy = sum(kinetic_energies)
    
    # For potential energy, we need to be careful with parallelization
    # because of the double loop. We'll use atomic operations for thread safety.
    potential_energy_chunks = zeros(Float64, nthreads())
    
    # Each thread handles a chunk of the outer loop
    @threads for i in 1:(n-1)
        thread_id = threadid()
        local_potential = 0.0
        
        for j in (i+1):n
            body_i = system.bodies[i]
            body_j = system.bodies[j]
            
            r_vec = body_j.position - body_i.position
            r = sqrt(sum(r_vec.^2) + ε^2)  # Softened distance
            
            local_potential -= G * body_i.mass * body_j.mass / r
        end
        
        potential_energy_chunks[thread_id] += local_potential
    end
    
    potential_energy = sum(potential_energy_chunks)
    
    return kinetic_energy + potential_energy
end

"""
    run_simulation!(system, dt, num_steps, theta)

Run the N-body simulation for a specified number of steps using
Barnes-Hut approximation with parameter theta.
Reports thread usage and performance metrics.
"""
function run_simulation!(system::System, dt::Float64, num_steps::Int, theta::Float64)
    println("Running simulation with $(nthreads()) threads")
    println("Using Barnes-Hut approximation with theta = $theta")
    
    initial_energy = calculate_energy(system)
    println("Initial total energy: ", initial_energy)
    
    # Start timing
    start_time = time()
    
    # Run the simulation for num_steps
    for step in 1:num_steps
        if step % (num_steps ÷ 10) == 0
            elapsed_so_far = time() - start_time
            steps_per_second = step / elapsed_so_far
            estimated_total = num_steps / steps_per_second
            
            println("Completed step $step of $num_steps ($(round(100.0 * step / num_steps, digits=1))%)")
            println("Current performance: $(round(steps_per_second, digits=2)) steps/second")
            println("Estimated total runtime: $(round(estimated_total / 60, digits=1)) minutes")
        end
        
        # Kick-drift-kick leapfrog integrator
        kick_step!(system, dt/2, theta)  # Half kick
        drift_step!(system, dt)         # Full drift
        kick_step!(system, dt/2, theta)  # Half kick
    end
    
    # End timing
    end_time = time()
    elapsed = end_time - start_time
    
    final_energy = calculate_energy(system)
    println("Final total energy: ", final_energy)
    println("Energy change: ", final_energy - initial_energy)
    println("Relative energy change: ", (final_energy - initial_energy) / initial_energy)
    println("Simulation completed in $(elapsed) seconds")
    println("Average performance: $(round(num_steps / elapsed, digits=2)) steps/second")
    
    return system
end

"""
    main(num_orbiting)

Create and run an N-body simulation with a central body and a specified number of orbiting bodies.
Uses Barnes-Hut algorithm with theta=0.3.
Shows thread information and performance metrics.
"""
function main(num_orbiting=10^6)
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
    num_steps = 1000           # Number of simulation steps
    
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
    reduced_main(num_orbiting=1000)

Run the simulation with a smaller number of particles for testing.
"""
function reduced_main(num_orbiting=1000)
    main(num_orbiting)
end

end  # module NBodySimulation

# Execute the simulation with reduced particle count for testing
# Comment out the next line to run the full simulation
#NBodySimulation.reduced_main()

# Uncomment the next line to run the full simulation with 1 million bodies
NBodySimulation.main()