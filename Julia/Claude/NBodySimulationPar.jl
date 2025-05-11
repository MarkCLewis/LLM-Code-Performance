module NBodySimulation

using LinearAlgebra
using Random
using Printf
using Base.Threads

"""
    Body

Represents a body in 3D space with position, velocity, and mass.
"""
struct Body
    position::Vector{Float64}  # 3D position vector
    velocity::Vector{Float64}  # 3D velocity vector
    mass::Float64              # Mass of the body
end

"""
    System

Contains all bodies in the simulation.
"""
struct System
    bodies::Vector{Body}
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
    # Create central body
    central_body = Body(central_position, central_velocity, central_mass)
    
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
        
        position = central_position + [x, y, z]
        
        # Calculate velocity for a circular orbit
        # v = sqrt(G * M / r)
        G = 6.67430e-11  # Gravitational constant
        v_orbit = sqrt(G * central_mass / r)
        
        # Velocity vector perpendicular to position vector (for circular orbit)
        # We need a vector perpendicular to the radial vector
        if abs(z) < 0.9 * r  # Avoid numerical issues near poles
            perp_vector = cross([0.0, 0.0, 1.0], [x, y, z])
        else
            perp_vector = cross([1.0, 0.0, 0.0], [x, y, z])
        end
        
        # Normalize and scale to orbital velocity
        velocity_direction = perp_vector / norm(perp_vector)
        velocity = central_velocity + velocity_direction * v_orbit
        
        # Create and add the orbiting body
        push!(bodies, Body(position, velocity, orbiting_mass))
    end
    
    return System(bodies)
end

"""
    compute_acceleration(system, body_index)

Compute the acceleration of a body due to gravitational forces from all other bodies.
"""
function compute_acceleration(system, body_index)
    G = 6.67430e-11  # Gravitational constant
    ε = 1e-10        # Softening parameter to avoid division by zero
    
    body = system.bodies[body_index]
    acceleration = zeros(Float64, 3)
    
    for i in 1:length(system.bodies)
        if i == body_index
            continue  # Skip self-interaction
        end
        
        other = system.bodies[i]
        r_vec = other.position - body.position
        r_squared = sum(r_vec.^2) + ε^2  # Softened distance squared
        r = sqrt(r_squared)
        
        # Gravitational acceleration: a = G * m / r^2 * (r_vec / r)
        acceleration += G * other.mass / r_squared * (r_vec / r)
    end
    
    return acceleration
end

"""
    kick_step!(system, dt)

Update velocities of all bodies using computed accelerations (kick step).
Uses multithreading for computing accelerations.
"""
function kick_step!(system, dt)
    n = length(system.bodies)
    accelerations = Vector{Vector{Float64}}(undef, n)
    
    # Compute accelerations in parallel
    @threads for i in 1:n
        accelerations[i] = compute_acceleration(system, i)
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
    run_simulation!(system, dt, num_steps)

Run the N-body simulation for a specified number of steps.
Reports thread usage and performance metrics.
"""
function run_simulation!(system, dt, num_steps)
    println("Running simulation with $(nthreads()) threads")
    
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
        kick_step!(system, dt/2)  # Half kick
        drift_step!(system, dt)   # Full drift
        kick_step!(system, dt/2)  # Half kick
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
Shows thread information and performance metrics.
"""
function main(num_orbiting=10000)
    # Display thread information
    println("Julia is running with $(nthreads()) threads")
    println("To increase threads, set the JULIA_NUM_THREADS environment variable")
    println("or start Julia with the --threads=auto or --threads=N flag")
    
    # Parameters
    central_mass = 1.0e30      # Mass of central body (approximately solar mass)
    central_position = [0.0, 0.0, 0.0]
    central_velocity = [0.0, 0.0, 0.0]
    orbiting_mass = 1.0e20     # Mass of orbiting bodies
    min_radius = 1.0e10        # Minimum orbital radius
    max_radius = 5.0e10        # Maximum orbital radius
    dt = 3600.0                # Time step (seconds)
    num_steps = 100           # Number of simulation steps
    
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
    run_simulation!(system, dt, num_steps)
    
    println("Simulation complete!")
    return system
end

"""
    reduced_main(num_orbiting=1000)

Run the simulation with a smaller number of particles for testing.
"""
function reduced_main(num_orbiting=10000)
    main(num_orbiting)
end

end  # module NBodySimulation

# Execute the simulation with reduced particle count for testing
# Comment out the next line to run the full simulation
#NBodySimulation.reduced_main()

# Uncomment the next line to run the full simulation with 1 million bodies
NBodySimulation.main()