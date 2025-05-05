using LinearAlgebra

const GRAVITATIONAL_CONSTANT = 6.67430e-11

mutable struct Body
    mass::Float64
    position::Vec3{Float64}
    velocity::Vec3{Float64}
end

struct System
    bodies::Vector{Body}
end

distance_squared(b1::Body, b2::Body) = sum((b1.position - b2.position).^2)
distance_squared(p1::Vec3{Float64}, p2::Vec3{Float64}) = sum((p1 - p2).^2)
distance(b1::Body, b2::Body) = sqrt(distance_squared(b1, b2))
distance(p1::Vec3{Float64}, p2::Vec3{Float64}) = sqrt(distance_squared(p1, p2))

function calculate_force!(body_i::Body, body_j::Body, force::Vec3{Float64})
    r_sq = distance_squared(body_i, body_j)
    if r_sq > 1e-9
        r = sqrt(r_sq)
        magnitude = (GRAVITATIONAL_CONSTANT * body_i.mass * body_j.mass) / r_sq
        force .+= magnitude * (body_j.position - body_i.position) / r
    end
    return force
end

function calculate_total_energy(system::System)
    kinetic_energy = 0.0
    potential_energy = 0.0
    num_bodies = length(system.bodies)

    for i in 1:num_bodies
        body_i = system.bodies[i]
        v_sq = sum(body_i.velocity.^2)
        kinetic_energy += 0.5 * body_i.mass * v_sq

        for j in i+1:num_bodies
            body_j = system.bodies[j]
            r = distance(body_i, body_j)
            potential_energy -= (GRAVITATIONAL_CONSTANT * body_i.mass * body_j.mass) / r
        end
    end
    return kinetic_energy + potential_energy
end

function initialize_circular_orbits(num_orbiting::Int, central_mass::Float64, orbit_radius::Float64, orbiting_mass::Float64)
    bodies = Vector{Body}()
    # Initialize the central body
    push!(bodies, Body(central_mass, zeros(Vec3{Float64}), zeros(Vec3{Float64})))

    # Initialize the orbiting bodies
    for i in 0:(num_orbiting - 1)
        angle = 2 * pi * i / num_orbiting
        x = orbit_radius * cos(angle)
        y = orbit_radius * sin(angle)
        z = 0.0 # Orbiting in the xy-plane
        position = Vec3(x, y, z)

        # Calculate the orbital velocity for a circular orbit
        orbital_speed = sqrt(GRAVITATIONAL_CONSTANT * central_mass / orbit_radius)
        vx = -orbital_speed * sin(angle)
        vy = orbital_speed * cos(angle)
        vz = 0.0
        velocity = Vec3(vx, vy, vz)

        push!(bodies, Body(orbiting_mass, position, velocity))
    end
    return System(bodies)
end

function kick_step!(system::System, dt::Float64)
    num_bodies = length(system.bodies)
    forces = [zeros(Vec3{Float64}) for _ in 1:num_bodies]

    # Calculate forces
    for i in 1:num_bodies
        body_i = system.bodies[i]
        force_on_i = forces[i]
        for j in 1:num_bodies
            if i != j
                body_j = system.bodies[j]
                calculate_force!(body_i, body_j, force_on_i)
            end
        end
    end

    # Update velocities (kick)
    for i in 1:num_bodies
        body = system.bodies[i]
        force = forces[i]
        body.velocity .+= (force / body.mass) * (dt / 2)
    end
end

function drift_step!(system::System, dt::Float64)
    for body in system.bodies
        body.position .+= body.velocity * dt
    end
end

function first_order_kick_step!(system::System, dt::Float64)
    kick_step!(system, dt)
    drift_step!(system, dt)
    kick_step!(system, dt)
end

function main()
    num_orbiting_bodies = 1_000_000
    central_mass = 1.989e30    # Mass of the Sun (kg)
    orbit_radius = 1.496e11    # 1 AU (m)
    orbiting_mass = 5.972e24   # Mass of the Earth (kg)
    num_steps = 1000
    time_step = 3600.0 * 24.0 * 7.0 # 1 week in seconds

    # Initialize the system
    system = initialize_circular_orbits(num_orbiting_bodies, central_mass, orbit_radius, orbiting_mass)
    initial_system = deepcopy(system)

    println("Initial number of bodies: $(length(system.bodies))")

    # Calculate initial energy
    initial_energy = calculate_total_energy(initial_system)
    println("Initial total energy: $(initial_energy) J")

    println("Running simulation for $num_steps steps...")
    start_time = time()
    for step in 1:num_steps
        first_order_kick_step!(system, time_step)
        if step % 100 == 0
            println("Step $step completed.")
        end
    end
    end_time = time()
    elapsed_time = end_time - start_time
    println("Simulation finished in $(elapsed_time) seconds.")

    # Calculate final energy
    final_energy = calculate_total_energy(system)
    println("Final total energy: $(final_energy) J")

    # Calculate the energy difference
    energy_difference = abs(final_energy - initial_energy)
    relative_energy_difference = energy_difference / abs(initial_energy)
    println("Absolute energy difference: $(energy_difference) J")
    println("Relative energy difference: $(relative_energy_difference)")
end

main()