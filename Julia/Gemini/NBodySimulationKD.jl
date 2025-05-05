using LinearAlgebra
using Base.Threads
using NearestNeighbors  # For kD-tree implementation

const GRAVITATIONAL_CONSTANT = 6.67430e-11
const THETA = 0.3

mutable struct Body
    mass::Float64
    position::Vec3{Float64}
    velocity::Vec3{Float64}
    force::Vec3{Float64} # To store calculated force
end

struct System
    bodies::Vector{Body}
end

distance_squared(p1::Vec3{Float64}, p2::Vec3{Float64}) = sum((p1 - p2).^2)
distance(p1::Vec3{Float64}, p2::Vec3{Float64}) = sqrt(distance_squared(p1, p2))

function calculate_force_element!(body_i::Body, mass_element_cm::Vec3{Float64}, mass_element_mass::Float64)
    r_sq = distance_squared(body_i.position, mass_element_cm)
    if r_sq > 1e-9
        r = sqrt(r_sq)
        magnitude = (GRAVITATIONAL_CONSTANT * body_i.mass * mass_element_mass) / r_sq
        body_i.force .+= magnitude * (mass_element_cm - body_i.position) / r
    end
end

function calculate_total_energy(system::System)
    kinetic_energy = Threads.Atomic{Float64}(0.0)
    potential_energy = Threads.Atomic{Float64}(0.0)
    num_bodies = length(system.bodies)

    @threads for i in 1:num_bodies
        body_i = system.bodies[i]
        v_sq = sum(body_i.velocity.^2)
        Threads.atomic_add!(kinetic_energy, 0.5 * body_i.mass * v_sq)

        for j in i+1:num_bodies
            body_j = system.bodies[j]
            r = distance(body_i.position, body_j.position)
            Threads.atomic_add!(potential_energy, -(GRAVITATIONAL_CONSTANT * body_i.mass * body_j.mass) / r)
        end
    end
    return Threads.load(kinetic_energy) + Threads.load(potential_energy)
end

function initialize_circular_orbits(num_orbiting::Int, central_mass::Float64, orbit_radius::Float64, orbiting_mass::Float64)
    bodies = Vector{Body}()
    # Initialize the central body
    push!(bodies, Body(central_mass, zeros(Vec3{Float64}), zeros(Vec3{Float64}), zeros(Vec3{Float64})))

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

        push!(bodies, Body(orbiting_mass, position, velocity, zeros(Vec3{Float64})))
    end
    return System(bodies)
end

struct KDNode
    body_index::Int # -1 for internal node
    center_of_mass::Vec3{Float64}
    total_mass::Float64
    min_bound::Vec3{Float64}
    max_bound::Vec3{Float64}
    left::Union{Nothing, KDNode}
    right::Union{Nothing, KDNode}

    KDNode(index::Int, cm::Vec3{Float64}, mass::Float64, minb::Vec3{Float64}, maxb::Vec3{Float64}) =
        new(index, cm, mass, minb, maxb, nothing, nothing)

    KDNode(minb::Vec3{Float64}, maxb::Vec3{Float64}) =
        new(-1, zeros(Vec3{Float64}), 0.0, minb, maxb, nothing, nothing)
end

function build_kd_tree(system::System)
    num_bodies = length(system.bodies)
    if num_bodies == 0
        return nothing
    end
    indices = 1:num_bodies
    min_bound = minimum(b.position for b in system.bodies)
    max_bound = maximum(b.position for b in system.bodies)
    return build_kd_tree_recursive(system, indices, min_bound, max_bound, 1)
end

function build_kd_tree_recursive(system::System, indices::UnitRange{Int}, min_bound::Vec3{Float64}, max_bound::Vec3{Float64}, depth::Int)
    num_indices = length(indices)
    if num_indices == 0
        return nothing
    elseif num_indices == 1
        index = first(indices)
        body = system.bodies[index]
        return KDNode(index, body.position, body.mass, min_bound, max_bound)
    else
        node = KDNode(min_bound, max_bound)
        node.total_mass = 0.0
        center_of_mass = zeros(Vec3{Float64})
        for i in indices
            body = system.bodies[i]
            node.total_mass += body.mass
            center_of_mass .+= body.mass * body.position
        end
        node.center_of_mass = center_of_mass / node.total_mass

        split_dim = (depth - 1) % 3 + 1
        median = (min_bound[split_dim] + max_bound[split_dim]) / 2

        left_indices = filter(i -> system.bodies[i].position[split_dim] <= median, indices)
        right_indices = filter(i -> system.bodies[i].position[split_dim] > median, indices)

        left_min_bound = copy(min_bound)
        left_max_bound = copy(max_bound)
        left_max_bound[split_dim] = median

        right_min_bound = copy(min_bound)
        right_max_bound = copy(max_bound)
        right_min_bound[split_dim] = median

        node.left = build_kd_tree_recursive(system, left_indices, left_min_bound, left_max_bound, depth + 1)
        node.right = build_kd_tree_recursive(system, right_indices, right_min_bound, right_max_bound, depth + 1)
        return node
    end
end

function calculate_force_kd_tree!(body::Body, node::Union{Nothing, KDNode})
    if isnothing(node)
        return
    elseif node.body_index > 0
        if node.body_index != findfirst(b -> b === body, sys.bodies) # Avoid self-interaction
            calculate_force_element!(body, node.center_of_mass, node.total_mass)
        end
    else
        s = node.max_bound - node.min_bound
        d = distance(body.position, node.center_of_mass)
        if d == 0.0 || maximum(s) / d < THETA
            calculate_force_element!(body, node.center_of_mass, node.total_mass)
        else
            calculate_force_kd_tree!(body, node.left)
            calculate_force_kd_tree!(body, node.right)
        end
    end
end

function kick_step!(system::System, dt::Float64)
    num_bodies = length(system.bodies)
    kd_tree = build_kd_tree(system)

    # Reset forces
    @threads for i in 1:num_bodies
        system.bodies[i].force .= 0.0
    end

    # Calculate forces using kD-tree in parallel
    @threads for i in 1:num_bodies
        calculate_force_kd_tree!(system.bodies[i], kd_tree)
    end

    # Update velocities (kick)
    @threads for i in 1:num_bodies
        body = system.bodies[i]
        body.velocity .+= (body.force / body.mass) * (dt / 2)
    end
end

function drift_step!(system::System, dt::Float64)
    @threads for body in system.bodies
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
    global sys = initialize_circular_orbits(num_orbiting_bodies, central_mass, orbit_radius, orbiting_mass)
    initial_system = deepcopy(sys)

    println("Initial number of bodies: $(length(sys.bodies))")

    # Calculate initial energy
    initial_energy = calculate_total_energy(initial_system)
    println("Initial total energy: $(initial_energy) J")

    println("Running simulation for $num_steps steps...")
    start_time = time()
    for step in 1:num_steps
        first_order_kick_step!(sys, time_step)
        if step % 100 == 0
            println("Step $step completed.")
        end
    end
    end_time = time()
    elapsed_time = end_time - start_time
    println("Simulation finished in $(elapsed_time) seconds.")

    # Calculate final energy
    final_energy = calculate_total_energy(sys)
    println("Final total energy: $(final_energy) J")

    # Calculate the energy difference
    energy_difference = abs(final_energy - initial_energy)
    relative_energy_difference = energy_difference / abs(initial_energy)
    println("Absolute energy difference: $(energy_difference) J")
    println("Relative energy difference: $(relative_energy_difference)")
end

main()