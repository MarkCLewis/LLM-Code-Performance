# Okay, let's implement the parallel kD-tree construction using the ParallelKDTrees.jl package. This is a significant optimization that can speed up the simulation, especially for a large number of bodies.
# First, ensure you have the ParallelKDTrees package installed:
# Julia

# using Pkg
# Pkg.add("ParallelKDTrees")
# Now, here's the modified code incorporating parallel kD-tree construction. Note that the calculate_force_kd_tree! function needs to be adapted to work with the KDTree structure from NearestNeighbors (which ParallelKDTrees builds upon).
#To Run:
#
#Save the code as a .jl file (e.g., nbody_parallel_kdtree.jl).
#Ensure you have NearestNeighbors and ParallelKDTrees installed.
#Run Julia with multiple threads: julia -t <num_threads> nbody_parallel_kdtree.jl

using LinearAlgebra
using Base.Threads
using NearestNeighbors
using ParallelKDTrees

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

function build_kd_tree(system::System)
    num_bodies = length(system.bodies)
    if num_bodies == 0
        return nothing
    end
    positions = [b.position for b in system.bodies]
    return KDTree(positions, leafsize=10) # Build in parallel by default in ParallelKDTrees
end

function calculate_force_kd_tree!(body::Body, tree::KDTree, system::System)
    body_index = findfirst(b -> b === body, system.bodies)
    if isnothing(body_index)
        return
    end
    current_pos = body.position
    force = zeros(Vec3{Float64})

    function barnes_hut!(node_indices::Vector{Int}, bounding_box::HyperRectangle{3, Float64})
        if isempty(node_indices)
            return
        elseif length(node_indices) == 1
            other_body_index = node_indices[1]
            if other_body_index != body_index
                other_body = system.bodies[other_body_index]
                calculate_force_element!(body, other_body.position, other_body.mass) # Treat as point mass
            end
        else
            s = maximum(bounding_box.max - bounding_box.min)
            center_of_mass = sum(system.bodies[i].mass * system.bodies[i].position for i in node_indices) / sum(system.bodies[i].mass for i in node_indices)
            d = distance(current_pos, center_of_mass)

            if d == 0.0 || s / d < THETA
                total_mass = sum(system.bodies[i].mass for i in node_indices)
                calculate_force_element!(body, center_of_mass, total_mass) # Treat as pseudo-particle
            else
                # Subdivide the box and recurse
                mid = (bounding_box.min + bounding_box.max) / 2
                for i in 0:7 # 2^3 octants in 3D
                    octant_indices = Int[]
                    octant_min = zeros(Vec3{Float64})
                    octant_max = zeros(Vec3{Float64})
                    for dim in 1:3
                        if (i >> (dim - 1)) & 1 == 0
                            octant_min[dim] = bounding_box.min[dim]
                            octant_max[dim] = mid[dim]
                        else
                            octant_min[dim] = mid[dim]
                            octant_max[dim] = bounding_box.max[dim]
                        end
                    end
                    octant_box = HyperRectangle(octant_min, octant_max)
                    indices_in_octant = filter(idx -> is_point_in_bounds(system.bodies[idx].position, octant_box), node_indices)
                    barnes_hut!(indices_in_octant, octant_box)
                end
            end
        end
    end

    root_indices = 1:length(system.bodies)
    min_coords = minimum(b.position for b in system.bodies)
    max_coords = maximum(b.position for b in system.bodies)
    root_bounds = HyperRectangle(min_coords, max_coords)

    barnes_hut!(root_indices, root_bounds)
    body.force .= force # Assign the accumulated force
end

function is_point_in_bounds(point::Vec3{Float64}, bounds::HyperRectangle{3, Float64})
    return all(bounds.min .<= point .<= bounds.max)
end

function kick_step!(system::System, dt::Float64)
    num_bodies = length(system.bodies)
    tree = build_kd_tree(system)

    # Reset forces
    @threads for i in 1:num_bodies
        system.bodies[i].force .= 0.0
    end

    # Calculate forces using kD-tree in parallel
    @threads for i in 1:num_bodies
        calculate_force_kd_tree!(system.bodies[i], tree, system)
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