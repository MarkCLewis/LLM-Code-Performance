using Random, Printf, Base.Threads

const G = 6.67430e-11  # Gravitational constant
const DT = 1.0         # Time step
const NUM_BODIES = 1_000_000  # Number of bodies
const STEPS = 1000     # Simulation steps
const THETA = 0.3      # Barnes-Hut approximation threshold

mutable struct Body
    position::Vector{Float64}
    velocity::Vector{Float64}
    mass::Float64
end

mutable struct KDTree
    center_of_mass::Vector{Float64}
    total_mass::Float64
    boundary_min::Vector{Float64}
    boundary_max::Vector{Float64}
    children::Vector{Union{Nothing, KDTree}}
    body::Union{Nothing, Body}
end

function build_kdtree(bodies::Vector{Body})::Union{Nothing, KDTree}
    if isempty(bodies)
        return nothing
    end
    if length(bodies) == 1
        return KDTree(bodies[1].position, bodies[1].mass, bodies[1].position, bodies[1].position, fill(nothing, 8), bodies[1])
    end
    
    total_mass = sum(b.mass for b in bodies)
    center_of_mass = sum(b.position .* b.mass for b in bodies) / total_mass
    boundary_min = minimum(reduce(hcat, (b.position for b in bodies)), dims=2)
    boundary_max = maximum(reduce(hcat, (b.position for b in bodies)), dims=2)
    
    node = KDTree(center_of_mass, total_mass, boundary_min, boundary_max, fill(nothing, 8), nothing)
    
    sorted_bodies = sort(bodies, by=b -> b.position)
    median_idx = length(sorted_bodies) ÷ 2
    left_bodies = sorted_bodies[1:median_idx]
    right_bodies = sorted_bodies[median_idx+1:end]
    
    node.children[1] = build_kdtree(left_bodies)
    node.children[2] = build_kdtree(right_bodies)
    
    return node
end

function compute_force(body::Body, tree::KDTree, theta::Float64)
    if tree === nothing || tree.body === body
        return zeros(3)
    end
    r_vec = tree.center_of_mass - body.position
    r = norm(r_vec)
    s = maximum(tree.boundary_max - tree.boundary_min)
    if s / r < theta || tree.body !== nothing
        return G * tree.total_mass / r^3 * r_vec
    else
        force = zeros(3)
        for child in tree.children
            if child !== nothing
                force .+= compute_force(body, child, theta)
            end
        end
        return force
    end
end

function update_velocities(bodies::Vector{Body}, tree::KDTree)
    @threads for body in bodies
        force = compute_force(body, tree, THETA)
        body.velocity .+= force .* DT
    end
end

function initialize_orbiting_bodies(num_bodies::Int, central_mass::Float64)
    bodies = Vector{Body}(undef, num_bodies)
    bodies[1] = Body([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], central_mass)
    radius = 1.0e9
    @threads for i in 2:num_bodies
        angle = 2π * rand()
        pos = [radius * cos(angle), radius * sin(angle), 0.0]
        speed = sqrt(G * central_mass / radius)
        vel = [-speed * sin(angle), speed * cos(angle), 0.0]
        bodies[i] = Body(pos, vel, 1.0)
    end
    return bodies
end

function main()
    Random.seed!(1234)
    bodies = initialize_orbiting_bodies(NUM_BODIES, 1.989e30)
    
    initial_energy = compute_energy(bodies)
    @printf("Initial Energy: %e\n", initial_energy)
    
    for step in 1:STEPS
        tree = build_kdtree(bodies)
        update_velocities(bodies, tree)
        update_positions(bodies)
    end
    
    final_energy = compute_energy(bodies)
    @printf("Final Energy: %e\n", final_energy)
end

main()
