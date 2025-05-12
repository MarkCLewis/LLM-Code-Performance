using Random
using Printf
using Base.Threads

const G = 6.67430e-11
const DT = 1.0
const N_BODIES = 100_000
const STEPS = 10
const THETA = 0.3

mutable struct Body
    x::Float64; y::Float64; z::Float64
    vx::Float64; vy::Float64; vz::Float64
    ax::Float64; ay::Float64; az::Float64
    mass::Float64
end

mutable struct KDNode
    min::NTuple{3, Float64}
    max::NTuple{3, Float64}
    cm::NTuple{3, Float64}
    mass::Float64
    body::Union{Body, Nothing}
    left::Union{KDNode, Nothing}
    right::Union{KDNode, Nothing}
end

function initialize_system(n)
    bodies = Vector{Body}(undef, n + 1)
    bodies[1] = Body(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1e20)
    rng = MersenneTwister(42)
    radius = 1e7

    for i in 2:n+1
        angle = 2π * (i - 1) / n
        r = radius * (1 + 0.1 * rand(rng))
        x, y, z = r * cos(angle), r * sin(angle), 0.0
        v = sqrt(G * bodies[1].mass / r)
        vx, vy, vz = -v * sin(angle), v * cos(angle), 0.0
        bodies[i] = Body(x, y, z, vx, vy, vz, 0.0, 0.0, 0.0, 1.0)
    end
    return bodies
end

function build_kdtree(bodies::Vector{Body}, depth::Int = 0)::Union{KDNode, Nothing}
    if isempty(bodies)
        return nothing
    end

    axis = depth % 3 + 1
    sort!(bodies, by = b -> getfield(b, (:x, :y, :z)[axis]))
    mid = length(bodies) ÷ 2
    left = build_kdtree(view(bodies, 1:mid), depth + 1)
    right = build_kdtree(view(bodies, mid+1:length(bodies)), depth + 1)

    min_bound = mapreduce(b -> (b.x, b.y, b.z), (a, b) -> (min(a[1], b[1]), min(a[2], b[2]), min(a[3], b[3])), bodies)
    max_bound = mapreduce(b -> (b.x, b.y, b.z), (a, b) -> (max(a[1], b[1]), max(a[2], b[2]), max(a[3], b[3])), bodies)

    total_mass = sum(b.mass for b in bodies)
    cm_x = sum(b.mass * b.x for b in bodies) / total_mass
    cm_y = sum(b.mass * b.y for b in bodies) / total_mass
    cm_z = sum(b.mass * b.z for b in bodies) / total_mass

    return KDNode(min_bound, max_bound, (cm_x, cm_y, cm_z), total_mass,
                  length(bodies) == 1 ? bodies[1] : nothing, left, right)
end

function compute_force_from_node!(b::Body, node::KDNode)
    if node === nothing || node.mass == 0.0 || node.body === b
        return
    end

    dx = node.cm[1] - b.x
    dy = node.cm[2] - b.y
    dz = node.cm[3] - b.z
    dist2 = dx^2 + dy^2 + dz^2 + 1e-10
    dist = sqrt(dist2)

    size = maximum(node.max[i] - node.min[i] for i in 1:3)
    if node.body !== nothing || (size / dist < THETA)
        force = G * node.mass / (dist2 * dist)
        b.ax += dx * force
        b.ay += dy * force
        b.az += dz * force
    else
        compute_force_from_node!(b, node.left)
        compute_force_from_node!(b, node.right)
    end
end

function compute_forces!(bodies, tree)
    @threads for i in eachindex(bodies)
        b = bodies[i]
        b.ax = b.ay = b.az = 0.0
        compute_force_from_node!(b, tree)
    end
end

function update_bodies!(bodies)
    @threads for b in bodies
        b.vx += b.ax * DT
        b.vy += b.ay * DT
        b.vz += b.az * DT
        b.x += b.vx * DT
        b.y += b.vy * DT
        b.z += b.vz * DT
    end
end

function compute_energy(bodies)
    kinetic = zeros(Float64, length(bodies))
    potential = zeros(Float64, length(bodies))

    @threads for i in 1:length(bodies)
        bi = bodies[i]
        kinetic[i] = 0.5 * bi.mass * (bi.vx^2 + bi.vy^2 + bi.vz^2)
        pot = 0.0
        for j in i+1:length(bodies)
            bj = bodies[j]
            dx = bi.x - bj.x
            dy = bi.y - bj.y
            dz = bi.z - bj.z
            dist = sqrt(dx^2 + dy^2 + dz^2 + 1e-10)
            pot -= G * bi.mass * bj.mass / dist
        end
        potential[i] = pot
    end
    return sum(kinetic) + sum(potential)
end

function main()
    println("Initializing bodies...")
    bodies = initialize_system(N_BODIES)

    println("Computing initial energy...")
    energy0 = compute_energy(bodies)
    @printf("Initial energy: %.6e\n", energy0)

    for step in 1:STEPS
        tree = build_kdtree(copy(bodies))
        compute_forces!(bodies, tree)
        update_bodies!(bodies)
        if step % 100 == 0
            println("Step $step")
        end
    end

    println("Computing final energy...")
    energy1 = compute_energy(bodies)
    @printf("Final energy: %.6e\n", energy1)
    @printf("Energy difference: %.6e\n", abs(energy1 - energy0))
end

main()
