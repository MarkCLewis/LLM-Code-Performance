using Base.Threads
using StaticArrays
using Printf

const G = 6.67430e-11
const DT = 1.0
const THETA = 0.3
const EPS2 = 1e-10
const N_BODIES = 100_000
const STEPS = 10

typealias Vec3 SVector{3, Float64}

mutable struct Body
    pos::Vec3
    vel::Vec3
    acc::Vec3
    mass::Float64
end

mutable struct KDNode
    min::Vec3
    max::Vec3
    cm::Vec3
    mass::Float64
    left::Union{KDNode, Nothing}
    right::Union{KDNode, Nothing}
    bodies::Vector{Body}
end

function initialize_system(n::Int)
    bodies = Vector{Body}(undef, n + 1)
    bodies[1] = Body(SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), 1e20)
    radius = 1e7
    rng = MersenneTwister(42)

    @threads for i in 2:n+1
        angle = 2 * pi * i / n
        r = radius * (1.0 + 0.1 * rand(rng))
        x = r * cos(angle)
        y = r * sin(angle)
        z = 0.0
        v = sqrt(G * bodies[1].mass / r)
        vx = -v * sin(angle)
        vy = v * cos(angle)
        vz = 0.0
        bodies[i] = Body(SVector(x, y, z), SVector(vx, vy, vz), SVector(0.0, 0.0, 0.0), 1.0)
    end
    return bodies
end

function build_kdtree(bodies::Vector{Body}, depth::Int = 0)::Union{KDNode, Nothing}
    if isempty(bodies)
        return nothing
    end

    axis = depth % 3 + 1
    sort!(bodies, by = b -> b.pos[axis])
    mid = length(bodies) ÷ 2
    node_bodies = bodies

    min_coords = minimum(reduce(hcat, [b.pos for b in node_bodies]), dims=2)[:, 1]
    max_coords = maximum(reduce(hcat, [b.pos for b in node_bodies]), dims=2)[:, 1]
    total_mass = sum(b.mass for b in node_bodies)
    cm = sum(b.mass .* b.pos for b in node_bodies) / total_mass

    left = build_kdtree(view(node_bodies, 1:mid-1), depth + 1)
    right = build_kdtree(view(node_bodies, mid:end), depth + 1)

    return KDNode(SVector(min_coords...), SVector(max_coords...), cm, total_mass, left, right, node_bodies)
end

function compute_force!(b::Body, node::Union{KDNode, Nothing})
    if node === nothing || node.mass == 0.0 || b in node.bodies
        return
    end

    dx = node.cm - b.pos
    dist2 = sum(dx .^ 2) + EPS2
    size = maximum(node.max .- node.min)

    if length(node.bodies) <= 1 || (size^2 / dist2 < THETA^2)
        dist = sqrt(dist2)
        force = G * node.mass / (dist2 * dist)
        b.acc += force * dx
    else
        compute_force!(b, node.left)
        compute_force!(b, node.right)
    end
end

function compute_forces!(bodies::Vector{Body}, tree::KDNode)
    @threads for i in eachindex(bodies)
        b = bodies[i]
        b.acc = SVector(0.0, 0.0, 0.0)
        compute_force!(b, tree)
    end
end

function update_bodies!(bodies::Vector{Body})
    @threads for i in eachindex(bodies)
        b = bodies[i]
        b.vel += b.acc * DT
        b.pos += b.vel * DT
    end
end

function compute_energy(bodies::Vector{Body})
    kinetic = ThreadsX.sum(b -> 0.5 * b.mass * sum(b.vel .^ 2), bodies)
    potential = ThreadsX.sum(1:length(bodies)) do i
        pot = 0.0
        bi = bodies[i]
        for j in i+1:length(bodies)
            bj = bodies[j]
            dx = bi.pos - bj.pos
            dist = sqrt(sum(dx .^ 2) + EPS2)
            pot -= G * bi.mass * bj.mass / dist
        end
        pot
    end
    return kinetic + potential
end

function main()
    println("Initializing...")
    bodies = initialize_system(N_BODIES)
    println("Calculating initial energy...")
    initial_energy = compute_energy(bodies)
    @printf("Initial energy: %.6e\n", initial_energy)

    for step in 1:STEPS
        tree = build_kdtree(bodies)
        compute_forces!(bodies, tree)
        update_bodies!(bodies)
        if step % 100 == 0
            println("Step $step")
        end
    end

    println("Calculating final energy...")
    final_energy = compute_energy(bodies)
    @printf("Final energy: %.6e\n", final_energy)
    @printf("Energy difference: %.6e\n", abs(final_energy - initial_energy))
end

main()
