using Base.Threads

const G = 6.67430e-11 # Gravitational constant
const DT = 1e-3 # Time step
const NUM_BODIES = 1000000 # Number of bodies
const THETA = 0.3 # Theta value for approximation

struct Body
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
    mass::Float64
end

struct KDNode
    body::Body
    left::Union{KDNode, Nothing}
    right::Union{KDNode, Nothing}
    min::Vector{Float64}
    max::Vector{Float64}
end

function initialize_bodies()
    bodies = Vector{Body}(undef, NUM_BODIES)
    bodies[1] = Body(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1e30) # Central body mass

    for i in 2:NUM_BODIES
        angle = 2 * π * (i - 1) / (NUM_BODIES - 1)
        bodies[i] = Body(
            cos(angle) * 1e11,
            sin(angle) * 1e11,
            0.0,
            -sin(angle) * sqrt(G * bodies[1].mass / 1e11),
            cos(angle) * sqrt(G * bodies[1].mass / 1e11),
            0.0,
            1e24 # Small body mass
        )
    end
    return bodies
end

function calculate_energy(bodies)
    energy = 0.0
    @threads for i in 1:length(bodies)
        body = bodies[i]
        kinetic = 0.5 * body.mass * (body.vx^2 + body.vy^2 + body.vz^2)
        potential = 0.0
        for j in 1:length(bodies)
            other = bodies[j]
            if body != other
                dx = body.x - other.x
                dy = body.y - other.y
                dz = body.z - other.z
                distance = sqrt(dx^2 + dy^2 + dz^2)
                potential -= G * body.mass * other.mass / distance
            end
        end
        @atomic energy += kinetic + 0.5 * potential
    end
    return energy
end

function build_kdtree(bodies, depth=0)
    if isempty(bodies)
        return nothing
    end

    axis = depth % 3 + 1
    sort!(bodies, by=b -> getfield(b, (:x, :y, :z)[axis]))

    median = length(bodies) ÷ 2
    node = KDNode(
        bodies[median],
        build_kdtree(bodies[1:median-1], depth + 1),
        build_kdtree(bodies[median+1:end], depth + 1),
        [bodies[median].x, bodies[median].y, bodies[median].z],
        [bodies[median].x, bodies[median].y, bodies[median].z]
    )

    for i in 1:3
        if node.left != nothing
            node.min[i] = min(node.min[i], node.left.min[i])
            node.max[i] = max(node.max[i], node.left.max[i])
        end
        if node.right != nothing
            node.min[i] = min(node.min[i], node.right.min[i])
            node.max[i] = max(node.max[i], node.right.max[i])
        end
    end

    return node
end

function calculate_force(node, body, ax, ay, az)
    if node == nothing || node.body == nothing
        return
    end

    dx = node.body.x - body.x
    dy = node.body.y - body.y
    dz = node.body.z - body.z
    distance = sqrt(dx^2 + dy^2 + dz^2)

    size = max(node.max[1] - node.min[1], max(node.max[2] - node.min[2], node.max[3] - node.min[3]))

    if size / distance < THETA || (node.left == nothing && node.right == nothing)
        force = G * node.body.mass / (distance^3)
        ax[] += force * dx
        ay[] += force * dy
        az[] += force * dz
    else
        calculate_force(node.left, body, ax, ay, az)
        calculate_force(node.right, body, ax, ay, az)
    end
end

function kick_step!(bodies, root)
    @threads for i in 1:length(bodies)
        body = bodies[i]
        ax, ay, az = Ref(0.0), Ref(0.0), Ref(0.0)
        calculate_force(root, body, ax, ay, az)
        body.vx += ax[] * DT
        body.vy += ay[] * DT
        body.vz += az[] * DT
    end

    @threads for i in 1:length(bodies)
        body = bodies[i]
        body.x += body.vx * DT
        body.y += body.vy * DT
        body.z += body.vz * DT
    end
end

function main()
    bodies = initialize_bodies()

    initial_energy = calculate_energy(bodies)
    println("Initial energy: $initial_energy")

    for step in 1:1000
        root = build_kdtree(bodies)
        kick_step!(bodies, root)
    end

    final_energy = calculate_energy(bodies)
    println("Final energy: $final_energy")
end

main()
