using Random
using Base.Threads

const G = 6.67430e-11  # Gravitational constant
const dt = 1e-3        # Time step
const numBodies = 1000000  # Number of small bodies
const theta = 0.3      # Barnes-Hut approximation parameter

struct Body
    mass::Float64
    pos::Vector{Float64}
    vel::Vector{Float64}
    force::Vector{Float64}
end

mutable struct Node
    body::Union{Body, Nothing}
    center::Vector{Float64}
    size::Float64
    mass::Float64
    com::Vector{Float64}
    children::Vector{Union{Node, Nothing}}

    Node(center, size) = new(nothing, center, size, 0.0, [0.0, 0.0, 0.0], Vector{Union{Node, Nothing}}(undef, 8))
end

mutable struct NodePool
    pool::Vector{Node}
    index::Int

    NodePool(size::Int) = new(Vector{Node}(undef, size), 1)

    function allocate!(pool::NodePool, center::Vector{Float64}, size::Float64)
        if pool.index > length(pool.pool)
            push!(pool.pool, Node(center, size))
        else
            pool.pool[pool.index] = Node(center, size)
        end
        pool.index += 1
        return pool.pool[pool.index - 1]
    end

    function reset!(pool::NodePool)
        pool.index = 1
    end
end

function initializeBodies(n::Int)
    bodies = Vector{Body}(undef, n + 1)

    # Central body
    bodies[1] = Body(1e30, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])

    # Small bodies in circular orbits
    for i in 2:n+1
        radius = 1e11 * rand()
        speed = sqrt(G * bodies[1].mass / radius)
        angle = 2 * π * rand()

        bodies[i] = Body(1e24,
                         [radius * cos(angle), radius * sin(angle), 0.0],
                         [-speed * sin(angle), speed * cos(angle), 0.0],
                         [0.0, 0.0, 0.0])
    end

    return bodies
end

function buildTree(bodies::Vector{Body}, pool::NodePool)
    root = allocate!(pool, [0.0, 0.0, 0.0], 2 * 1e11)
    @threads for body in bodies
        insertBody!(root, body, pool)
    end
    return root
end

function insertBody!(node::Node, body::Body, pool::NodePool)
    if node.body === nothing && node.mass == 0.0
        node.body = body
        node.mass = body.mass
        node.com = copy(body.pos)
        return
    end

    if node.body !== nothing
        oldBody = node.body
        node.body = nothing
        insertBody!(node, oldBody, pool)
    end

    node.mass += body.mass
    for i in 1:3
        node.com[i] = (node.com[i] * (node.mass - body.mass) + body.pos[i] * body.mass) / node.mass
    end

    index = 0
    for i in 1:3
        if body.pos[i] > node.center[i]
            index |= 1 << (i - 1)
        end
    end

    if node.children[index] === nothing
        newCenter = copy(node.center)
        for i in 1:3
            newCenter[i] += (index & (1 << (i - 1)) != 0 ? 1 : -1) * node.size / 4
        end
        node.children[index] = allocate!(pool, newCenter, node.size / 2)
    end

    insertBody!(node.children[index], body, pool)
end

function calculateForce!(body::Body, node::Node)
    if node.body !== nothing && node.body !== body
        dx = node.body.pos[1] - body.pos[1]
        dy = node.body.pos[2] - body.pos[2]
        dz = node.body.pos[3] - body.pos[3]
        dist = sqrt(dx*dx + dy*dy + dz*dz)
        force = G * body.mass * node.body.mass / (dist * dist * dist)

        body.force[1] += force * dx
        body.force[2] += force * dy
        body.force[3] += force * dz
        return
    end

    dx = node.com[1] - body.pos[1]
    dy = node.com[2] - body.pos[2]
    dz = node.com[3] - body.pos[3]
    dist = sqrt(dx*dx + dy*dy + dz*dz)

    if node.size / dist < theta
        force = G * body.mass * node.mass / (dist * dist * dist)

        body.force[1] += force * dx
        body.force[2] += force * dy
        body.force[3] += force * dz
        return
    end

    for child in node.children
        if child !== nothing
            calculateForce!(body, child)
        end
    end
end

function simulateStep!(bodies::Vector{Body}, pool::NodePool)
    # Reset forces
    @threads for i in 1:length(bodies)
        bodies[i].force .= 0.0
    end

    # Build kD-tree
    reset!(pool)
    root = buildTree(bodies, pool)

    # Calculate forces using kD-tree
    @threads for i in 1:length(bodies)
        calculateForce!(bodies[i], root)
    end

    # Update velocities and positions
    @threads for i in 1:length(bodies)
        bodies[i].vel .+= bodies[i].force ./ bodies[i].mass .* dt
        bodies[i].pos .+= bodies[i].vel .* dt
    end
end

function calculateEnergy(bodies::Vector{Body})
    kineticEnergy = 0.0
    potentialEnergy = 0.0

    @threads for i in 1:length(bodies)
        kineticEnergy += 0.5 * bodies[i].mass * sum(bodies[i].vel .^ 2)

        for j in i+1:length(bodies)
            dx = bodies[j].pos[1] - bodies[i].pos[1]
            dy = bodies[j].pos[2] - bodies[i].pos[2]
            dz = bodies[j].pos[3] - bodies[i].pos[3]
            dist = sqrt(dx*dx + dy*dy + dz*dz)
            potentialEnergy -= G * bodies[i].mass * bodies[j].mass / dist
        end
    end

    return kineticEnergy + potentialEnergy
end

function main()
    Random.seed!(1234)

    # Initialize bodies
    bodies = initializeBodies(numBodies)

    # Create a node pool
    pool = NodePool(numBodies)

    # Calculate initial energy
    initialEnergy = calculateEnergy(bodies)
    println("Initial Energy: $initialEnergy")

    # Run simulation for 1000 steps
    for step in 1:1000
        simulateStep!(bodies, pool)
    end

    # Calculate final energy
    finalEnergy = calculateEnergy(bodies)
    println("Final Energy: $finalEnergy")
end

main()