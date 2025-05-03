using Random

const G = 6.67430e-11  # Gravitational constant
const dt = 1e-3        # Time step
const numBodies = 1000000  # Number of small bodies

struct Body
    mass::Float64
    pos::Vector{Float64}
    vel::Vector{Float64}
    force::Vector{Float64}
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

function simulateStep!(bodies::Vector{Body})
    # Reset forces
    for body in bodies
        body.force .= 0.0
    end

    # Calculate forces
    for i in 1:length(bodies)
        for j in i+1:length(bodies)
            dx = bodies[j].pos[1] - bodies[i].pos[1]
            dy = bodies[j].pos[2] - bodies[i].pos[2]
            dz = bodies[j].pos[3] - bodies[i].pos[3]
            dist = sqrt(dx*dx + dy*dy + dz*dz)
            force = G * bodies[i].mass * bodies[j].mass / (dist * dist * dist)

            bodies[i].force[1] += force * dx
            bodies[i].force[2] += force * dy
            bodies[i].force[3] += force * dz

            bodies[j].force[1] -= force * dx
            bodies[j].force[2] -= force * dy
            bodies[j].force[3] -= force * dz
        end
    end

    # Update velocities and positions
    for body in bodies
        body.vel .+= body.force ./ body.mass .* dt
        body.pos .+= body.vel .* dt
    end
end

function calculateEnergy(bodies::Vector{Body})
    kineticEnergy = 0.0
    potentialEnergy = 0.0

    for i in 1:length(bodies)
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

    # Calculate initial energy
    initialEnergy = calculateEnergy(bodies)
    println("Initial Energy: $initialEnergy")

    # Run simulation for 1000 steps
    for step in 1:1000
        simulateStep!(bodies)
    end

    # Calculate final energy
    finalEnergy = calculateEnergy(bodies)
    println("Final Energy: $finalEnergy")
end

main()