const G = 6.67430e-11 # Gravitational constant
const DT = 1e-3 # Time step
const NUM_BODIES = 1000000 # Number of bodies

struct Body
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
    mass::Float64
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
    Threads.@threads for i in 1:length(bodies)
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
        Threads.atomic_add!(energy, kinetic + 0.5 * potential)
    end
    return energy
end

function kick_step!(bodies)
    Threads.@threads for i in 1:length(bodies)
        body = bodies[i]
        ax, ay, az = 0.0, 0.0, 0.0
        for j in 1:length(bodies)
            other = bodies[j]
            if body != other
                dx = other.x - body.x
                dy = other.y - body.y
                dz = other.z - body.z
                distance = sqrt(dx^2 + dy^2 + dz^2)
                force = G * other.mass / (distance^3)
                ax += force * dx
                ay += force * dy
                az += force * dz
            end
        end
        body.vx += ax * DT
        body.vy += ay * DT
        body.vz += az * DT
    end

    Threads.@threads for i in 1:length(bodies)
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
        kick_step!(bodies)
    end

    final_energy = calculate_energy(bodies)
    println("Final energy: $final_energy")
end

main()
