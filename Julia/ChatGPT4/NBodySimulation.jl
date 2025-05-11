using Random
using Printf

const G = 6.67430e-11
const DT = 3600.0 * 24.0
const N_BODIES = 10_000
const STEPS = 100

mutable struct Body
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
    ax::Float64
    ay::Float64
    az::Float64
    mass::Float64
end

function initialize_system(n)
    bodies = Vector{Body}(undef, n + 1)

    # Central body
    bodies[1] = Body(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1e20)

    radius = 1e7
    rng = MersenneTwister(42)

    for i in 2:n+1
        angle = 2π * (i-1) / n
        r = radius * (1.0 + 0.1 * rand(rng))
        x = r * cos(angle)
        y = r * sin(angle)
        z = 0.0
        v = sqrt(G * bodies[1].mass / r)
        vx = -v * sin(angle)
        vy = v * cos(angle)
        vz = 0.0
        bodies[i] = Body(x, y, z, vx, vy, vz, 0.0, 0.0, 0.0, 1.0)
    end

    return bodies
end

function compute_forces!(bodies)
    n = length(bodies)
    for i in 1:n
        bi = bodies[i]
        bi.ax = bi.ay = bi.az = 0.0
        for j in 1:n
            if i == j
                continue
            end
            bj = bodies[j]
            dx = bj.x - bi.x
            dy = bj.y - bi.y
            dz = bj.z - bi.z
            dist2 = dx^2 + dy^2 + dz^2 + 1e-10
            dist = sqrt(dist2)
            force = G * bj.mass / (dist2 * dist)
            bi.ax += dx * force
            bi.ay += dy * force
            bi.az += dz * force
        end
    end
end

function update_bodies!(bodies)
    for b in bodies
        b.vx += b.ax * DT
        b.vy += b.ay * DT
        b.vz += b.az * DT
        b.x += b.vx * DT
        b.y += b.vy * DT
        b.z += b.vz * DT
    end
end

function compute_energy(bodies)
    kinetic = 0.0
    potential = 0.0
    n = length(bodies)

    for i in 1:n
        b = bodies[i]
        v2 = b.vx^2 + b.vy^2 + b.vz^2
        kinetic += 0.5 * b.mass * v2

        for j in i+1:n
            b2 = bodies[j]
            dx = b.x - b2.x
            dy = b.y - b2.y
            dz = b.z - b2.z
            dist = sqrt(dx^2 + dy^2 + dz^2 + 1e-10)
            potential -= G * b.mass * b2.mass / dist
        end
    end

    return kinetic + potential
end

function main()
    println("Initializing bodies...")
    bodies = initialize_system(N_BODIES)

    println("Computing initial energy...")
    energy0 = compute_energy(bodies)
    @printf("Initial energy: %.6e\n", energy0)

    for step in 1:STEPS
        compute_forces!(bodies)
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
