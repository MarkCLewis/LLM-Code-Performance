using Random, Printf

const G = 6.67430e-11  # Gravitational constant
const DT = 1.0         # Time step
const NUM_BODIES = 1_000_000  # Number of bodies
const STEPS = 1000     # Simulation steps

mutable struct Body
    position::Vector{Float64}
    velocity::Vector{Float64}
    mass::Float64
end

function compute_energy(bodies::Vector{Body})
    kinetic = 0.0
    potential = 0.0
    for i in eachindex(bodies)
        kinetic += 0.5 * bodies[i].mass * sum(bodies[i].velocity .^ 2)
        for j in i+1:length(bodies)
            r = norm(bodies[i].position - bodies[j].position)
            potential -= G * bodies[i].mass * bodies[j].mass / r
        end
    end
    return kinetic + potential
end

function update_positions(bodies::Vector{Body})
    for body in bodies
        body.position .+= body.velocity .* DT
    end
end

function update_velocities(bodies::Vector{Body})
    for i in eachindex(bodies)
        acceleration = zeros(3)
        for j in eachindex(bodies)
            if i != j
                r_vec = bodies[j].position - bodies[i].position
                r = norm(r_vec)
                if r > 0.0
                    acceleration .+= G * bodies[j].mass / r^3 * r_vec
                end
            end
        end
        bodies[i].velocity .+= acceleration .* DT
    end
end

function initialize_orbiting_bodies(num_bodies::Int, central_mass::Float64)
    bodies = Vector{Body}(undef, num_bodies)
    bodies[1] = Body([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], central_mass)
    radius = 1.0e9
    for i in 2:num_bodies
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
        update_velocities(bodies)
        update_positions(bodies)
    end
    
    final_energy = compute_energy(bodies)
    @printf("Final Energy: %e\n", final_energy)
end

main()
