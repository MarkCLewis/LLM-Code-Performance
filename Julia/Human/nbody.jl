# The Computer Language Benchmarks Game
# https://salsa.debian.org/benchmarksgame-team/benchmarksgame/
#
# Contributed by Andrei Fomiga, Stefan Karpinski, Viral B. Shah, Jeff
# Bezanson, and Adam Beckmeyer.
# Based on Mark C. Lewis's Java version.

using Printf

# Constants
const SOLAR_MASS = 4 * pi * pi
const DAYS_PER_YEAR = 365.24

# Use a struct instead of mutable struct since a struct can be stored
# inline in an array avoiding the overhead of following a pointer
struct Body
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
    m::Float64
end

function init_sun(bodies)
    px = py = pz = 0.0
    for b in bodies
        px -= b.vx * b.m
        py -= b.vy * b.m
        pz -= b.vz * b.m
    end
    Body(0.0, 0.0, 0.0, px / SOLAR_MASS, py / SOLAR_MASS, pz / SOLAR_MASS, SOLAR_MASS)
end

function advance!(bodies, dt)
    n = length(bodies)
    @inbounds for i=1:n-1
        bi = bodies[i]

        # Since the fields of bi aren't mutable, we track the changing
        # value of bi's velocity outside of the Body struct
        ivx = bi.vx
        ivy = bi.vy
        ivz = bi.vz

        for j=i+1:n
            bj = bodies[j]
            
            dx = bi.x - bj.x
            dy = bi.y - bj.y
            dz = bi.z - bj.z

            dsq = dx^2 + dy^2 + dz^2
            mag = dt / (dsq * √dsq)

            ivx -= dx * bj.m * mag
            ivy -= dy * bj.m * mag
            ivz -= dz * bj.m * mag

            bodies[j] = Body(
                bj.x, bj.y, bj.z,
                bj.vx + dx * bi.m * mag,
                bj.vy + dy * bi.m * mag,
                bj.vz + dz * bi.m * mag,
                bj.m
            )
        end

        bodies[i] = Body(
            bi.x, bi.y, bi.z,
            ivx, ivy, ivz,
            bi.m
        )
    end

    @inbounds for i=1:n
        bi = bodies[i]
        bodies[i] = Body(
            bi.x + dt * bi.vx, bi.y + dt * bi.vy, bi.z + dt * bi.vz,
            bi.vx, bi.vy, bi.vz,
            bi.m
        )
    end
end

function energy(bodies)
    n = length(bodies)
    e = 0.0
    @inbounds for i=1:n
        bi = bodies[i]

        e += 0.5 * bi.m * (bi.vx^2 + bi.vy^2 + bi.vz^2)
        for j=i+1:n
            bj = bodies[j]
            
            d = √((bi.x - bj.x)^2 + (bi.y - bj.y)^2 + (bi.z - bj.z)^2)
            e -= bi.m * bodies[j].m / d
        end
    end
    e
end

function circular_orbits(n::Int64)::Array{Body}
  first = Body(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)
  bods = [first]
  for i in 1:n
      d = .1 + (i * 5.0 / n)
      v = sqrt(1.0 / d)
      theta = rand(Float64)*2*pi
      x = d * cos(theta)
      y = d * sin(theta)
      vx = -v * sin(theta)
      vy = v * cos(theta)
      temp = Body(x, y, 0.0, vx, vy, 0, 1.0e-7)
      push!(bods, temp) #is this pushing a reference to temp in which case is everything going to explode and how do i fix that
  end
  bods
  #print(bods)
end

function nbody(steps::Int64, numBodies::Int64)
    bodies = circular_orbits(numBodies)


    # @printf("%.9f\n", energy(bodies))
    for i = 1:steps
        advance!(bodies, 0.01)
    end
    # @printf("%.9f\n", energy(bodies))
end

isinteractive() ||  nbody(parse(Int, ARGS[1]), parse(Int, ARGS[2]))