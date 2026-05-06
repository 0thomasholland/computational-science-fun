module nBodySolver

using LinearAlgebra
using StaticArrays

export State, SimState, compute_forces, total_energy, step!

const G = 6.674e-11

struct State
    positions::Vector{SVector{3,Float64}}
    velocities::Vector{SVector{3,Float64}}
    masses::Vector{Float64}
end

mutable struct SimState
    bodies::State
    time::Float64
end

function compute_forces(bodies::State)
    N = length(bodies.masses)
    forces = [SVector(0.0, 0.0, 0.0) for _ in 1:N]
    for i in 1:N
        for j in (i+1):N
            r = bodies.positions[j] - bodies.positions[i]
            dist = norm(r)
            f = G * bodies.masses[i] * bodies.masses[j] / dist^2 * (r / dist)
            forces[i] += f
            forces[j] -= f
        end
    end
    return forces
end

function total_energy(bodies::State)
    N = length(bodies.masses)
    ke = sum(0.5 * bodies.masses[i] * dot(bodies.velocities[i], bodies.velocities[i]) for i in 1:N)
    pe = 0.0
    for i in 1:N
        for j in (i+1):N
            r = norm(bodies.positions[j] - bodies.positions[i])
            pe -= G * bodies.masses[i] * bodies.masses[j] / r
        end
    end
    return ke + pe
end

function step!(sim::SimState, dt::Float64)
    bodies = sim.bodies
    N = length(bodies.masses)

    forces = compute_forces(bodies)

    half_vel = [bodies.velocities[i] + 0.5 * dt * forces[i] / bodies.masses[i] for i in 1:N]

    new_pos = [bodies.positions[i] + dt * half_vel[i] for i in 1:N]

    new_bodies_mid = State(new_pos, half_vel, bodies.masses)
    forces2 = compute_forces(new_bodies_mid)

    new_vel = [half_vel[i] + 0.5 * dt * forces2[i] / bodies.masses[i] for i in 1:N]

    sim.bodies = State(new_pos, new_vel, bodies.masses)
    sim.time += dt
    return sim
end

end # module
