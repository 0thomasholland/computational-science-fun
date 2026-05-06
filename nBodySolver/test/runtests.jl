using Test
using nBodySolver
using StaticArrays

@testset "nBodySolver" begin

    @testset "Energy conservation (many steps)" begin
        m = 1e30
        r = 1.5e11
        bodies = State(
            [SVector(-r/2, 0.0, 0.0), SVector(r/2, 0.0, 0.0)],
            [SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0)],
            [m, m]
        )
        sim = SimState(bodies, 0.0)
        E0 = total_energy(sim.bodies)

        dt = 1e5
        nsteps = 10_000
        for _ in 1:nsteps
            step!(sim, dt)
        end

        E1 = total_energy(sim.bodies)
        @test abs((E1 - E0) / E0) < 1e-4
    end

    @testset "Kepler circular orbit energy conservation" begin
        G_val = nBodySolver.G
        m = 1e30
        r = 1.5e11
        v_circ = sqrt(G_val * m / (2 * r))

        bodies = State(
            [SVector(-r/2, 0.0, 0.0), SVector(r/2, 0.0, 0.0)],
            [SVector(0.0, -v_circ, 0.0), SVector(0.0, v_circ, 0.0)],
            [m, m]
        )
        sim = SimState(bodies, 0.0)

        T = 2π * sqrt(r^3 / (2 * G_val * m))
        dt = T / 10_000
        nsteps = 10_000

        E0 = total_energy(sim.bodies)
        for _ in 1:nsteps
            step!(sim, dt)
        end
        E1 = total_energy(sim.bodies)

        @test abs((E1 - E0) / E0) < 1e-4

        dx = sim.bodies.positions[2] - sim.bodies.positions[1]
        separation = sqrt(sum(dx .^ 2))
        @test abs(separation / r - 1.0) < 0.01
    end

end
