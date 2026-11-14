# Whole simulations: a bus fault, the CSV output, and a step without memory
# allocation (paper, Table I).

using Test
using NodalEMT

@testset "bus fault: the voltage collapses and recovers" begin
    case = test_case(:two_area)
    sim = build_simulation(case)
    k = bus_index(case.pfd, 8); nb = nbus(case.pfd)
    schedule_bus_fault!(sim; bus=8, t_on=0.02, t_off=0.08, Zflt=1e-3)
    res = run!(sim, 0.3; record_every=1, progress=false)
    t = res.t
    pre = maximum(vmag(res, k, nb, j) for j in eachindex(t) if t[j] < 0.02)
    during = maximum(vmag(res, k, nb, j) for j in eachindex(t) if 0.03 < t[j] < 0.08)
    after = sum(vmag(res, k, nb, j) for j in eachindex(t) if t[j] > 0.25) / count(>(0.25), t)
    @test pre > 0.9
    @test during < 0.05 * pre
    @test after > 0.8 * pre
    @test maximum(abs.(res.omega[:, end] .- case.pfd.ws)) / case.pfd.ws < 0.02
end

@testset "run! and save_results (WECC)" begin
    case = test_case(:wecc)
    sim = build_simulation(case)
    schedule_reference_step!(sim; t=0.002, gen=1, kind=:exciter, dsp=0.05)
    res = run!(sim, 200 * 50e-6; record_every=10, progress=false)
    @test length(res.t) == 21 && size(res.V) == (sim.net.N, 21) && all(isfinite, res.V)
    dir = mktempdir()
    save_results(res, sim, dir)
    @test isfile(joinpath(dir, "voltages.csv")) && isfile(joinpath(dir, "generators.csv"))
    @test isfile(joinpath(dir, "buses.csv")) && isfile(joinpath(dir, "ibr.csv"))
end

@testset "a time step allocates no memory" begin
    case = test_case(:wecc)
    sim = build_simulation(case)
    schedule_bus_fault!(sim; bus=2404, t_on=0.001, t_off=1.0, Zflt=1e-3)   # compensation in effect
    for _ in 1:50; step!(sim); end
    step!(sim)
    @test (@allocated step!(sim)) == 0
end
