# Scenario batches (paper, Section III): sharing the constants, the
# factorization and the worker must change nothing in any scenario.

using Test
using NodalEMT
using NodalEMT: chunk_ranges

@testset "chunks" begin
    @test chunk_ranges(10, 3) == [1:4, 5:7, 8:10]
    @test chunk_ranges(3, 8) == [1:1, 2:2, 3:3]           # never more chunks than scenarios
    for (n, k) in ((7, 3), (100, 7), (500, 16), (1, 1))
        r = chunk_ranges(n, k)
        @test vcat(collect.(r)...) == collect(1:n)         # contiguous, in order
        @test maximum(length, r) - minimum(length, r) <= 1  # equal sizes
    end
end

"Scenario s: a different disturbance in each, the last one undisturbed."
function configure!(sim, s)
    if s == 1
        schedule_bus_fault!(sim; bus=8, t_on=0.01, t_off=0.04, Zflt=1e-3)
    elseif s == 2
        schedule_generator_trip!(sim; t=0.02, gen=2)
    elseif s == 3
        schedule_reference_step!(sim; t=0.01, gen=1, kind=:governor, dsp=-0.02)
    elseif s == 4
        schedule_bus_fault!(sim; bus=7, t_on=0.01, t_off=0.03, Zflt=5e-2, phases=[1])
        schedule_generator_trip!(sim; t=0.025, gen=4)
        schedule_reference_step!(sim; t=0.005, gen=3, kind=:exciter, dsp=0.05)
    elseif s == 5
        schedule_bus_fault!(sim; bus=8, t_on=0.01, t_off=0.02, Zflt=1e-3)
        schedule_bus_fault!(sim; bus=8, t_on=0.03, t_off=0.04, Zflt=1e-3)
    elseif s == 6
        schedule_generator_trip!(sim; t=0.015, gen=1)
        schedule_generator_trip!(sim; t=0.035, gen=3)
    end
    return sim
end

@testset "every scenario equals an independent simulation (factorization = $fact)" for fact in
        (:per_worker, :per_scenario, :shared)
    case = test_case(:two_area)
    S = 7; nsteps = 1200
    batch = build_batch(case, S; nworkers=3, factorization=fact)
    @test nscenarios(batch) == S && nworkers(batch) == 3
    @test length(batch.factors) == (fact === :per_worker ? 3 : fact === :per_scenario ? S : 1)
    for s in 1:S; configure!(batch[s], s); end
    run_steps!(batch, nsteps)
    for s in 1:S
        ref = configure!(build_simulation(case), s)
        for _ in 1:nsteps; step!(ref); end
        @test batch[s].net.V == ref.net.V                                   # bit for bit
        @test generator_state_vector(batch[s]) == generator_state_vector(ref)
        @test bus_state_vector(batch[s]) == bus_state_vector(ref)
        @test batch[s].net.ber.Hs == ref.net.ber.Hs
    end
end

@testset "constants are shared, states are not" begin
    case = test_case(:wecc)
    batch = build_batch(case, 4; nworkers=2)
    a, b = batch[1], batch[3]                    # different workers
    for f in (:F, :T, :Req, :icf, :Gv1, :G, :Vt0, :It0, :ZbaseA)
        @test getfield(a.net, f) === getfield(b.net, f)
    end
    @test a.mc === b.mc && a.ibrc === b.ibrc && a.case === b.case
    @test a.net.ber.zc === b.net.ber.zc          # Bergeron parameters shared ...
    @test a.net.ber.Hs !== b.net.ber.Hs          # ... ring buffers not
    for f in (:brch_Ihis, :brch_Ipre, :node_Ihis, :V, :V_prev, :rhs, :ber_node_Ihis)
        @test getfield(a.net, f) !== getfield(b.net, f)
    end
    @test a.x !== b.x && a.gref !== b.gref
    @test a.net.lu !== b.net.lu                  # one factorization per worker ...
    @test batch[1].net.lu === batch[2].net.lu    # ... shared by the scenarios of that worker
    schedule_bus_fault!(a; bus=2404, t_on=0.001, t_off=0.01, Zflt=1e-3)
    run_steps!(batch, 400)
    @test batch[3].net.V == batch[4].net.V       # both undisturbed: identical
    @test maxabs(a.net.V, batch[3].net.V) > 1e-3 # the faulted one is not
end

@testset "run! records every scenario; the worker count does not matter" begin
    case = test_case(:two_area)
    ref = nothing
    for nw in (1, 2, 4)
        batch = build_batch(case, 4; nworkers=nw)
        schedule_bus_fault!(batch[2]; bus=8, t_on=0.01, t_off=0.03, Zflt=1e-3)
        schedule_generator_trip!(batch[3]; t=0.02, gen=1)
        res = run!(batch, 800 * 50e-6; record_every=10, progress=false)
        @test length(res) == 4
        @test all(r -> length(r.t) == 81 && all(isfinite, r.V), res)
        k = bus_index(case.pfd, 8); nb = nbus(case.pfd)
        @test minimum(vmag(res[2], k, nb, j) for j in eachindex(res[2].t)) < 0.1     # faulted
        @test minimum(vmag(res[1], k, nb, j) for j in eachindex(res[1].t)) > 0.9     # undisturbed
        v = vcat([copy(batch[s].net.V) for s in 1:4]...)
        ref === nothing ? (ref = v) : (@test v == ref)
    end
end
