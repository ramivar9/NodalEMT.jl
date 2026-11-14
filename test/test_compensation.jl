# Compensation (paper, Section III): keeping G and correcting the solution
# must give the result of the modified network.

using Test
using LinearAlgebra
using NodalEMT
using NodalEMT: fault_change, trip_change, build_compensation, set_active!, apply_compensation!, phase_node

@testset "one solution: compensation equals the modified network" begin
    case = test_case(:wecc)
    sim = build_simulation(case)
    net = sim.net; nb = net.nbus
    fnodes = [phase_node(bus_index(case.pfd, 2404), ph, nb) for ph in 1:3]
    gnodes = [phase_node(sim.mc.bus[5], ph, nb) for ph in 1:3]
    changes = [fault_change(fnodes, 1e-3), trip_change(gnodes, sim.mc.Gequiv[5])]
    comp = build_compensation(net, changes)
    b = randn(net.N)
    for on in ((true, false), (false, true), (true, true), (false, false))
        set_active!(comp, 1, on[1]); set_active!(comp, 2, on[2])
        G = copy(net.G)
        for (k, ch) in enumerate(changes)
            on[k] || continue
            for (a, i) in enumerate(ch.nodes), (c, j) in enumerate(ch.nodes)
                G[i, j] += ch.dY[a, c]
            end
        end
        Vdirect = G \ b
        Vcomp = apply_compensation!(net.lu \ b, comp)
        @test maxabs(Vcomp, Vdirect) / maximum(abs, Vdirect) < 1e-9
    end
end

@testset "a whole run equals refactorization (two-area)" begin
    fault(s) = schedule_bus_fault!(s; bus=8, t_on=0.05, t_off=0.11, Zflt=1e-3)
    trip(s) = schedule_generator_trip!(s; t=0.05, gen=2)
    both(s) = (schedule_bus_fault!(s; bus=8, t_on=0.03, t_off=0.09, Zflt=1e-2);
               schedule_generator_trip!(s; t=0.06, gen=1))
    twotrips(s) = (schedule_generator_trip!(s; t=0.04, gen=1); schedule_generator_trip!(s; t=0.09, gen=3))
    onephase(s) = schedule_bus_fault!(s; bus=7, t_on=0.02, t_off=0.05, Zflt=5e-2, phases=[1])
    reclose(s) = (schedule_bus_fault!(s; bus=8, t_on=0.02, t_off=0.04, Zflt=1e-3);
                  schedule_bus_fault!(s; bus=8, t_on=0.08, t_off=0.10, Zflt=1e-3))
    case = test_case(:two_area)
    worst = 0.0
    for thr in (Inf, 1.0), cfg in (identity, fault, trip, both, twotrips, onephase, reclose)
        a = build_simulation(case; network_events=:refactorize, bergeron_threshold=thr)
        b = build_simulation(case; network_events=:compensate, bergeron_threshold=thr)
        cfg(a); cfg(b)
        for _ in 1:3000
            worst = max(worst, maxabs(step!(a), step!(b)),
                        maxabs(generator_state_vector(a), generator_state_vector(b)),
                        maxabs(bus_state_vector(a), bus_state_vector(b)))
        end
        @test all(isfinite, b.net.V)
    end
    @info "compensation vs refactorization (two-area)" worst_difference = worst
    @test worst < 1e-9
end

@testset "a whole run equals refactorization (WECC)" begin
    case = test_case(:wecc)
    worst = 0.0
    for thr in (Inf, 1.0)
        a = build_simulation(case; network_events=:refactorize, bergeron_threshold=thr)
        b = build_simulation(case; network_events=:compensate, bergeron_threshold=thr)
        for s in (a, b)
            schedule_bus_fault!(s; bus=2404, t_on=0.02, t_off=0.06, Zflt=1e-2)
            schedule_generator_trip!(s; t=0.04, gen=5)
        end
        m = dynamics_mask(a.mc.n, length(generator_state_vector(a)))
        for _ in 1:2000
            worst = max(worst, maxabs(step!(a), step!(b)),
                        maxabs(generator_state_vector(a)[m], generator_state_vector(b)[m]))
        end
    end
    @info "compensation vs refactorization (WECC)" worst_difference = worst
    @test worst < 1e-9
end

@testset "G and its factors are never changed" begin
    case = test_case(:two_area)
    sim = build_simulation(case)
    schedule_bus_fault!(sim; bus=8, t_on=0.005, t_off=0.02, Zflt=1e-3)
    schedule_generator_trip!(sim; t=0.01, gen=1)
    G_before = copy(sim.net.G)
    F_before = sim.net.lu
    for _ in 1:600; step!(sim); end
    @test sim.net.G == G_before
    @test sim.net.lu === F_before
    @test sim.comp !== nothing && sim.gen_tripped[1]
end

@testset "a disturbance scheduled during the run" begin
    case = test_case(:two_area)
    a = build_simulation(case; network_events=:compensate)
    b = build_simulation(case; network_events=:refactorize)
    for _ in 1:200; step!(a); step!(b); end
    for s in (a, b)
        schedule_bus_fault!(s; bus=8, t_on=0.02, t_off=0.05, Zflt=1e-3)
    end
    @test a.comp === nothing                   # discarded, rebuilt on the next step
    worst = 0.0
    for _ in 1:2000
        worst = max(worst, maxabs(step!(a), step!(b)))
    end
    @test a.comp !== nothing
    @test worst < 1e-9
end

@testset "invalid settings are rejected" begin
    case = test_case(:two_area)
    @test_throws ErrorException build_simulation(case; network_events=:nonsense)
    @test_throws ErrorException build_simulation(case; variant=:nonsense)
    @test_throws ErrorException build_batch(case, 2; network_events=:refactorize)
end
