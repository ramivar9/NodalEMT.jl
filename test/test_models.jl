# The updated model set (paper, Section II) against the ParaEMT models,
# each change checked where it acts.

using Test
using SparseArrays
using NodalEMT
using NodalEMT: damping_constants, rl_series_companion, sexs_step, CaseData, PowerFlowData

@testset "numerical damping" begin
    par = damping_constants(:paraemt)
    upd = damping_constants(:updated)
    @test par.line_damptrap == 1.0 && upd.line_damptrap == 10.0       # ten times lighter for lines
    @test upd.chg_rs_factor == par.chg_rs_factor / 10
    @test par.load_damping && !upd.load_damping                        # none for loads ...
    @test upd.shunt_rs_factor == 0.0 && upd.load_c_rs_factor == 0.0    # ... and shunts
    @test_throws ErrorException damping_constants(:nonsense)

    ws = 2pi * 60; ts = 50e-6; R = 0.01; X = 0.2
    r1, i1, g1 = rl_series_companion(R, X, ws, ts; damptrap=1.0)
    r10, i10, g10 = rl_series_companion(R, X, ws, ts; damptrap=10.0)
    r0, i0, g0 = rl_series_companion(R, X, ws, ts; damptrap=0.0)
    @test r1 < r10 < r0                         # lighter damping: closer to the plain trapezoidal model
    a = ts / (2X / ws)
    @test i0 ≈ (1 - R * a) / (1 + R * a) && r0 ≈ R + 1 / a && g0 ≈ a / (1 + R * a)

    case = test_case(:wecc)
    p = build_simulation(case; variant=:paraemt, bergeron_threshold=Inf)
    u = build_simulation(case; variant=:updated, bergeron_threshold=Inf)
    @test p.net.F == u.net.F && p.net.Req != u.net.Req && nnz(p.net.G) == nnz(u.net.G)
end

@testset "SEXS: the stabilizer derivative enters the lead term" begin
    args = (0.3, 1.5, 1.02, 1.02, 1.0, 0.0, 0.01)           # v1, EFD, vref, vref_prev, vt, dvt, vs
    p = (0.1, 1.0, 100.0, 0.05, -5.0, 5.0, 50e-6)           # TA, TB, K, TE, Emin, Emax, ts
    v_paraemt, _ = sexs_step(args..., 0.0, p...)
    v_updated, _ = sexs_step(args..., 200.0, p...)          # dvs = 200
    @test v_updated - v_paraemt ≈ p[1] / p[2] * 200.0 * p[7]
end

@testset "motoring machines start in equilibrium" begin
    case = test_case(:two_area)
    pfd = case.pfd
    mw = copy(pfd.gen_MW); mw[2] = -abs(mw[2])              # generator 2 becomes a motor
    args = [getfield(pfd, f) for f in fieldnames(PowerFlowData)]
    args[findfirst(==(:gen_MW), fieldnames(PowerFlowData))] = mw
    motor = CaseData(case.name, PowerFlowData(args...), case.dyd)
    up = build_simulation(motor; variant=:updated)
    pa = @test_logs (:warn, r"initialization residual") build_simulation(motor; variant=:paraemt)
    pe(s, i) = s.x.mac.psid[i] * s.x.mac.iq[i] - s.x.mac.psiq[i] * s.x.mac.id[i]
    @test pe(up, 2) < 0                          # absorbs power, as the power flow says
    @test pe(pa, 2) > 0                          # ParaEMT starts it as a generator
    @test up.x.mac.id[1] == pa.x.mac.id[1]       # generating machines are unchanged
end

@testset "REECB and REPCA start in equilibrium" begin
    case = test_case(:wecc)
    up = build_simulation(case; variant=:updated)
    pa = build_simulation(case; variant=:paraemt)
    # REECB: ParaEMT starts the reactive-power PI integrators at zero
    @test all(iszero, pa.x_ibr.reecb.s2) && all(iszero, pa.x_ibr.reecb.s3)
    @test any(!iszero, up.x_ibr.reecb.s2)
    @test up.x_ibr.reecb.v2iPI ≈ up.x_ibr.reecb.s2 .- up.x_ibr.reecb.s0
    # REPCA without frequency regulation (FFlag = 0; every bundled IBR has FFlag = 1, so it is forced here)
    flat = test_case(:wecc)
    fill!(flat.dyd.repca.FFlag, 0)
    fup = build_simulation(flat; variant=:updated)
    fpa = build_simulation(flat; variant=:paraemt)
    @test fup.ibrc.Pref_out == fup.x_ibr.reecb.Pref     # dispatch kept
    @test all(iszero, fpa.ibrc.Pref_out)                # reset to zero by ParaEMT
    P0 = copy(fup.x_ibr.reecb.s5)
    for _ in 1:2000; step!(fup); step!(fpa); end
    @test maximum(abs, fpa.x_ibr.reecb.s5 .- P0) > 10 * maximum(abs, fup.x_ibr.reecb.s5 .- P0)
end
