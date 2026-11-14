# Bergeron lines (paper, Section IV).
#
# The decisive test is the first one: driven by sinusoidal terminal voltages,
# a conductor must carry the analytic steady-state current of the same
# discrete model (the two-port of eq. (8)) from the very first step.  That
# holds to machine precision only if the ring buffers were filled with the
# exact periodic solution and with the newest sample at t = 0.

using Test
using SparseArrays
using NodalEMT
using NodalEMT: PowerFlowData, build_bergeron, bergeron_coefficients, bergeron_delay_phasor, bergeron_twoport,
             bergeron_injection!, update_bergeron!, select_bergeron_lines, phase_node, damping_constants

const WS = 2pi * 60
const TS = 50e-6

"A two-bus power flow with one line, from which a single Bergeron line is built."
function two_bus_pfd(R, X, B; Vm=(1.02, 0.99), Va=(0.1, -0.05))
    return PowerFlowData(100.0, WS,
        [1, 2], ["A", "B"], [3, 1], collect(Vm), collect(Va), [230.0, 230.0],
        Int[], String[], Float64[], Float64[], Float64[],
        Int[], String[], Float64[], Float64[], Float64[],
        [1], [2], ["1"], [complex(R, X)], [B],
        Int[], Int[], String[], ComplexF64[], Float64[],
        Int[], String[], Float64[], Float64[],
        Int[], String[], ComplexF64[],
        Int[], ComplexF64[])
end

"Two-port (Y11, Y12) of a symmetric two-port with chain matrix [A B; C A]."
chain_twoport(M) = (M[1, 1] / M[1, 2], -1 / M[1, 2])
lossless_chain(zc, th) = [cos(th) im * zc * sin(th); im * sin(th) / zc cos(th)]
series_chain(r) = [1 r; 0 1]

@testset "exact initialization: a conductor starts in its discrete steady state" begin
    worst = 0.0
    for (R, X, B) in ((0.0, 1.0, 0.2), (0.01, 1.0, 0.2), (0.02, 0.35, 0.06), (0.001, 0.08, 0.012), (0.05, 2.5, 0.5))
        pfd = two_bus_pfd(R, X, B)
        Vt0 = zeros(ComplexF64, 6)
        for k in 1:2
            v = pfd.bus_Vm[k] * cis(pfd.bus_Va[k])
            Vt0[k] = v; Vt0[k+2] = v * cis(-2pi / 3); Vt0[k+4] = v * cis(2pi / 3)
        end
        ber = build_bergeron(pfd, TS, [1], Vt0, 2)
        zc = sqrt(X / B)
        @test ber.n == 3
        @test ber.zc[1] ≈ zc
        @test ber.delay[1] ≈ sqrt(X * B) / WS / TS
        @test ber.z[1] ≈ zc + R / 4
        @test ber.h[1] ≈ (zc - R / 4) / (zc + R / 4)
        @test ber.len[1] == ceil(Int, ber.delay[1])
        Y11, Y12 = bergeron_twoport(ber.z[1], ber.h[1], ber.wt[1], ber.wr[1],
                                    bergeron_delay_phasor(ber.delay[1], TS, WS))
        Is = Y11 * Vt0[ber.F[1]] + Y12 * Vt0[ber.T[1]]
        Ir = Y12 * Vt0[ber.F[1]] + Y11 * Vt0[ber.T[1]]
        V = zeros(6); node_I = zeros(6)
        for n in 1:400                          # prescribed sinusoidal terminal voltages
            t = n * TS
            V .= real.(Vt0 .* cis(WS * t))
            bergeron_injection!(node_I, ber)
            update_bergeron!(ber, V)
            worst = max(worst, abs(ber.i_from[1] - real(Is * cis(WS * t))), abs(ber.i_to[1] - real(Ir * cis(WS * t))))
        end
    end
    @test worst < 1e-12
    @info "Bergeron steady state from the first step" max_current_error = worst
end

@testset "two-port: eq. (8) of the paper and the lumped-resistance circuit" begin
    for (R, X, B) in ((0.011, 0.11, 0.1925), (0.0025, 0.025, 0.04375), (0.05, 2.5, 0.5))
        zc = sqrt(X / B); th = sqrt(X * B)          # surge impedance, electrical length at 60 Hz
        z, h, wt, wr = bergeron_coefficients(R, zc)
        @test wt ≈ (1 + h) / 2 && wr ≈ (1 - h) / 2
        for D in (cis(-th), bergeron_delay_phasor(th / WS / TS, TS, WS))
            Y11, Y12 = bergeron_twoport(z, h, wt, wr, D)
            ND = 2z * (1 + D * h) * (1 - D * h^2)
            @test Y11 ≈ (2 - D * (1 - h)^2 + 2D^2 * h^2) / ND rtol = 1e-12
            @test Y12 ≈ -D * (1 + h)^2 / ND rtol = 1e-12
        end
        # with the continuous delay, the model is exactly R/4 -- lossless half -- R/2 -- lossless half -- R/4
        Y11, Y12 = bergeron_twoport(z, h, wt, wr, cis(-th))
        M = series_chain(R / 4) * lossless_chain(zc, th / 2) * series_chain(R / 2) *
            lossless_chain(zc, th / 2) * series_chain(R / 4)
        A11, A12 = chain_twoport(M)
        @test Y11 ≈ A11 rtol = 1e-12
        @test Y12 ≈ A12 rtol = 1e-12
    end
    # the interpolated delay never amplifies: |D| <= 1
    @test all(abs(bergeron_delay_phasor(d, TS, WS)) <= 1 + 1e-15 for d in 1.01:0.37:40)
end

@testset "line selection and the entries of G" begin
    case = test_case(:wecc)
    pfd = case.pfd
    d = line_travel_times(pfd, TS)
    @test length(d) == length(pfd.line_from) && all(d .>= 0)
    for thr in (1.0, 3.0, 10.0)
        b, p = select_bergeron_lines(pfd, TS, thr)
        @test length(b) + length(p) == length(pfd.line_from)
        @test all(d[b] .> thr) && all(d[p] .<= thr)
    end
    @test isempty(select_bergeron_lines(pfd, TS, Inf)[1])
    @test_throws ErrorException select_bergeron_lines(pfd, TS, 0.5)

    # at 50 us, 276 of the 329 WECC lines become Bergeron lines (paper, Section V)
    @test length(select_bergeron_lines(pfd, TS, 1.0)[1]) == 276

    # a Bergeron line replaces the pi stamp of its line by 1/Z at both terminals
    npi = build_simulation(case; bergeron_threshold=Inf)
    nhy = build_simulation(case; bergeron_threshold=3.0)
    ber_idx, _ = select_bergeron_lines(pfd, TS, 3.0)
    @test nhy.net.ber.n == 3length(ber_idx)
    @test length(nhy.net.F) == length(npi.net.F) - 9length(ber_idx)   # 3 series + 6 charging branches per line
    dmp = damping_constants(:updated)
    nb = nbus(pfd)
    I_ = Int[]; J_ = Int[]; V_ = Float64[]
    add(i, j, v) = (push!(I_, i); push!(J_, j); push!(V_, v))
    for i in ber_idx
        R = real(pfd.line_RX[i]); X = imag(pfd.line_RX[i])
        req, = NodalEMT.rl_series_companion(R, X, pfd.ws, TS; damptrap=dmp.line_damptrap)
        reqc, = NodalEMT.c_shunt_companion(pfd.line_chg[i] / 2 / pfd.ws, TS; rs_factor=dmp.chg_rs_factor)
        f = bus_index(pfd, pfd.line_from[i]); t = bus_index(pfd, pfd.line_to[i])
        for ph in 1:3
            fn = phase_node(f, ph, nb); tn = phase_node(t, ph, nb)
            add(fn, fn, -1 / req); add(tn, tn, -1 / req); add(fn, tn, 1 / req); add(tn, fn, 1 / req)
            add(fn, fn, -1 / reqc); add(tn, tn, -1 / reqc)
        end
    end
    for k in 1:nhy.net.ber.n
        add(nhy.net.ber.F[k], nhy.net.ber.F[k], 1 / nhy.net.ber.z[k])
        add(nhy.net.ber.T[k], nhy.net.ber.T[k], 1 / nhy.net.ber.z[k])
    end
    expected = sparse(I_, J_, V_, npi.net.N, npi.net.N)
    @test maximum(abs, nhy.net.G - npi.net.G - expected) < 1e-9 * maximum(abs, npi.net.G)
    # and the Bergeron lines remove off-diagonal entries of G
    @test nnz(build_simulation(case).net.G) < nnz(npi.net.G)
end

@testset "pi+Bergeron network: undisturbed run stays at the operating point" begin
    case = test_case(:two_area)
    function drift(; kw...)
        sim = build_simulation(case; kw...)
        Vm0 = copy(sim.x_bus.vtm)
        worst = 0.0
        for _ in 1:4000
            step!(sim)
            worst = max(worst, maximum(abs, sim.x_bus.vtm .- Vm0))
        end
        @test all(isfinite, sim.net.V)
        return worst
    end
    dev_pi = drift(bergeron_threshold=Inf)
    dev_20 = drift(bergeron_threshold=20.0)      # no two-area line is that long: identical to pi
    dev_hy = drift(bergeron_threshold=1.0)       # six Bergeron lines (Fig. 3 of the paper)
    @test dev_20 == dev_pi
    @test dev_hy < 5e-2
    @info "bus voltage drift over 0.2 s" pi = dev_pi pi_bergeron = dev_hy
end

@testset "Bergeron lines through a fault and a trip" begin
    case = test_case(:two_area)
    sim = build_simulation(case)
    @test sim.net.ber.n == 18                    # six lines, three conductors each
    schedule_bus_fault!(sim; bus=8, t_on=0.05, t_off=0.11, Zflt=1e-3)
    schedule_generator_trip!(sim; t=0.3, gen=1)
    res = run!(sim, 0.6; record_every=10, progress=false)
    @test all(isfinite, res.V) && all(isfinite, res.omega)
    k = bus_index(case.pfd, 8); nb = nbus(case.pfd)
    during = maximum(vmag(res, k, nb, j) for j in eachindex(res.t) if 0.06 < res.t[j] < 0.11)
    before = maximum(vmag(res, k, nb, j) for j in eachindex(res.t) if res.t[j] < 0.05)
    @test during < 0.05 * before
    @test maximum(abs, res.omega[:, end] .- case.pfd.ws) / case.pfd.ws < 0.05
end
