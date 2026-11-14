#=
make_figures.jl -- the figures and numbers of docs/verification.md.

    julia --project=docs -e "using Pkg; Pkg.develop(path=\".\"); Pkg.instantiate()"   # once
    julia --project=docs --threads=16 docs/make_figures.jl                            # all
    julia --project=docs docs/make_figures.jl --skip-timing                           # without fig. 6

Writes docs/figures/*.png and prints the numbers quoted in the text.  The
timing figure (batch speedup) takes about ten minutes on 16 threads; its
numbers depend on the machine.
=#
using NodalEMT
using NodalEMT: build_bergeron, bergeron_delay_phasor, bergeron_twoport, bergeron_injection!, update_bergeron!,
                power_flow_mismatch
using Plots
using Printf
using Base.Threads
using LinearAlgebra

BLAS.set_num_threads(1)
const PSSE = joinpath(@__DIR__, "..", "data", "psse")
const FIG = mkpath(joinpath(@__DIR__, "figures"))
const TS = 50e-6

# ---- style: three validated categorical slots, light surface, recessive grid ----
const C1, C2, C3 = colorant"#2a78d6", colorant"#eb6834", colorant"#1baf7a"    # blue, orange, aqua
const INK, INK2 = colorant"#0b0b0b", colorant"#52514e"
const SURFACE, GRID, REF = colorant"#fcfcfb", colorant"#d9d8d3", colorant"#8a8984"
gr()
default(; size=(760, 360), dpi=150, fontfamily="Helvetica", background_color=SURFACE,
        foreground_color_axis=INK2, foreground_color_border=INK2, foreground_color_text=INK2,
        tickfontcolor=INK2, guidefontcolor=INK2, titlefontcolor=INK, titlefontsize=10, titlelocation=:left,
        guidefontsize=9, tickfontsize=8, legendfontsize=8, legend_background_color=SURFACE,
        legend_foreground_color=SURFACE, gridcolor=GRID, gridalpha=1.0, gridlinewidth=0.5,
        framestyle=:axes, linewidth=2, left_margin=4Plots.mm, bottom_margin=4Plots.mm)
save(p, name) = (savefig(p, joinpath(FIG, name)); println("  wrote docs/figures/", name))

two_area(; kw...) = load_psse(joinpath(PSSE, "twoarea.raw"), joinpath(PSSE, "twoarea.dyr"); vm_te=0.1, kw...)
wecc(; kw...) = load_psse(joinpath(PSSE, "WECC240.raw"), joinpath(PSSE, "WECC240.dyr"); kw...)

# ==========================================================================
println("Table 1: the power flow of PowerFlows.jl against PSS/E's stored solution")
# ==========================================================================
TA, W = two_area(), wecc()
for (name, solved, stored) in (("two-area", TA, two_area(solve_power_flow=false)),
                               ("WECC 240", W, wecc(solve_power_flow=false)))
    p, q = solved.pfd, stored.pfd
    @printf("  %-9s max|dVm| %.1e pu  max|dVa| %.1e rad  max|dP| %.1e MW  balance %.1e pu (PSS/E's: %.1e pu)\n",
            name, maximum(abs, p.bus_Vm .- q.bus_Vm), maximum(abs, p.bus_Va .- q.bus_Va),
            maximum(abs, vcat(p.gen_MW .- q.gen_MW, p.ibr_MW .- q.ibr_MW)),
            maximum(abs, power_flow_mismatch(p)), maximum(abs, power_flow_mismatch(q)))
end

# ==========================================================================
println("Figure 1: Bergeron line 7-8 driven at its power-flow voltages, three initializations")
# ==========================================================================
let case = TA, pfd = TA.pfd
    ws = pfd.ws; nb = nbus(pfd)
    Vt0 = build_simulation(case).net.Vt0
    i = findfirst(k -> pfd.line_from[k] == 7 && pfd.line_to[k] == 8, eachindex(pfd.line_from))
    "Bergeron conductors of line i with the ring buffers filled for the delay phasor D and the newest sample at t = shift."
    function rings(D_of, shift)
        b = build_bergeron(pfd, TS, [i], Vt0, nb)
        for k in 1:b.n
            L = b.len[k]; o = b.off[k]
            Y11, Y12 = bergeron_twoport(b.z[k], b.h[k], b.wt[k], b.wr[k], D_of(b.delay[k]))
            Vs = Vt0[b.F[k]]; Vr = Vt0[b.T[k]]
            Hs = Vs / b.z[k] + b.h[k] * (Y11 * Vs + Y12 * Vr); Hr = Vr / b.z[k] + b.h[k] * (Y12 * Vs + Y11 * Vr)
            for j in 0:L-1
                t = shift - (L - 1 - j) * TS
                b.Hs[o+j+1] = real(Hs * cis(ws * t)); b.Hr[o+j+1] = real(Hr * cis(ws * t))
            end
        end
        return b
    end
    discrete(d) = bergeron_delay_phasor(d, TS, ws)
    continuous(d) = cis(-ws * d * TS)
    variants = [("exact discrete initialization", rings(discrete, 0.0)),
                ("continuous delay e^(-jωτ)", rings(continuous, 0.0)),
                ("newest sample at -Δt", rings(discrete, -TS))]
    ref = build_bergeron(pfd, TS, [i], Vt0, nb)
    Y11, Y12 = bergeron_twoport(ref.z[1], ref.h[1], ref.wt[1], ref.wr[1], discrete(ref.delay[1]))
    Is = [Y11 * Vt0[ref.F[k]] + Y12 * Vt0[ref.T[k]] for k in 1:ref.n]       # steady state of the discrete line
    n = 2000; t = (1:n) .* TS
    p = plot(; yscale=:log10, xlabel="time (ms)", ylabel="current error / current amplitude", legend=:right,
             ylims=(1e-16, 1e-1), yticks=10.0 .^ (-16:2:-2),
             title="Fig. 1  Bergeron line 7-8 (τ/Δt = $(round(ref.delay[1], digits=2))) at its steady-state terminal voltages")
    node_I = zeros(3nb); V = zeros(3nb)
    for ((label, b), c) in zip(variants, (C1, C2, C3))
        err = zeros(n)
        for s in 1:n
            V .= real.(Vt0 .* cis(ws * t[s]))
            bergeron_injection!(node_I, b); update_bergeron!(b, V)
            err[s] = maximum(abs(b.i_from[k] - real(Is[k] * cis(ws * t[s]))) / abs(Is[k]) for k in 1:b.n)
        end
        plot!(p, t .* 1e3, max.(err, 1e-16); label=label, color=c)
        @printf("  %-31s largest error %.1e of the amplitude (first 100 ms)\n", label, maximum(err))
    end
    save(p, "fig1_bergeron_initialization.png")
end

# ==========================================================================
println("Figure 2: three-phase fault at bus 8 of the two-area system (the PSCAD/EMTDC case of the paper)")
# ==========================================================================
let pfd = TA.pfd
    k = bus_index(pfd, 8); nb = nbus(pfd)
    vpk = pfd.bus_basekV[k] * sqrt(2 / 3); ipk = 2 * pfd.basemva / (3 * vpk)
    function run_fault(t_on, threshold)
        sim = build_simulation(TA; bergeron_threshold=threshold)
        zf = 1.0 / sim.net.ZbaseA[k]                                         # 1 ohm
        schedule_bus_fault!(sim; bus=8, t_on=t_on, t_off=t_on + 0.1, Zflt=zf)
        ti = Float64[]; ia = Float64[]; tv = Float64[]; vb = Float64[]
        while current_time(sim) < t_on + 0.1 + 0.0105
            Vn = step!(sim); t = current_time(sim)
            t_on <= t <= t_on + 0.006 && (push!(ti, (t - t_on) * 1e3); push!(ia, Vn[k] / zf * ipk))
            t_on + 0.1 <= t && (push!(tv, (t - t_on - 0.1) * 1e3); push!(vb, Vn[k+nb] * vpk))
        end
        return ti, ia, tv, vb
    end
    # inception at the point on wave of the PSCAD run, as in examples/two_area_fault.jl
    bt, bi, btv, bv = run_fault(5.00185, 1.0)
    pt, pi_, ptv, pv = run_fault(5.00175, Inf)
    @printf("  peak fault current: %.2f kA with Bergeron lines, %.2f kA with pi lines (PSCAD: 8.06 kA)\n",
            maximum(abs, bi), maximum(abs, pi_))
    a = plot(bt, bi; label="π+Bergeron lines", color=C1, xlabel="time after fault inception (ms)",
             ylabel="phase-a fault current (kA)", legend=:bottomright, title="(a) fault current")
    plot!(a, pt, pi_; label="single π section per line", color=C2)
    hline!(a, [8.06]; color=REF, linestyle=:dash, linewidth=1, label="PSCAD/EMTDC peak, 8.06 kA")
    b = plot(btv, bv; label="π+Bergeron lines", color=C1, xlabel="time after fault clearing (ms)",
             ylabel="bus-8 phase-b voltage (kV)", legend=:bottomright, title="(b) voltage after clearing")
    plot!(b, ptv, pv; label="single π section per line", color=C2)
    p = plot(a, b; layout=(1, 2), size=(900, 360),
             plot_title="Fig. 2  Three-phase fault at bus 8 of the two-area system (1 Ω, 100 ms)",
             plot_titlefontsize=10, plot_titlefontcolor=INK)
    save(p, "fig2_two_area_fault.png")
end

# ==========================================================================
println("Figure 3: compensation against refactorization, WECC 240-bus system")
# ==========================================================================
let
    sims = [build_simulation(W; network_events=e) for e in (:refactorize, :compensate)]
    for s in sims
        schedule_bus_fault!(s; bus=2404, t_on=0.02, t_off=0.08, Zflt=1e-3)
        schedule_generator_trip!(s; t=0.05, gen=5)
    end
    n = round(Int, 0.2 / TS); d = zeros(n)
    for j in 1:n
        Va = step!(sims[1]); Vb = step!(sims[2])
        d[j] = maximum(abs, Va .- Vb)
    end
    @printf("  largest node-voltage difference over 0.2 s: %.1e pu\n", maximum(d))
    p = plot((1:n) .* TS .* 1e3, max.(d, 1e-17); color=C1, yscale=:log10, ylims=(1e-17, 1e-8),
             yticks=10.0 .^ (-16:2:-8), legend=false, xlabel="time (ms)",
             ylabel="largest node-voltage difference (pu)",
             title="Fig. 3  Fault and trip applied by compensation vs. by refactorizing G (WECC 240-bus)")
    for (t, label) in ((20, "fault"), (50, "trip"), (80, "clearing"))
        vline!(p, [t]; color=REF, linestyle=:dash, linewidth=1)
        annotate!(p, t + 1, 2e-9, text(label, 8, INK2, :left))
    end
    save(p, "fig3_compensation.png")
end

# ==========================================================================
println("Figure 4: undisturbed start from the power flow, WECC 240-bus system")
# ==========================================================================
let
    a = plot(; ylabel="largest deviation (pu)", legend=:topright, title="(a) bus voltage magnitudes")
    b = plot(; xlabel="time (s)", ylabel="largest deviation (Hz)", legend=:topright, title="(b) machine speeds")
    ws = W.pfd.ws
    for (label, kw, c) in (("π lines", (bergeron_threshold=Inf,), C2), ("π+Bergeron lines", (), C1))
        sim = build_simulation(W; kw...)
        v0 = copy(sim.x_bus.vtm)
        n = round(Int, 10.0 / TS); tt = Float64[]; dv = Float64[]; dw = Float64[]
        for j in 1:n
            step!(sim)
            if j % 100 == 0
                push!(tt, j * TS); push!(dv, maximum(abs, sim.x_bus.vtm .- v0))
                push!(dw, maximum(abs, sim.x.mac.omega ./ ws .- 1) * 60)
            end
        end
        plot!(a, tt, dv; label=label, color=c); plot!(b, tt, dw; label=label, color=c)
        late = tt .> 5
        @printf("  %-17s voltage: largest %.1e pu, after 5 s %.1e pu;  speed: largest %.1e Hz, at 10 s %.1e Hz\n",
                label, maximum(dv), maximum(dv[late]), maximum(dw), dw[end])
    end
    p = plot(a, b; layout=(2, 1), size=(760, 600), link=:x,
             plot_title="Fig. 4  Undisturbed run from the power-flow solution (WECC 240-bus)",
             plot_titlefontsize=10, plot_titlefontcolor=INK)
    save(p, "fig4_steady_start.png")
end

# ==========================================================================
println("Figure 5: example output, fault and governor step on the two-area system")
# ==========================================================================
let pfd = TA.pfd
    sim = build_simulation(TA)
    schedule_bus_fault!(sim; bus=8, t_on=1.0, t_off=1.1, Zflt=0.01)
    schedule_reference_step!(sim; t=4.0, gen=1, kind=:governor, dsp=0.002)
    r = run!(sim, 8.0; record_every=20, progress=false)
    a = plot(; ylabel="voltage magnitude (pu)", legend=:bottomright, title="(a) bus voltages (filtered measurement)")
    for (b, c) in zip((7, 8, 9), (C1, C2, C3))
        plot!(a, r.t, r.vtm[bus_index(pfd, b), :]; label="bus $b", color=c)
    end
    s = plot(; xlabel="time (s)", ylabel="speed (Hz)", legend=:topright, title="(b) generator speeds")
    for (g, c) in zip((1, 3), (C1, C2))
        plot!(s, r.t, r.omega[g, :] ./ 2pi; label="generator $g (area $(g == 1 ? 1 : 2))", color=c)
    end
    p = plot(a, s; layout=(2, 1), size=(760, 600), link=:x,
             plot_title="Fig. 5  examples/basic_simulation.jl: fault at bus 8 (1.0-1.1 s), governor step of G1 (4 s)",
             plot_titlefontsize=10, plot_titlefontcolor=INK)
    save(p, "fig5_example_two_area.png")
end

# ==========================================================================
"--skip-timing" in ARGS && exit()
println("Figure 6: batch speedup, WECC 240-bus system, $(nthreads()) threads")
# ==========================================================================
let S = 64, nsteps = round(Int, 0.6 / TS)
    ongrid(t) = round(Int, t / TS) * TS
    function contingency!(sim, s)                                     # as examples/wecc_batch.jl
        pfd = sim.case.pfd
        schedule_bus_fault!(sim; bus=pfd.bus_num[mod1(3s + 1, nbus(pfd))], t_on=ongrid(0.20), t_off=ongrid(0.26),
                            Zflt=0.01, phases=iseven(s) ? [1] : [1, 2, 3])
        schedule_reference_step!(sim; t=ongrid(0.40), gen=mod1(s, ngen(pfd)), kind=:governor, dsp=-0.01)
        s % 3 == 0 && schedule_generator_trip!(sim; t=ongrid(0.50), gen=mod1(s + 1, ngen(pfd)))
    end
    function timed(n, workers, fact, steps)
        b = build_batch(W, n; nworkers=workers, factorization=fact)
        foreach(s -> contingency!(b[s], s), 1:n)
        GC.gc()
        return @elapsed run_steps!(b, steps)
    end
    foreach(f -> timed(2, min(2, nthreads()), f, 500), (:per_worker, :per_scenario, :shared))    # compile
    t1 = timed(S, 1, :per_worker, nsteps)
    @printf("  1 worker: %.1f s for %d scenarios of 0.6 s\n", t1, S)
    P = filter(<=(nthreads()), [1, 2, 4, 8, 16])
    p = plot(P, P; color=REF, linestyle=:dash, linewidth=1, label="ideal", xscale=:log2, yscale=:log2,
             xticks=(P, string.(P)), yticks=(P, string.(P)), xlabel="worker threads", ylabel="speedup over one worker",
             legend=:topleft, title="Fig. 6  $S WECC contingencies of 0.6 s, three factorization policies")
    for (f, c, label) in ((:per_worker, C1, "one LU per worker"), (:per_scenario, C2, "one LU per scenario"),
                          (:shared, C3, "one shared LU"))
        sp = [t1 / timed(S, w, f, nsteps) for w in P]
        plot!(p, P, sp; label=label, color=c, marker=:circle, markersize=4, markerstrokewidth=0)
        @printf("  %-20s speedup %s\n", label, join((@sprintf("%.1f", x) for x in sp), " / "))
    end
    save(p, "fig6_batch_speedup.png")
end
