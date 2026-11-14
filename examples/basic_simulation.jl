#=
basic_simulation.jl -- one simulation from PSS/E files to results.

    julia --project=. examples/basic_simulation.jl

Reads the Kundur two-area system from its .raw and .dyr files (the power flow
is solved on the way), applies a fault and a governor reference step, runs
8 s and shows how to read the results.  CSV files are written to
examples/results_basic/.
=#
using NodalEMT
using Printf

const PSSE = joinpath(@__DIR__, "..", "data", "psse")

# --- 1. The case: network, power-flow solution and dynamic models ----------
# vm_te = 0.1 s is ParaEMT's bus voltage filter for this system (the default,
# 0.02 s, is that of the WECC system).
case = load_psse(joinpath(PSSE, "twoarea.raw"), joinpath(PSSE, "twoarea.dyr"); vm_te=0.1)
pfd = case.pfd
println(case)
for i in 1:ngen(pfd)
    @printf("  generator %d at bus %d: %6.1f MW %6.1f Mvar\n", i, pfd.gen_bus[i], pfd.gen_MW[i], pfd.gen_Mvar[i])
end

# --- 2. The simulation: 50 us step, Bergeron lines, the updated model set --
sim = build_simulation(case)

# --- 3. Disturbances ----------------------------------------------------------
# Times in s; buses by their number in the case; generators by their position
# in pfd.gen_bus.  A governor reference step dsp changes the mechanical power
# by dsp/R (TGOV1 droop R = 0.04 here): 0.002 -> +0.05 pu on 900 MVA.
schedule_bus_fault!(sim; bus=8, t_on=1.0, t_off=1.1, Zflt=0.01)       # three-phase, 100 ms
schedule_reference_step!(sim; t=4.0, gen=1, kind=:governor, dsp=0.002)
# also available:  schedule_bus_fault!(...; phases=[1])  (single-phase),
#                  schedule_generator_trip!(sim; t, gen),
#                  schedule_reference_step!(sim; t, gen, kind=:exciter, dsp)

# --- 4. Run -------------------------------------------------------------------
res = run!(sim, 8.0; record_every=20)          # store every 20 steps = 1 ms

# --- 5. Results: column j of every matrix belongs to time res.t[j] ------------
# res.vtm  bus voltage magnitudes [pu]     res.freq  bus frequencies [pu]
# res.omega generator speeds [rad/s]       res.pe    electrical powers [pu]
# res.V    instantaneous node voltages (phase a of all buses, then b, then c)
k = bus_index(pfd, 8)
during = 1.02 .< res.t .< 1.1
@printf("\nbus 8 voltage: %.3f pu before the fault, %.3f pu during it, %.3f pu at %.1f s\n",
        res.vtm[k, 1], minimum(res.vtm[k, during]), res.vtm[k, end], res.t[end])
after = res.t .> 1.2                           # after the fault (the PLLs swing during it)
f = res.freq[:, after] .* 60                   # Hz
@printf("lowest bus frequency after the fault: %.3f Hz\n", minimum(f))
j = findlast(<(4.0), res.t)                    # the last instant before the governor step
for i in 1:ngen(pfd)
    @printf("  generator %d: %.3f Hz at the end; electrical power %.3f pu at %.1f s, %.3f pu at %.1f s (machine base)\n",
            i, res.omega[i, end] / 2pi, res.pe[i, j], res.t[j], res.pe[i, end], res.t[end])
end

# --- 6. Save everything as CSV ------------------------------------------------
out = save_results(res, sim, joinpath(@__DIR__, "results_basic"))
println("\nCSV files in ", out, ": ", join(readdir(out), ", "))
