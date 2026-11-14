#=
own_system.jl -- a template for your own PSS/E case.

    julia --project=. examples/own_system.jl my_case.raw my_case.dyr
    julia --project=. examples/own_system.jl                  # the bundled WECC 240-bus system

1. reads the .raw and .dyr files and solves the power flow (with a more
   robust solver if Newton-Raphson fails);
2. reports the power flow and checks that the EMT network reproduces it;
3. saves the case as JSON, which load_case reads in about a second, without
   solving the power flow again;
4. runs an undisturbed second to check that the simulation starts in steady
   state, then a fault at the best-connected bus of the highest voltage level.
Outputs go to examples/results_own_system/.
=#
using NodalEMT
using Printf

const PSSE = joinpath(@__DIR__, "..", "data", "psse")
raw = length(ARGS) >= 1 ? ARGS[1] : joinpath(PSSE, "WECC240.raw")
dyr = length(ARGS) >= 2 ? ARGS[2] : joinpath(PSSE, "WECC240.dyr")
out = mkpath(joinpath(@__DIR__, "results_own_system"))
json = joinpath(out, splitext(basename(raw))[1] * ".json")

# --- 1. read and solve -----------------------------------------------------------
# A .dyr model that NodalEMT does not support, or a device it cannot simulate
# (HVDC, FACTS, phase shifters), stops here with a message naming it.
case = try
    load_psse(raw, dyr; json=json)
catch err
    occursin("did not converge", sprint(showerror, err)) || rethrow()
    @warn "Newton-Raphson did not converge; trying the trust-region solver"
    load_psse(raw, dyr; json=json, solver=:trust_region)
end
# other keywords: reactive_power_limits=true, solve_power_flow=false (keep the
# solution stored in the .raw file), pll_ke, pll_te, vm_te (bus measurements)
pfd = case.pfd
println(case)

# --- 2. the power flow -------------------------------------------------------------
@printf("generation %.0f MW (%d machines, %d IBRs), load %.0f MW\n",
        sum(pfd.gen_MW) + sum(pfd.ibr_MW; init=0.0), ngen(pfd), nibr(pfd), sum(pfd.load_MW))
@printf("bus voltages %.3f ... %.3f pu\n", extrema(pfd.bus_Vm)...)
@printf("power-flow balance of the EMT network: %.1e pu (PowerFlows.jl solves in single precision)\n",
        maximum(abs, NodalEMT.power_flow_mismatch(pfd)))
@printf("%d of %d lines are Bergeron lines at 50 us\n",
        count(>(1), line_travel_times(pfd, 50e-6)), length(pfd.line_from))

# --- 3. the saved case reads back without the power-flow packages ------------------
t = @elapsed again = load_case(json)
@printf("saved to %s; load_case reads it in %.2f s\n", json, t)

# --- 4. steady state, then a fault ----------------------------------------------------
sim = build_simulation(again)
v0 = copy(sim.x_bus.vtm)
for _ in 1:round(Int, 1.0 / 50e-6)
    step!(sim)
end
@printf("after 1 s without disturbance the bus voltages moved by at most %.1e pu\n", maximum(abs, sim.x_bus.vtm .- v0))

# the fault: at the bus of the highest voltage level with the most lines
top = maximum(pfd.bus_basekV)
nlines(b) = count(==(b), pfd.line_from) + count(==(b), pfd.line_to)
fbus = argmax(nlines, [b for (b, kv) in zip(pfd.bus_num, pfd.bus_basekV) if kv == top])
sim = build_simulation(again)
schedule_bus_fault!(sim; bus=fbus, t_on=0.1, t_off=0.16, Zflt=1e-3)
res = run!(sim, 1.0; record_every=20, progress=false)
k = bus_index(pfd, fbus); nb = nbus(pfd)
vmag(j) = sqrt(2 / 3 * (res.V[k, j]^2 + res.V[k+nb, j]^2 + res.V[k+2nb, j]^2))   # instantaneous magnitude
@printf("fault at bus %d (%.0f kV, %d lines): voltage %.3f pu during the fault, %.3f pu at 1 s\n",
        fbus, top, nlines(fbus), minimum(vmag(j) for j in eachindex(res.t) if 0.11 < res.t[j] < 0.16), vmag(length(res.t)))
save_results(res, sim, out)
println("results in ", out)
