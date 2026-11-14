#=
contingency_screening.jl -- a batch of contingencies on all threads, ranked.

    julia --project=. --threads=8 examples/contingency_screening.jl

Builds one scenario per 500-kV bus of the WECC 240-bus system: a three-phase
fault through 0.001 pu at that bus at 0.1 s, cleared after 60 ms, followed by
the trip of the largest unit at that bus, if any.  The batch runs on all
threads (one LU factorization per worker; the faults and trips are applied by
compensation, so G never changes), and every scenario is ranked by the lowest
bus voltage 0.2 s after clearing; the lowest speed of the machines in service
is listed too.  The ranking is written to examples/results_screening/ranking.csv.
=#
using NodalEMT
using Printf
using Base.Threads

const PSSE = joinpath(@__DIR__, "..", "data", "psse")
const T_ON, T_OFF, T_END = 0.1, 0.16, 1.0

case = load_psse(joinpath(PSSE, "WECC240.raw"), joinpath(PSSE, "WECC240.dyr"))
pfd = case.pfd
buses = [b for (b, kv) in zip(pfd.bus_num, pfd.bus_basekV) if kv >= 500]
S = length(buses)
println(case, "\n$S scenarios (one per 500-kV bus) on $(nthreads()) threads")

# --- the scenarios: same system, different disturbances ----------------------
batch = build_batch(case, S)                   # nworkers = nthreads(), factorization = :per_worker
tripped = zeros(Int, S)                        # the unit tripped in each scenario (0: none)
for (s, b) in enumerate(buses)
    schedule_bus_fault!(batch[s]; bus=b, t_on=T_ON, t_off=T_OFF, Zflt=1e-3)
    units = findall(==(b), pfd.gen_bus)
    if !isempty(units)
        g = units[argmax(pfd.gen_MW[units])]
        schedule_generator_trip!(batch[s]; t=T_OFF, gen=g)
        tripped[s] = g
    end
end

# --- run all of them -----------------------------------------------------------
t = @elapsed results = run!(batch, T_END; record_every=40)   # one Results per scenario, every 2 ms

# --- one row per scenario --------------------------------------------------------
rows = map(1:S) do s
    r = results[s]
    late = r.t .>= T_OFF + 0.2                                 # 0.2 s after clearing
    vmin = minimum(r.vtm[:, late])                             # lowest bus voltage [pu]
    online = setdiff(1:ngen(pfd), tripped[s])                  # machine speeds, not the bus PLLs,
    nadir = minimum(r.omega[online, :]) / 2pi                  # which swing at the fault [Hz]
    (bus=buses[s], trip=tripped[s] == 0 ? "-" : "$(pfd.gen_bus[tripped[s]]) $(pfd.gen_id[tripped[s]])",
     trip_MW=tripped[s] == 0 ? 0.0 : pfd.gen_MW[tripped[s]], vmin=vmin, nadir=nadir)
end
sort!(rows, by=r -> r.vmin)

@printf("\nsimulated %d x %.1f s in %.1f s; the 10 most severe by post-fault voltage:\n", S, T_END, t)
println("  fault bus  tripped unit   MW   lowest V 0.2 s after clearing   lowest machine speed")
for r in rows[1:min(10, end)]
    @printf("  %9d  %12s %6.0f   %10.3f pu                    %8.3f Hz\n", r.bus, r.trip, r.trip_MW, r.vmin, r.nadir)
end

out = joinpath(mkpath(joinpath(@__DIR__, "results_screening")), "ranking.csv")
open(out, "w") do io
    println(io, "fault_bus,tripped_unit,tripped_MW,vmin_pu,min_speed_Hz")
    for r in rows
        @printf(io, "%d,%s,%.1f,%.5f,%.5f\n", r.bus, r.trip, r.trip_MW, r.vmin, r.nadir)
    end
end
println("ranking written to ", out)
