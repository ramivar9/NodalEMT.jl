#=
wecc_batch.jl -- a batch of WECC 240-bus contingencies on all threads, with
the three factorization policies of the paper (Table III).

    julia --project=. --threads=16 examples/wecc_batch.jl [S]

S scenarios (default 32; the paper uses 500) of 0.6 s at 50 us, as in the
paper: scenario s applies a 0.01-pu fault at one of 81 buses at 0.2 s for
60 ms (single-phase in even and three-phase in odd scenarios), a governor
reference step at 0.4 s, and every third scenario trips a generator at
0.5 s.  The batch is timed on one worker and on all threads with

    :per_worker    one LU factorization of G per worker (NodalEMT.jl)
    :per_scenario  one LU factorization per scenario
    :shared        one LU factorization for the whole batch

The time of the time loop is measured; building the batch is not.  Every
policy gives the same result, which is checked at the end.

The system is read from data/psse/.
=#
using NodalEMT
using LinearAlgebra
using Printf
using Base.Threads

BLAS.set_num_threads(1)                     # one thread per worker, none inside BLAS

const PSSE = joinpath(@__DIR__, "..", "data", "psse")
const S = isempty(ARGS) ? 32 : parse(Int, ARGS[1])
const TS = 50e-6
const NSTEPS = round(Int, 0.6 / TS)

ongrid(t) = round(Int, t / TS) * TS

"The disturbances of scenario `s` (paper, Section V-D)."
function contingency!(sim, s)
    pfd = sim.case.pfd
    schedule_bus_fault!(sim; bus=pfd.bus_num[mod1(3s + 1, nbus(pfd))], t_on=ongrid(0.20), t_off=ongrid(0.26),
                        Zflt=0.01, phases=iseven(s) ? [1] : [1, 2, 3])
    schedule_reference_step!(sim; t=ongrid(0.40), gen=mod1(s, ngen(pfd)), kind=:governor, dsp=-0.01)
    s % 3 == 0 && schedule_generator_trip!(sim; t=ongrid(0.50), gen=mod1(s + 1, ngen(pfd)))
    return sim
end

function timed_batch(case, nscen, workers, factorization, nsteps)
    batch = build_batch(case, nscen; ts=TS, nworkers=workers, factorization)
    for s in 1:nscen
        contingency!(batch[s], s)
    end
    GC.gc()
    t = @elapsed run_steps!(batch, nsteps)
    return t, batch
end

case = load_psse(joinpath(PSSE, "WECC240.raw"), joinpath(PSSE, "WECC240.dyr"))
println(case)
P = nthreads()
println("$S scenarios of 0.6 s ($NSTEPS steps) on $P threads")
# compile everything once on a small batch
for f in (:per_worker, :per_scenario, :shared)
    timed_batch(case, 2, min(2, P), f, 5500)
end

t1, ref = timed_batch(case, S, 1, :per_worker, NSTEPS)
@printf("  1 worker:                       %8.1f s\n", t1)
for f in (:per_worker, :per_scenario, :shared)
    t, batch = timed_batch(case, S, P, f, NSTEPS)
    same = all(batch[s].net.V == ref[s].net.V for s in 1:S)
    @printf("  %2d workers, %-13s %3d LU: %8.1f s  speedup %5.2f  identical results: %s\n",
            P, string(f), length(batch.factors), t, t1 / t, same)
end
