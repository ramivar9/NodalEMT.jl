# ==========================================================================
#  batch.jl -- a batch of scenarios on worker threads (paper, Section III)
#
#  The scenarios of a batch are the same power system with different
#  disturbances, so they never interact.  With compensation (compensation.jl)
#  G never changes, so G, its factors, the companion-model table, the
#  machine and IBR constants and the Bergeron-line parameters are stored once
#  for all scenarios.  Each scenario owns its states, histories, ring buffers
#  and compensation matrices.
#
#  The S scenarios are split into P equal contiguous chunks, one per worker
#  thread.  Each worker factorizes G once and simulates the scenarios of its
#  chunk one after the other, each from the initial state to the end.  No
#  worker waits for another, and every scenario gives exactly the result of a
#  single-threaded run.
#
#  Why one factorization per worker: Julia's UMFPACK factorization object
#  serializes its substitutions (it takes a lock), so one factorization
#  shared by all threads makes them wait for each other.  On the 500-scenario
#  WECC batch of the paper, 16 workers run 12.2 times faster than one with a
#  factorization per worker, and 7.6 times faster with a single shared one.
#  A factorization per scenario is as fast as one per worker but computes and
#  stores S factorizations instead of P.  The three choices are available
#  through the `factorization` keyword of `build_batch`.
# ==========================================================================

using Base.Threads

"""
    ScenarioBatch

`S` scenarios of one system that share every constant.  Build it with
[`build_batch`](@ref), schedule the disturbances of scenario `s` on
`batch[s]` with the usual `schedule_*!` functions, and advance with
[`run!`](@ref) or [`run_steps!`](@ref).
"""
struct ScenarioBatch
    sims::Vector{EmtSimulation}
    chunks::Vector{UnitRange{Int}}     # the scenarios simulated by each worker thread
    factors::Vector{Any}               # the distinct LU factorizations of G
    template::EmtSimulation            # the initialized system all scenarios were copied from
end

nscenarios(b::ScenarioBatch) = length(b.sims)
nworkers(b::ScenarioBatch) = length(b.chunks)
Base.getindex(b::ScenarioBatch, s::Int) = b.sims[s]
Base.length(b::ScenarioBatch) = length(b.sims)

"Splits `n` items into `k` contiguous chunks whose lengths differ by at most one."
function chunk_ranges(n::Int, k::Int)
    k = max(1, min(k, n))
    base, rem = divrem(n, k)
    out = UnitRange{Int}[]
    lo = 1
    for i in 1:k
        len = base + (i <= rem ? 1 : 0)
        push!(out, lo:(lo+len-1))
        lo += len
    end
    return out
end

"""
    new_scenario(src, lu) -> EmtSimulation

A scenario in the initial state of `src` that shares its constants and uses
the factorization `lu`, with its own states and no disturbances.
"""
function new_scenario(src::EmtSimulation, lu_shared)
    dyd = src.case.dyd
    n = src.mc.n
    x = GeneratorStates(dyd, n); copyto!(x, src.x)
    x1 = GeneratorStates(dyd, n); copyto!(x1, src.x_h1)
    x2 = GeneratorStates(dyd, n); copyto!(x2, src.x_h2)
    xi = IbrState(src.ibrc.n); copyto!(xi, src.x_ibr)
    xb = BusState(src.net.nbus); copyto!(xb, src.x_bus)
    return EmtSimulation(src.case, src.settings, copy_state(src.net, lu_shared),
                         src.mc, src.ibrc,                                  # shared constants
                         x, x1, x2, xi, xb, copy(src.gref), copy(src.vref),
                         MachineWork(n), zeros(src.net.N), zeros(3src.ibrc.n),
                         src.tn, copy(src.gen_tripped), src.reinit_pending,
                         ReferenceStep[], GeneratorTrip[], BusFault[], nothing)
end

"""
    build_batch(case, S; nworkers=Threads.nthreads(), factorization=:per_worker, kwargs...) -> ScenarioBatch

Initializes the system once and copies it into `S` scenarios, split into
`nworkers` chunks.  `factorization` selects who owns the LU factors of G:

* `:per_worker`   -- one factorization per worker (NodalEMT.jl's choice)
* `:per_scenario` -- one factorization per scenario
* `:shared`       -- one factorization for the whole batch

`kwargs` are the settings of [`build_simulation`](@ref); faults and trips must
be applied by compensation (`network_events = :compensate`, the default),
since a shared G cannot be modified.
"""
function build_batch(case::CaseData, S::Int; nworkers::Int=Threads.nthreads(),
                     factorization::Symbol=:per_worker, kwargs...)
    S >= 1 || error("a batch needs at least one scenario")
    factorization in (:per_worker, :per_scenario, :shared) ||
        error("factorization must be :per_worker, :per_scenario or :shared")
    template = build_simulation(case; kwargs...)
    template.settings.network_events === :compensate ||
        error("the scenarios of a batch share G, so network_events must be :compensate")
    chunks = chunk_ranges(S, nworkers)
    G = template.net.G
    if factorization === :per_worker
        factors = Any[k == 1 ? template.net.lu : lu(G) for k in eachindex(chunks)]
        owner = zeros(Int, S)
        for (k, r) in enumerate(chunks), s in r
            owner[s] = k
        end
    elseif factorization === :per_scenario
        factors = Any[s == 1 ? template.net.lu : lu(G) for s in 1:S]
        owner = collect(1:S)
    else
        factors = Any[template.net.lu]
        owner = ones(Int, S)
    end
    sims = [new_scenario(template, factors[owner[s]]) for s in 1:S]
    return ScenarioBatch(sims, chunks, factors, template)
end

"""
    run_steps!(batch, nsteps) -> batch

Advances every scenario by `nsteps` steps, without recording.  Each worker
runs the scenarios of its chunk one after the other.
"""
function run_steps!(b::ScenarioBatch, nsteps::Int)
    ch = b.chunks
    @threads for k in eachindex(ch)
        @inbounds for s in ch[k]
            sim = b.sims[s]
            for _ in 1:nsteps
                step!(sim)
            end
        end
    end
    return b
end

"""
    run!(batch, tend; record_every=10, progress=true) -> Vector{Results}

Runs every scenario to `tend` [s] and returns one `Results` per scenario.
`record_every = 0` records nothing.
"""
function run!(b::ScenarioBatch, tend::Real; record_every::Int=10, progress::Bool=true)
    ts = b.template.settings.ts
    pfd = b.template.case.pfd
    nsteps = round(Int, tend / ts) - b.sims[1].tn
    rec = record_every > 0
    nt = rec ? nsteps ÷ record_every + 1 : 1
    res = [Results(b.template.net.N, nbus(pfd), ngen(pfd), nibr(pfd), nt, max(record_every, 1))
           for _ in eachindex(b.sims)]
    t0 = time()
    ch = b.chunks
    @threads for k in eachindex(ch)
        @inbounds for s in ch[k]
            sim = b.sims[s]
            rec && record!(res[s], sim)
            for _ in 1:nsteps
                step!(sim)
                rec && sim.tn % record_every == 0 && record!(res[s], sim)
            end
            rec && finalize!(res[s])
        end
    end
    if progress
        el = time() - t0
        @printf("  %d scenarios x %d steps in %.2f s on %d workers\n", length(b.sims), nsteps, el, length(ch))
    end
    return res
end
