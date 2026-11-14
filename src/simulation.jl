# ==========================================================================
#  simulation.jl -- one EMT simulation: initialization, disturbances, time
#  step and time loop
#
#  build_simulation(case)  initializes the network, machines, controllers,
#                          IBRs and bus measurements at the power flow
#  schedule_*!(sim; ...)   bus faults, generator trips, reference steps
#  step!(sim)              advances the whole system by one time step
#  run!(sim, tend)         time loop with result recording
#
#  One time step (the order of ParaEMT's main loop):
#     1. apply the disturbances that fall on this step
#     2. extrapolate the machine quantities needed before the solution
#     3. machine Norton currents and IBR currents       -> i_s
#     4. network solution  G v = i_s + i_h, corrected by compensation for
#        the faults and trips in effect
#     5. bus measurements (PLL, voltage filter) from the new voltages
#     6. machines: algebraic solution, swing equation, stabilizer, exciter,
#        governor
#     7. IBR controllers
#     8. history currents of the branches and of the Bergeron lines for the
#        next step
#  Every array is allocated by build_simulation; a step allocates no memory.
# ==========================================================================

using Printf

"""
    SimSettings(; ts=50e-6, variant=:updated, bergeron_threshold=1.0, network_events=:compensate)

* `ts`                 -- time step [s]
* `variant`            -- model set: `:updated` (the paper's) or `:paraemt`
                          (ParaEMT), see [`MODEL_VARIANTS`](@ref)
* `bergeron_threshold` -- lines whose travel time exceeds this many time steps
                          are Bergeron lines; `Inf` keeps every line a pi section
* `network_events`     -- how faults and trips change the network:
                          `:compensate` keeps G and its factorization and
                          corrects the solution (paper, Section III);
                          `:refactorize` adds the change to G and factorizes
                          again, as ParaEMT does for a trip.  Both
                          give the same result up to round-off.
"""
Base.@kwdef struct SimSettings
    ts::Float64 = 50e-6
    variant::Symbol = :updated
    bergeron_threshold::Float64 = 1.0
    network_events::Symbol = :compensate
end

"""
The machine predictor uses one past point for the first 10 steps, two for the
next 10 and three afterwards (ParaEMT ties this to its output down-sampling
rate of 10).
"""
const PREDICTOR_RAMP = 10

# --------------------------------------------------------------------------
# The states of all generators (machine and its controllers)
# --------------------------------------------------------------------------

"States of all synchronous generators and their controllers."
mutable struct GeneratorStates
    mac::MachineState
    exc::SexsState
    tgov1::Tgov1State
    hygov::HygovState
    gast::GastState
    pss::IeeestState
end
function GeneratorStates(dyd::DynamicData, n::Int)
    return GeneratorStates(MachineState(n), SexsState(n), Tgov1State(length(dyd.tgov1.gen)),
                           HygovState(length(dyd.hygov.gen)), GastState(length(dyd.gast.gen)),
                           IeeestState(length(dyd.ieeest.gen)))
end
function Base.copyto!(dst::GeneratorStates, src::GeneratorStates)
    for f in fieldnames(GeneratorStates)
        copyto!(getfield(dst, f), getfield(src, f))
    end
    return dst
end

"Mechanical power of generator `i` from the governor that drives it."
@inline function mechanical_power(x::GeneratorStates, dyd::DynamicData, i::Int)
    k = dyd.gov_instance[i]
    T = dyd.gov_type[i]
    return T == TGOV1 ? x.tgov1.pm[k] : (T == HYGOV ? x.hygov.pm[k] : x.gast.pm[k])
end

# --------------------------------------------------------------------------
# Disturbances
# --------------------------------------------------------------------------

"Step `dsp` of the governor (`:governor`) or exciter (`:exciter`) reference of generator `gen` at time `t`."
mutable struct ReferenceStep
    t::Float64
    gen::Int
    kind::Symbol
    dsp::Float64
    done::Bool
end

"Generator `gen` is disconnected at time `t`."
mutable struct GeneratorTrip
    t::Float64
    gen::Int
    done::Bool
end

"Fault through `Zflt` [pu] from the given phases of a bus to ground, from `t_on` to `t_off`."
mutable struct BusFault
    t_on::Float64
    t_off::Float64
    bus::Int              # bus index (position in pfd.bus_num)
    phases::Vector{Int}   # subset of (1, 2, 3) = (a, b, c)
    Zflt::Float64
    active::Bool
    done::Bool
end

# --------------------------------------------------------------------------
# The simulation
# --------------------------------------------------------------------------

"""
    EmtSimulation

The case, the constants, the states and the disturbances of one simulation
(one scenario).  Created by [`build_simulation`](@ref), advanced by
[`step!`](@ref) and [`run!`](@ref).
"""
mutable struct EmtSimulation
    case::CaseData
    settings::SimSettings
    net::Network
    mc::MachineConstants
    ibrc::IbrConstants
    # states
    x::GeneratorStates          # generator states of the current step
    x_h1::GeneratorStates       # one step earlier
    x_h2::GeneratorStates       # two steps earlier
    x_ibr::IbrState
    x_bus::BusState
    # controller references
    gref::Vector{Float64}       # governor references
    vref::Vector{Float64}       # exciter voltage references
    # work arrays of a step
    work::MachineWork
    I_inj::Vector{Float64}      # machine and IBR source currents i_s [network pu]
    Iibr::Vector{Float64}       # three-phase current of every IBR
    # bookkeeping
    tn::Int                     # number of completed steps
    gen_tripped::Vector{Bool}
    reinit_pending::Bool        # a trip happened on this step: rebuild the histories
    ref_steps::Vector{ReferenceStep}
    trips::Vector{GeneratorTrip}
    faults::Vector{BusFault}
    # compensation matrices of this scenario (network_events = :compensate);
    # built on the first step, once the disturbances are known
    comp::Union{Nothing,Compensation}
end

"Simulated time [s]."
current_time(sim::EmtSimulation) = sim.tn * sim.settings.ts

"""
    build_simulation(case; ts=50e-6, variant=:updated, bergeron_threshold=1.0,
                     network_events=:compensate) -> EmtSimulation

Initializes every part of the system at the power-flow operating point.  The
keywords are those of [`SimSettings`](@ref).
"""
function build_simulation(case::CaseData; kwargs...)
    settings = SimSettings(; kwargs...)
    settings.variant in MODEL_VARIANTS || error("variant must be one of $MODEL_VARIANTS")
    settings.network_events in (:compensate, :refactorize) ||
        error("network_events must be :compensate or :refactorize")
    pfd = case.pfd; dyd = case.dyd
    ts = settings.ts
    n = ngen(pfd)

    # network: companion models, Bergeron lines, G of the passive network
    net = build_network(pfd, ts; variant=settings.variant, bergeron_threshold=settings.bergeron_threshold)
    # machines: equivalent circuits and their Norton conductances in G; then factorize G
    mc = build_machine_constants(pfd, dyd, net, ts)
    stamp_generators!(net, mc)
    refactorize!(net)

    # generator states in equilibrium with the power flow
    x = GeneratorStates(dyd, n)
    EFD0, pm0 = init_machines!(x.mac, mc, pfd, net; variant=settings.variant)
    resid = check_machine_equilibrium(x.mac, mc)
    resid > 1e-8 && @warn "machine initialization residual $resid"
    Vt0 = [pfd.bus_Vm[mc.bus[i]] for i in 1:n]
    vref = init_sexs!(x.exc, dyd.sexs, EFD0, Vt0)
    gref = init_governors!(x.tgov1, x.hygov, x.gast, dyd, pm0)
    # the IEEEST states start at zero (no speed deviation)

    ibrc = build_ibr_constants(pfd, dyd, net, ts)
    x_ibr = IbrState(ibrc.n)
    init_ibr!(x_ibr, ibrc, pfd; variant=settings.variant)

    x_bus = BusState(nbus(pfd))
    init_bus_measurement!(x_bus, pfd)

    x_h1 = GeneratorStates(dyd, n); copyto!(x_h1, x)
    x_h2 = GeneratorStates(dyd, n); copyto!(x_h2, x)

    return EmtSimulation(case, settings, net, mc, ibrc, x, x_h1, x_h2, x_ibr, x_bus, gref, vref,
                         MachineWork(n), zeros(net.N), zeros(3ibrc.n), 0, falses(n), false,
                         ReferenceStep[], GeneratorTrip[], BusFault[], nothing)
end

# --------------------------------------------------------------------------
# Scheduling disturbances
# --------------------------------------------------------------------------

"""
    schedule_reference_step!(sim; t, gen, kind=:governor, dsp)

Adds `dsp` to the governor reference (`kind = :governor`) or to the exciter
voltage reference (`kind = :exciter`) of generator `gen` at time `t`.
"""
function schedule_reference_step!(sim::EmtSimulation; t::Real, gen::Int, kind::Symbol=:governor, dsp::Real)
    kind in (:governor, :exciter) || error("kind must be :governor or :exciter")
    push!(sim.ref_steps, ReferenceStep(t, gen, kind, dsp, false))
    return sim
end

"""
    schedule_generator_trip!(sim; t, gen)

Disconnects generator `gen` (index into `pfd.gen_bus`) at time `t`.
"""
function schedule_generator_trip!(sim::EmtSimulation; t::Real, gen::Int)
    push!(sim.trips, GeneratorTrip(t, gen, false))
    sim.comp = nothing                 # rebuilt on the next step to include this event
    return sim
end

"""
    schedule_bus_fault!(sim; bus, t_on, t_off=Inf, Zflt=1e-3, phases=[1, 2, 3])

Fault of impedance `Zflt` [pu] from the given phases of bus number `bus` to
ground, from `t_on` to `t_off`.
"""
function schedule_bus_fault!(sim::EmtSimulation; bus::Int, t_on::Real, t_off::Real=Inf, Zflt::Real=1e-3,
                             phases::AbstractVector{Int}=[1, 2, 3])
    k = bus_index(sim.case.pfd, bus)
    push!(sim.faults, BusFault(t_on, t_off, k, collect(phases), Zflt, false, false))
    sim.comp = nothing                 # rebuilt on the next step to include this event
    return sim
end

fault_nodes(sim::EmtSimulation, f::BusFault) = [phase_node(f.bus, ph, sim.net.nbus) for ph in f.phases]
machine_nodes(sim::EmtSimulation, gen::Int) = [phase_node(sim.mc.bus[gen], ph, sim.net.nbus) for ph in 1:3]

# position of fault k and of trip k among the changes of the compensation (faults first)
@inline fault_slot(sim::EmtSimulation, k::Int) = k
@inline trip_slot(sim::EmtSimulation, k::Int) = length(sim.faults) + k

"""
    ensure_compensation!(sim) -> Compensation or nothing

With `network_events = :compensate`: builds the compensation matrices of all
scheduled faults and trips the first time they are needed (paper,
Section III) and puts into effect whatever has already happened.
"""
function ensure_compensation!(sim::EmtSimulation)
    sim.settings.network_events === :compensate || return nothing
    sim.comp === nothing || return sim.comp
    changes = NetworkChange[]
    for f in sim.faults
        push!(changes, fault_change(fault_nodes(sim, f), f.Zflt))
    end
    for e in sim.trips
        push!(changes, trip_change(machine_nodes(sim, e.gen), sim.mc.Gequiv[e.gen]))
    end
    sim.comp = build_compensation(sim.net, changes)
    for (k, f) in enumerate(sim.faults)
        f.active && set_active!(sim.comp, fault_slot(sim, k), true)
    end
    for (k, e) in enumerate(sim.trips)
        e.done && set_active!(sim.comp, trip_slot(sim, k), true)
    end
    return sim.comp
end

"""
    apply_events!(sim)

Applies the reference steps, trips and fault switchings that fall on the step
about to be computed.  With `:compensate` an event only changes which
compensation changes are in effect; with `:refactorize` it is added to G,
which is factorized again.
"""
function apply_events!(sim::EmtSimulation)
    t = current_time(sim)
    comp = ensure_compensation!(sim)
    for e in sim.ref_steps
        (e.done || t < e.t) && continue
        e.kind == :governor ? (sim.gref[e.gen] += e.dsp) : (sim.vref[e.gen] += e.dsp)
        e.done = true
    end
    changed = false
    for (k, e) in enumerate(sim.trips)
        (e.done || t < e.t) && continue
        if comp === nothing
            add_conductance_block!(sim.net, machine_nodes(sim, e.gen), -sim.mc.Gequiv[e.gen])
            changed = true
        else
            set_active!(comp, trip_slot(sim, k), true)
        end
        sim.gen_tripped[e.gen] = true         # its source current is no longer injected
        sim.reinit_pending = true             # and the histories are rebuilt at the end of the step
        e.done = true
    end
    for (k, f) in enumerate(sim.faults)
        if !f.active && !f.done && t >= f.t_on                # inception
            if comp === nothing
                add_conductance_block!(sim.net, fault_nodes(sim, f), fault_change(fault_nodes(sim, f), f.Zflt).dY)
                changed = true
            else
                set_active!(comp, fault_slot(sim, k), true)
            end
            f.active = true
        elseif f.active && t >= f.t_off                      # clearing
            if comp === nothing
                add_conductance_block!(sim.net, fault_nodes(sim, f), -fault_change(fault_nodes(sim, f), f.Zflt).dY)
                changed = true
            else
                set_active!(comp, fault_slot(sim, k), false)
            end
            f.active = false; f.done = true
        end
    end
    changed && refactorize!(sim.net)
    return sim
end

# --------------------------------------------------------------------------
# One time step
# --------------------------------------------------------------------------

"Number of past points of the machine predictor (1, 2 or 3) for the step about to be computed."
@inline predictor_order(sim::EmtSimulation) = min(sim.tn ÷ PREDICTOR_RAMP + 1, 3)

"""
    step!(sim) -> V

Advances the simulation by one time step and returns the node voltages.
"""
function step!(sim::EmtSimulation)
    st = sim.settings; ts = st.ts
    net = sim.net; mc = sim.mc; dyd = sim.case.dyd
    nb = net.nbus
    x = sim.x; w = sim.work

    # 1. disturbances (an event at time t applies to the step that reaches t)
    order = predictor_order(sim)
    sim.tn += 1
    apply_events!(sim)

    # 2. machine quantities extrapolated from the last three steps
    predict!(w, mc, x.mac, sim.x_h1.mac, sim.x_h2.mac, x.exc.EFD, sim.x_h1.exc.EFD, sim.x_h2.exc.EFD, order)

    # 3. source currents i_s (tripped machines inject nothing)
    fill!(sim.I_inj, 0.0)
    norton_current!(sim.I_inj, w, mc, x.mac, x.exc.EFD, nb, sim.gen_tripped)
    ibr_current!(sim.I_inj, sim.Iibr, sim.x_ibr, sim.ibrc, sim.x_bus, nb)

    # 4. network solution, eq. (1), with the compensation of eq. (3)
    V = solve_network!(net, sim.I_inj, sim.comp)

    # 5. bus measurements from the new voltages
    update_bus_measurement!(sim.x_bus, dyd.meas, V, nb, ts, mc.ws)

    # 6. generators: keep the two previous states, then update in place
    copyto!(sim.x_h2, sim.x_h1)
    copyto!(sim.x_h1, x)
    pss_derivative = st.variant === :updated       # the updated SEXS lead term sees d(vs)/dt
    @inbounds for i in 1:mc.n
        pm = mechanical_power(x, dyd, i)
        pe = machine_algebraic!(x.mac, w, mc, V, nb, i)
        x.mac.omega[i], x.mac.delta[i] = swing_step(x.mac.omega[i], x.mac.delta[i], pe, pm,
                                                    mc.H[i], mc.D[i], mc.ws, ts)
        dw = x.mac.omega[i] / mc.ws - 1                             # speed deviation [pu]
        # stabilizer (if any) -> exciter -> governor
        kp = dyd.pss_instance[i]
        vs_1 = kp == 0 ? 0.0 : sim.x_h1.pss.vs[kp]
        vs = kp == 0 ? 0.0 : update_ieeest!(x.pss, dyd.ieeest, kp, dw, ts)
        dvs = pss_derivative ? (vs - vs_1) / ts : 0.0
        k = mc.bus[i]
        update_sexs!(x.exc, dyd.sexs, i, sim.vref[i], sim.x_bus.vtm[k], sim.x_bus.dvtm[k], vs, dvs, ts)
        kg = dyd.gov_instance[i]
        T = dyd.gov_type[i]
        if T == TGOV1
            update_tgov1!(x.tgov1, dyd.tgov1, kg, sim.gref[i], dw, ts)
        elseif T == HYGOV
            update_hygov!(x.hygov, dyd.hygov, kg, sim.gref[i], dw, ts)
        else
            update_gast!(x.gast, dyd.gast, kg, sim.gref[i], dw, ts)
        end
    end

    # 7. IBR controllers
    update_ibr!(sim.x_ibr, sim.ibrc, V, sim.Iibr, sim.x_bus, nb)

    # 8. history currents for the next step
    if sim.reinit_pending
        reinitialize_history!(sim)
        sim.reinit_pending = false
    else
        update_history!(net)
    end
    if !isempty(net.ber)
        update_bergeron!(net.ber, V)
        bergeron_injection!(net.ber_node_Ihis, net.ber)
    end
    return V
end

"""
    reinitialize_history!(sim)

ParaEMT's `Re_Init`, applied on the step on which a generator trips: the
branch history currents are rebuilt so that the network memory is consistent
with the modified network.
  1. intermediate histories from the previous-step voltages (per branch the
     current term or the voltage term of the history source, chosen by the
     sign of `icf`, as in ParaEMT);
  2. machine sources from ParaEMT's two-point extrapolation of the generator
     states of the two previous steps, 2 x(t_{n-2}) - x(t_{n-1});
  3. a solution of the MODIFIED network with these sources (the Bergeron
     history currents are added unchanged);
  4. the final histories from that voltage and the previous branch currents.
"""
function reinitialize_history!(sim::EmtSimulation)
    net = sim.net; mc = sim.mc
    nb = net.nbus; N = net.N; ws = mc.ws; ts = mc.ts
    Vp = net.V_prev; Ipre = net.brch_Ipre
    xa = sim.x_h2; xb = sim.x_h1               # states two steps and one step before the current one

    # 1. intermediate histories
    node_mid = zeros(N)
    for k in eachindex(net.F)
        f = net.F[k]; t = net.T[k]; r = net.tr[k]
        c_i = net.icf[k] > 0 ? 1.0 : 0.0
        c_v = 1.0 - c_i
        dV = t == 0 ? Vp[f] : Vp[f] * r - Vp[t]
        h = c_i * net.icf[k] * Ipre[k] + c_v * net.Gv1[k] * dV
        t == 0 || (node_mid[t] += h)
        node_mid[f] -= h * r
    end
    # 2. machine sources from the extrapolated states
    w = MachineWork(mc.n)
    for i in 1:mc.n
        w.pd_omega[i] = 2 * xa.mac.omega[i] - xb.mac.omega[i]
        w.pd_id[i] = 2 * xa.mac.id[i] - xb.mac.id[i]
        w.pd_iq[i] = 2 * xa.mac.iq[i] - xb.mac.iq[i]
        w.pd_EFD[i] = 2 * xa.exc.EFD[i] - xb.exc.EFD[i]
        w.pd_ud[i] = 2 * speed_voltage_d(xa.mac, i, ws) - speed_voltage_d(xb.mac, i, ws)
        w.pd_uq[i] = 2 * speed_voltage_q(xa.mac, i, ws) - speed_voltage_q(xb.mac, i, ws)
        w.pd_delta[i] = xa.mac.delta[i] + ts / 2 * (w.pd_omega[i] + xa.mac.omega[i]) / 2
        w.theta[i] = w.pd_delta[i] - pi / 2
    end
    Igs = zeros(N)
    norton_current!(Igs, w, mc, xa.mac, xa.exc.EFD, nb, sim.gen_tripped)
    # 3. solution of the modified network
    rhs = Igs .+ node_mid
    isempty(net.ber) || (rhs .+= net.ber_node_Ihis)
    Vr = net.lu \ rhs
    sim.comp === nothing || apply_compensation!(Vr, sim.comp)
    # 4. final histories
    history_from!(net.brch_Ihis, net.node_Ihis, net, Ipre, Vr)
    return sim
end

# --------------------------------------------------------------------------
# Time loop
# --------------------------------------------------------------------------

"""
    run!(sim, tend; record_every=10, progress=true) -> Results

Runs the simulation from the current time to `tend` [s] and records the
results every `record_every` steps.
"""
function run!(sim::EmtSimulation, tend::Real; record_every::Int=10, progress::Bool=true)
    ts = sim.settings.ts
    nsteps = round(Int, tend / ts) - sim.tn
    pfd = sim.case.pfd
    nt = nsteps ÷ record_every + 1
    res = Results(sim.net.N, nbus(pfd), ngen(pfd), nibr(pfd), nt, record_every)
    record!(res, sim)
    t0 = time()
    for k in 1:nsteps
        step!(sim)
        sim.tn % record_every == 0 && record!(res, sim)
        if progress && (k % max(1, nsteps ÷ 20) == 0 || k == nsteps)
            @printf("  t = %7.4f s  (%3.0f%%)  elapsed %.1f s\n", current_time(sim), 100k / nsteps, time() - t0)
        end
    end
    finalize!(res)
    return res
end

"Stores the current state of `sim` in the results."
function record!(r::Results, sim::EmtSimulation)
    j = r.n + 1
    j > length(r.t) && return r
    r.t[j] = current_time(sim)
    r.V[:, j] .= sim.net.V
    r.vtm[:, j] .= sim.x_bus.vtm
    r.freq[:, j] .= sim.x_bus.we
    x = sim.x
    r.omega[:, j] .= x.mac.omega
    r.delta[:, j] .= x.mac.delta
    r.pe[:, j] .= x.mac.psid .* x.mac.iq .- x.mac.psiq .* x.mac.id
    for i in 1:sim.mc.n
        r.pm[i, j] = mechanical_power(x, sim.case.dyd, i)
    end
    r.EFD[:, j] .= x.exc.EFD
    r.ibr_P[:, j] .= sim.x_ibr.repca.Pe
    r.ibr_Q[:, j] .= sim.x_ibr.repca.Qe
    r.n = j
    return r
end
