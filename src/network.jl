# ==========================================================================
#  network.jl -- the electrical network (Dommel's nodal formulation)
#
#  Every passive element (line R-L, line charging C, transformer R-L,
#  constant-impedance load, shunt capacitor) is replaced by its trapezoidal
#  "companion model": a conductance 1/Req in parallel with a history current
#  source that depends on the previous time step.  Collecting the
#  conductances gives the nodal matrix G, which is constant for a fixed time
#  step, so the whole network is solved at every step from (paper, eq. (1))
#
#          G v(t) = i_s(t) + i_h(t)
#
#  where i_s are the machine and IBR source currents and i_h the history
#  currents (of the companion models and of the Bergeron lines).
#
#  Node numbering: three nodes per bus,
#      phase a of bus k = k,   phase b = k + nbus,   phase c = k + 2 nbus.
#  Branch k of the companion table connects node F[k] to node T[k] (T[k] = 0
#  is ground) and has the constants
#      Req : companion resistance (stamped into G as 1/Req)
#      icf : history coefficient of the branch current
#      Gv1 : history coefficient of the branch voltage
#      tr  : 1/t for a transformer with the off-nominal ratio t at its F end,
#            1 for every other branch
#  so that after each solution
#      dV       = tr V[F] - V[T]
#      I_branch = dV/Req + I_his             (current of this step)
#      I_his    = icf I_branch + Gv1 dV      (history source of the next step)
#  The branch current enters node T and leaves node F as tr I_branch.  This is
#  the transformer model of the power flow (tap at the from end, impedance at
#  the to end): G gets tr^2/Req at (F,F), 1/Req at (T,T) and -tr/Req off the
#  diagonal.
#
#  Lines whose travel time exceeds `bergeron_threshold` time steps are not in
#  this table; they are Bergeron lines (bergeron.jl).
# ==========================================================================

using SparseArrays
using LinearAlgebra

# --------------------------------------------------------------------------
# Model sets: the ParaEMT models and the updated ones
# --------------------------------------------------------------------------

"""
    MODEL_VARIANTS = (:updated, :paraemt)

`:updated` (the default) is the model set of the paper (Section II);
`:paraemt` reproduces the models of ParaEMT (public release, commit d79d735).
The variant selects

* the numerical damping of the companion models ([`damping_constants`](@ref));
* the PSS derivative in the SEXS lead term (`sexs_step` in controls.jl);
* the initialization of REECB, REPCA and motoring machines (ibr.jl, machine.jl);
* the companion model of series capacitors (lines with X < 0): `:paraemt`
  reproduces ParaEMT, which uses for their voltage coefficient Gv1 the
  value of the previously processed inductive line.
"""
const MODEL_VARIANTS = (:updated, :paraemt)

"""
    DampingConstants

Numerical damping of the trapezoidal rule.  ParaEMT places a large resistance
Rp in parallel with every inductor and a small resistance Rs in series with
every capacitor:

* `line_damptrap`    -- Rp = line_damptrap * (20/3) * (2L/ts) for line inductors
* `load_damping`     -- inductive loads get the same parallel resistance
* `chg_rs_factor`    -- Rs = chg_rs_factor * ts/(2C) for the line-charging capacitors
* `shunt_rs_factor`  -- the same factor for fixed and switched shunt capacitors
* `load_c_rs_factor` -- the same factor for the capacitor of a capacitive load
"""
struct DampingConstants
    line_damptrap::Float64
    load_damping::Bool
    chg_rs_factor::Float64
    shunt_rs_factor::Float64
    load_c_rs_factor::Float64
end

"""
    damping_constants(variant) -> DampingConstants

ParaEMT damps every inductor and capacitor.  The updated model set
weakens the damping of the lines tenfold and removes that of loads and
shunts, whose damping resistors carry power-frequency currents that the power
flow does not contain.
"""
function damping_constants(variant::Symbol)
    variant === :paraemt && return DampingConstants(1.0, true, 0.15, 0.15, 0.15e-3)
    variant === :updated && return DampingConstants(10.0, false, 0.015, 0.0, 0.0)
    error("variant must be one of $MODEL_VARIANTS, got :$variant")
end

# --------------------------------------------------------------------------
# The network
# --------------------------------------------------------------------------

"""
    Network

The companion-model table, G and its LU factorization (constants, shared by
all scenarios of a batch), and the history currents and node voltages (the
state of one scenario).
"""
mutable struct Network
    nbus::Int
    N::Int                              # number of nodes = 3 nbus
    ZbaseA::Vector{Float64}             # per-bus impedance base [ohm]
    IbaseA::Vector{Float64}             # per-bus current base [kA]
    # companion-model table (one row per branch and phase)
    F::Vector{Int}
    T::Vector{Int}                      # 0 = ground
    Req::Vector{Float64}
    icf::Vector{Float64}
    Gv1::Vector{Float64}
    tr::Vector{Float64}                 # 1/t of a transformer's F end, 1 otherwise
    # nodal matrix (passive network + machine Norton conductances) and its UMFPACK factorization
    G::SparseMatrixCSC{Float64,Int}
    lu::Any
    # state of the scenario
    brch_Ihis::Vector{Float64}          # history current of each branch
    brch_Ipre::Vector{Float64}          # branch current of the previous step
    node_Ihis::Vector{Float64}          # history currents summed per node
    V::Vector{Float64}                  # node voltages of the current step [pu]
    V_prev::Vector{Float64}             # node voltages of the previous step
    rhs::Vector{Float64}                # work vector: right-hand side of eq. (1)
    ber::BergeronLines                  # Bergeron lines (empty if every line is a pi section)
    ber_node_Ihis::Vector{Float64}      # Bergeron history currents per node
    # power-flow phasors (used to initialize the devices)
    Vt0::Vector{ComplexF64}             # node voltage phasors at t = 0
    It0::Vector{ComplexF64}             # generator current phasors at t = 0 (system base)
end

"Node index of phase `ph` (1, 2, 3 = a, b, c) of bus `k`."
@inline phase_node(k::Int, ph::Int, nbus::Int) = k + (ph - 1) * nbus

# --------------------------------------------------------------------------
# Companion models
# --------------------------------------------------------------------------

"""
    rl_series_companion(R, X, ws, ts; damptrap) -> (Req, icf, Gv1)

Series R-L branch (X = ws L) with the damping resistance
Rp = damptrap (20/3) (2L/ts) in parallel with the inductor; `damptrap = 0`
means no damping resistor (transformers).
"""
function rl_series_companion(R::Float64, X::Float64, ws::Float64, ts::Float64; damptrap::Float64)
    L = X / ws
    a = ts / (2L)                      # trapezoidal factor of the inductor
    Rp_inv = damptrap == 0 ? 0.0 : 1.0 / (damptrap * 20.0 / 3.0 * 2.0 * L / ts)
    Req = (1 + R * (a + Rp_inv)) / (a + Rp_inv)
    icf = (1 - R * (a - Rp_inv)) / (1 + R * (a + Rp_inv))
    Gv1 = (a - Rp_inv) / (1 + R * (a + Rp_inv))
    return Req, icf, Gv1
end

"Series R-C branch (X = -1/(ws C) < 0, a series-compensated line)."
function rc_series_companion(R::Float64, X::Float64, ws::Float64, ts::Float64)
    C = -1 / (X * ws)
    Req = R + ts / (2C)
    icf = (2R * C - ts) / (2R * C + ts)
    Gv1 = -1 / Req
    return Req, icf, Gv1
end

"Capacitor C to ground with the damping resistance Rs = rs_factor ts/(2C) in series."
function c_shunt_companion(C::Float64, ts::Float64; rs_factor::Float64)
    Rc = ts / (2C)
    Rs = rs_factor * Rc
    Req = Rs + Rc
    icf = -(Rc - Rs) / (Rc + Rs)
    Gv1 = -1.0 / (Rc + Rs)
    return Req, icf, Gv1
end

# --------------------------------------------------------------------------
# Building the network
# --------------------------------------------------------------------------

"""
    build_network(pfd, ts; variant, bergeron_threshold) -> Network

Builds the companion-model table, the Bergeron lines and G of the passive
network, and the history currents of the power-flow steady state, so that the
first step starts from the operating point.  The machine conductances are
added to G later (`stamp_generators!`).
"""
function build_network(pfd::PowerFlowData, ts::Float64; variant::Symbol, bergeron_threshold::Real)
    dmp = damping_constants(variant)
    nb = nbus(pfd)
    N = 3nb
    ws = pfd.ws

    # ---- per-bus bases and the three-phase phasor voltages of the power flow ----
    VbaseA = pfd.bus_basekV ./ sqrt(3.0)
    ZbaseA = VbaseA .^ 2 .* 3 ./ pfd.basemva
    IbaseA = pfd.basemva ./ (3 .* VbaseA)
    Vt0 = zeros(ComplexF64, N)
    for k in 1:nb
        va = pfd.bus_Vm[k] * cis(pfd.bus_Va[k])
        Vt0[phase_node(k, 1, nb)] = va
        Vt0[phase_node(k, 2, nb)] = va * cis(-2pi / 3)
        Vt0[phase_node(k, 3, nb)] = va * cis(2pi / 3)
    end
    It0 = zeros(ComplexF64, ngen(pfd))
    for i in 1:ngen(pfd)
        k = bus_index(pfd, pfd.gen_bus[i])
        S = complex(pfd.gen_MW[i], pfd.gen_Mvar[i]) / pfd.basemva
        It0[i] = conj(S) / conj(Vt0[k])
    end

    # ---- companion table and the entries of G ----
    F = Int[]; T = Int[]; Req = Float64[]; icf = Float64[]; Gv1 = Float64[]; tr = Float64[]
    Ipre0 = ComplexF64[]                      # phasor current of each branch at t = 0
    rows = Int[]; cols = Int[]; vals = Float64[]
    stamp!(r, c, v) = (push!(rows, r); push!(cols, c); push!(vals, v))
    function add_branch!(f, t, req, ic, gv, ipre, r=1.0)
        push!(F, f); push!(T, t); push!(Req, req); push!(icf, ic); push!(Gv1, gv); push!(Ipre0, ipre); push!(tr, r)
    end
    function stamp_series!(fn, tn, req, r=1.0) # a branch between two nodes, ratio 1/r at fn
        g = 1 / req
        stamp!(fn, fn, g * r * r); stamp!(tn, tn, g)
        stamp!(fn, tn, -g * r); stamp!(tn, fn, -g * r)
    end
    # a shunt Y = G + jB [pu] at a node: capacitor (B > 0) or reactor (B < 0), and a resistor (G != 0)
    function add_shunt!(node, Y, rs_factor)
        B = imag(Y); Gs = real(Y)
        if B > 0
            C = B / ws
            reqc, icc, gvc = c_shunt_companion(C, ts; rs_factor=rs_factor)
            stamp!(node, node, 1 / reqc)
            add_branch!(node, 0, reqc, icc, gvc, Vt0[node] * im * ws * C)
        elseif B < 0
            reql, icl, gvl = rl_series_companion(0.0, -1 / B, ws, ts; damptrap=dmp.load_damping ? dmp.line_damptrap : 0.0)
            stamp!(node, node, 1 / reql)
            add_branch!(node, 0, reql, icl, gvl, Vt0[node] * im * B)
        end
        if Gs != 0
            stamp!(node, node, Gs)
            add_branch!(node, 0, 1 / Gs, 0.0, 0.0, Vt0[node] * Gs)
        end
    end

    # long lines become Bergeron lines: 1/Z at both terminals of each conductor
    ber_idx, pi_idx = select_bergeron_lines(pfd, ts, bergeron_threshold)
    ber = build_bergeron(pfd, ts, ber_idx, Vt0, nb)
    for k in 1:ber.n
        stamp!(ber.F[k], ber.F[k], 1 / ber.z[k])
        stamp!(ber.T[k], ber.T[k], 1 / ber.z[k])
    end

    # the other lines are pi sections: series R-L (or R-C), half the charging at each end
    last_Gv1 = 0.0
    for i in pi_idx
        f = bus_index(pfd, pfd.line_from[i]); t = bus_index(pfd, pfd.line_to[i])
        R = real(pfd.line_RX[i]); X = imag(pfd.line_RX[i])
        if X > 0
            req, ic, gv = rl_series_companion(R, X, ws, ts; damptrap=dmp.line_damptrap)
            last_Gv1 = gv
        elseif X < 0
            req, ic, gv = rc_series_companion(R, X, ws, ts)
            variant === :paraemt && (gv = last_Gv1)        # reproduces ParaEMT (see MODEL_VARIANTS)
        else
            error("line $i has X = 0")
        end
        C = pfd.line_chg[i] / 2 / ws
        for ph in 1:3
            fn = phase_node(f, ph, nb); tn = phase_node(t, ph, nb)
            stamp_series!(fn, tn, req)
            add_branch!(fn, tn, req, ic, gv, (Vt0[fn] - Vt0[tn]) / complex(R, X))
        end
        for endbus in (f, t), ph in 1:3          # three sending-end capacitors, then three receiving-end ones
            node = phase_node(endbus, ph, nb)
            if C > 0
                reqc, icc, gvc = c_shunt_companion(C, ts; rs_factor=dmp.chg_rs_factor)
                stamp!(node, node, 1 / reqc)
                add_branch!(node, 0, reqc, icc, gvc, Vt0[node] * im * ws * C)
            else
                # no charging: ParaEMT keeps a table entry with Req = 1e10 that is not in G
                # and carries a negligible current; it is kept to reproduce ParaEMT exactly
                add_branch!(node, 0, 1e10, -1.0, 0.0, 0.0 + 0im)
            end
        end
    end

    # transformers: series R-L without damping resistor, off-nominal ratio k at the from end
    for i in eachindex(pfd.xfmr_from)
        f = bus_index(pfd, pfd.xfmr_from[i]); t = bus_index(pfd, pfd.xfmr_to[i])
        R = real(pfd.xfmr_RX[i]); X = imag(pfd.xfmr_RX[i])
        k = pfd.xfmr_k[i]
        k > 0 || error("transformer $i has the ratio $k")
        r = 1 / k
        req, ic, gv = rl_series_companion(R, X, ws, ts; damptrap=0.0)
        for ph in 1:3
            fn = phase_node(f, ph, nb); tn = phase_node(t, ph, nb)
            stamp_series!(fn, tn, req, r)
            add_branch!(fn, tn, req, ic, gv, (Vt0[fn] * r - Vt0[tn]) / complex(R, X), r)
        end
    end

    # loads: constant impedance Z = |V|^2 / conj(S) at the power-flow voltage
    for i in eachindex(pfd.load_bus)
        k = bus_index(pfd, pfd.load_bus[i])
        Vm = abs(Vt0[k])
        Z = Vm^2 / (complex(pfd.load_MW[i], -pfd.load_Mvar[i]) / pfd.basemva)
        R = real(Z); X = imag(Z)
        if X > 0                                   # inductive: series R-L
            req, ic, gv = rl_series_companion(R, X, ws, ts; damptrap=dmp.load_damping ? dmp.line_damptrap : 0.0)
        elseif X < 0                               # capacitive: R in parallel with a damped C
            Y = 1 / Z
            Rl = 1 / real(Y); C = imag(Y) / ws
            Rc = ts / (2C); Rs = dmp.load_c_rs_factor * Rc
            req = 1 / (1 / Rl + 1 / (Rs + Rc))
            ic = -(Rc - Rs) / (Rc + Rs)
            gv = 1 / Rl - 1 / (Rc + Rs)
        else                                       # purely resistive
            req = R; ic = 0.0; gv = 0.0
        end
        for ph in 1:3
            node = phase_node(k, ph, nb)
            stamp!(node, node, 1 / req)
            add_branch!(node, 0, req, ic, gv, Vt0[node] / Z)
        end
    end

    # fixed shunts: G + jB in MW and Mvar at 1 pu (B > 0: capacitor)
    for i in eachindex(pfd.shnt_bus)
        k = bus_index(pfd, pfd.shnt_bus[i])
        Y = pfd.shnt_gb[i] / pfd.basemva
        for ph in 1:3
            add_shunt!(phase_node(k, ph, nb), Y, dmp.shunt_rs_factor)
        end
    end

    # switched shunts: MW and -Mvar drawn at the power-flow voltage (PSS/E's YSWACT),
    # converted to an admittance
    for i in eachindex(pfd.shnt_sw_bus)
        abs(pfd.shnt_sw_gb[i]) < 1e-15 && continue
        k = bus_index(pfd, pfd.shnt_sw_bus[i])
        Ysw = conj(-pfd.shnt_sw_gb[i]) / pfd.basemva / pfd.bus_Vm[k]^2
        Y = complex(real(Ysw), -imag(Ysw))
        for ph in 1:3
            add_shunt!(phase_node(k, ph, nb), Y, dmp.shunt_rs_factor)
        end
    end

    G = sparse(rows, cols, vals, N, N)      # repeated (row, col) entries are summed

    # ---- initial history currents from the phasor solution ----
    # the instantaneous value of a phasor at t = 0 is its real part
    V0 = real.(Vt0)
    brch_Ipre = real.(Ipre0)
    brch_Ihis = zeros(length(F))
    node_Ihis = zeros(N)
    net = Network(nb, N, ZbaseA, IbaseA, F, T, Req, icf, Gv1, tr, G, nothing,
                  brch_Ihis, brch_Ipre, node_Ihis, copy(V0), copy(V0), zeros(N),
                  ber, zeros(N), Vt0, It0)
    history_from!(net.brch_Ihis, net.node_Ihis, net, brch_Ipre, V0)
    isempty(ber) || bergeron_injection!(net.ber_node_Ihis, ber)
    return net
end

"""
    history_from!(brch_Ihis, node_Ihis, net, Ibr, V)

History source of every branch, I_his = icf Ibr + Gv1 dV, from the branch
currents `Ibr` and the node voltages `V`, and its sum per node.
"""
function history_from!(brch_Ihis::AbstractVector, node_Ihis::AbstractVector, net::Network,
                       Ibr::AbstractVector, V::AbstractVector)
    fill!(node_Ihis, 0.0)
    for k in eachindex(net.F)
        f = net.F[k]; t = net.T[k]; r = net.tr[k]
        if t == 0
            h = net.icf[k] * Ibr[k] + net.Gv1[k] * V[f]
        else
            h = net.icf[k] * Ibr[k] + net.Gv1[k] * (V[f] * r - V[t])
            node_Ihis[t] += h
        end
        brch_Ihis[k] = h
        node_Ihis[f] -= h * r
    end
    return node_Ihis
end

"""
    add_conductance_block!(net, nodes, Gb)

Adds the small dense matrix `Gb` to G at `nodes` (a machine's Norton
conductance, or a fault admittance with `network_events = :refactorize`).
Call [`refactorize!`](@ref) afterwards.
"""
function add_conductance_block!(net::Network, nodes::AbstractVector{Int}, Gb::AbstractMatrix)
    for (a, na) in enumerate(nodes), (b, nb_) in enumerate(nodes)
        net.G[na, nb_] += Gb[a, b]
    end
    return net
end

"Computes the sparse LU factorization (UMFPACK) of G."
function refactorize!(net::Network)
    net.lu = lu(net.G)
    return net
end

"""
    solve_network!(net, I_inj, comp) -> V

Solves eq. (1), G v = i_s + i_h, for the node voltages of the current step.
`I_inj` holds the machine and IBR source currents; the history currents of
the companion models and of the Bergeron lines are added here.  With a
[`Compensation`](@ref), G is the unmodified network and the solution is
corrected for the faults and trips in effect (eq. (3)).
"""
function solve_network!(net::Network, I_inj::AbstractVector, comp)
    copyto!(net.V_prev, net.V)
    if isempty(net.ber)
        net.rhs .= I_inj .+ net.node_Ihis
    else
        net.rhs .= I_inj .+ net.node_Ihis .+ net.ber_node_Ihis
    end
    ldiv!(net.V, net.lu, net.rhs)
    comp === nothing || apply_compensation!(net.V, comp)
    return net.V
end

"""
    update_history!(net)

After the solution: the current of every branch from the solved voltages,
and its history source for the next step, summed per node.
"""
function update_history!(net::Network)
    V = net.V
    fill!(net.node_Ihis, 0.0)
    @inbounds for k in eachindex(net.F)
        f = net.F[k]; t = net.T[k]; r = net.tr[k]
        dV = t == 0 ? V[f] : V[f] * r - V[t]
        Ibr = dV / net.Req[k] + net.brch_Ihis[k]
        Ihis = net.icf[k] * Ibr + net.Gv1[k] * dV
        net.brch_Ipre[k] = Ibr
        net.brch_Ihis[k] = Ihis
        net.node_Ihis[f] -= Ihis * r
        t != 0 && (net.node_Ihis[t] += Ihis)
    end
    return net
end

"""
    network_admittance(pfd) -> SparseMatrixCSC{ComplexF64}

Bus admittance matrix [pu] of the network that `build_network` models: pi
lines, transformers with their ratios at the from end, fixed and switched
shunts (loads and units excluded), buses in the order of `pfd.bus_num`.
"""
function network_admittance(pfd::PowerFlowData)
    nb = nbus(pfd)
    I_ = Int[]; J_ = Int[]; Y_ = ComplexF64[]
    add(i, j, y) = (push!(I_, i); push!(J_, j); push!(Y_, y))
    for i in eachindex(pfd.line_from)
        f = bus_index(pfd, pfd.line_from[i]); t = bus_index(pfd, pfd.line_to[i])
        y = 1 / pfd.line_RX[i]; b = im * pfd.line_chg[i] / 2
        add(f, f, y + b); add(t, t, y + b); add(f, t, -y); add(t, f, -y)
    end
    for i in eachindex(pfd.xfmr_from)
        f = bus_index(pfd, pfd.xfmr_from[i]); t = bus_index(pfd, pfd.xfmr_to[i])
        y = 1 / pfd.xfmr_RX[i]; r = 1 / pfd.xfmr_k[i]
        add(f, f, y * r^2); add(t, t, y); add(f, t, -y * r); add(t, f, -y * r)
    end
    for i in eachindex(pfd.shnt_bus)
        k = bus_index(pfd, pfd.shnt_bus[i]); add(k, k, pfd.shnt_gb[i] / pfd.basemva)
    end
    for i in eachindex(pfd.shnt_sw_bus)
        k = bus_index(pfd, pfd.shnt_sw_bus[i])
        Ysw = conj(-pfd.shnt_sw_gb[i]) / pfd.basemva / pfd.bus_Vm[k]^2
        add(k, k, complex(real(Ysw), -imag(Ysw)))
    end
    return sparse(I_, J_, Y_, nb, nb)
end

"""
    power_flow_mismatch(pfd) -> Vector{ComplexF64}

Complex power mismatch [pu] at every bus of the network that `build_network`
models ([`network_admittance`](@ref), loads at their power-flow power, units
at their outputs) at the power-flow voltages.  It vanishes when the
power-flow solution belongs to the same network, so that the EMT simulation
starts from that operating point (Bergeron lines aside, whose steady state is
that of the distributed line).
"""
function power_flow_mismatch(pfd::PowerFlowData)
    nb = nbus(pfd)
    V = pfd.bus_Vm .* cis.(pfd.bus_Va)
    Y = network_admittance(pfd)
    S = zeros(ComplexF64, nb)
    for i in eachindex(pfd.gen_bus)
        S[bus_index(pfd, pfd.gen_bus[i])] += complex(pfd.gen_MW[i], pfd.gen_Mvar[i]) / pfd.basemva
    end
    for i in eachindex(pfd.ibr_bus)
        S[bus_index(pfd, pfd.ibr_bus[i])] += complex(pfd.ibr_MW[i], pfd.ibr_Mvar[i]) / pfd.basemva
    end
    for i in eachindex(pfd.load_bus)
        S[bus_index(pfd, pfd.load_bus[i])] -= complex(pfd.load_MW[i], pfd.load_Mvar[i]) / pfd.basemva
    end
    return S .- V .* conj.(Y * V)
end

"A copy of the network for a new scenario: shares every constant and `lu`, owns its own state."
function copy_state(src::Network, lu_shared)
    return Network(src.nbus, src.N, src.ZbaseA, src.IbaseA,
                   src.F, src.T, src.Req, src.icf, src.Gv1, src.tr, src.G, lu_shared,   # shared constants
                   copy(src.brch_Ihis), copy(src.brch_Ipre), copy(src.node_Ihis),
                   copy(src.V), copy(src.V_prev), similar(src.rhs),
                   copy_state(src.ber), copy(src.ber_node_Ihis),
                   src.Vt0, src.It0)
end
