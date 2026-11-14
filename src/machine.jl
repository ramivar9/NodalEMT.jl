# ==========================================================================
#  machine.jl -- GENROU synchronous machine
#
#  The machine is modelled in its rotor (d, q) frame with three d-axis
#  circuits (stator d, field fd, damper 1d) and three q-axis circuits (stator
#  q, dampers 1q and 2q).  The voltage equations of one axis are
#
#        e = R i + dpsi/dt / ws + (speed voltage),        psi = L i
#
#  and are discretized with a damped trapezoidal rule (alpha = 99/101):
#
#        Rd1 i(t) = e(t) - his(t),      his(t) = -alpha e(t-ts) + Rd2 i(t-ts)
#
#  Eliminating the two rotor circuits gives the stator current as a Norton
#  source seen from the network,  i_stator = (e_mod - e_stator) / Rav,  where
#  e_mod collects the history terms and Rav = (Rd + Rq)/2.  The constant 3x3
#  conductance Gequiv of every machine is part of G; only the source current
#  changes from step to step.  It is computed BEFORE the network solution
#  from quantities extrapolated from the previous steps (`predict!`).
# ==========================================================================

"""
    MachineConstants

Per-generator equivalent-circuit constants and discretization matrices.
"""
struct MachineConstants
    n::Int
    bus::Vector{Int}                # bus index of each generator
    ws::Float64
    ts::Float64
    H::Vector{Float64}              # inertia constant [s]
    D::Vector{Float64}              # damping
    # equivalent circuit (machine base, pu)
    Lad::Vector{Float64}; Laq::Vector{Float64}; Ll::Vector{Float64}
    Ld::Vector{Float64};  Lq::Vector{Float64}
    Lffd::Vector{Float64}; Lf1d::Vector{Float64}; L11d::Vector{Float64}
    L11q::Vector{Float64}; L22q::Vector{Float64}
    Ra::Vector{Float64}; Rfd::Vector{Float64}
    # bases: peak phase current and impedance on the machine's own MVA base
    base_Is::Vector{Float64}; base_Zs::Vector{Float64}
    # discretization
    alpha::Vector{Float64}
    Rd1::Vector{Matrix{Float64}}; Rd2::Vector{Matrix{Float64}}       # d axis, 3x3
    Rq1::Vector{Matrix{Float64}}; Rq2::Vector{Matrix{Float64}}       # q axis, 3x3
    Rd1inv::Vector{Matrix{Float64}}; Rq1inv::Vector{Matrix{Float64}} # inverses of the 2x2 rotor blocks
    Rd_coe::Vector{Vector{Float64}}; Rq_coe::Vector{Vector{Float64}} # elimination coefficients
    Rd::Vector{Float64}; Rq::Vector{Float64}; Rav::Vector{Float64}   # stator Thevenin resistances
    Gequiv::Vector{Matrix{Float64}}  # 3x3 Norton conductance [network pu]
    EFD2efd::Vector{Float64}         # exciter output EFD -> field voltage efd = Rfd/Lad EFD
    inj_scale::Vector{Float64}       # machine-base current -> network-pu current
end

"""
    MachineState

Dynamic and algebraic variables of every machine.  `psi` are flux linkages,
`te` is ParaEMT's electrical torque (ed id + eq iq) and `qe` the reactive power.
"""
mutable struct MachineState
    delta::Vector{Float64}   # rotor angle [rad]
    omega::Vector{Float64}   # rotor speed [rad/s]
    id::Vector{Float64};  iq::Vector{Float64}
    ifd::Vector{Float64}; i1d::Vector{Float64}; i1q::Vector{Float64}; i2q::Vector{Float64}
    ed::Vector{Float64};  eq::Vector{Float64}
    psid::Vector{Float64}; psiq::Vector{Float64}
    psifd::Vector{Float64}; psi1q::Vector{Float64}; psi1d::Vector{Float64}; psi2q::Vector{Float64}
    te::Vector{Float64};  qe::Vector{Float64}
end
MachineState(n::Int) = MachineState((zeros(n) for _ in 1:18)...)

function Base.copyto!(dst::MachineState, src::MachineState)
    for f in fieldnames(MachineState)
        copyto!(getfield(dst, f), getfield(src, f))
    end
    return dst
end

# --------------------------------------------------------------------------
# Constants from the data sheet
# --------------------------------------------------------------------------

"""
    build_machine_constants(pfd, dyd, net, ts) -> MachineConstants

Converts the data-sheet parameters (Xd, Xd', Xd'', Td0', ...) into the
equivalent circuit and builds the discretization matrices.
"""
function build_machine_constants(pfd::PowerFlowData, dyd::DynamicData, net::Network, ts::Float64)
    g = dyd.genrou
    n = ngen(pfd)
    ws = pfd.ws
    bus = [bus_index(pfd, pfd.gen_bus[i]) for i in 1:n]

    # equivalent circuit (ParaEMT's ToEquiCirData)
    Lad = g.Xd .- g.Xl
    Laq = g.Xq .- g.Xl
    Ll = copy(g.Xl)
    Ld = Lad .+ Ll
    Lq = Laq .+ Ll
    Lfd = (g.Xdp .- g.Xl) .* Lad ./ (Lad .- (g.Xdp .- g.Xl))
    L1q = (g.Xqp .- g.Xl) .* Laq ./ (Laq .- (g.Xqp .- g.Xl))
    R1d = zeros(n); R2q = zeros(n)
    for i in 1:n
        z = g.Xdpp[i] - g.Xl[i]
        y = Lad[i] * Lfd[i] / (Lad[i] + Lfd[i])
        L1d = y * z / (y - z)
        R1d[i] = (y + L1d) / g.Td0pp[i]
        y = Laq[i] * L1q[i] / (Laq[i] + L1q[i])
        L2q = y * z / (y - z)
        R2q[i] = (y + L2q) / g.Tq0pp[i]
    end
    Rfd = (Lad .+ Lfd) ./ g.Td0p
    R1q = (Laq .+ L1q) ./ g.Tq0p
    Ra = copy(g.Ra)
    Lf1d = copy(Lad)
    Lffd = Lad .^ 2 ./ (g.Xd .- g.Xdp)
    L11d = Lad .^ 2 ./ (g.Xd .- g.Xdpp)
    L11q = Laq .^ 2 ./ (g.Xq .- g.Xqp)
    L22q = Laq .^ 2 ./ (g.Xq .- g.Xdpp)

    # bases: peak phase voltage and current on the machine's own MVA base
    base_es = [pfd.bus_basekV[bus[i]] * sqrt(2 / 3) * 1000 for i in 1:n]
    base_is = [pfd.gen_MVA_base[i] * 1e6 / (base_es[i] * 3 / 2) for i in 1:n]
    base_Is = base_is ./ sqrt(2)
    base_Zs = base_es ./ base_is

    # damped trapezoidal discretization
    alpha = fill(99.0 / 101.0, n)
    Rd1 = Vector{Matrix{Float64}}(undef, n); Rd2 = similar(Rd1); Rq1 = similar(Rd1); Rq2 = similar(Rd1)
    Rd1inv = similar(Rd1); Rq1inv = similar(Rd1)
    Rd_coe = Vector{Vector{Float64}}(undef, n); Rq_coe = similar(Rd_coe)
    Rd = zeros(n); Rq = zeros(n); Rav = zeros(n)
    Gequiv = similar(Rd1)
    for i in 1:n
        a = alpha[i]
        Ldm = [Lad[i]+Ll[i]  Lad[i]   Lad[i];
               Lad[i]        Lffd[i]  Lf1d[i];
               Lad[i]        Lf1d[i]  L11d[i]]
        Lqm = [Laq[i]+Ll[i]  Laq[i]   Laq[i];
               Laq[i]        L11q[i]  Laq[i];
               Laq[i]        Laq[i]   L22q[i]]
        Rdm = Diagonal([Ra[i], Rfd[i], R1d[i]])
        Rqm = Diagonal([Ra[i], R1q[i], R2q[i]])
        # the stator row is written with time in seconds (factor 1/ws on dpsi/dt),
        # the rotor rows with time in radians: hence the row scaling
        rowscale = [1 / (ts * ws), 1 / ts, 1 / ts]
        Rd1[i] = Rdm + (1 + a) .* rowscale .* Ldm
        Rd2[i] = a .* Rdm - (1 + a) .* rowscale .* Ldm
        Rq1[i] = Rqm + (1 + a) .* rowscale .* Lqm
        Rq2[i] = a .* Rqm - (1 + a) .* rowscale .* Lqm
        # eliminate the two rotor circuits: stator Thevenin resistance of each axis
        Rd1inv[i] = inv(Rd1[i][2:3, 2:3])
        Rd_coe[i] = vec(Rd1[i][1:1, 2:3] * Rd1inv[i])
        Rd[i] = Rd1[i][1, 1] - dot(Rd_coe[i], Rd1[i][2:3, 1])
        Rq1inv[i] = inv(Rq1[i][2:3, 2:3])
        Rq_coe[i] = vec(Rq1[i][1:1, 2:3] * Rq1inv[i])
        Rq[i] = Rq1[i][1, 1] - dot(Rq_coe[i], Rq1[i][2:3, 1])
        Rav[i] = (Rd[i] + Rq[i]) / 2
        # phase-domain Norton conductance: Rs on the diagonal, Rm off the diagonal
        R0 = Ra[i] + (1 + a) / (ts * ws) * g.X0[i]
        Rs = (R0 + 2Rav[i]) / 3
        Rm = (R0 - Rav[i]) / 3
        Requiv = [Rs Rm Rm; Rm Rs Rm; Rm Rm Rs]
        Gequiv[i] = inv(Requiv) ./ base_Zs[i] .* net.ZbaseA[bus[i]]   # machine ohms -> network pu
    end
    EFD2efd = Rfd ./ Lad
    inj_scale = [base_Is[i] / (net.IbaseA[bus[i]] * 1000) for i in 1:n]

    return MachineConstants(n, bus, ws, ts, copy(g.H), copy(g.D),
        Lad, Laq, Ll, Ld, Lq, Lffd, Lf1d, L11d, L11q, L22q, Ra, Rfd, base_Is, base_Zs,
        alpha, Rd1, Rd2, Rq1, Rq2, Rd1inv, Rq1inv, Rd_coe, Rq_coe, Rd, Rq, Rav, Gequiv, EFD2efd, inj_scale)
end

"Adds the Norton conductance of every machine to G."
function stamp_generators!(net::Network, mc::MachineConstants)
    for i in 1:mc.n
        nodes = [phase_node(mc.bus[i], ph, net.nbus) for ph in 1:3]
        add_conductance_block!(net, nodes, mc.Gequiv[i])
    end
    return net
end

# --------------------------------------------------------------------------
# Initialization from the power flow
# --------------------------------------------------------------------------

"""
    init_machines!(x, mc, pfd, net; variant) -> (EFD0, pm0)

Steady-state initialization of every machine at its power-flow dispatch.
Returns the field voltage `EFD0` and the mechanical power `pm0` with which
the exciters and governors are initialized.

A motoring machine (`gen_MW < 0`) has its rotor angle reflected through
`sign(gen_MW)`, and its stator current must be reflected too.  The updated
model set does that; `variant = :paraemt` reproduces ParaEMT, which
starts motoring machines with the stator current of a generator.
"""
function init_machines!(x::MachineState, mc::MachineConstants, pfd::PowerFlowData, net::Network; variant::Symbol)
    n = mc.n
    EFD0 = zeros(n); pm0 = zeros(n)
    for i in 1:n
        k = mc.bus[i]
        S = hypot(pfd.gen_MW[i], pfd.gen_Mvar[i])
        phi = S > 1e-10 ? asin(pfd.gen_Mvar[i] / S) : 0.0          # power-factor angle
        Ig = net.It0[i] * net.IbaseA[k] / (mc.base_Is[i] / 1000)   # terminal current, machine base
        Im = abs(Ig)
        Vm = abs(net.Vt0[k])
        # rotor angle relative to the terminal voltage (Xq behind Ra), then absolute
        dt = sign(pfd.gen_MW[i]) * atan((mc.Lq[i] * Im * cos(phi) - mc.Ra[i] * Im * sin(phi)) /
                                        (Vm + mc.Ra[i] * Im * cos(phi) + mc.Lq[i] * Im * sin(phi)))
        x.delta[i] = dt + pfd.bus_Va[k]
        x.omega[i] = mc.ws
        x.ed[i] = Vm * sin(dt)
        x.eq[i] = Vm * cos(dt)
        if variant === :updated && pfd.gen_MW[i] < 0
            x.id[i] = -Im * sin(dt - phi)
            x.iq[i] = -Im * cos(dt - phi)
        else
            x.id[i] = Im * sin(dt + phi)
            x.iq[i] = Im * cos(dt + phi)
        end
        x.i1d[i] = 0.0; x.i1q[i] = 0.0; x.i2q[i] = 0.0
        x.psid[i] = x.eq[i] + mc.Ra[i] * x.iq[i]
        x.psiq[i] = -(x.ed[i] + mc.Ra[i] * x.id[i])
        x.ifd[i] = (x.eq[i] + mc.Ld[i] * x.id[i] + mc.Ra[i] * x.iq[i]) / mc.Lad[i]
        x.psifd[i] = mc.Lffd[i] * x.ifd[i] - mc.Lad[i] * x.id[i]
        x.psi1d[i] = mc.Lad[i] * (x.ifd[i] - x.id[i])
        x.psi1q[i] = -mc.Laq[i] * x.iq[i]
        x.psi2q[i] = -mc.Laq[i] * x.iq[i]
        pm0[i] = x.psid[i] * x.iq[i] - x.psiq[i] * x.id[i]          # electrical = mechanical power
        x.te[i] = x.ed[i] * x.id[i] + x.eq[i] * x.iq[i]
        x.qe[i] = x.psid[i] * x.id[i] + x.psiq[i] * x.iq[i]
        EFD0[i] = mc.Rfd[i] * x.ifd[i] * mc.Lad[i] / mc.Rfd[i]     # EFD = Lad ifd (written as in ParaEMT)
    end
    return EFD0, pm0
end

"Largest residual of the steady-state voltage and flux equations after initialization (ParaEMT's CheckMacEq)."
function check_machine_equilibrium(x::MachineState, mc::MachineConstants)
    worst = 0.0
    for i in 1:mc.n
        r = (x.ed[i] + x.psiq[i] + mc.Ra[i] * x.id[i],
             x.eq[i] - x.psid[i] + mc.Ra[i] * x.iq[i],
             -(mc.Lad[i] + mc.Ll[i]) * x.id[i] + mc.Lad[i] * x.ifd[i] + mc.Lad[i] * x.i1d[i] - x.psid[i],
             -(mc.Laq[i] + mc.Ll[i]) * x.iq[i] + mc.Laq[i] * x.i1q[i] + mc.Laq[i] * x.i2q[i] - x.psiq[i],
             mc.Lffd[i] * x.ifd[i] + mc.Lf1d[i] * x.i1d[i] - mc.Lad[i] * x.id[i] - x.psifd[i],
             mc.Lf1d[i] * x.ifd[i] + mc.L11d[i] * x.i1d[i] - mc.Lad[i] * x.id[i] - x.psi1d[i],
             mc.L11q[i] * x.i1q[i] + mc.Laq[i] * x.i2q[i] - mc.Laq[i] * x.iq[i] - x.psi1q[i],
             mc.Laq[i] * x.i1q[i] + mc.L22q[i] * x.i2q[i] - mc.Laq[i] * x.iq[i] - x.psi2q[i])
        worst = max(worst, maximum(abs, r))
    end
    return worst
end

# --------------------------------------------------------------------------
# One step: predictor, Norton source, update after the network solution
# --------------------------------------------------------------------------

"""
    MachineWork

Values of the current step: the predicted quantities, the history terms of
the Norton source and the Park-transform angle, kept from the part of the
step before the network solution to the part after it.
"""
struct MachineWork
    pd_delta::Vector{Float64}; pd_omega::Vector{Float64}
    pd_id::Vector{Float64};    pd_iq::Vector{Float64}
    pd_EFD::Vector{Float64};   pd_ud::Vector{Float64}; pd_uq::Vector{Float64}
    his_fd::Vector{Float64};   his_1d::Vector{Float64}
    his_1q::Vector{Float64};   his_2q::Vector{Float64}
    ed_mod::Vector{Float64};   eq_mod::Vector{Float64}
    theta::Vector{Float64}
end
MachineWork(n::Int) = MachineWork((zeros(n) for _ in 1:fieldcount(MachineWork))...)

# speed voltages u_d = -psi_q omega/ws and u_q = psi_d omega/ws
@inline speed_voltage_d(x::MachineState, i, ws) = -x.psiq[i] * x.omega[i] / ws
@inline speed_voltage_q(x::MachineState, i, ws) = x.psid[i] * x.omega[i] / ws

"""
    predict!(w, mc, x1, x2, x3, EFD1, EFD2, EFD3, order)

Extrapolates the quantities the Norton source needs before the network
solution (speed, stator currents, field voltage, speed voltages) from the
last one, two or three steps (`order` = 1, 2, 3):
    order 1:  y = y1
    order 2:  y = 2 y1 - y2
    order 3:  y = 1.25 y1 + 0.5 y2 - 0.75 y3
The rotor angle is advanced from the predicted speed by the trapezoidal rule.
"""
function predict!(w::MachineWork, mc::MachineConstants, x1::MachineState, x2::MachineState, x3::MachineState,
                  EFD1::Vector{Float64}, EFD2::Vector{Float64}, EFD3::Vector{Float64}, order::Int)
    ws = mc.ws
    @inbounds for i in 1:mc.n
        if order == 1
            w.pd_omega[i] = x1.omega[i]; w.pd_id[i] = x1.id[i]; w.pd_iq[i] = x1.iq[i]; w.pd_EFD[i] = EFD1[i]
            w.pd_ud[i] = speed_voltage_d(x1, i, ws); w.pd_uq[i] = speed_voltage_q(x1, i, ws)
        elseif order == 2
            w.pd_omega[i] = 2 * x1.omega[i] - x2.omega[i]
            w.pd_id[i] = 2 * x1.id[i] - x2.id[i]
            w.pd_iq[i] = 2 * x1.iq[i] - x2.iq[i]
            w.pd_EFD[i] = 2 * EFD1[i] - EFD2[i]
            w.pd_ud[i] = 2 * speed_voltage_d(x1, i, ws) - speed_voltage_d(x2, i, ws)
            w.pd_uq[i] = 2 * speed_voltage_q(x1, i, ws) - speed_voltage_q(x2, i, ws)
        else
            w.pd_omega[i] = 1.25 * x1.omega[i] + 0.5 * x2.omega[i] - 0.75 * x3.omega[i]
            w.pd_id[i] = 1.25 * x1.id[i] + 0.5 * x2.id[i] - 0.75 * x3.id[i]
            w.pd_iq[i] = 1.25 * x1.iq[i] + 0.5 * x2.iq[i] - 0.75 * x3.iq[i]
            w.pd_EFD[i] = 1.25 * EFD1[i] + 0.5 * EFD2[i] - 0.75 * EFD3[i]
            w.pd_ud[i] = 1.25 * speed_voltage_d(x1, i, ws) + 0.5 * speed_voltage_d(x2, i, ws) - 0.75 * speed_voltage_d(x3, i, ws)
            w.pd_uq[i] = 1.25 * speed_voltage_q(x1, i, ws) + 0.5 * speed_voltage_q(x2, i, ws) - 0.75 * speed_voltage_q(x3, i, ws)
        end
        w.pd_delta[i] = x1.delta[i] + mc.ts * 0.5 * (x1.omega[i] + w.pd_omega[i])
        w.theta[i] = w.pd_delta[i] - pi / 2          # Park-transform angle of the d axis
    end
    return w
end

"""
    norton_current!(I_inj, w, mc, x, EFD, nbus, tripped)

History terms of the discretized machine equations and the resulting Norton
source current of every machine, transformed to the three phases and added
to `I_inj` (network pu).  Tripped machines inject nothing.
"""
function norton_current!(I_inj::Vector{Float64}, w::MachineWork, mc::MachineConstants, x::MachineState,
                         EFD::Vector{Float64}, nbus::Int, tripped::AbstractVector{Bool})
    ws = mc.ws
    @inbounds for i in 1:mc.n
        tripped[i] && continue
        a = mc.alpha[i]
        Rd2 = mc.Rd2[i]; Rq2 = mc.Rq2[i]
        efd = EFD[i] * mc.EFD2efd[i]
        ud = speed_voltage_d(x, i, ws); uq = speed_voltage_q(x, i, ws)
        # previous-step currents of the d and q axes (the stator current enters as -id, -iq)
        id_v = (-x.id[i], x.ifd[i], x.i1d[i])
        iq_v = (-x.iq[i], x.i1q[i], x.i2q[i])
        # history terms his = -alpha e_old + Rd2 i_old of the three rows of each axis
        his_d = -a * x.ed[i] + a * ud + Rd2[1, 1] * id_v[1] + Rd2[1, 2] * id_v[2] + Rd2[1, 3] * id_v[3]
        his_fd = -a * efd + Rd2[2, 1] * id_v[1] + Rd2[2, 2] * id_v[2] + Rd2[2, 3] * id_v[3]
        his_1d = Rd2[3, 1] * id_v[1] + Rd2[3, 2] * id_v[2] + Rd2[3, 3] * id_v[3]
        his_q = -a * x.eq[i] + a * uq + Rq2[1, 1] * iq_v[1] + Rq2[1, 2] * iq_v[2] + Rq2[1, 3] * iq_v[3]
        his_1q = Rq2[2, 1] * iq_v[1] + Rq2[2, 2] * iq_v[2] + Rq2[2, 3] * iq_v[3]
        his_2q = Rq2[3, 1] * iq_v[1] + Rq2[3, 2] * iq_v[2] + Rq2[3, 3] * iq_v[3]
        # reduce the rotor rows into the stator row (Thevenin form)
        his_red_d = his_d - (mc.Rd_coe[i][1] * (his_fd - w.pd_EFD[i] * mc.EFD2efd[i]) + mc.Rd_coe[i][2] * his_1d)
        his_red_q = his_q - (mc.Rq_coe[i][1] * his_1q + mc.Rq_coe[i][2] * his_2q)
        # Norton source; the d/q asymmetry (Rd - Rq)/2 is moved into the source with the predicted currents
        ed_mod = w.pd_ud[i] + his_red_d - (mc.Rd[i] - mc.Rq[i]) / 2 * w.pd_id[i]
        eq_mod = w.pd_uq[i] + his_red_q + (mc.Rd[i] - mc.Rq[i]) / 2 * w.pd_iq[i]
        w.his_fd[i] = his_fd; w.his_1d[i] = his_1d; w.his_1q[i] = his_1q; w.his_2q[i] = his_2q
        w.ed_mod[i] = ed_mod; w.eq_mod[i] = eq_mod
        id_src = ed_mod / mc.Rav[i]
        iq_src = eq_mod / mc.Rav[i]
        # inverse Park transform (d, q) -> (a, b, c) with the predicted rotor angle
        th = w.theta[i]
        k = mc.bus[i]
        for (ph, shift) in ((1, 0.0), (2, -2pi / 3), (3, 2pi / 3))
            i_ph = cos(th + shift) * id_src - sin(th + shift) * iq_src
            I_inj[phase_node(k, ph, nbus)] += i_ph * mc.inj_scale[i]
        end
    end
    return I_inj
end

"""
    machine_algebraic!(x, w, mc, V, nbus, i) -> pe

After the network solution: transforms the terminal voltage of machine `i`
to the (d, q) frame, solves the discretized machine equations for all
currents and flux linkages, stores them in `x` and returns the electrical
power.
"""
@inline function machine_algebraic!(x::MachineState, w::MachineWork, mc::MachineConstants, V::Vector{Float64},
                                    nbus::Int, i::Int)
    th = w.theta[i]
    k = mc.bus[i]
    va = V[phase_node(k, 1, nbus)]; vb = V[phase_node(k, 2, nbus)]; vc = V[phase_node(k, 3, nbus)]
    ed = 2 / 3 * (cos(th) * va + cos(th - 2pi / 3) * vb + cos(th + 2pi / 3) * vc)
    eq = -2 / 3 * (sin(th) * va + sin(th - 2pi / 3) * vb + sin(th + 2pi / 3) * vc)
    id = (w.ed_mod[i] - ed) / mc.Rav[i]
    iq = (w.eq_mod[i] - eq) / mc.Rav[i]
    # rotor currents from the two remaining rows of Rd1 i = e - his
    Rd1 = mc.Rd1[i]; Rq1 = mc.Rq1[i]; Rd1inv = mc.Rd1inv[i]; Rq1inv = mc.Rq1inv[i]
    v1 = w.pd_EFD[i] * mc.EFD2efd[i] - w.his_fd[i] + Rd1[2, 1] * id
    v2 = -w.his_1d[i] + Rd1[3, 1] * id
    ifd = Rd1inv[1, 1] * v1 + Rd1inv[1, 2] * v2
    i1d = Rd1inv[2, 1] * v1 + Rd1inv[2, 2] * v2
    v1q = -w.his_1q[i] + Rq1[2, 1] * iq
    v2q = -w.his_2q[i] + Rq1[3, 1] * iq
    i1q = Rq1inv[1, 1] * v1q + Rq1inv[1, 2] * v2q
    i2q = Rq1inv[2, 1] * v1q + Rq1inv[2, 2] * v2q
    # flux linkages psi = L i
    psid = -(mc.Lad[i] + mc.Ll[i]) * id + mc.Lad[i] * ifd + mc.Lad[i] * i1d
    psiq = -(mc.Laq[i] + mc.Ll[i]) * iq + mc.Laq[i] * i1q + mc.Laq[i] * i2q
    psifd = -mc.Lad[i] * id + mc.Lffd[i] * ifd + mc.Lf1d[i] * i1d
    psi1q = -mc.Laq[i] * iq + mc.L11q[i] * i1q + mc.Laq[i] * i2q
    psi1d = -mc.Lad[i] * id + mc.Lf1d[i] * ifd + mc.L11d[i] * i1d
    psi2q = -mc.Laq[i] * iq + mc.Laq[i] * i1q + mc.L22q[i] * i2q
    pe = psid * iq - psiq * id
    x.id[i] = id; x.iq[i] = iq; x.ifd[i] = ifd; x.i1d[i] = i1d; x.i1q[i] = i1q; x.i2q[i] = i2q
    x.ed[i] = ed; x.eq[i] = eq
    x.psid[i] = psid; x.psiq[i] = psiq; x.psifd[i] = psifd; x.psi1q[i] = psi1q; x.psi1d[i] = psi1d; x.psi2q[i] = psi2q
    x.te[i] = ed * id + eq * iq
    x.qe[i] = psid * id + psiq * iq
    return pe
end
