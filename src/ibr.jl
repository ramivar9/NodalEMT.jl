# ==========================================================================
#  ibr.jl -- inverter-based resources (WECC generic renewable models)
#
#  An IBR is a controlled current source.  Three cascaded models are used:
#     REPCA  plant controller    : P/Q (or V/f) references -> P and Q commands
#     REECB  electrical control  : P and Q commands -> current commands Ipcmd, Iqcmd
#     REGCA  generator/converter : current commands -> injected current (with
#                                  rate limits and low-voltage logic)
#  The injected current is placed on the d/q axes of the bus PLL angle and
#  transformed to the three phases.  All models are integrated with the
#  forward-Euler rule.  Within one step, each model reads the PREVIOUS-step
#  outputs of the others (a one-step delay), as in ParaEMT.
#
#  Each model is a scalar function (`regca_step`, `reecb_step`,
#  `repca_step`); the parameters of one unit are passed as small immutable
#  structs (`RegcaP`, `ReecbP`, `RepcaP`).
# ==========================================================================

"Per-IBR constants."
struct IbrConstants
    n::Int
    bus::Vector{Int}                 # bus index
    ts::Float64
    basemva::Float64
    MVA_base::Vector{Float64}        # power-flow MVA base (used in the power measurement)
    inj_scale::Vector{Float64}       # converter-base current -> network-pu current
    Qgen0::Vector{Float64}           # initial reactive power (selects the sign of the Iq ramp limit)
    pfaref::Vector{Float64}          # power-factor angle reference (REECB, PFFLAG = 1)
    Vref0::Vector{Float64}           # REECB voltage reference (parameter, or the initial voltage)
    Pref_out::Vector{Float64}        # REECB power reference when REPCA does not regulate frequency (FFlag = 0)
    regca::RegcaParams
    reecb::ReecbParams
    repca::RepcaParams
end

"REGCA states."
mutable struct RegcaState
    s0::Vector{Float64}     # active current
    s1::Vector{Float64}     # reactive current (sign: -Iq)
    s2::Vector{Float64}     # filtered voltage for the low-voltage power logic
    Vmp::Vector{Float64}; Vap::Vector{Float64}   # measured voltage magnitude and angle
    i1::Vector{Float64}     # high-voltage reactive-current reduction
    i2::Vector{Float64}     # low-voltage active-current scaling
    ip2rr::Vector{Float64}  # memory of the active-current rate limiter
end
RegcaState(n::Int) = RegcaState((zeros(n) for _ in 1:8)...)

"REECB states."
mutable struct ReecbState
    s0::Vector{Float64}     # filtered terminal voltage
    s1::Vector{Float64}     # filtered active power
    s2::Vector{Float64}     # Q -> V PI output
    s3::Vector{Float64}     # V -> Iq PI output
    s4::Vector{Float64}     # Q/V -> Iq lag output
    s5::Vector{Float64}     # P-order lag output
    Ipcmd::Vector{Float64}; Iqcmd::Vector{Float64}
    Pref::Vector{Float64}; Qext::Vector{Float64}
    q2vPI::Vector{Float64}; v2iPI::Vector{Float64}   # PI error memories
end
ReecbState(n::Int) = ReecbState((zeros(n) for _ in 1:12)...)

"REPCA states."
mutable struct RepcaState
    s0::Vector{Float64}     # filtered regulated voltage
    s1::Vector{Float64}     # filtered reactive power
    s2::Vector{Float64}     # Q/V PI output
    s3::Vector{Float64}     # lead-lag output -> Qext of REECB
    s4::Vector{Float64}     # filtered active power
    s5::Vector{Float64}     # frequency-droop PI output
    s6::Vector{Float64}     # P-order lag output -> Pref of REECB
    Vref::Vector{Float64}; Qref::Vector{Float64}; Freq_ref::Vector{Float64}; Plant_pref::Vector{Float64}
    LineMW::Vector{Float64}; LineMvar::Vector{Float64}; LineMVA::Vector{Float64}
    QVdbout::Vector{Float64}; fdbout::Vector{Float64}
    vq2qPI::Vector{Float64}; p2pPI::Vector{Float64}
    Vf::Vector{Float64}; Pe::Vector{Float64}; Qe::Vector{Float64}   # measured frequency and power
end
RepcaState(n::Int) = RepcaState((zeros(n) for _ in 1:21)...)

"All IBR states."
mutable struct IbrState
    regca::RegcaState
    reecb::ReecbState
    repca::RepcaState
end
IbrState(n::Int) = IbrState(RegcaState(n), ReecbState(n), RepcaState(n))

for T in (RegcaState, ReecbState, RepcaState)
    @eval function Base.copyto!(dst::$T, src::$T)
        for f in fieldnames($T)
            copyto!(getfield(dst, f), getfield(src, f))
        end
        return dst
    end
end
Base.copyto!(dst::IbrState, src::IbrState) = (copyto!(dst.regca, src.regca); copyto!(dst.reecb, src.reecb);
                                             copyto!(dst.repca, src.repca); dst)

# --------------------------------------------------------------------------
# Constants and initialization
# --------------------------------------------------------------------------

function build_ibr_constants(pfd::PowerFlowData, dyd::DynamicData, net::Network, ts::Float64)
    n = nibr(pfd)
    bus = [bus_index(pfd, pfd.ibr_bus[i]) for i in 1:n]
    inj_scale = zeros(n)
    for i in 1:n
        base_es = pfd.bus_basekV[bus[i]] * sqrt(2 / 3) * 1000
        base_is = pfd.ibr_MVA_base[i] * 1e6 / (base_es * 3 / 2)
        Ibase = base_is / sqrt(2)
        inj_scale[i] = Ibase / (net.IbaseA[bus[i]] * 1000)
    end
    return IbrConstants(n, bus, ts, pfd.basemva, copy(pfd.ibr_MVA_base), inj_scale, zeros(n), zeros(n), zeros(n),
                        zeros(n), dyd.regca, dyd.reecb, dyd.repca)
end

"""
    init_ibr!(x, c, pfd; variant)

Steady-state initialization of REGCA, REECB and REPCA at the power-flow
dispatch of every IBR.  The updated model set (`variant = :updated`)
starts two quantities in equilibrium that ParaEMT does not:

1. the reactive-power PI integrators of REECB (`s2`, `s3`) and the memory
   `v2iPI` start at their steady-state values instead of zero;
2. a REPCA plant controller without frequency regulation (FFlag = 0) hands
   REECB the power-flow dispatch as power reference; ParaEMT resets
   it to zero by a coding error, so such an IBR ramps its power down.

`variant = :paraemt` reproduces ParaEMT.
"""
function init_ibr!(x::IbrState, c::IbrConstants, pfd::PowerFlowData; variant::Symbol)
    rg = c.regca; rb = c.reecb; rp = c.repca
    for i in 1:c.n
        k = c.bus[i]
        P = pfd.ibr_MW[i] / rg.MVAbase[i]
        Q = pfd.ibr_Mvar[i] / rg.MVAbase[i]
        Vm = pfd.bus_Vm[k]; Va = pfd.bus_Va[k]
        Vt = Vm * cis(Va)
        It = conj(complex(P, Q) / Vt)
        # REGCA
        Ip_out = real(It) * cos(Va) + imag(It) * sin(Va)
        Iq_out = imag(It) * cos(Va) - real(It) * sin(Va)
        i1 = max(0.0, (Vm - rg.Volim[i]) * rg.Khv[i])
        Iq = Iq_out + i1
        i2 = if Vm > rg.Lvpnt1[i]
            1.0
        elseif Vm < rg.Lvpnt0[i]
            error("voltage at bus $(pfd.bus_num[k]) too low to initialize IBR $i")
        else
            (Vm - rg.Lvpnt0[i]) / (rg.Lvpnt1[i] - rg.Lvpnt0[i])
        end
        Ip = Ip_out / i2
        x.regca.s0[i] = Ip; x.regca.s1[i] = -Iq; x.regca.s2[i] = Vm
        x.regca.Vmp[i] = Vm; x.regca.Vap[i] = Va; x.regca.i1[i] = i1; x.regca.i2[i] = i2; x.regca.ip2rr[i] = 0.0
        c.Qgen0[i] = Q
        # REECB
        x.reecb.s0[i] = Vm; x.reecb.s1[i] = P
        x.reecb.s4[i] = Q / Vm; x.reecb.s5[i] = P
        c.Vref0[i] = rb.Vref0[i] == 0.0 ? Vm : rb.Vref0[i]
        c.pfaref[i] = atan(Q / P)
        x.reecb.Ipcmd[i] = Ip; x.reecb.Iqcmd[i] = -Iq
        x.reecb.Pref[i] = Ip * Vm; x.reecb.Qext[i] = Q
        x.reecb.q2vPI[i] = 0.0
        if variant === :updated
            # back-solve the cascaded Q -> V -> Iq PI outputs at the power-flow point
            s2 = clamp(rb.VFLAG[i] == 1 ? Vm : Q, rb.Vmin[i], rb.Vmax[i])
            v2i = s2 - Vm
            Iqmax = rb.PQFLAG[i] == 0 ? rb.Imax[i] : sqrt(max(0.0, rb.Imax[i]^2 - Ip^2))
            s3 = v2i < 0 ? -rb.Imax[i] : (v2i > 0 ? Iqmax : (rb.QFLAG[i] == 1 ? -Iq : 0.0))
            x.reecb.s2[i] = s2; x.reecb.s3[i] = s3; x.reecb.v2iPI[i] = v2i
        else
            x.reecb.s2[i] = 0.0; x.reecb.s3[i] = 0.0; x.reecb.v2iPI[i] = 0.0
        end
        # REPCA (local measurement only)
        (rp.remote_bus[i] == 0 && rp.branch_From_bus[i] == 0 && rp.branch_To_bus[i] == 0) ||
            error("REPCA remote-bus or branch metering is not supported (IBR $i)")
        Pbranch = P; Qbranch = Q; Ibranch = It
        V1_in1 = abs(Vt + complex(rp.Rc[i], rp.Xc[i]) * Ibranch)
        V1_in0 = Qbranch * rp.Kc[i] + Vm
        V1 = rp.VCFlag[i] == 0 ? V1_in0 : V1_in1
        x.repca.s0[i] = V1; x.repca.s1[i] = Qbranch
        if rp.FFlag[i] == 0
            x.repca.s4[i] = 0.0; x.repca.s5[i] = 0.0; x.repca.s6[i] = 0.0
            c.Pref_out[i] = variant === :updated ? x.reecb.Pref[i] : 0.0
        else
            x.repca.s4[i] = Pbranch; x.repca.s5[i] = x.reecb.Pref[i]; x.repca.s6[i] = x.reecb.Pref[i]
            c.Pref_out[i] = 0.0
        end
        x.repca.s2[i] = Q; x.repca.s3[i] = Q
        x.repca.Vref[i] = V1; x.repca.Qref[i] = Qbranch; x.repca.Freq_ref[i] = 1.0; x.repca.Plant_pref[i] = Pbranch
        x.repca.LineMW[i] = Pbranch; x.repca.LineMvar[i] = Qbranch; x.repca.LineMVA[i] = hypot(P, Q)
        x.repca.QVdbout[i] = 0.0; x.repca.fdbout[i] = 0.0; x.repca.vq2qPI[i] = 0.0; x.repca.p2pPI[i] = 0.0
        x.repca.Vf[i] = 1.0; x.repca.Pe[i] = Pbranch; x.repca.Qe[i] = Qbranch
    end
    return x
end

# --------------------------------------------------------------------------
# One step of each model (scalar functions)
# --------------------------------------------------------------------------

"Dead band [d1, d2] as in the WECC models: zero inside, shifted outside."
@inline function deadband(v, d1, d2)
    if d1 <= v <= d2
        return zero(v)
    end
    v > d2 && (v -= d2)
    v < d1 && (v -= d1)
    return v
end

"""
    ibr_phase_currents(de, we, s0, s1, i1, i2, scale, ts) -> (ia, ib, ic)

Three-phase current of one IBR from its REGCA currents (ip = s0 i2,
iq = -s1 - i1), placed on the PLL angle advanced by one step.
"""
@inline function ibr_phase_currents(de, we, s0, s1, i1, i2, scale, ts)
    th = de + ts * we * (2pi * 60)
    ip = s0 * i2
    iq = -s1 - i1
    ia = (cos(th) * ip - sin(th) * iq) * scale
    ib = (cos(th - 2pi / 3) * ip - sin(th - 2pi / 3) * iq) * scale
    ic = (cos(th + 2pi / 3) * ip - sin(th + 2pi / 3) * iq) * scale
    return ia, ib, ic
end

"Active and reactive power from instantaneous phase voltages and currents (`sb` = base ratio)."
@inline function ibr_power(va, vb, vc, ia, ib, ic, sb)
    pe = (va * ia + vb * ib + vc * ic) * (2 / 3) * sb
    qe = ((vb - vc) * ia + (vc - va) * ib + (va - vb) * ic) / sqrt(3.0) * (2 / 3) * sb
    return pe, qe
end

"Parameters of one REGCA unit."
struct RegcaP
    Tg::Float64; Rrpwr::Float64; Brkpt::Float64; Zerox::Float64; Volim::Float64
    Lvpnt1::Float64; Lvpnt0::Float64; Tfltr::Float64; Khv::Float64; Iqrmax::Float64; Iqrmin::Float64
    Qgen0::Float64
end
@inline RegcaP(p::RegcaParams, Qgen0, i) = RegcaP(p.Tg[i], p.Rrpwr[i], p.Brkpt[i], p.Zerox[i], p.Volim[i],
    p.Lvpnt1[i], p.Lvpnt0[i], p.Tfltr[i], p.Khv[i], p.Iqrmax[i], p.Iqrmin[i], Qgen0[i])

"Parameters of one REECB unit."
struct ReecbP
    PFFLAG::Int; VFLAG::Int; QFLAG::Int; PQFLAG::Int
    Vdip::Float64; Vup::Float64; Trv::Float64; dbd1::Float64; dbd2::Float64; Kqv::Float64
    Iqhl::Float64; Iqll::Float64; Vref0::Float64; Tp::Float64; Qmax::Float64; Qmin::Float64
    Vmax::Float64; Vmin::Float64; Kqp::Float64; Kqi::Float64; Kvp::Float64; Kvi::Float64
    Tiq::Float64; dPmax::Float64; dPmin::Float64; Pmax::Float64; Pmin::Float64; Imax::Float64; Tpord::Float64
    pfaref::Float64; FFlag::Int; Pref_out::Float64
end
@inline ReecbP(p::ReecbParams, Vref0, pfaref, FFlag, Pref_out, i) = ReecbP(p.PFFLAG[i], p.VFLAG[i], p.QFLAG[i],
    p.PQFLAG[i], p.Vdip[i], p.Vup[i], p.Trv[i], p.dbd1[i], p.dbd2[i], p.Kqv[i], p.Iqhl[i], p.Iqll[i], Vref0[i],
    p.Tp[i], p.Qmax[i], p.Qmin[i], p.Vmax[i], p.Vmin[i], p.Kqp[i], p.Kqi[i], p.Kvp[i], p.Kvi[i], p.Tiq[i],
    p.dPmax[i], p.dPmin[i], p.Pmax[i], p.Pmin[i], p.Imax[i], p.Tpord[i], pfaref[i], FFlag[i], Pref_out[i])

"Parameters of one REPCA unit."
struct RepcaP
    VCFlag::Int; RefFlag::Int
    Tfltr::Float64; Kp::Float64; Ki::Float64; Tft::Float64; Tfv::Float64; Vfrz::Float64
    Rc::Float64; Xc::Float64; Kc::Float64; emax::Float64; emin::Float64; dbd1::Float64; dbd2::Float64
    Qmax::Float64; Qmin::Float64; Kpg::Float64; Kig::Float64; Tp::Float64; fdbd1::Float64; fdbd2::Float64
    femax::Float64; femin::Float64; Pmax::Float64; Pmin::Float64; Tg::Float64; Ddn::Float64; Dup::Float64
end
@inline RepcaP(p::RepcaParams, i) = RepcaP(p.VCFlag[i], p.RefFlag[i], p.Tfltr[i], p.Kp[i], p.Ki[i], p.Tft[i],
    p.Tfv[i], p.Vfrz[i], p.Rc[i], p.Xc[i], p.Kc[i], p.emax[i], p.emin[i], p.dbd1[i], p.dbd2[i], p.Qmax[i],
    p.Qmin[i], p.Kpg[i], p.Kig[i], p.Tp[i], p.fdbd1[i], p.fdbd2[i], p.femax[i], p.femin[i], p.Pmax[i],
    p.Pmin[i], p.Tg[i], p.Ddn[i], p.Dup[i])

"""
    regca_step(s0, s1, s2, Vmp, ip2rr, Ipcmd, Iqcmd, Vm, p::RegcaP, ts)
        -> (s0_new, s1_new, s2_new, i1, i2, ip2rr_new)

REGCA: high-voltage reactive-current reduction `i1`, low-voltage active-current
scaling `i2`, reactive-current lag with ramp limit (sign from the initial Q),
low-voltage power logic on the filtered voltage, and active-current lag with
the rate limit Rrpwr.
"""
@inline function regca_step(s0, s1, s2, Vmp, ip2rr, Ipcmd, Iqcmd, Vm, p::RegcaP, ts)
    i1 = max(0.0, (Vm - p.Volim) * p.Khv)
    i2 = Vm < p.Lvpnt0 ? 0.0 : (Vm > p.Lvpnt1 ? 1.0 : (Vm - p.Lvpnt0) / (p.Lvpnt1 - p.Lvpnt0))
    s1n = s1 + (Iqcmd - s1) / p.Tg * ts
    if p.Qgen0 > 0
        s1n > s1 + p.Iqrmax * ts && (s1n = s1 + p.Iqrmax * ts)
    else
        s1n < s1 + p.Iqrmin * ts && (s1n = s1 + p.Iqrmin * ts)
    end
    s2n = s2 + (Vmp - s2) / p.Tfltr * ts
    lvpl = s2n < p.Zerox ? 0.0 : (s2n > p.Brkpt ? 1.0 : (s2n - p.Zerox) / (p.Brkpt - p.Zerox))
    tin = Ipcmd - s0
    if s0 > 0
        (tin - ip2rr) / ts > p.Rrpwr && (tin = p.Rrpwr * ts + ip2rr)
    else
        (tin - ip2rr) / ts < -p.Rrpwr && (tin = -p.Rrpwr * ts + ip2rr)
    end
    s0n = s0 + tin / p.Tg * ts
    s0n > lvpl && (s0n = lvpl)
    return s0n, s1n, s2n, i1, i2, tin
end

"""
    reecb_step(s0..s5, Pref, q2vPI, v2iPI, Ipcmd, Iqcmd, repca_s3, repca_s6, Vm, pe, qe, p::ReecbP, ts)
        -> (s0..s5, Ipcmd, Iqcmd, Pref, Qext, q2vPI, v2iPI)

REECB: current limits (P or Q priority), freeze of the integrators during a
voltage dip, reactive path (voltage dead band -> Iq injection; Q or power
factor -> V PI -> Iq PI, or Q/V lag), active path (rate-limited Pref ->
P-order lag -> Ipcmd = Pord/V).
"""
@inline function reecb_step(s0, s1, s2, s3, s4, s5, Pref_1, q2vPI_1, v2iPI_1, Ipcmd_1, Iqcmd_1,
                            repca_s3_1, repca_s6_1, Vm, pe, qe, p::ReecbP, ts)
    if p.PQFLAG == 0            # Q priority
        Ipmax = sqrt(p.Imax^2 - Iqcmd_1^2); Ipmin = 0.0; Iqmax = p.Imax; Iqmin = -p.Imax
    else                        # P priority
        Ipmax = p.Imax; Ipmin = 0.0; Iqmax = sqrt(p.Imax^2 - Ipcmd_1^2); Iqmin = -p.Imax
    end
    dip = (Vm < p.Vdip) || (Vm > p.Vup)
    s0n = s0 + (Vm - s0) / p.Trv * ts
    v2 = deadband(p.Vref0 - s0n, p.dbd1, p.dbd2)
    Iqinj = clamp(v2 * p.Kqv, p.Iqll, p.Iqhl)
    s1n = s1 + (pe - s1) / p.Tp * ts
    Q0 = p.PFFLAG == 1 ? s1n * tan(p.pfaref) : repca_s3_1
    Q1 = Q0 - qe
    Q0 > p.Qmax && (Q1 = p.Qmax - qe)
    Q0 < p.Qmin && (Q1 = p.Qmin - qe)
    dQ1dt = (Q1 - q2vPI_1) / ts
    s2n = dip ? s2 : clamp(s2 + (dQ1dt * p.Kqp + Q1 * p.Kqi) * ts, p.Vmin, p.Vmax)
    V0 = clamp(p.VFLAG == 1 ? s2n : Q0, p.Vmin, p.Vmax)
    V1 = V0 - s0n
    dV1dt = (V1 - v2iPI_1) / ts
    s3n = dip ? s3 : clamp(s3 + (dV1dt * p.Kvp + V1 * p.Kvi) * ts, Iqmin, Iqmax)
    V3 = max(s0n, 0.01)
    s4n = dip ? s4 : s4 + (Q0 / V3 - s4) / p.Tiq * ts
    Iqcmd = clamp((p.QFLAG == 1 ? s3n : s4n) + Iqinj, Iqmin, Iqmax)
    Pref_in = p.FFlag == 1 ? repca_s6_1 : p.Pref_out   # FFlag = 0: fixed reference from the initialization
    dPdt = (Pref_in - Pref_1) / ts
    Pref = Pref_in
    dPdt > p.dPmax && (Pref = p.dPmax * ts + Pref_1)
    dPdt < p.dPmin && (Pref = p.dPmin * ts + Pref_1)
    s5n = dip ? s5 : clamp(s5 + (Pref - s5) / p.Tpord * ts, p.Pmin, p.Pmax)
    Ipcmd = clamp(s5n / V3, Ipmin, Ipmax)
    return s0n, s1n, s2n, s3n, s4n, s5n, Ipcmd, Iqcmd, Pref, repca_s3_1, Q1, V1
end

"""
    repca_step(s0..s6, Vref, Qref, Freq_ref, Plant_pref, LineMW, LineMvar, vq2qPI, p2pPI, Vm, Va, Vf, p::RepcaP, ts)
        -> (s0..s6, QVdbout, fdbout, vq2qPI, p2pPI)

REPCA: frequency-droop path (dead band, asymmetric droop, filtered line MW,
PI, P-order lag -> s6) and voltage/reactive path (line-drop compensation or
Q droop, filters, dead band, PI with voltage freeze, lead-lag -> s3).
"""
@inline function repca_step(s0, s1, s2, s3, s4, s5, s6, Vref, Qref, Freq_ref, Plant_pref, LineMW, LineMvar,
                            vq2qPI_1, p2pPI_1, Vm, Va, Vf, p::RepcaP, ts)
    f_err = deadband(Freq_ref - Vf, p.fdbd1, p.fdbd2)
    pfdroop = min(f_err * p.Ddn, 0.0) + max(f_err * p.Dup, 0.0)
    s4n = s4 + (LineMW - s4) / p.Tp * ts
    P1 = clamp(Plant_pref - s4n + pfdroop, p.femin, p.femax)
    dP1dt = (P1 - p2pPI_1) / ts
    s5n = clamp(s5 + (dP1dt * p.Kpg + P1 * p.Kig) * ts, p.Pmin, p.Pmax)
    s6n = s6 + (s5n - s6) / p.Tg * ts
    # voltage / reactive-power path (complex arithmetic written out)
    Vre = Vm * cos(Va); Vim = Vm * sin(Va)
    den = Vre * Vre + Vim * Vim
    Ire = (LineMW * Vre + LineMvar * Vim) / den          # Ibranch = conj(S / V)
    Iim = (LineMW * Vim - LineMvar * Vre) / den
    dre = p.Rc * Ire - p.Xc * Iim; dim = p.Rc * Iim + p.Xc * Ire
    V1_in1 = sqrt((Vre - dre)^2 + (Vim - dim)^2)
    V1_in0 = Vm + LineMvar * p.Kc
    V1 = p.VCFlag == 0 ? V1_in0 : V1_in1
    s0n = s0 + (V1 - s0) / p.Tfltr * ts
    s1n = s1 + (LineMvar - s1) / p.Tfltr * ts
    Q1_in = deadband(p.RefFlag == 0 ? Qref - s1n : Vref - s0n, p.dbd1, p.dbd2)
    Q1 = clamp(Q1_in, p.emin, p.emax)
    dQ1dt = (Q1 - vq2qPI_1) / ts
    s2n = Vm < p.Vfrz ? s2 : clamp(s2 + (dQ1dt * p.Kp + Q1 * p.Ki) * ts, p.Qmin, p.Qmax)
    ds2dt = (s2n - s2) / ts
    s3n = s3 + (p.Tft * ds2dt + s2 - s3) / p.Tfv * ts
    return s0n, s1n, s2n, s3n, s4n, s5n, s6n, Q1_in, f_err, Q1, P1
end

# --------------------------------------------------------------------------
# All IBRs: current injection before the solution, controls after it
# --------------------------------------------------------------------------

"""
    ibr_current!(I_inj, Iibr, x, c, bus, nbus)

Three-phase current of every IBR from its previous-step REGCA currents, on
the previous-step PLL angle, added to `I_inj`.  `Iibr` (3 per IBR) keeps the
currents for the power measurement after the solution.
"""
function ibr_current!(I_inj::Vector{Float64}, Iibr::Vector{Float64}, x::IbrState, c::IbrConstants,
                      bus::BusState, nbus::Int)
    g = x.regca
    @inbounds for i in 1:c.n
        k = c.bus[i]
        ia, ib, ic = ibr_phase_currents(bus.de[k], bus.we[k], g.s0[i], g.s1[i], g.i1[i], g.i2[i], c.inj_scale[i], c.ts)
        Iibr[3i - 2] = ia; Iibr[3i - 1] = ib; Iibr[3i] = ic
        I_inj[phase_node(k, 1, nbus)] += ia
        I_inj[phase_node(k, 2, nbus)] += ib
        I_inj[phase_node(k, 3, nbus)] += ic
    end
    return I_inj
end

"""
    update_ibr!(x, c, V, Iibr, bus, nbus)

One step of REGCA, REECB and REPCA for every IBR, with this step's node
voltages `V`, bus measurements `bus` and injected currents `Iibr`.
"""
function update_ibr!(x::IbrState, c::IbrConstants, V::Vector{Float64}, Iibr::Vector{Float64},
                     bus::BusState, nbus::Int)
    ts = c.ts
    g = x.regca; e = x.reecb; p = x.repca
    @inbounds for i in 1:c.n
        k = c.bus[i]
        Vm = bus.vtm[k]; Va = bus.de[k]; Vf = bus.we[k]
        pe, qe = ibr_power(V[k], V[k + nbus], V[k + 2nbus], Iibr[3i - 2], Iibr[3i - 1], Iibr[3i],
                           c.basemva / c.MVA_base[i])
        # previous-step outputs that couple the three models
        Ipcmd_1 = e.Ipcmd[i]; Iqcmd_1 = e.Iqcmd[i]; repca_s3_1 = p.s3[i]; repca_s6_1 = p.s6[i]

        g.s0[i], g.s1[i], g.s2[i], g.i1[i], g.i2[i], g.ip2rr[i] =
            regca_step(g.s0[i], g.s1[i], g.s2[i], g.Vmp[i], g.ip2rr[i], Ipcmd_1, Iqcmd_1, Vm,
                       RegcaP(c.regca, c.Qgen0, i), ts)
        g.Vmp[i] = Vm; g.Vap[i] = Va

        e.s0[i], e.s1[i], e.s2[i], e.s3[i], e.s4[i], e.s5[i], e.Ipcmd[i], e.Iqcmd[i], e.Pref[i], e.Qext[i],
        e.q2vPI[i], e.v2iPI[i] =
            reecb_step(e.s0[i], e.s1[i], e.s2[i], e.s3[i], e.s4[i], e.s5[i], e.Pref[i], e.q2vPI[i], e.v2iPI[i],
                       Ipcmd_1, Iqcmd_1, repca_s3_1, repca_s6_1, Vm, pe, qe,
                       ReecbP(c.reecb, c.Vref0, c.pfaref, c.repca.FFlag, c.Pref_out, i), ts)

        p.s0[i], p.s1[i], p.s2[i], p.s3[i], p.s4[i], p.s5[i], p.s6[i], p.QVdbout[i], p.fdbout[i], p.vq2qPI[i],
        p.p2pPI[i] =
            repca_step(p.s0[i], p.s1[i], p.s2[i], p.s3[i], p.s4[i], p.s5[i], p.s6[i], p.Vref[i], p.Qref[i],
                       p.Freq_ref[i], p.Plant_pref[i], p.LineMW[i], p.LineMvar[i], p.vq2qPI[i], p.p2pPI[i],
                       Vm, Va, Vf, RepcaP(c.repca, i), ts)
        p.LineMW[i] = pe; p.LineMvar[i] = qe; p.LineMVA[i] = hypot(pe, qe)
        p.Vf[i] = Vf; p.Pe[i] = pe; p.Qe[i] = qe
    end
    return x
end
