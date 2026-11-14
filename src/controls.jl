# ==========================================================================
#  controls.jl -- swing equation, SEXS exciter, TGOV1/HYGOV/GAST governors
#                 and IEEEST power system stabilizer
#
#  All controllers are integrated with the forward-Euler rule
#        x(t + ts) = x(t) + ts f(x(t), u(t)),
#  which is accurate because the EMT time step (tens of microseconds) is far
#  smaller than any controller time constant.  Limits are applied after the
#  update.
#
#  Every model is a scalar function `*_step(...)` that maps the states of one
#  unit at time t, its inputs and its parameters to the states at t + ts.
#  The `update_*!` wrappers at the end apply it to one unit of the state
#  vectors.
# ==========================================================================

# --------------------------------------------------------------------------
# States (one vector entry per unit)
# --------------------------------------------------------------------------

"SEXS exciter: lead-lag output `v1`, field voltage `EFD`, and the voltage reference of the previous step."
mutable struct SexsState
    v1::Vector{Float64}
    EFD::Vector{Float64}
    vref_prev::Vector{Float64}
end
SexsState(n::Int) = SexsState(zeros(n), zeros(n), zeros(n))

"TGOV1: valve position `p1`, turbine lead-lag output `p2`, mechanical power `pm`."
mutable struct Tgov1State
    p1::Vector{Float64}; p2::Vector{Float64}; pm::Vector{Float64}
end
Tgov1State(n::Int) = Tgov1State(zeros(n), zeros(n), zeros(n))

"HYGOV: filter `xe`, gate command `xc`, gate `xg`, water flow `xq`, mechanical power `pm`."
mutable struct HygovState
    xe::Vector{Float64}; xc::Vector{Float64}; xg::Vector{Float64}; xq::Vector{Float64}; pm::Vector{Float64}
end
HygovState(n::Int) = HygovState(zeros(n), zeros(n), zeros(n), zeros(n), zeros(n))

"GAST: fuel valve `p1`, fuel system `p2`, exhaust temperature `p3`, mechanical power `pm`."
mutable struct GastState
    p1::Vector{Float64}; p2::Vector{Float64}; p3::Vector{Float64}; pm::Vector{Float64}
end
GastState(n::Int) = GastState(zeros(n), zeros(n), zeros(n), zeros(n))

"IEEEST: filter stages `y1`..`y7`, the last two inputs `x1`, `x2`, and the output `vs`."
mutable struct IeeestState
    y1::Vector{Float64}; y2::Vector{Float64}; y3::Vector{Float64}; y4::Vector{Float64}
    y5::Vector{Float64}; y6::Vector{Float64}; y7::Vector{Float64}
    x1::Vector{Float64}; x2::Vector{Float64}; vs::Vector{Float64}
end
IeeestState(n::Int) = IeeestState((zeros(n) for _ in 1:10)...)

for T in (SexsState, Tgov1State, HygovState, GastState, IeeestState)
    @eval function Base.copyto!(dst::$T, src::$T)
        for f in fieldnames($T)
            copyto!(getfield(dst, f), getfield(src, f))
        end
        return dst
    end
end

# --------------------------------------------------------------------------
# Initialization (steady state: every derivative zero)
# --------------------------------------------------------------------------

"SEXS: with zero error, v1 = EFD/K and vref = v1 + Vt.  Returns the voltage references."
function init_sexs!(s::SexsState, p::SexsParams, EFD0::Vector{Float64}, Vt0::Vector{Float64})
    for i in eachindex(EFD0)
        s.EFD[i] = EFD0[i]
        s.v1[i] = EFD0[i] / p.K[i]
        s.vref_prev[i] = s.v1[i] + Vt0[i]
    end
    return s.v1 .+ Vt0
end

"""
    init_governors!(t, h, g, dyd, pm0) -> gref

Initializes every governor at the mechanical power `pm0` of its generator and
returns the reference `gref` that holds it there.  TGOV1 and HYGOV store the
reference multiplied by their droop R; GAST stores the power reference
(ParaEMT's convention).
"""
function init_governors!(t::Tgov1State, h::HygovState, g::GastState, dyd::DynamicData, pm0::Vector{Float64})
    gref = zeros(length(pm0))
    for (k, i) in enumerate(dyd.tgov1.gen)
        t.p1[k] = pm0[i]; t.p2[k] = pm0[i]; t.pm[k] = pm0[i]
        gref[i] = pm0[i] * dyd.tgov1.R[k]
    end
    for (k, i) in enumerate(dyd.hygov.gen)
        q0 = pm0[i] / dyd.hygov.At[k] + dyd.hygov.qNL[k]
        h.xe[k] = 0.0; h.xc[k] = q0; h.xg[k] = q0; h.xq[k] = q0; h.pm[k] = pm0[i]
        gref[i] = q0 * dyd.hygov.R[k]
    end
    for (k, i) in enumerate(dyd.gast.gen)
        g.p1[k] = pm0[i]; g.p2[k] = pm0[i]; g.p3[k] = pm0[i]; g.pm[k] = pm0[i]
        gref[i] = pm0[i]
    end
    return gref
end

# --------------------------------------------------------------------------
# One step of each model (scalar functions)
# --------------------------------------------------------------------------

"""
    swing_step(w, delta, pe, pm, H, D, ws, ts) -> (w_new, delta_new)

Swing equation 2H dw/dt = ws (pm/(w/ws) - pe) - D (w - ws): forward Euler for
the speed, trapezoidal rule for the angle.  Speeds in rad/s.
"""
@inline function swing_step(w, delta, pe, pm, H, D, ws, ts)
    w_new = w + (ws * (pm / (w / ws) - pe) - D * (w - ws)) / H / 2 * ts
    delta_new = delta + (w_new + w) / 2 * ts
    return w_new, delta_new
end

"""
    sexs_step(v1, EFD, vref, vref_prev, vt, dvt, vs, dvs, TA, TB, K, TE, Emin, Emax, ts) -> (v1_new, EFD_new)

SEXS: error vref - vt + vs -> lead-lag (1 + s TA)/(1 + s TB) -> K/(1 + s TE)
-> EFD limited to [Emin, Emax].  `vt`, `dvt` are the measured terminal
voltage and its derivative, `vs`, `dvs` the stabilizer output and its
derivative.  The lead term acts on the derivative of the whole error,
dvref + dvs - dvt.  ParaEMT omits dvs, so that the stabilizer
signal loses most of its phase advance; `dvs = 0` reproduces it.
"""
@inline function sexs_step(v1, EFD, vref, vref_prev, vt, dvt, vs, dvs, TA, TB, K, TE, Emin, Emax, ts)
    dvref = (vref - vref_prev) / ts
    v1_new = v1 + ((vref - vt - v1 + vs) / TB + TA / TB * (dvref + dvs - dvt)) * ts
    EFD_new = TE == 0 ? K * v1 : EFD + (K * v1 - EFD) * ts / TE
    return v1_new, clamp(EFD_new, Emin, Emax)
end

"""
    tgov1_step(p1, p2, gref, dw, R, T1, T2, T3, Vmin, Vmax, Dt, ts) -> (p1_new, p2_new, pm)

TGOV1: droop 1/R on (gref - dw), valve 1/(1 + s T1) limited to [Vmin, Vmax],
turbine lead-lag (1 + s T2)/(1 + s T3), pm = p2 - Dt dw.  `dw` is the speed
deviation omega/ws - 1.
"""
@inline function tgov1_step(p1, p2, gref, dw, R, T1, T2, T3, Vmin, Vmax, Dt, ts)
    p1_new = clamp(p1 + ((gref - dw) / R - p1) / T1 * ts, Vmin, Vmax)
    dp1dt = (p1_new - p1) / ts
    p2_new = p2 + (p1 - p2 + dp1dt * T2) / T3 * ts
    return p1_new, p2_new, p2_new - Dt * dw
end

"""
    hygov_step(xe, xc, xg, xq, gref, dw, R, r, Tr, Tf, Tg, VELM, GMAX, GMIN, TW, At, Dturb, qNL, ts)
        -> (xe_new, xc_new, xg_new, xq_new, pm)

HYGOV: filter Tf on the droop error, transient-droop PI (r, Tr) producing the
gate command with rate limit VELM and limits [GMIN, GMAX], gate servo Tg,
penstock water time constant TW, turbine power At (q - qNL)(q/g)^2 - Dturb dw g.
"""
@inline function hygov_step(xe, xc, xg, xq, gref, dw, R, r, Tr, Tf, Tg, VELM, GMAX, GMIN, TW, At, Dturb, qNL, ts)
    n1 = gref - (xc * R + dw)
    xe_new = xe + (n1 - xe) / Tf * ts
    xc_new = xc + ((xe_new - xe) / ts * Tr + xe) / r / Tr * ts
    dxcdt = (xc_new - xc) / ts
    dxcdt > VELM && (xc_new = VELM * ts + xc)
    dxcdt < -VELM && (xc_new = -VELM * ts + xc)
    xc_new = clamp(xc_new, GMIN, GMAX)
    xg_new = xg + (xc - xg) / Tg * ts
    ratio = xq / xg
    xq_new = xq + (1 - ratio * ratio) / TW * ts
    pm = (xq - qNL) * ratio * ratio * At - Dturb * dw * xg
    return xe_new, xc_new, xg_new, xq_new, pm
end

"""
    gast_step(p1, p2, p3, gref, dw, R, T1, T2, T3, LdLmt, KT, VMAX, VMIN, Dturb, ts) -> (p1_new, p2_new, p3_new, pm)

GAST: fuel demand min(gref - dw/R, temperature limit), valve 1/(1 + s T1)
limited to [VMIN, VMAX], fuel system 1/(1 + s T2), exhaust-temperature
measurement 1/(1 + s T3) feeding the limit, pm = p2 - Dturb dw.
"""
@inline function gast_step(p1, p2, p3, gref, dw, R, T1, T2, T3, LdLmt, KT, VMAX, VMIN, Dturb, ts)
    p3_new = p3 + (p2 - p3) / T3 * ts
    p2_new = p2 + (p1 - p2) / T2 * ts
    pl = min(gref - dw / R, (LdLmt - p3_new) * KT + LdLmt)
    p1_new = clamp(p1 + (pl - p1) / T1 * ts, VMIN, VMAX)
    return p1_new, p2_new, p3_new, p2_new - Dturb * dw
end

"""
    ieeest_step(y1..y7, x1, x2, dw, A1..A6, T1..T6, KS, LSMAX, LSMIN, VCU, ts) -> (y1..y7, x1, x2, vs)

IEEEST stabilizer with the speed deviation `dw` as input: two second-order
filters (A1..A6), two lead-lag stages (T1/T2, T3/T4), the washout and gain
stage (KS, T5, T6), the output limits [LSMIN, LSMAX] and the VCU cut-off.
Returns the new states and the signal `vs` added to the exciter error.
"""
@inline function ieeest_step(y1, y2, y3, y4, y5, y6, y7, x1, x2, dw,
                             A1, A2, A3, A4, A5, A6, T1, T2, T3, T4, T5, T6, KS, LSMAX, LSMIN, VCU, ts)
    dxdt = (dw - x1) / ts
    dxdt_1 = (x1 - x2) / ts
    dx2dt2 = (dxdt - dxdt_1) / ts
    u = dw + A5 * dxdt + A6 * dx2dt2
    # stage 1: (1 + A5 s + A6 s^2) / (1 + A1 s + A2 s^2)
    if A2 == 0
        y1n = 0.0
        y2n = A1 == 0 ? u : y2 + (u - y2) / A1 * ts
    else
        y2n = y2 + y1 * ts
        y1n = y1 + (u - y2 - A1 * y1) / A2 * ts
    end
    # stage 2: 1 / (1 + A3 s + A4 s^2)
    if A4 == 0
        y3n = 0.0
        y4n = A3 == 0 ? y2n : y4 + (y2 - y4) / A3 * ts
    else
        y4n = y4 + y3 * ts
        y3n = y3 + (y2 - y4 - A3 * y3) / A4 * ts
    end
    # stage 3: (1 + s T1) / (1 + s T2)
    dy4dt = (y4n - y4) / ts
    y5n = T2 == 0 ? y4n + dy4dt * T1 : y5 + (y4 - y5) / T2 * ts + (y4n - y4) * T1 / T2
    # stage 4: (1 + s T3) / (1 + s T4)
    dy5dt = (y5n - y5) / ts
    y6n = T4 == 0 ? y5n - dy5dt * T3 : y6 + (y5 - y6 + T3 * dy5dt) / T4 * ts
    # stage 5: KS s T5 / (1 + s T6) (washout), or KS / (1 + s T6) when T5 = 0
    y7n = T5 == 0 ? y7 + (KS * y6 - y7) / T6 * ts : y7 + (KS * y6 - y7 / T5) / T6 * T5 * ts
    # output: limited derivative of y7, then the VCU cut-off
    dy7dt = (y7n - y7) / ts
    vss = clamp(dy7dt, LSMIN, LSMAX)
    vs = (VCU == 0 || vss <= VCU) ? vss : 0.0
    return y1n, y2n, y3n, y4n, y5n, y6n, y7n, dw, x1, vs
end

# --------------------------------------------------------------------------
# Update one unit in the state vectors
# --------------------------------------------------------------------------

@inline function update_sexs!(s::SexsState, p::SexsParams, i::Int, vref, vt, dvt, vs, dvs, ts)
    v1n, EFDn = sexs_step(s.v1[i], s.EFD[i], vref, s.vref_prev[i], vt, dvt, vs, dvs,
                          p.TA[i], p.TB[i], p.K[i], p.TE[i], p.Emin[i], p.Emax[i], ts)
    s.vref_prev[i] = vref
    s.v1[i] = v1n; s.EFD[i] = EFDn
    return EFDn
end

@inline function update_tgov1!(t::Tgov1State, p::Tgov1Params, k::Int, gref, dw, ts)
    t.p1[k], t.p2[k], t.pm[k] = tgov1_step(t.p1[k], t.p2[k], gref, dw, p.R[k], p.T1[k], p.T2[k], p.T3[k],
                                           p.Vmin[k], p.Vmax[k], p.Dt[k], ts)
    return t.pm[k]
end

@inline function update_hygov!(h::HygovState, p::HygovParams, k::Int, gref, dw, ts)
    h.xe[k], h.xc[k], h.xg[k], h.xq[k], h.pm[k] = hygov_step(h.xe[k], h.xc[k], h.xg[k], h.xq[k], gref, dw,
        p.R[k], p.r[k], p.Tr[k], p.Tf[k], p.Tg[k], p.VELM[k], p.GMAX[k], p.GMIN[k], p.TW[k], p.At[k],
        p.Dturb[k], p.qNL[k], ts)
    return h.pm[k]
end

@inline function update_gast!(g::GastState, p::GastParams, k::Int, gref, dw, ts)
    g.p1[k], g.p2[k], g.p3[k], g.pm[k] = gast_step(g.p1[k], g.p2[k], g.p3[k], gref, dw, p.R[k], p.T1[k],
        p.T2[k], p.T3[k], p.LdLmt[k], p.KT[k], p.VMAX[k], p.VMIN[k], p.Dturb[k], ts)
    return g.pm[k]
end

@inline function update_ieeest!(s::IeeestState, p::IeeestParams, k::Int, dw, ts)
    s.y1[k], s.y2[k], s.y3[k], s.y4[k], s.y5[k], s.y6[k], s.y7[k], s.x1[k], s.x2[k], s.vs[k] =
        ieeest_step(s.y1[k], s.y2[k], s.y3[k], s.y4[k], s.y5[k], s.y6[k], s.y7[k], s.x1[k], s.x2[k], dw,
                    p.A1[k], p.A2[k], p.A3[k], p.A4[k], p.A5[k], p.A6[k], p.T1[k], p.T2[k], p.T3[k],
                    p.T4[k], p.T5[k], p.T6[k], p.KS[k], p.LSMAX[k], p.LSMIN[k], p.VCU[k], ts)
    return s.vs[k]
end
