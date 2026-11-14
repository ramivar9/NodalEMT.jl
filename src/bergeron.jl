# ==========================================================================
#  bergeron.jl -- Bergeron travelling-wave lines (paper, Section IV)
#
#  A pi section has no propagation delay.  A line whose travel time
#
#        tau = sqrt(X B) / ws                                        [s]
#
#  exceeds the time step is therefore modelled as the lossy Bergeron line of
#  EMTP: two lossless halves of surge impedance Zc = sqrt(X/B), with R/4 at
#  both terminals and R/2 at the midpoint.  Eliminating the midpoint gives,
#  at the sending terminal s (paper, eq. (4)),
#
#        i_s(t)   = v_s(t)/Z - Ih_s(t)
#        Ih_s(t)  = (1+h)/2 * H_r(t - tau)  +  (1-h)/2 * H_s(t - tau)
#        H_s(t)   = v_s(t)/Z + h i_s(t)
#
#  with Z = Zc + R/4 and h = (Zc - R/4)/Z; the receiving terminal r follows
#  by exchanging s and r.  (1+h)/2 is the part of a wave transmitted through
#  the midpoint resistance, (1-h)/2 the part reflected by it.  Each terminal
#  is a resistance Z in parallel with a history current that is known before
#  the network is solved, so a Bergeron line adds only 1/Z to two diagonal
#  entries of G and decouples the buses at its two ends.
#
#  Ring buffers and fractional delay (eq. (5)).  With tau/ts = m + phi,
#  m = floor(tau/ts) and 0 <= phi < 1, the delayed history is interpolated
#  between the samples at t_{n-m} and t_{n-m-1}:
#
#        H(t_n - tau) ~ (1 - phi) H(t_{n-m}) + phi H(t_{n-m-1})
#
#  Every conductor end keeps its histories in a ring buffer of length
#  L = ceil(tau/ts).  `head` points at the newest sample, that of t_{n-1}, so
#  the two samples above lie m-1 and m positions behind it.  After the
#  network solution the new H(t_n) overwrites the oldest sample and `head`
#  advances.  The rings of all conductors are packed into one vector per
#  line end (Hs for the sending ends, Hr for the receiving ends) with
#  per-conductor offsets, so that lines of different lengths are advanced by
#  one loop.
#
#  Exact periodic initialization (eqs. (6)-(8)).  At t = 0 the rings must
#  already hold one travel time of history.  The interpolation realizes the
#  delay D = exp(-j ws m ts) [(1 - phi) + phi exp(-j ws ts)] instead of
#  exp(-j ws tau).  Using D in the phasor form of eq. (4) gives the two-port
#  of the DISCRETE line, hence the steady-state terminal currents and history
#  phasors, and every ring is filled with Re{H exp(j ws t)} at
#  t = -(L-1) ts, ..., 0.  The newest sample must belong to t = 0: the first
#  solution, at t_1 = ts, reads the samples m-1 and m positions behind it,
#  which bracket t_1 - tau.  The line then starts in its exact steady state.
#
#  Each phase of a line is an independent conductor with the positive-
#  sequence data, so `n` is three times the number of Bergeron lines.
# ==========================================================================

"""
    BergeronLines

The Bergeron conductors of one system.  The fields up to `len` are constants
shared by all scenarios; `head`, `Hs`, `Hr`, `i_from` and `i_to` are the state
of one scenario.
"""
struct BergeronLines
    n::Int                             # number of conductors (3 per line)
    line::Vector{Int}                  # parent line (index into pfd.line_*)
    F::Vector{Int}                     # network node of the sending end s
    T::Vector{Int}                     # network node of the receiving end r
    zc::Vector{Float64}                # surge impedance Zc [pu]
    z::Vector{Float64}                 # terminal resistance Z = Zc + R/4; 1/Z is stamped into G
    h::Vector{Float64}                 # h = (Zc - R/4)/Z
    wt::Vector{Float64}                # (1+h)/2: weight of the wave transmitted from the other end
    wr::Vector{Float64}                # (1-h)/2: weight of the wave reflected at the midpoint
    delay::Vector{Float64}             # tau/ts (not an integer in general)
    off::Vector{Int}                   # start of this conductor's ring inside Hs/Hr (0-based)
    len::Vector{Int}                   # ring length L = ceil(tau/ts)
    head::Vector{Int}                  # position of the newest sample in the ring (0-based)
    Hs::Vector{Float64}                # history H at the sending ends, all rings back to back
    Hr::Vector{Float64}                # history H at the receiving ends
    i_from::Vector{Float64}            # terminal current into the line at s (for output)
    i_to::Vector{Float64}              # terminal current into the line at r
end

"No Bergeron lines (every line is a pi section)."
BergeronLines() = BergeronLines(0, Int[], Int[], Int[], Float64[], Float64[], Float64[], Float64[], Float64[],
                                Float64[], Int[], Int[], Int[], Float64[], Float64[], Float64[], Float64[])

Base.isempty(b::BergeronLines) = b.n == 0

"""
    line_travel_times(pfd, ts) -> Vector{Float64}

Travel time of every line in time steps, `sqrt(X B)/(ws ts)`.  Lines without
a travelling-wave model (X <= 0 or B <= 0) get 0 and stay pi sections.
"""
function line_travel_times(pfd::PowerFlowData, ts::Float64)
    d = zeros(length(pfd.line_from))
    for i in eachindex(pfd.line_from)
        X = imag(pfd.line_RX[i]); B = pfd.line_chg[i]
        (X > 0 && B > 0) && (d[i] = sqrt(X * B) / pfd.ws / ts)
    end
    return d
end

"""
    select_bergeron_lines(pfd, ts, threshold) -> (bergeron, pi)

Indices of the lines whose travel time exceeds `threshold` time steps
(Bergeron lines) and of the others (pi sections).  `threshold = Inf` keeps
every line a pi section.  The threshold must be at least one step, because
both interpolated samples must be known before the solution.
"""
function select_bergeron_lines(pfd::PowerFlowData, ts::Float64, threshold::Real)
    isfinite(threshold) && threshold < 1 &&
        error("bergeron_threshold must be >= 1 (a Bergeron line needs at least one step of delay)")
    d = line_travel_times(pfd, ts)
    mask = d .> threshold
    return findall(mask), findall(.!mask)
end

"Terminal resistance and history weights `(Z, h, (1+h)/2, (1-h)/2)` of a conductor with resistance `R` and surge impedance `zc`."
function bergeron_coefficients(R::Real, zc::Real)
    z = zc + R / 4
    h = (zc - R / 4) / z
    return Float64(z), Float64(h), Float64((1 + h) / 2), Float64((1 - h) / 2)
end

"""
    bergeron_delay_phasor(delay, ts, ws) -> D

Phasor at `ws` of the interpolated delay of `delay` = tau/ts steps (paper,
eq. (6)): D = (1 - phi) exp(-j ws m ts) + phi exp(-j ws (m+1) ts).
"""
function bergeron_delay_phasor(delay::Real, ts::Real, ws::Real)
    m = floor(Int, delay); phi = delay - m
    return (1 - phi) * cis(-ws * m * ts) + phi * cis(-ws * (m + 1) * ts)
end

"""
    bergeron_twoport(z, h, wt, wr, D) -> (Y11, Y12)

Steady-state two-port of a conductor whose delay has the phasor `D`: the
currents into the line are I_s = Y11 V_s + Y12 V_r and I_r = Y12 V_s + Y11 V_r.
Obtained from I_s = V_s/Z - D (wt H_r + wr H_s) and H = V/Z + h I at both
ends; with wt = (1+h)/2 and wr = (1-h)/2 it equals eq. (8) of the paper.
"""
function bergeron_twoport(z::Real, h::Real, wt::Real, wr::Real, D::Number)
    p = 1 + D * wr * h; q = D * wt * h
    u = (1 - D * wr) / z; v = -D * wt / z
    det = p^2 - q^2
    return (p * u - q * v) / det, (p * v - q * u) / det
end

"""
    build_bergeron(pfd, ts, idx, Vt0, nbus) -> BergeronLines

Builds the Bergeron conductors of the lines `idx` and fills their ring
buffers with the exact steady state of the power-flow phasors `Vt0` (the
four steps of the paper's Section IV-B).
"""
function build_bergeron(pfd::PowerFlowData, ts::Float64, idx::Vector{Int}, Vt0::Vector{ComplexF64}, nbus::Int)
    isempty(idx) && return BergeronLines()
    n = 3length(idx)
    line = zeros(Int, n); F = zeros(Int, n); T = zeros(Int, n)
    zc = zeros(n); z = zeros(n); h = zeros(n); wt = zeros(n); wr = zeros(n); delay = zeros(n)
    off = zeros(Int, n); len = zeros(Int, n); head = zeros(Int, n)
    Hs = Float64[]; Hr = Float64[]
    ws = pfd.ws
    for (m, i) in enumerate(idx)
        f = bus_index(pfd, pfd.line_from[i]); t = bus_index(pfd, pfd.line_to[i])
        R = real(pfd.line_RX[i]); X = imag(pfd.line_RX[i]); B = pfd.line_chg[i]
        zci = sqrt(X / B)                                   # surge impedance
        di = sqrt(X * B) / ws / ts                          # tau/ts
        zi, hi, wti, wri = bergeron_coefficients(R, zci)
        L = ceil(Int, di)                                   # ring length
        # steps 1-2: delay phasor D (eq. (6)) and the discrete two-port (eq. (8))
        Y11, Y12 = bergeron_twoport(zi, hi, wti, wri, bergeron_delay_phasor(di, ts, ws))
        # sample times of the ring, oldest first; the newest sample (head) is t = 0
        times = [-(L - 1 - j) * ts for j in 0:(L-1)]
        for ph in 1:3
            k = 3(m - 1) + ph
            line[k] = i
            F[k] = phase_node(f, ph, nbus); T[k] = phase_node(t, ph, nbus)
            zc[k] = zci; z[k] = zi; h[k] = hi; wt[k] = wti; wr[k] = wri; delay[k] = di
            off[k] = length(Hs); len[k] = L; head[k] = L - 1
            Vs = Vt0[F[k]]; Vr = Vt0[T[k]]
            Is = Y11 * Vs + Y12 * Vr                         # terminal currents into the line
            Ir = Y12 * Vs + Y11 * Vr
            # step 3: history phasors H = V/Z + h I (eq. (7)); step 4: fill the rings
            append!(Hs, real.((Vs / zi + hi * Is) .* cis.(ws .* times)))
            append!(Hr, real.((Vr / zi + hi * Ir) .* cis.(ws .* times)))
        end
    end
    return BergeronLines(n, line, F, T, zc, z, h, wt, wr, delay, off, len, head, Hs, Hr, zeros(n), zeros(n))
end

"""
    delayed_samples(head, len, delay) -> (i_new, i_old, phi)

Ring positions (0-based) of the two samples that bracket t_n - tau, and the
interpolation weight phi of the older one (eq. (5)).  `head` holds the
newest sample, that of t_{n-1}, so the samples of t_{n-m} and t_{n-m-1} lie
m-1 and m positions behind it.
"""
@inline function delayed_samples(head::Integer, len::Integer, delay::Real)
    m = unsafe_trunc(Int, floor(delay))
    phi = delay - m
    return mod(head - (m - 1), len), mod(head - m, len), phi
end

"""
    bergeron_injection!(node_I, ber)

History current Ih of every Bergeron terminal for the step about to be
solved (eq. (4)), added to the right-hand side of the network equation.
"""
function bergeron_injection!(node_I::Vector{Float64}, ber::BergeronLines)
    fill!(node_I, 0.0)
    @inbounds for k in 1:ber.n
        o = ber.off[k]
        i0, i1, phi = delayed_samples(ber.head[k], ber.len[k], ber.delay[k])
        hs = (1 - phi) * ber.Hs[o+i0+1] + phi * ber.Hs[o+i1+1]     # H_s(t - tau)
        hr = (1 - phi) * ber.Hr[o+i0+1] + phi * ber.Hr[o+i1+1]     # H_r(t - tau)
        node_I[ber.F[k]] += ber.wt[k] * hr + ber.wr[k] * hs
        node_I[ber.T[k]] += ber.wt[k] * hs + ber.wr[k] * hr
    end
    return node_I
end

"""
    update_bergeron!(ber, V)

After the network solution: terminal currents from eq. (4), new histories
H(t_n) = v/Z + h i written over the oldest ring sample, and `head` advanced.
"""
function update_bergeron!(ber::BergeronLines, V::AbstractVector)
    @inbounds for k in 1:ber.n
        o = ber.off[k]; L = ber.len[k]
        old = ber.head[k]
        new = mod(old + 1, L)                                       # the oldest slot
        ber.head[k] = new
        i0, i1, phi = delayed_samples(old, L, ber.delay[k])
        hs = (1 - phi) * ber.Hs[o+i0+1] + phi * ber.Hs[o+i1+1]     # histories injected this step
        hr = (1 - phi) * ber.Hr[o+i0+1] + phi * ber.Hr[o+i1+1]
        vs = V[ber.F[k]]; vr = V[ber.T[k]]
        z = ber.z[k]
        i_s = vs / z - (ber.wt[k] * hr + ber.wr[k] * hs)
        i_r = vr / z - (ber.wt[k] * hs + ber.wr[k] * hr)
        ber.Hs[o+new+1] = vs / z + ber.h[k] * i_s
        ber.Hr[o+new+1] = vr / z + ber.h[k] * i_r
        ber.i_from[k] = i_s
        ber.i_to[k] = i_r
    end
    return ber
end

"A copy of the Bergeron lines that shares the constants and owns its own state (for a new scenario)."
function copy_state(b::BergeronLines)
    isempty(b) && return BergeronLines()
    return BergeronLines(b.n, b.line, b.F, b.T, b.zc, b.z, b.h, b.wt, b.wr, b.delay, b.off, b.len,
                         copy(b.head), copy(b.Hs), copy(b.Hr), copy(b.i_from), copy(b.i_to))
end
