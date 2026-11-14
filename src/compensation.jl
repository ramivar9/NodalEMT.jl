# ==========================================================================
#  compensation.jl -- faults and generator trips without changing G
#                     (paper, Section III)
#
#  A bus fault connects the faulted phase nodes to ground through the fault
#  impedance Zf, and a generator trip removes the 3x3 Norton conductance Geq
#  of the machine.  Both change the network matrix only at a few "port"
#  nodes, from G to
#
#        G + E dY E'
#
#  where E (N x m) selects the m port nodes and dY = I/Zf or dY = -Geq.  For
#  the base solution v0 = G^-1 (i_s + i_h), the Sherman-Morrison-Woodbury
#  identity gives the solution of the modified network (eq. (2)):
#
#        v = v0 - ZT (dY^-1 + ZP)^-1 E' v0,
#        ZT = G^-1 E   (N x m transfer impedances),
#        ZP = E' ZT     (m x m port impedances).
#
#  At the first step of a scenario, the ports of all of its events are
#  collected in E and the m columns of ZT are computed by one forward and
#  backward substitution each.  The events in effect at a step select a
#  subset a of the ports, and (eq. (3))
#
#        v = v0 - ZT_a M_a E_a' v0,      M_a = (dY_a^-1 + ZP_a)^-1.
#
#  M_a is recomputed only when the subset changes (fault inception, fault
#  clearing, trip); every other step costs the m_a port voltages of v0, one
#  product with M_a and a rank-m_a update.  G and its factors are never
#  touched, so they can be shared by every scenario of a batch (batch.jl).
#  As eq. (2) is an identity, the result equals that of the modified matrix
#  up to round-off.
# ==========================================================================

using LinearAlgebra

"""
    NetworkChange(nodes, dY)

One modification of the network: the admittance block `dY` (m x m) added
between the m network `nodes`.
"""
struct NetworkChange
    nodes::Vector{Int}
    dY::Matrix{Float64}
end

"A fault of impedance `Zf` [pu] from each of the given phase nodes to ground: dY = I/Zf."
fault_change(nodes::Vector{Int}, Zf::Real) =
    NetworkChange(nodes, Matrix{Float64}((1 / Zf) * I, length(nodes), length(nodes)))

"A tripped machine: its Norton conductance leaves the three terminal nodes, dY = -Geq."
trip_change(nodes::Vector{Int}, Geq::AbstractMatrix) = NetworkChange(nodes, -Matrix{Float64}(Geq))

"""
    Compensation

The compensation matrices of one scenario: ZT and ZP for the ports of all of
its events, and M_a for the events in effect at the current step.
"""
mutable struct Compensation
    ports::Vector{Int}                 # port nodes of every change, back to back (the columns of E)
    slots::Vector{UnitRange{Int}}      # the ports of each change
    dYinv::Matrix{Float64}             # block-diagonal dY^-1
    ZT::Matrix{Float64}                # N x m transfer impedances G^-1 E
    ZP::Matrix{Float64}                # m x m port impedances E' ZT
    active::Vector{Bool}               # which changes are in effect
    cur::Vector{Int}                   # the ports of the changes in effect (the subset a)
    Ma::Matrix{Float64}                # M_a = (dY_a^-1 + ZP_a)^-1
    ZTa::Matrix{Float64}               # the columns ZT_a, stored contiguously
    lam::Vector{Float64}               # work: M_a E_a' v0
    vp::Vector{Float64}                # work: E_a' v0
end

"""
    build_compensation(net, changes) -> Compensation

Computes ZT by one substitution with the factors of G per port.  Nothing
is in effect until [`set_active!`](@ref) is called.
"""
function build_compensation(net::Network, changes::Vector{NetworkChange})
    N = net.N
    ports = Int[]
    slots = UnitRange{Int}[]
    for ch in changes
        push!(slots, (length(ports)+1):(length(ports)+length(ch.nodes)))
        append!(ports, ch.nodes)
    end
    m = length(ports)
    dYinv = zeros(m, m)
    for (k, ch) in enumerate(changes)
        dYinv[slots[k], slots[k]] = inv(ch.dY)
    end
    ZT = zeros(N, m)
    e = zeros(N)
    for j in 1:m
        fill!(e, 0.0); e[ports[j]] = 1.0
        ldiv!(view(ZT, :, j), net.lu, e)
    end
    return Compensation(ports, slots, dYinv, ZT, ZT[ports, :], falses(length(changes)),
                        Int[], zeros(0, 0), zeros(N, m), zeros(m), zeros(m))
end

"""
    set_active!(c, k, on)

Puts change `k` into effect (`on = true`) or removes it, and recomputes M_a
for the new subset.  This is the whole cost of a fault inception, a fault
clearing or a trip.
"""
function set_active!(c::Compensation, k::Int, on::Bool)
    c.active[k] == on && return c
    c.active[k] = on
    empty!(c.cur)
    for j in eachindex(c.active)
        c.active[j] && append!(c.cur, c.slots[j])
    end
    if isempty(c.cur)
        c.Ma = zeros(0, 0)
    else
        c.Ma = inv(c.dYinv[c.cur, c.cur] + c.ZP[c.cur, c.cur])
        @inbounds for (i, j) in enumerate(c.cur)
            copyto!(view(c.ZTa, :, i), view(c.ZT, :, j))
        end
    end
    return c
end

"""
    apply_compensation!(v, c)

Corrects the base solution in place for the changes in effect (eq. (3)):
v <- v - ZT_a M_a E_a' v.
"""
function apply_compensation!(v::AbstractVector, c::Compensation)
    m = length(c.cur)
    m == 0 && return v
    vp = view(c.vp, 1:m); lam = view(c.lam, 1:m)
    @inbounds for i in 1:m
        vp[i] = v[c.ports[c.cur[i]]]                 # E_a' v0
    end
    mul!(lam, c.Ma, vp)                              # M_a E_a' v0
    mul!(v, view(c.ZTa, :, 1:m), lam, -1.0, 1.0)     # v0 - ZT_a (M_a E_a' v0)
    return v
end
