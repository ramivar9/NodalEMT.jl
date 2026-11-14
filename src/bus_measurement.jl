# ==========================================================================
#  bus_measurement.jl -- phase-locked loop (PLL) and voltage-magnitude
#  filter of every bus
#
#  The controllers need the voltage magnitude, angle and frequency of a bus,
#  which an EMT simulation does not have directly: only the instantaneous
#  phase voltages va, vb, vc are known.  Every bus therefore carries
#
#    * a PLL that tracks the angle `de` and the frequency `we` (pu) of the
#      voltage by driving its q-axis component vq to zero:
#          vq  = -2/3 [ sin(th) va + sin(th - 2pi/3) vb + sin(th + 2pi/3) vc ]
#          ze' = (ke/te) vq          (integral part)
#          we  = 1 + ke vq + ze      (frequency, pu)
#          de' = we ws               (angle)
#    * a first-order filter with time constant vm_te on the instantaneous
#      magnitude vt = sqrt(2/3 (va^2 + vb^2 + vc^2)):
#          dvtm = (vt - vtm) / vm_te,      vtm' = dvtm
# ==========================================================================

"Measurement states of every bus."
mutable struct BusState
    ze::Vector{Float64}     # PLL integrator
    de::Vector{Float64}     # PLL angle [rad]
    we::Vector{Float64}     # PLL frequency [pu]
    vt::Vector{Float64}     # instantaneous voltage magnitude [pu]
    vtm::Vector{Float64}    # filtered voltage magnitude [pu]
    dvtm::Vector{Float64}   # derivative of vtm
end
BusState(n::Int) = BusState((zeros(n) for _ in 1:6)...)

function Base.copyto!(dst::BusState, src::BusState)
    for f in fieldnames(BusState)
        copyto!(getfield(dst, f), getfield(src, f))
    end
    return dst
end

"Starts every PLL at the power-flow angle and every filter at the power-flow magnitude."
function init_bus_measurement!(b::BusState, pfd::PowerFlowData)
    for k in 1:nbus(pfd)
        b.ze[k] = 0.0
        b.de[k] = pfd.bus_Va[k]
        b.we[k] = 1.0
        b.vt[k] = pfd.bus_Vm[k]
        b.vtm[k] = pfd.bus_Vm[k]
        b.dvtm[k] = 0.0
    end
    return b
end

"""
    bus_measurement_step(ze, de, we, vtm, va, vb, vc, ke, te, vm_te, ws, ts)
        -> (ze_new, de_new, we_new, vt, vtm_new, dvtm)

One PLL and filter step of one bus.
"""
@inline function bus_measurement_step(ze, de, we, vtm, va, vb, vc, ke, te, vm_te, ws, ts)
    vt = sqrt((va * va + vb * vb + vc * vc) * 2 / 3)
    dvtm = (vt - vtm) / vm_te
    th = de + ts * we * ws
    vq = -2 / 3 * (sin(th) * va + sin(th - 2pi / 3) * vb + sin(th + 2pi / 3) * vc)
    ze_new = ze + ke / te * vq * ts
    de_new = de + we * ws * ts
    we_new = 1 + ke * vq + ze
    return ze_new, de_new, we_new, vt, vtm + dvtm * ts, dvtm
end

"Advances the PLL and the filter of every bus with the node voltages `V` of this step."
function update_bus_measurement!(b::BusState, meas::MeasurementParams, V::Vector{Float64}, nbus::Int,
                                 ts::Float64, ws::Float64)
    @inbounds for k in 1:nbus
        b.ze[k], b.de[k], b.we[k], b.vt[k], b.vtm[k], b.dvtm[k] =
            bus_measurement_step(b.ze[k], b.de[k], b.we[k], b.vtm[k], V[k], V[k + nbus], V[k + 2nbus],
                                 meas.pll_ke[k], meas.pll_te[k], meas.vm_te[k], ws, ts)
    end
    return b
end
