#=
two_area_fault.jl -- three-phase fault at bus 8 of the Kundur two-area system,
with pi+Bergeron lines and with a single pi section per line.

    julia --project=. examples/two_area_fault.jl

This is the disturbance of the PSCAD/EMTDC comparison in the paper (Fig. 6):
a three-phase fault to ground through 1 ohm per phase at bus 8 (230 kV),
cleared after 100 ms by opening all phases at once, at a 50 us step.  The
fault is applied at the same point on the wave as in the PSCAD run, which
puts it at 5.00185 s with the pi+Bergeron lines and at 5.00175 s with pi lines.

Printed: the peak phase-a fault current in the first 6 ms and the peak
bus-8 phase-b voltage in the first 20 ms after clearing.  PSCAD gives
8.06 kA and 818 kV.  Waveforms around the fault are written to
examples/results_two_area/.

The system is read from data/psse/.  vm_te = 0.1 s is ParaEMT's bus voltage
filter for this system.
=#
using NodalEMT
using Printf

const PSSE = joinpath(@__DIR__, "..", "data", "psse")
const OUT = joinpath(@__DIR__, "results_two_area")
const FAULT_BUS = 8
const R_FAULT = 1.0          # fault resistance per phase [ohm]
const DURATION = 0.1         # fault duration [s]

const CASE = load_psse(joinpath(PSSE, "twoarea.raw"), joinpath(PSSE, "twoarea.dyr"); vm_te=0.1)

function fault_run(label, t_on; bergeron_threshold)
    case = CASE
    pfd = case.pfd
    sim = build_simulation(case; ts=50e-6, bergeron_threshold)
    nb = nbus(pfd)
    k = bus_index(pfd, FAULT_BUS)
    zflt = R_FAULT / sim.net.ZbaseA[k]                     # ohm -> pu
    vpk = pfd.bus_basekV[k] * sqrt(2 / 3)                  # peak phase voltage base [kV]
    ipk = 2 * pfd.basemva / (3 * vpk)                      # peak phase current base [kA]
    t_off = t_on + DURATION
    schedule_bus_fault!(sim; bus=FAULT_BUS, t_on, t_off, Zflt=zflt, phases=[1, 2, 3])

    println(label, ": ", sim.net.ber.n ÷ 3, " Bergeron lines")
    peak_i = 0.0; peak_v = 0.0
    mkpath(OUT)
    open(joinpath(OUT, "$(label).csv"), "w") do io
        println(io, "t_s,va_kV,vb_kV,vc_kV,ia_fault_kA")
        while current_time(sim) < t_off + 0.15
            V = step!(sim)
            t = current_time(sim)
            va, vb, vc = V[k] * vpk, V[k+nb] * vpk, V[k+2nb] * vpk
            ia = any(f -> f.active, sim.faults) ? V[k] / zflt * ipk : 0.0
            t_on <= t <= t_on + 0.006 && (peak_i = max(peak_i, abs(ia)))
            t_off <= t <= t_off + 0.020 && (peak_v = max(peak_v, abs(vb)))
            t >= t_on - 0.01 && @printf(io, "%.6f,%.6f,%.6f,%.6f,%.6f\n", t, va, vb, vc, ia)
        end
    end
    @printf("  peak fault current (phase a, first 6 ms):          %7.2f kA\n", peak_i)
    @printf("  peak bus-8 voltage (phase b, 20 ms after clearing): %7.1f kV\n", peak_v)
    return peak_i, peak_v
end

fault_run("pi_bergeron", 5.00185; bergeron_threshold=1.0)       # lines longer than one step are Bergeron lines
fault_run("single_pi", 5.00175; bergeron_threshold=Inf)    # every line a single pi section
println("waveforms written to ", OUT)
