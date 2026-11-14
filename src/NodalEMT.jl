"""
    NodalEMT

NodalEMT.jl: an open-source Julia simulator of
electromagnetic transients (EMT) for batches of fault and outage scenarios
on a multicore CPU.  It implements the network, machine, exciter, governor,
stabilizer and IBR models of NREL's ParaEMT with

* faults and generator trips applied by compensation, so that the network
  matrix G and its factorization never change;
* scenario batches on worker threads, each worker owning one factorization;
* pi+Bergeron lines: long lines are lossy Bergeron models with fractional-
  delay interpolation in packed ring buffers and an exact periodic
  initialization.

Quick start:

    using NodalEMT
    psse = joinpath(pkgdir(NodalEMT), "data", "psse")
    case = load_psse(joinpath(psse, "twoarea.raw"), joinpath(psse, "twoarea.dyr"); vm_te=0.1)
    sim  = build_simulation(case)                  # 50 us, Bergeron lines, updated models
    schedule_bus_fault!(sim; bus=8, t_on=1.0, t_off=1.1)
    res  = run!(sim, 2.0)
    save_results(res, sim, "results_2area")

Source files, in reading order:

    case_data.jl        power flow and dynamic parameters from a JSON case file
    bergeron.jl         lossy Bergeron model: ring buffers, fractional delay, exact initialization
    network.jl          companion models, the network matrix G, history currents
    compensation.jl     faults and generator trips without changing G
    machine.jl          GENROU synchronous machine as a Norton source
    controls.jl         swing equation, SEXS, TGOV1/HYGOV/GAST, IEEEST
    bus_measurement.jl  PLL and voltage filter of every bus
    ibr.jl              inverter-based resources (REGCA, REECB, REPCA)
    results.jl          recording and CSV output
    simulation.jl       initialization, disturbances, one time step, time loop
    batch.jl            scenario batches on worker threads
    psse.jl             the .dyr reader and the assembly of a PSS/E case
    power_flow.jl       load_psse: the .raw file and its power flow (PowerSystems.jl, PowerFlows.jl)
    precompile.jl       compiling the common path when the package is installed

A case comes from PSS/E .raw and .dyr files with the power flow solved by
PowerFlows.jl (`load_psse`), or from a JSON file that `load_psse` saved
(`load_case`).
"""
module NodalEMT

using LinearAlgebra
using SparseArrays
import PowerSystems
import PowerFlows

include("case_data.jl")
include("bergeron.jl")
include("network.jl")
include("compensation.jl")
include("machine.jl")
include("controls.jl")
include("bus_measurement.jl")
include("ibr.jl")
include("results.jl")
include("simulation.jl")
include("batch.jl")
include("psse.jl")
include("power_flow.jl")

export CaseData, load_case, load_psse, nbus, ngen, nibr, bus_index, line_travel_times
export EmtSimulation, build_simulation, step!, run!, current_time
export schedule_bus_fault!, schedule_generator_trip!, schedule_reference_step!
export ScenarioBatch, build_batch, run_steps!, nscenarios, nworkers
export Results, save_results

include("precompile.jl")

end # module
