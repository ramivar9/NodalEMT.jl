# ==========================================================================
#  precompile.jl -- compile the common path when the package is installed
#
#  Reading a .raw file with PowerSystems.jl and solving its first power flow
#  with PowerFlows.jl compile for about 100 s.  Running that path here, on the
#  two bundled systems, stores the compiled code with the package, so the
#  first load_psse of a Julia session takes seconds.  The workload also
#  builds and steps a simulation and a small batch with a fault, a generator
#  trip and a reference step.
# ==========================================================================

using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    psse = joinpath(@__DIR__, "..", "data", "psse")
    @compile_workload begin
        with_log_level(Base.CoreLogging.Error) do
            for name in ("twoarea", "WECC240")
                case = load_psse(joinpath(psse, name * ".raw"), joinpath(psse, name * ".dyr"))
                bus = case.pfd.bus_num[end]
                sim = build_simulation(case)
                schedule_bus_fault!(sim; bus=bus, t_on=0.0005, t_off=0.001)
                schedule_generator_trip!(sim; t=0.0015, gen=1)
                schedule_reference_step!(sim; t=0.0005, gen=1, dsp=0.001)
                run!(sim, 0.002; progress=false)
                batch = build_batch(case, 2; nworkers=1)
                schedule_bus_fault!(batch[2]; bus=bus, t_on=0.0005, t_off=0.001)
                run!(batch, 0.002; progress=false)
            end
        end
    end
end
