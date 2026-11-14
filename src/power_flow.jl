# ==========================================================================
#  power_flow.jl -- PSS/E cases with the power flow of PowerFlows.jl
#
#  `load_psse` reads the .raw file with PowerSystems.jl, solves the AC power
#  flow with PowerFlows.jl and hands the solved network to psse_case
#  (psse.jl), which adds the .dyr models.  The network is converted component
#  by component into the models of NodalEMT's network (network.jl), with the
#  conventions of the bus admittance matrix that PowerFlows.jl solved
#  (PowerNetworkMatrices.jl):
#
#    Line            pi section; unequal end charging and line conductances
#                    become fixed shunts at the line ends
#    2-winding xfmr  series R + jX with the ratio t at the from end
#                    (Y_ff = y/t^2, Y_ft = -y/t, Y_tt = y), magnetizing
#                    admittance as a fixed shunt at the from bus
#    3-winding xfmr  three such transformers from the windings to the star bus
#    FixedAdmittance fixed shunt
#    SwitchedAdmittance  switched shunt at Y + sum(status .* Y_increase)
#    StandardLoad / PowerLoad   load drawing its ZIP power at the solved voltage
#
#  power_flow_mismatch (network.jl) then checks that this network reproduces the
#  power-flow solution.  PowerNetworkMatrices.jl stores the admittance matrix
#  in single precision (ComplexF32), so the solution balances the exact
#  network only to about 1e-7 of the largest admittance (1e-5 pu for the
#  two-area system, 7e-4 pu for the WECC system with its 1e4-pu jumpers),
#  not to the solver tolerance; a warning is given above 1e-2 pu, which a
#  wrongly modelled branch or shunt exceeds.  The solution stored in a .raw
#  file, rounded to five decimals, can miss by more (5e-2 pu for WECC).
# ==========================================================================

const PSY = PowerSystems
const PF = PowerFlows

const SOLVERS = Dict(:newton_raphson => PF.NewtonRaphsonACPowerFlow, :trust_region => PF.TrustRegionACPowerFlow,
                     :levenberg_marquardt => PF.LevenbergMarquardtACPowerFlow)

"""
    load_psse(raw, dyr; kwargs...) -> CaseData

Builds a case from a PSS/E power-flow file `raw` and a dynamic-data file `dyr`:
reads the network with PowerSystems.jl, solves the AC power flow with
PowerFlows.jl, reads the dynamic models and initializes nothing yet (pass the
result to [`build_simulation`](@ref) or [`build_batch`](@ref)).

Keywords:

* `solve_power_flow = true` -- solve the power flow; `false` keeps the
  solution stored in the .raw file
* `reactive_power_limits = false` -- enforce the reactive limits of the units
  (PV buses at a limit become PQ buses)
* `solver = :newton_raphson` -- or `:trust_region` or `:levenberg_marquardt`,
  PowerFlows.jl's more robust solvers
* `tol = 1e-9`, `maxiter = 50` -- convergence tolerance [pu] and iteration limit
* `pll_ke = 0.1`, `pll_te = 0.05`, `vm_te = 0.02` -- bus PLL gain and time
  constant and voltage-filter time constant (not part of PSS/E data)
* `json = nothing` -- a path: also write the case as a JSON file that
  [`load_case`](@ref) reads without solving the power flow again
* `verbose = false` -- show the messages of PowerSystems.jl and PowerFlows.jl

Supported data: buses, lines, two- and three-winding transformers with
off-nominal ratios (no phase shift), fixed and switched shunts, loads (any
ZIP mix), and units modelled by GENROU + SEXS + TGOV1/HYGOV/GAST (+ IEEEST)
or by REGCA1 + REECB1 + REPCA1.  RAW versions 33, 34 and 35.

The power flow follows PSS/E where they differ: a PV bus without an
in-service unit is solved as PQ, vector groups (VECGRP) add no phase shift,
and every unit keeps its scheduled output except for its share (by MVA base)
of the change the power flow makes at its bus.
"""
function load_psse(raw::AbstractString, dyr::AbstractString; solve_power_flow::Bool=true,
                   reactive_power_limits::Bool=false, solver::Symbol=:newton_raphson,
                   tol::Real=1e-9, maxiter::Integer=50,
                   pll_ke::Real=0.1, pll_te::Real=0.05, vm_te::Real=0.02,
                   json::Union{Nothing,AbstractString}=nothing, verbose::Bool=false,
                   name::AbstractString=splitext(basename(raw))[1])
    isfile(raw) || error("no such file: $raw")
    isfile(dyr) || error("no such file: $dyr")
    records = read_dyr(dyr)
    sys = with_log_level(verbose ? Base.CoreLogging.Info : Base.CoreLogging.Error) do
        read_raw(raw)
    end
    PSY.set_units_base_system!(sys, "SYSTEM_BASE")
    ignore_vector_groups!(sys)
    if solve_power_flow
        pq_without_units!(sys)
        scheduled = Dict(PSY.get_name(g) => (unit_P(g), PSY.get_reactive_power(g)) for g in units_of(sys))
        haskey(SOLVERS, solver) || error("solver must be one of $(sort(collect(keys(SOLVERS))))")
        pf = PF.ACPowerFlow{SOLVERS[solver]}(; check_reactive_power_limits=reactive_power_limits,
                                             correct_bustypes=true)
        converged = with_log_level(verbose ? Base.CoreLogging.Info : Base.CoreLogging.Error) do
            PF.solve_and_store_power_flow!(pf, sys; tol=tol, maxIterations=maxiter)
        end
        converged || error("the AC power flow of $(basename(raw)) did not converge with the $solver solver" *
                           (solver === :newton_raphson ? "; try solver = :trust_region" : ""))
        share_bus_injections!(sys, scheduled)
    end
    network, units = network_and_units(sys)
    pfd_d, dyd_d = psse_case(network, units, records; pll_ke=pll_ke, pll_te=pll_te, vm_te=vm_te)
    pfd = PowerFlowData(pfd_d)
    case = CaseData(name, pfd, DynamicData(dyd_d, pfd))
    mis = maximum(abs, power_flow_mismatch(pfd))
    mis > (solve_power_flow ? 1e-2 : 1e-1) &&
        @warn "the network as NodalEMT models it misses the power-flow balance by up to $(mis) pu"
    json === nothing || write_case_json(json, name, pfd_d, dyd_d)
    return case
end

"Runs `f` with only the messages of level `level` and above."
function with_log_level(f, level)
    return Base.CoreLogging.with_logger(f, Base.CoreLogging.SimpleLogger(stderr, level))
end

"Reads a .raw file (v33, v34 or v35) into a PowerSystems.System."
function read_raw(raw::AbstractString)
    path, _ = raw_for_parser(raw)
    try
        return PSY.System(path; runchecks=false,
                          bus_name_formatter=d -> strip(string(d["name"])) * "-" * string(d["index"]),
                          gen_name_formatter=d -> join(strip.(string.(d["source_id"][2:3])), "|"))
    finally
        path == raw || rm(path; force=true)
    end
end

"""
PowerSystems.jl turns the PSS/E vector group (VECGRP) into a phase shift of
30 degrees per clock number; PSS/E's power flow ignores it, and so does the
per-phase EMT network.  Phase shifts given as winding angles are not supported.
"""
function ignore_vector_groups!(sys)
    for t in PSY.get_components(Union{PSY.Transformer2W, PSY.TapTransformer}, sys)
        PSY.get_winding_group_number(t) in (PSY.WindingGroupNumber.UNDEFINED, PSY.WindingGroupNumber.GROUP_0) ||
            PSY.set_winding_group_number!(t, PSY.WindingGroupNumber.GROUP_0)
    end
    for t in PSY.get_components(PSY.Transformer3W, sys)
        PSY.set_primary_group_number!(t, PSY.WindingGroupNumber.GROUP_0)
        PSY.set_secondary_group_number!(t, PSY.WindingGroupNumber.GROUP_0)
        PSY.set_tertiary_group_number!(t, PSY.WindingGroupNumber.GROUP_0)
    end
    return sys
end

"""
The in-service units: generators, and synchronous condensers (PowerSystems.jl
makes a condenser of a PSS/E machine with PG = PT = PB = 0).
"""
units_of(sys) = [g for g in PSY.get_components(Union{PSY.Generator, PSY.SynchronousCondenser}, sys) if PSY.get_available(g)]
unit_P(g) = g isa PSY.SynchronousCondenser ? 0.0 : PSY.get_active_power(g)

"""
A PV bus without an in-service unit -- for instance one whose voltage a
switched shunt regulated in PSS/E -- is solved as a PQ bus, as PSS/E does;
its switched shunt keeps the admittance of the .raw file.
"""
function pq_without_units!(sys)
    units = Set(PSY.get_number(PSY.get_bus(g)) for g in units_of(sys))
    for b in PSY.get_components(PSY.ACBus, sys)
        bt = PSY.get_bustype(b)
        PSY.get_number(b) in units && continue
        bt == PSY.ACBusTypes.REF && error("the slack bus $(PSY.get_number(b)) has no in-service unit")
        bt in (PSY.ACBusTypes.PV, PSY.ACBusTypes.SLACK) && PSY.set_bustype!(b, PSY.ACBusTypes.PQ)
    end
    return sys
end

"""
    share_bus_injections!(sys, scheduled)

PowerFlows.jl solves for the total injection of each bus and shares it among
the bus's units in its own way (in proportion to active power, clipped at the
unit limits), which can move a unit's dispatch.  Here every unit keeps its
scheduled P and Q of the .raw file, and only the change that the power flow
made to a bus total (P at the slack bus, Q at PV buses) is shared among the
bus's units in proportion to their MVA bases.  For a .raw file that PSS/E has
solved, this returns PSS/E's unit outputs.
"""
function share_bus_injections!(sys, scheduled)
    atbus = Dict{Int,Vector{Any}}()
    for g in units_of(sys)
        push!(get!(atbus, PSY.get_number(PSY.get_bus(g)), Any[]), g)
    end
    for units in values(atbus)
        gens = [g for g in units if !(g isa PSY.SynchronousCondenser)]    # condensers stay at P = 0
        if !isempty(gens)
            Psol = sum(unit_P, gens); Psch = sum(g -> scheduled[PSY.get_name(g)][1], gens)
            mb = sum(PSY.get_base_power, gens)
            for g in gens
                PSY.set_active_power!(g, scheduled[PSY.get_name(g)][1] + (Psol - Psch) * PSY.get_base_power(g) / mb)
            end
        end
        Qsol = sum(PSY.get_reactive_power, units); Qsch = sum(g -> scheduled[PSY.get_name(g)][2], units)
        mb = sum(PSY.get_base_power, units)
        for g in units
            PSY.set_reactive_power!(g, scheduled[PSY.get_name(g)][2] + (Qsol - Qsch) * PSY.get_base_power(g) / mb)
        end
    end
    return sys
end

"PSS/E bus type codes; SLACK marks PSS/E's area-interchange bus, which the power flow solves as PV."
const BUS_TYPE = Dict(PSY.ACBusTypes.PQ => 1, PSY.ACBusTypes.PV => 2, PSY.ACBusTypes.REF => 3,
                      PSY.ACBusTypes.SLACK => 2)

"Circuit id from a PowerSystems branch name (`...-i_<ckt>`)."
circuit_id(name) = (m = match(r"-i_([^-]*)$", name); m === nothing ? "1" : String(m.captures[1]))

"""
    network_and_units(sys) -> (network, units)

The "pfd" entries of the solved system (everything but the units) and its
in-service units, sorted by bus number and id.
"""
function network_and_units(sys)
    PSY.set_units_base_system!(sys, "SYSTEM_BASE")
    base = PSY.get_base_power(sys)
    buses = sort([b for b in PSY.get_components(PSY.ACBus, sys) if PSY.get_bustype(b) != PSY.ACBusTypes.ISOLATED],
                 by=PSY.get_number)
    inservice = Set(PSY.get_number(b) for b in buses)
    num(b) = (n = PSY.get_number(b); n in inservice || error("an in-service device is connected to isolated bus $n"); n)
    Vm = Dict(PSY.get_number(b) => PSY.get_magnitude(b) for b in buses)

    net = Dict{String,Any}(
        "basemva" => base, "ws" => 2pi * PSY.get_frequency(sys),
        "bus_num" => PSY.get_number.(buses), "bus_name" => PSY.get_name.(buses),
        "bus_type" => [BUS_TYPE[PSY.get_bustype(b)] for b in buses],
        "bus_Vm" => PSY.get_magnitude.(buses), "bus_Va" => PSY.get_angle.(buses),
        "bus_basekV" => [something(PSY.get_base_voltage(b), 0.0) for b in buses])

    lf = Int[]; lt = Int[]; lid = String[]; lrx = ComplexF64[]; lb = Float64[]
    xf = Int[]; xt = Int[]; xid = String[]; xrx = ComplexF64[]; xk = Float64[]
    sb = Int[]; sid = String[]; sgb = ComplexF64[]                 # fixed shunts [MW, Mvar at 1 pu]
    addshunt!(bus, y, id) = abs(y) > 0 && (push!(sb, bus); push!(sid, id); push!(sgb, y * base))
    unsupported = String[]

    for br in PSY.get_components(PSY.Branch, sys)
        PSY.get_available(br) || continue
        if br isa Union{PSY.Line, PSY.MonitoredLine}
            f = num(PSY.get_from(PSY.get_arc(br))); t = num(PSY.get_to(PSY.get_arc(br)))
            b = PSY.get_b(br); g = PSY.get_g(br)
            bsym = min(b.from, b.to)
            push!(lf, f); push!(lt, t); push!(lid, circuit_id(PSY.get_name(br)))
            push!(lrx, complex(PSY.get_r(br), PSY.get_x(br))); push!(lb, 2bsym)
            addshunt!(f, complex(g.from, b.from - bsym), "line"); addshunt!(t, complex(g.to, b.to - bsym), "line")
        elseif br isa Union{PSY.Transformer2W, PSY.TapTransformer, PSY.PhaseShiftingTransformer}
            abs(PSY.get_α(br)) < 1e-12 || (push!(unsupported, "phase-shifting transformer $(PSY.get_name(br))"); continue)
            f = num(PSY.get_from(PSY.get_arc(br))); t = num(PSY.get_to(PSY.get_arc(br)))
            push!(xf, f); push!(xt, t); push!(xid, circuit_id(PSY.get_name(br)))
            push!(xrx, complex(PSY.get_r(br), PSY.get_x(br)))
            push!(xk, br isa PSY.Transformer2W ? 1.0 : PSY.get_tap(br))
            addshunt!(f, PSY.get_primary_shunt(br), "magnetizing")
        elseif br isa PSY.Transformer3W
            star = num(PSY.get_star_bus(br))
            windings = ((PSY.get_available_primary(br), PSY.get_primary_star_arc(br), PSY.get_r_primary(br),
                         PSY.get_x_primary(br), PSY.get_primary_turns_ratio(br)),
                        (PSY.get_available_secondary(br), PSY.get_secondary_star_arc(br), PSY.get_r_secondary(br),
                         PSY.get_x_secondary(br), PSY.get_secondary_turns_ratio(br)),
                        (PSY.get_available_tertiary(br), PSY.get_tertiary_star_arc(br), PSY.get_r_tertiary(br),
                         PSY.get_x_tertiary(br), PSY.get_tertiary_turns_ratio(br)))
            for (k, (on, arc, r, x, ratio)) in enumerate(windings)
                on || continue
                w = num(PSY.get_from(arc))
                push!(xf, w); push!(xt, star); push!(xid, "$(circuit_id(PSY.get_name(br)))w$k")
                push!(xrx, complex(r, x)); push!(xk, ratio)
            end
            addshunt!(num(PSY.get_from(PSY.get_primary_star_arc(br))), complex(PSY.get_g(br), PSY.get_b(br)), "magnetizing")
        else
            push!(unsupported, "$(nameof(typeof(br))) $(PSY.get_name(br))")
        end
    end

    # loads: the power drawn at the solved voltage
    ldb = Int[]; ldid = String[]; lP = Float64[]; lQ = Float64[]
    for ld in PSY.get_components(PSY.ElectricLoad, sys)
        (PSY.get_available(ld) && !(ld isa Union{PSY.FixedAdmittance, PSY.SwitchedAdmittance})) || continue
        k = num(PSY.get_bus(ld)); v = Vm[k]
        if ld isa PSY.StandardLoad
            P = PSY.get_constant_active_power(ld) + PSY.get_current_active_power(ld) * v +
                PSY.get_impedance_active_power(ld) * v^2
            Q = PSY.get_constant_reactive_power(ld) + PSY.get_current_reactive_power(ld) * v +
                PSY.get_impedance_reactive_power(ld) * v^2
        elseif ld isa PSY.PowerLoad
            P = PSY.get_active_power(ld); Q = PSY.get_reactive_power(ld)
        else
            push!(unsupported, "$(nameof(typeof(ld))) $(PSY.get_name(ld))"); continue
        end
        push!(ldb, k); push!(ldid, "1"); push!(lP, P * base); push!(lQ, Q * base)
    end
    for s in PSY.get_components(PSY.FixedAdmittance, sys)
        PSY.get_available(s) && addshunt!(num(PSY.get_bus(s)), PSY.get_Y(s), "1")
    end
    # switched shunts: PSS/E's YSWACT convention, MW and -Mvar at the solved voltage
    swb = Int[]; swgb = ComplexF64[]
    for s in PSY.get_components(PSY.SwitchedAdmittance, sys)
        PSY.get_available(s) || continue
        k = num(PSY.get_bus(s))
        Y = PSY.get_Y(s) + sum(PSY.get_initial_status(s) .* PSY.get_Y_increase(s); init=0.0im)
        push!(swb, k); push!(swgb, complex(-real(Y), -imag(Y)) * Vm[k]^2 * base)
    end

    # units
    units = PsseUnit[]
    for g in PSY.get_components(PSY.StaticInjection, sys)
        PSY.get_available(g) || continue
        g isa Union{PSY.ElectricLoad} && continue
        if g isa Union{PSY.Generator, PSY.SynchronousCondenser}
            bus = num(PSY.get_bus(g))
            id = String(last(split(PSY.get_name(g), "|"; limit=2)))
            push!(units, PsseUnit(bus, id, unit_P(g) * base, PSY.get_reactive_power(g) * base,
                                  PSY.get_base_power(g), Float64(get(PSY.get_ext(g), "r", 0.0))))
        else
            push!(unsupported, "$(nameof(typeof(g))) $(PSY.get_name(g))")
        end
    end
    sort!(units, by=u -> (u.bus, u.id))
    isempty(unsupported) || error("the .raw data contain in-service devices that NodalEMT does not model:\n  " *
                                  join(unsupported, "\n  "))
    # loads and shunts in bus order, as PSS/E lists them
    p = sortperm(ldb); ldb, ldid, lP, lQ = ldb[p], ldid[p], lP[p], lQ[p]
    p = sortperm(sb); sb, sid, sgb = sb[p], sid[p], sgb[p]
    p = sortperm(swb); swb, swgb = swb[p], swgb[p]

    merge!(net, Dict{String,Any}(
        "line_from" => lf, "line_to" => lt, "line_id" => lid, "line_RX" => complex_entry(lrx), "line_chg" => lb,
        "xfmr_from" => xf, "xfmr_to" => xt, "xfmr_id" => xid, "xfmr_RX" => complex_entry(xrx), "xfmr_k" => xk,
        "load_bus" => ldb, "load_id" => ldid, "load_MW" => lP, "load_Mvar" => lQ,
        "shnt_bus" => sb, "shnt_id" => sid, "shnt_gb" => complex_entry(sgb),
        "shnt_sw_bus" => swb, "shnt_sw_gb" => complex_entry(swgb)))
    return net, units
end
