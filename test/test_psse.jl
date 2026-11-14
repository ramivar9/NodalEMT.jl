# Cases from PSS/E .raw and .dyr files (load_psse): the .dyr reader, the RAW
# v34 conversion, the network handed to the EMT model, and the end-to-end
# path .raw/.dyr -> power flow (PowerFlows.jl) -> EMT, checked against the
# power-flow solution that PSS/E stored in the bundled .raw files.

using Test
using NodalEMT
using PowerSystems, PowerFlows
using NodalEMT: read_dyr, dyr_parameters, DYR_PARAMETERS, psse_case, PsseUnit, raw_for_parser, split_record,
                power_flow_mismatch, network_admittance, DynamicData

const PNM = PowerFlows.PNM

"Compares two dynamic-data structs field by field (ids without blanks)."
function same_dynamics(a::DynamicData, b::DynamicData)
    ok = true
    for f in fieldnames(DynamicData)
        x = getfield(a, f); y = getfield(b, f)
        if x isa AbstractVector
            ok &= x == y
        else
            for g in fieldnames(typeof(x))
                u = getfield(x, g); v = getfield(y, g)
                same = eltype(u) <: AbstractString ? strip.(u) == strip.(v) :
                       eltype(u) <: AbstractFloat ? isapprox(u, v; rtol=1e-6) : u == v
                same || @info "dynamic data differ" f g
                ok &= same
            end
        end
    end
    return ok
end

"Largest deviation [pu] of the bus voltage magnitudes from their initial values in an undisturbed run."
function drift(case; nsteps=2000, kw...)
    sim = build_simulation(case; kw...)
    v0 = copy(sim.x_bus.vtm); w = 0.0
    for _ in 1:nsteps
        step!(sim); w = max(w, maximum(abs, sim.x_bus.vtm .- v0))
    end
    return w
end

@testset ".dyr reader: free format, comments, errors" begin
    path = tempname() * ".dyr"
    write(path, """
    @! bus  model  id  parameters
       1 'GENROU' 1  8.0 0.03 0.4 0.05
          6.5 0 1.8 1.7 0.3 0.55 0.25 0.2 0 0 /  text after the slash is a comment
    2,'SEXS','G1',0.1,10,20,0.1D0,0,3/
    /
    """)
    r = read_dyr(path)
    @test length(r) == 2
    @test (r[1].bus, r[1].model, r[1].id, r[1].line) == (1, "GENROU", "1", 2)
    p = dyr_parameters(r[1])
    @test p["Td0p"] == 8.0 && p["H"] == 6.5 && p["Xl"] == 0.2 && p["S12"] == 0.0
    @test (r[2].bus, r[2].model, r[2].id) == (2, "SEXS", "G1") && dyr_parameters(r[2])["TE"] == 0.1
    write(path, "1 'SEXS' 1 0.1 10 /\n")
    @test_throws ErrorException dyr_parameters(only(read_dyr(path)))        # too few values
    write(path, "1 'SEXS' 1 0.1 10 20 0.1 0 3\n")
    @test_throws ErrorException read_dyr(path)                              # no '/'
    rm(path)
    # the bundled files hold one record per model and unit
    count_models(f) = (c = Dict{String,Int}(); foreach(x -> c[x.model] = get(c, x.model, 0) + 1, read_dyr(f)); c)
    c = count_models(joinpath(PSSEDIR, "WECC240.dyr"))
    @test c["GENROU"] == c["SEXS"] == 111 && c["TGOV1"] + c["HYGOV"] + c["GAST"] == 111
    @test c["REGCA1"] == c["REECB1"] == c["REPCA1"] == 91 && c["IEEEST"] == 10
end

@testset ".dyr parameter order agrees with PowerSystems.jl's PSS/E mapping" begin
    # PowerSystems.jl maps the .dyr values to its own model structs by position; translate its
    # positions back to NodalEMT's names and compare.  N: a value that is not read from the record.
    pm = PowerSystems.YAML.load_file(PowerSystems.PSSE_DYR_MAPPING_FILE)["parameter_mapping"][1]
    N = nothing
    table = [
        ("GENROU", pm["RoundRotorQuadratic"], [N, "Td0p", "Td0pp", "Tq0p", "Tq0pp", "Xd", "Xq", "Xdp", "Xqp", "Xdpp", "Xl", "S10", "S12"]),
        ("GENROU", pm["RoundRotorShaft"], ["H", "D"]),
        ("SEXS", pm["SEXS"], ["TA_o_TB", "TB", "K", "TE", "Emin", "Emax"]),
        ("TGOV1", pm["SteamTurbineGov1"], ["R", "T1", "Vmin", "Vmax", "T2", "T3", "Dt", N, N, N]),
        ("HYGOV", pm["HydroTurbineGov"], ["R", "r", "Tr", "Tf", "Tg", "VELM", "GMIN", "GMAX", "TW", "At", "Dturb", "qNL"]),
        ("GAST", pm["GasTG"], ["R", "T1", "T2", "T3", "LdLmt", "KT", "VMIN", "VMAX", "Dturb"]),
        ("IEEEST", pm["IEEEST"], ["ICS", "IB", "A1", "A2", "A3", "A4", "A5", "A6", "T1", "T2", "T3", "T4", "T5", "T6",
                                  "KS", "LSMIN", "LSMAX", "VCU", "VCL"]),
        ("REGCA1", pm["RenewableEnergyConverterTypeA"], ["Tg", "Rrpwr", "Brkpt", "Zerox", "Lvpl1", "Volim", "Lvpnt0",
                                                        "Lvpnt1", "Iolim", "Tfltr", "Khv", "Iqrmin", "Iqrmax", "Accel", "LVPLsw"]),
        ("REECB1", pm["RECurrentControlB"]["REECB1"], ["QFLAG", "PQFLAG", "Vdip", "Vup", "Trv", "dbd1", "dbd2", "Kqv",
                                                       "Iqll", "Iqhl", "Vref0", "Kvp", "Kvi", "Tiq", "Imax"]),
        ("REECB1", pm["ActiveRenewableControllerAB"]["REECB1"], vcat(fill(N, 14), ["dPmin", "dPmax", "Pmin", "Pmax", "Tpord"])),
        ("REECB1", pm["ReactiveRenewableControllerAB"]["REECB1"], vcat(fill(N, 6), ["PFFLAG", "VFLAG"], fill(N, 12),
                                                                       ["Tp", "Qmin", "Qmax", "Vmin", "Vmax", "Kqp", "Kqi"])),
        ("REPCA1", pm["ActiveRenewableControllerAB"]["REPCA1"], ["remote_bus", "branch_From_bus", "branch_To_bus", "branch_id",
            "FFlag", "Kpg", "Kig", "Tp", "fdbd1", "fdbd2", "femin", "femax", "Pmin", "Pmax", "Tg", "Ddn", "Dup", N, N, N]),
        ("REPCA1", pm["ReactiveRenewableControllerAB"]["REPCA1"], vcat(["remote_bus", "branch_From_bus", "branch_To_bus",
            "branch_id", "VCFlag", "RefFlag", N, N, "Tfltr", "Kp", "Ki", "Tft", "Tfv", "Vfrz", "Rc", "Xc", "Kc", "emin", "emax",
            "dbd1", "dbd2", "Qmin", "Qmax"], fill(N, 5))),
    ]
    # an entry is a position, a tuple "(a, b)" of positions, or "NaN"/a number for a value not read from the record
    positions(e) = e isa Integer ? Any[e] : e isa AbstractString && startswith(e, "(") ?
                   Any[parse(Int, strip(s)) for s in split(strip(e, ['(', ')']), ",")] : Any[nothing]
    checked = 0
    for (model, entries, names) in table
        pos = reduce(vcat, positions.(entries))
        @test length(pos) == length(names)
        for (p, n) in zip(pos, names)
            (p === nothing || n === nothing) && continue
            @test DYR_PARAMETERS[model][p] == n
            checked += 1
        end
    end
    @test checked >= 130
end

@testset "RAW v34 is rewritten in the v35 layout" begin
    raw = joinpath(PSSEDIR, "WECC240.raw")
    path, rev = raw_for_parser(raw)
    @test rev == 34 && path != raw
    new = readlines(path); old = readlines(raw); rm(path)
    @test split_record(new[2])[1][3] == "35"
    # first data line of a section
    data(L, key) = findnext(l -> !startswith(lstrip(l), "@!"), L, findfirst(l -> occursin(key, l), L) + 1)
    # generators: NREG (absent here, default 0) after IREG, BASLOD (default 0) after PB
    f = split_record(new[data(new, "BEGIN GENERATOR")])[1]
    g = split_record(old[data(old, "BEGIN GENERATOR")])[1]
    @test length(g) == 28 && length(f) == 30
    @test f[1:8] == g[1:8] && f[9] == "0" && f[10:19] == g[9:18] && f[20] == "0" && f[21:30] == g[19:28]
    # transformer windings: NOD from the end of the line to after CONT
    i = data(new, "BEGIN TRANSFORMER"); j = data(old, "BEGIN TRANSFORMER")
    wnew = split_record(new[i+2])[1]; wold = split_record(old[j+2])[1]
    @test length(wnew) == length(wold) == 27
    @test wnew[1:17] == wold[1:17] && wnew[18] == wold[27] && wnew[19:27] == wold[18:26]
    s = split_record(new[data(new, "BEGIN SWITCHED SHUNT")])[1]
    @test s[1] == "4001" && s[2] == "'1 '" && s[12] == "600.00" && s[13:15] == ["1", "5", "200.00"]
    # versions 33 and 35 are read as they are
    @test raw_for_parser(joinpath(PSSEDIR, "twoarea.raw")) == (joinpath(PSSEDIR, "twoarea.raw"), 33)
end

@testset ".dyr data that do not fit NodalEMT's models are rejected" begin
    recs = read_dyr(joinpath(PSSEDIR, "twoarea.dyr"))
    units = [PsseUnit(b, "1", 700.0, 100.0, 900.0, 0.0) for b in 1:4]
    net = Dict{String,Any}("bus_num" => collect(1:4), "ws" => 2pi * 60)
    @test psse_case(net, units, recs) isa Tuple                                       # as bundled
    @test_throws r"exactly one governor" psse_case(net, units, filter(r -> !(r.bus == 2 && r.model == "TGOV1"), recs))
    @test_throws r"no SEXS" psse_case(net, units, filter(r -> !(r.bus == 3 && r.model == "SEXS"), recs))
    esst = NodalEMT.DyrRecord(4, "ESST1A", "1", Any[], 0)
    @test_throws r"unsupported models ESST1A" psse_case(net, units, vcat(recs, esst))
    @test_throws r"no GENROU or REGCA1" psse_case(net, vcat(units, PsseUnit(4, "2", 0.0, 0.0, 100.0, 0.0)), recs)
    # records of devices that are not in-service units are ignored with a warning
    stray = NodalEMT.DyrRecord(9, "CLOD", "1", Any[], 0)
    @test_logs (:warn, r"ignored .dyr records.*CLOD") psse_case(net, units, vcat(recs, stray))
end

@testset "the power flow of PowerFlows.jl reproduces PSS/E's solution: $name" for (name, raw, dyr, sizes, vtol, qtol) in (
        ("two-area, RAW v33", "twoarea.raw", "twoarea.dyr", (11, 4, 0, 8, 4), 1e-5, 1e-2),
        # WECC: up to 5.4 Mvar between units at buses joined by near-zero-impedance branches, which
        # PSS/E's stored, rounded solution does not balance (see the mismatch below)
        ("WECC, RAW v34", "WECC240.raw", "WECC240.dyr", (243, 111, 91, 329, 122), 2e-5, 10.0))
    solved = load_psse(joinpath(PSSEDIR, raw), joinpath(PSSEDIR, dyr))
    stored = load_psse(joinpath(PSSEDIR, raw), joinpath(PSSEDIR, dyr); solve_power_flow=false)
    p, q = solved.pfd, stored.pfd
    @test (nbus(p), ngen(p), nibr(p), length(p.line_from), length(p.xfmr_from)) == sizes
    for f in (:bus_num, :bus_basekV, :line_from, :line_to, :line_RX, :line_chg, :xfmr_from, :xfmr_to, :xfmr_RX,
              :xfmr_k, :gen_bus, :gen_id, :gen_MVA_base, :ibr_bus, :ibr_id, :load_bus, :shnt_bus, :shnt_gb)
        @test getfield(p, f) == getfield(q, f)                                # the same network and units
    end
    @test maximum(abs, p.bus_Vm .- q.bus_Vm) < vtol && maximum(abs, p.bus_Va .- q.bus_Va) < vtol
    @test maximum(abs, vcat(p.gen_MW .- q.gen_MW, p.ibr_MW .- q.ibr_MW)) < 1e-2
    @test maximum(abs, vcat(p.gen_Mvar .- q.gen_Mvar, p.ibr_Mvar .- q.ibr_Mvar)) < qtol
    @test same_dynamics(solved.dyd, stored.dyd)
    # the re-solved power flow balances the network better than PSS/E's rounded solution
    mis_solved = maximum(abs, power_flow_mismatch(p)); mis_stored = maximum(abs, power_flow_mismatch(q))
    @test mis_solved < 2e-3 && mis_solved < mis_stored
    # and an undisturbed EMT simulation stays at the operating point
    d = drift(solved)
    @test d < 1e-2
    @info "$name: power flow against PSS/E" max_dVm = maximum(abs, p.bus_Vm .- q.bus_Vm) mismatch_solved = mis_solved mismatch_stored = mis_stored emt_drift = d
end

"""
The two-area system with an off-nominal tap and a magnetizing branch, a line-end
shunt, a shunt reactor with conductance, a three-winding transformer to two new
buses and a synchronous condenser (PG = PT = PB = 0) on a new PV bus: writes
`path.raw` and `path.dyr`.
"""
function modified_twoarea(path)
    L = readlines(joinpath(PSSEDIR, "twoarea.raw"))
    before!(key, lines) = (i = findfirst(l -> occursin(key, l), L); foreach(l -> insert!(L, i, l), reverse(lines)))
    before!("End of Bus data", ["    12,'SEC         ',  20.0000,1,   1,   1,   1,1.00000, -25.0000,1.10000,0.90000,1.10000,0.90000",
                                "    13,'TERT        ',  20.0000,1,   1,   1,   1,1.00000, -25.0000,1.10000,0.90000,1.10000,0.90000",
                                "    14,'SC          ', 230.0000,2,   1,   1,   1,1.00000, -38.0000,1.10000,0.90000,1.10000,0.90000"])
    before!("End of Generator data", ["    14,'1 ',     0.000,     0.000,   100.000,  -100.000,1.00000,     0,   200.000, " *
                                      "0.00000E+0, 2.50000E-1, 0.00000E+0, 0.00000E+0,1.00000,1,  100.0,     0.000,     0.000,   1,1.0000"])
    before!("End of Branch data", ["     8,    14,'1 ', 1.00000E-3, 1.00000E-2,   0.00000,  400.00,  400.00,    0.00,  " *
                                   "0.00000,  0.00000,  0.00000,  0.00000,1,1,  10.00,   1,1.0000"])
    before!("End of Load data", ["    12,'1 ',1,   1,   1,    50.000,    10.000,     0.000,     0.000,     0.000,     0.000,   1,1,0"])
    before!("End of Fixed shunt data", ["     8,'R ',1,     2.000,   -40.000", "    13,'1 ',1,     0.000,   -10.000"])
    l = findfirst(x -> startswith(lstrip(x), "5,     6,'1 '"), L)
    L[l] = replace(L[l], "0.00,  0.00000,  0.00000,  0.00000,  0.00000" => "0.00,  0.00000,  0.02000,  0.00000,  0.00000")
    t = findfirst(x -> startswith(lstrip(x), "5,     1,     0"), L)
    L[t] = replace(L[t], "0.00000E+0, 0.00000E+0,1,'G1" => "1.00000E-3,-1.00000E-2,1,'G1")    # MAG1, MAG2
    L[t+2] = replace(L[t+2], r"^1\.00000" => "1.04000")                                       # WINDV1
    w(v, kv) = "$v, $kv,   0.000,   100.00,     0.00,     0.00, 0,      0, 1.10000, 0.90000, 1.10000, 0.90000,  33, 0, 0.00000, 0.00000,  0.000"
    before!("End of Transformer data", [
        "     7,    12,    13,'1 ',1,1,1, 0.00000E+0, 0.00000E+0,2,'T3W         ',1,   1,1.0000,   0,1.0000,   0,1.0000,   0,1.0000,'            '",
        " 1.00000E-3, 2.00000E-2,   100.00, 1.00000E-3, 3.00000E-2,   100.00, 1.00000E-3, 2.50000E-2,   100.00,1.00000, -25.0000",
        w("1.02000", "230.000"), w("1.00000", " 20.000"), w("0.97000", " 20.000")])
    write(path * ".raw", join(L, "\n") * "\n")
    # the condenser gets the machine models of generator 1
    D = readlines(joinpath(PSSEDIR, "twoarea.dyr"))
    first_unit = D[1:findfirst(l -> startswith(lstrip(l), "2 "), D)-1]
    write(path * ".dyr", join(vcat(D, [replace(l, r"^\s*1 '" => "    14 '") for l in first_unit]), "\n") * "\n")
    return path * ".raw", path * ".dyr"
end

@testset "taps, magnetizing branches, reactors, 3-winding transformers and condensers reach the EMT network" begin
    raw, dyr = modified_twoarea(tempname())
    @test length(collect(get_components(SynchronousCondenser, System(raw; runchecks=false)))) == 1
    c = load_psse(raw, dyr; vm_te=0.1)
    p = c.pfd
    k = findfirst(==(14), p.gen_bus)                                          # the condenser is a machine at P = 0
    @test k !== nothing && p.gen_MW[k] == 0 && abs(p.gen_Mvar[k]) > 1
    star = setdiff(p.bus_num, 1:14)
    @test length(star) == 1                                                    # the star bus of the 3-winding transformer
    @test sort([x for x in zip(p.xfmr_from, p.xfmr_to, p.xfmr_k) if x[2] == star[1]]) ==
          [(7, star[1], 1.02), (12, star[1], 1.0), (13, star[1], 0.97)]
    @test (5, 1, 1.04) in collect(zip(p.xfmr_from, p.xfmr_to, p.xfmr_k))
    @test (8, 2.0 - 40.0im) in collect(zip(p.shnt_bus, p.shnt_gb))             # reactor with conductance
    @test (5, 0.1 - 1.0im) in collect(zip(p.shnt_bus, p.shnt_gb))              # magnetizing branch at the from bus
    @test any(b == 5 && isapprox(y, 2.0im) for (b, y) in zip(p.shnt_bus, p.shnt_gb))   # line-end susceptance BI
    # NodalEMT's network is the one PowerFlows.jl solved: same admittance matrix up to its single precision
    Yp = PNM.Ybus(System(raw; runchecks=false))
    perm = [findfirst(==(b), p.bus_num) for b in PNM.get_bus_axis(Yp)]
    Y = network_admittance(p)[perm, perm]
    @test maximum(abs, Matrix(Y) .- ComplexF64.(Matrix(Yp.data))) < 1e-6 * maximum(abs, Yp.data)
    @test maximum(abs, power_flow_mismatch(p)) < 1e-4
    rm(raw); rm(dyr)
    # an undisturbed EMT run with pi lines stays at the operating point; without the taps it would not
    ok = drift(c; bergeron_threshold=Inf)
    bad = deepcopy(c); fill!(bad.pfd.xfmr_k, 1.0)
    notaps = drift(bad; bergeron_threshold=Inf)
    @test ok < 1e-3
    @test notaps > 100 * ok
    @info "taps, reactors, 3-winding transformer, condenser" drift = ok drift_without_taps = notaps
end

@testset "a case written to JSON reads back identically" begin
    json = tempname() * ".json"
    c = load_psse(joinpath(PSSEDIR, "twoarea.raw"), joinpath(PSSEDIR, "twoarea.dyr"); json=json)
    r = load_case(json); rm(json)
    for f in fieldnames(typeof(c.pfd))
        @test getfield(c.pfd, f) == getfield(r.pfd, f)
    end
    @test same_dynamics(c.dyd, r.dyd)
    a = build_simulation(c); b = build_simulation(r)
    for _ in 1:50; step!(a); step!(b); end
    @test a.net.V == b.net.V
end
