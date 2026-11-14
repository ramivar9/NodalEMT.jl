# ==========================================================================
#  psse.jl -- cases from PSS/E .raw and .dyr files
#
#  load_psse(raw, dyr) builds a case end to end:
#    1. the .raw file is read by PowerSystems.jl and its power flow solved by
#       PowerFlows.jl (power_flow.jl);
#    2. the .dyr file is read here (`read_dyr`);
#    3. both are turned into the "pfd" and "dyd" dictionaries of a JSON case
#       file (`psse_case`), so that a PSS/E case and a JSON case are built by
#       the same code (PowerFlowData and DynamicData in case_data.jl).
#
#  This file holds the parts that need no other package: the .dyr reader, the
#  assembly of the case and the conversion of PSS/E v34 .raw files to the v35
#  layout read by PowerSystems.jl.
# ==========================================================================

# --------------------------------------------------------------------------
# Free-format records (shared by .raw and .dyr)
# --------------------------------------------------------------------------

"""
    split_record(line) -> (fields, comment)

Splits one line of a PSS/E file into its fields and the comment after `/`.
Fields are separated by commas or blanks; text in single quotes is one field
(quotes kept); two commas in a row give an empty field.
"""
function split_record(line::AbstractString)
    pieces = String[]; buf = IOBuffer(); inq = false; comment = ""
    chars = collect(line)
    i = 1
    while i <= length(chars)
        c = chars[i]
        if c == '\''
            inq = !inq; write(buf, c)
        elseif !inq && c == '/'
            comment = String(chars[i:end]); break
        elseif !inq && c == ','
            push!(pieces, String(take!(buf)))
        else
            write(buf, c)
        end
        i += 1
    end
    last = String(take!(buf))
    (!isempty(strip(last)) || !isempty(pieces)) && push!(pieces, last)
    fields = String[]
    for p in pieces
        toks = blank_split(p)
        isempty(toks) ? push!(fields, "") : append!(fields, toks)
    end
    # a trailing empty field after the last comma carries no data
    !isempty(fields) && isempty(fields[end]) && !isempty(pieces) && isempty(strip(pieces[end])) && pop!(fields)
    return fields, comment
end

"Splits at blanks outside single quotes."
function blank_split(s::AbstractString)
    toks = String[]; buf = IOBuffer(); inq = false
    for c in s
        if c == '\''
            inq = !inq; write(buf, c)
        elseif !inq && isspace(c)
            t = String(take!(buf)); isempty(t) || push!(toks, t)
        else
            write(buf, c)
        end
    end
    t = String(take!(buf)); isempty(t) || push!(toks, t)
    return toks
end

unquote(s::AbstractString) = (length(s) >= 2 && s[1] == '\'' && s[end] == '\'') ? String(s[2:end-1]) : String(s)

"A number in PSS/E notation (Fortran `D` exponents allowed), or `nothing`."
parse_number(s::AbstractString) = tryparse(Float64, replace(s, r"[dD]" => "e"))

# --------------------------------------------------------------------------
# .dyr files
# --------------------------------------------------------------------------

"""
Parameter names of the supported models, in the order of their ICONs and
CONs in a .dyr record (PSS/E model library).  The names are those of the
"dyd" section of a case file.
"""
const DYR_PARAMETERS = Dict{String,Vector{String}}(
    "GENROU" => ["Td0p", "Td0pp", "Tq0p", "Tq0pp", "H", "D", "Xd", "Xq", "Xdp", "Xqp", "Xdpp", "Xl", "S10", "S12"],
    "SEXS"   => ["TA_o_TB", "TB", "K", "TE", "Emin", "Emax"],
    "TGOV1"  => ["R", "T1", "Vmax", "Vmin", "T2", "T3", "Dt"],
    "HYGOV"  => ["R", "r", "Tr", "Tf", "Tg", "VELM", "GMAX", "GMIN", "TW", "At", "Dturb", "qNL"],
    "GAST"   => ["R", "T1", "T2", "T3", "LdLmt", "KT", "VMAX", "VMIN", "Dturb"],
    "IEEEST" => ["ICS", "IB", "A1", "A2", "A3", "A4", "A5", "A6", "T1", "T2", "T3", "T4", "T5", "T6",
                 "KS", "LSMAX", "LSMIN", "VCU", "VCL"],
    "REGCA1" => ["LVPLsw", "Tg", "Rrpwr", "Brkpt", "Zerox", "Lvpl1", "Volim", "Lvpnt1", "Lvpnt0", "Iolim",
                 "Tfltr", "Khv", "Iqrmax", "Iqrmin", "Accel"],
    "REECB1" => ["IREG", "PFFLAG", "VFLAG", "QFLAG", "PQFLAG", "Vdip", "Vup", "Trv", "dbd1", "dbd2", "Kqv",
                 "Iqhl", "Iqll", "Vref0", "Tp", "Qmax", "Qmin", "Vmax", "Vmin", "Kqp", "Kqi", "Kvp", "Kvi",
                 "Tiq", "dPmax", "dPmin", "Pmax", "Pmin", "Imax", "Tpord"],
    "REPCA1" => ["remote_bus", "branch_From_bus", "branch_To_bus", "branch_id", "VCFlag", "RefFlag", "FFlag",
                 "Tfltr", "Kp", "Ki", "Tft", "Tfv", "Vfrz", "Rc", "Xc", "Kc", "emax", "emin", "dbd1", "dbd2",
                 "Qmax", "Qmin", "Kpg", "Kig", "Tp", "fdbd1", "fdbd2", "femax", "femin", "Pmax", "Pmin",
                 "Tg", "Ddn", "Dup"],
)

const GOVERNOR_MODELS = ("TGOV1", "HYGOV", "GAST")

"""
    DyrRecord

One record of a .dyr file: `bus`, `model` (upper case), unit `id` (without
quotes and blanks), the `values` (`Float64`, or `String` for quoted text) and
the `line` on which it starts.
"""
struct DyrRecord
    bus::Int
    model::String
    id::String
    values::Vector{Any}
    line::Int
end

"""
    read_dyr(path) -> Vector{DyrRecord}

Reads every record `IBUS 'MODEL' ID values... /` of a PSS/E .dyr file.  A
record may span several lines; `@!` lines and the text after `/` are comments.
"""
function read_dyr(path::AbstractString)
    recs = DyrRecord[]
    toks = String[]; start = 0
    for (ln, line) in enumerate(eachline(path))
        startswith(lstrip(line), "@!") && continue
        fields, comment = split_record(line)
        filter!(!isempty, fields)
        !isempty(fields) && isempty(toks) && (start = ln)
        append!(toks, fields)
        if !isempty(comment)                     # '/' ends the record
            isempty(toks) || push!(recs, dyr_record(toks, start, path))
            empty!(toks)
        end
    end
    isempty(toks) || error("$path: the record starting on line $start is not terminated by '/'")
    return recs
end

function dyr_record(toks::Vector{String}, line::Int, path)
    length(toks) >= 3 || error("$path, line $line: a record needs a bus, a model name and an id")
    bus = parse_number(toks[1])
    bus === nothing && error("$path, line $line: '$(toks[1])' is not a bus number")
    vals = Any[]
    for t in toks[4:end]
        if startswith(t, "'")
            push!(vals, unquote(t))
        else
            x = parse_number(t)
            x === nothing && error("$path, line $line: '$t' is not a number")
            push!(vals, x)
        end
    end
    return DyrRecord(round(Int, bus), uppercase(strip(unquote(toks[2]))), strip(unquote(toks[3])), vals, line)
end

"Parameters of a supported record by name; checks their number."
function dyr_parameters(r::DyrRecord)
    names = DYR_PARAMETERS[r.model]
    length(r.values) == length(names) ||
        error(".dyr line $(r.line): $(r.model) of unit ($(r.bus), '$(r.id)') has $(length(r.values)) values, " *
              "expected $(length(names))")
    return Dict(zip(names, r.values))
end

# --------------------------------------------------------------------------
# Assembling a case
# --------------------------------------------------------------------------

"""
    PsseUnit(bus, id, P, Q, mbase, ra)

A generating unit of the power flow: bus number, machine id, output [MW,
Mvar], MVA base and the resistance of its source impedance (ZSORCE, machine
base).  Whether it is a synchronous machine or an IBR follows from its .dyr
models.
"""
struct PsseUnit
    bus::Int
    id::String
    P::Float64
    Q::Float64
    mbase::Float64
    ra::Float64
end

"""
    psse_case(network, units, records; pll_ke, pll_te, vm_te) -> (pfd, dyd)

The "pfd" and "dyd" dictionaries of a case: `network` holds every "pfd" entry
except the units, `units` are the in-service units in the order in which
they are to appear, `records` the .dyr records.  A unit with a GENROU record
is a synchronous machine, a unit with a REGCA1 record an IBR.
"""
function psse_case(network::AbstractDict, units::Vector{PsseUnit}, records::Vector{DyrRecord};
                   pll_ke::Real=0.1, pll_te::Real=0.05, vm_te::Real=0.02)
    # records by unit and model
    byunit = Dict{Tuple{Int,String},Dict{String,DyrRecord}}()
    for r in records
        d = get!(byunit, (r.bus, r.id), Dict{String,DyrRecord}())
        haskey(d, r.model) && error(".dyr lines $(d[r.model].line) and $(r.line): two $(r.model) records for unit ($(r.bus), '$(r.id)')")
        d[r.model] = r
    end
    unitkeys = Set((u.bus, u.id) for u in units)
    ignored = Dict{String,Int}()
    for (k, d) in byunit, m in keys(d)
        k in unitkeys || (ignored[m] = get(ignored, m, 0) + 1)
    end

    machines = PsseUnit[]; ibrs = PsseUnit[]; problems = String[]
    found(d) = join(sort(collect(keys(d))), ", ")
    for u in units
        d = get(byunit, (u.bus, u.id), Dict{String,DyrRecord}())
        name = "unit ($(u.bus), '$(u.id)')"
        unsupported = [m for m in keys(d) if !haskey(DYR_PARAMETERS, m)]
        isempty(unsupported) || push!(problems, "$name: unsupported models $(join(sort(unsupported), ", "))")
        if haskey(d, "GENROU") && haskey(d, "REGCA1")
            push!(problems, "$name has both GENROU and REGCA1")
        elseif haskey(d, "GENROU")
            push!(machines, u)
            haskey(d, "SEXS") || push!(problems, "$name has no SEXS exciter (found: $(found(d)))")
            ng = count(m -> haskey(d, m), GOVERNOR_MODELS)
            ng == 1 || push!(problems, "$name needs exactly one governor (TGOV1, HYGOV or GAST), found: $(found(d))")
            any(haskey(d, m) for m in ("REECB1", "REPCA1")) && push!(problems, "$name: IBR models on a GENROU unit")
        elseif haskey(d, "REGCA1")
            push!(ibrs, u)
            for m in ("REECB1", "REPCA1")
                haskey(d, m) || push!(problems, "$name has no $m (found: $(found(d)))")
            end
            any(haskey(d, m) for m in ("SEXS", "IEEEST", GOVERNOR_MODELS...)) &&
                push!(problems, "$name: machine controls on a REGCA1 unit")
        else
            push!(problems, "$name has no GENROU or REGCA1 model" * (isempty(d) ? "" : " (found: $(found(d)))"))
        end
    end
    isempty(problems) || error("the .dyr data do not fit NodalEMT's models:\n  " * join(problems, "\n  "))
    isempty(ignored) || @warn "ignored .dyr records of devices that are not in-service units: " *
        join(("$m ($n)" for (m, n) in sort(collect(ignored))), ", ")

    pfd = Dict{String,Any}(network)
    pfd["gen_bus"] = [u.bus for u in machines]; pfd["gen_id"] = [u.id for u in machines]
    pfd["gen_MW"] = [u.P for u in machines]; pfd["gen_Mvar"] = [u.Q for u in machines]
    pfd["gen_MVA_base"] = [u.mbase for u in machines]
    pfd["ibr_bus"] = [u.bus for u in ibrs]; pfd["ibr_id"] = [u.id for u in ibrs]
    pfd["ibr_MW"] = [u.P for u in ibrs]; pfd["ibr_Mvar"] = [u.Q for u in ibrs]
    pfd["ibr_MVA_base"] = [u.mbase for u in ibrs]

    dyd = Dict{String,Any}()
    rec(u, m) = byunit[(u.bus, u.id)][m]
    function put!(prefix, m, list)
        dyd["$(prefix)_bus"] = [u.bus for u in list]
        dyd["$(prefix)_id"] = [u.id for u in list]
        ps = [dyr_parameters(rec(u, m)) for u in list]
        for p in DYR_PARAMETERS[m]
            dyd["$(prefix)_$p"] = [x[p] for x in ps]
        end
        return ps
    end
    # machines: one GENROU and one SEXS per machine, in the machine order
    g = put!("gen_genrou", "GENROU", machines)
    dyd["gen_H"] = dyd["gen_genrou_H"]; dyd["gen_D"] = dyd["gen_genrou_D"]
    dyd["gen_Ra"] = [u.ra for u in machines]
    dyd["gen_X0"] = copy(dyd["gen_genrou_Xl"])          # as ParaEMT
    nsat = count(x -> x["S10"] != 0 || x["S12"] != 0, g)
    nsat > 0 && @warn "GENROU saturation is not modelled; $nsat machines have nonzero S(1.0) or S(1.2)"
    put!("exc_sexs", "SEXS", machines)
    put!("gov_tgov1", "TGOV1", [u for u in machines if haskey(byunit[(u.bus, u.id)], "TGOV1")])
    put!("gov_hygov", "HYGOV", [u for u in machines if haskey(byunit[(u.bus, u.id)], "HYGOV")])
    put!("gov_gast", "GAST", [u for u in machines if haskey(byunit[(u.bus, u.id)], "GAST")])
    pss = put!("pss_ieeest", "IEEEST", [u for u in machines if haskey(byunit[(u.bus, u.id)], "IEEEST")])
    for (u, p) in zip([u for u in machines if haskey(byunit[(u.bus, u.id)], "IEEEST")], pss)
        (p["ICS"] == 1 && p["IB"] == 0) ||
            error("IEEEST of unit ($(u.bus), '$(u.id)'): only the local rotor speed is supported as input (ICS = 1, IB = 0)")
    end
    # IBRs: one REGCA1, REECB1 and REPCA1 per IBR, in the IBR order
    put!("ibr_regca", "REGCA1", ibrs)
    reecb = put!("ibr_reecb", "REECB1", ibrs)
    for (u, p) in zip(ibrs, reecb)
        p["IREG"] == 0 || error("REECB1 of unit ($(u.bus), '$(u.id)'): remote voltage regulation (IREG = $(Int(p["IREG"]))) is not supported")
    end
    put!("ibr_repca", "REPCA1", ibrs)
    dyd["ibr_MVAbase"] = [u.mbase for u in ibrs]
    dyd["ibr_fbase"] = fill(network["ws"] / 2pi, length(ibrs))
    # PLL and voltage filter of every bus
    nb = length(network["bus_num"])
    dyd["pll_bus"] = copy(network["bus_num"]); dyd["vm_bus"] = copy(network["bus_num"])
    dyd["pll_ke"] = fill(Float64(pll_ke), nb); dyd["pll_te"] = fill(Float64(pll_te), nb)
    dyd["vm_te"] = fill(Float64(vm_te), nb)
    return pfd, dyd
end

"""
    write_case_json(path, name, pfd, dyd)

Writes a case in the JSON format read by [`load_case`](@ref).
"""
function write_case_json(path::AbstractString, name::AbstractString, pfd::AbstractDict, dyd::AbstractDict)
    open(path, "w") do io
        JSON.print(io, Dict("meta" => Dict("name" => name, "source" => "NodalEMT.load_psse"),
                            "pfd" => pfd, "dyd" => dyd))
    end
    return path
end

"Complex vector in the case-file form `{\"re\": [...], \"im\": [...]}`."
complex_entry(v::AbstractVector) = Dict("re" => real.(v), "im" => imag.(v))

# --------------------------------------------------------------------------
# .raw files: version 34 to the layout read by PowerSystems.jl
# --------------------------------------------------------------------------
#
# PowerSystems.jl reads RAW versions 33 and 35.  Version 34 has the sections
# of version 35, but a few record types differ: generators (NREG at the end,
# no BASLOD), loads (no LOADTYPE), switched shunts (no ID and no block
# status), FACTS devices (NREG after MNAME) and transformer windings (NOD at
# the end of the line).  `raw_for_parser` rewrites those records by field
# name into the v35 layout; every other line is copied unchanged.

const RAW34_FIELDS = Dict(
    "LOAD" => split("I ID STAT AREA ZONE PL QL IP IQ YP YQ OWNER SCALE INTRPT DGENP DGENQ DGENM"),
    "GENERATOR" => split("I ID PG QG QT QB VS IREG MBASE ZR ZX RT XT GTAP STAT RMPCT PT PB O1 F1 O2 F2 O3 F3 O4 F4 WMOD WPF NREG"),
    "FACTS DEVICE" => split("NAME I J MODE PDES QDES VSET SHMX TRMX VTMN VTMX VSMX IMX LINX RMPCT OWNER SET1 SET2 VSREF FCREG MNAME NREG"),
    "SWITCHED SHUNT" => split("I MODSW ADJM STAT VSWHI VSWLO SWREM RMPCT RMIDNT BINIT N1 B1 N2 B2 N3 B3 N4 B4 N5 B5 N6 B6 N7 B7 N8 B8 NREG"),
)
const RAW35_FIELDS = Dict(
    "LOAD" => split("I ID STAT AREA ZONE PL QL IP IQ YP YQ OWNER SCALE INTRPT DGENP DGENQ DGENM LOADTYPE"),
    "GENERATOR" => split("I ID PG QG QT QB VS IREG NREG MBASE ZR ZX RT XT GTAP STAT RMPCT PT PB BASLOD O1 F1 O2 F2 O3 F3 O4 F4 WMOD WPF"),
    "FACTS DEVICE" => split("NAME I J MODE PDES QDES VSET SHMX TRMX VTMN VTMX VSMX IMX LINX RMPCT OWNER SET1 SET2 VSREF FCREG NREG MNAME"),
    "SWITCHED SHUNT" => split("I ID MODSW ADJM STAT VSWHI VSWLO SWREM NREG RMPCT RMIDNT BINIT " *
                              join(("S$k N$k B$k" for k in 1:8), " ")),
)
"Values of the v35 fields that a v34 record lacks (PSS/E defaults; PowerSystems.jl fills no default for them)."
const RAW35_NEW_DEFAULTS = Dict("ID" => "'1 '", "NREG" => "0", "BASLOD" => "0")

"Sections of a v34 file after the system-wide data, in file order."
const RAW34_SECTIONS = ["BUS", "LOAD", "FIXED SHUNT", "GENERATOR", "BRANCH", "SYSTEM SWITCHING DEVICE",
    "TRANSFORMER", "AREA", "TWO-TERMINAL DC", "VSC DC LINE", "IMPEDANCE CORRECTION", "MULTI-TERMINAL DC",
    "MULTI-SECTION LINE", "ZONE", "INTER-AREA TRANSFER", "OWNER", "FACTS DEVICE", "SWITCHED SHUNT", "GNE",
    "INDUCTION MACHINE", "SUBSTATION"]

"The RAW version (REV) of a .raw file."
function raw_version(lines::Vector{String})
    for l in lines
        startswith(lstrip(l), "@!") && continue
        f, _ = split_record(l)
        return length(f) >= 3 && parse_number(f[3]) !== nothing ? round(Int, parse_number(f[3])) : 33
    end
    error("empty .raw file")
end

"""
    raw_for_parser(path) -> (path_to_read, version)

The file to hand to PowerSystems.jl: `path` itself for versions 33 and 35, a
rewritten temporary copy for version 34 (and for a v33 file with `@!`
comments, which PowerSystems.jl would take for v35).
"""
function raw_for_parser(path::AbstractString)
    lines = readlines(path)
    rev = raw_version(lines)
    if rev <= 33
        any(l -> startswith(lstrip(l), "@!"), lines) || return path, rev
        out = [l for l in lines if !startswith(lstrip(l), "@!")]
    elseif rev == 34
        out = raw34_to_35(lines)
    else
        rev == 35 || @warn "RAW version $rev: read as version 35"
        return path, rev
    end
    tmp = tempname() * ".raw"
    write(tmp, join(out, "\n") * "\n")
    return tmp, rev
end

"""
Rewrites a record from the field order `from` to `to` by field name.  A
missing field is empty (the parser's default), except the new v35 fields in
RAW35_NEW_DEFAULTS and the status S_k of each switched-shunt block present.
"""
function reorder(fields, from, to)
    d = Dict(zip(from, fields))
    for k in 1:8
        isempty(strip(get(d, "N$k", ""))) || (d["S$k"] = "1")
    end
    return [get(d, n, get(RAW35_NEW_DEFAULTS, n, "")) for n in to]
end

"Moves NOD of a v34 transformer winding line from the end to after CONT."
function winding34_to_35(fields)
    f = copy(fields)
    if length(f) >= 27
        nod = f[27]; deleteat!(f, 27); insert!(f, 18, nod)
    elseif length(f) >= 18
        insert!(f, 18, "")
    end
    return f
end

function raw34_to_35(lines::Vector{String})
    out = String[]
    i = 1
    # comments, the case identification (REV -> 35) and two title lines
    while i <= length(lines) && startswith(lstrip(lines[i]), "@!")
        push!(out, lines[i]); i += 1
    end
    isempty(out) && push!(out, "@!IC,SBASE,REV,XFRRAT,NXFRAT,BASFRQ")
    f, c = split_record(lines[i]); f[3] = "35"
    push!(out, join(f, ", ") * (isempty(c) ? "" : " " * c)); i += 1
    for _ in 1:2
        i <= length(lines) && (push!(out, lines[i]); i += 1)
    end
    # system-wide data, if present, up to its terminator
    first_data = findnext(l -> !startswith(lstrip(l), "@!"), lines, i)
    if first_data !== nothing && occursin(r"^\s*(GENERAL|GAUSS|NEWTON|ADJUST|TYSL|SOLVER|RATING)", lines[first_data])
        while i <= length(lines)
            push!(out, lines[i]); i += 1
            occursin("END OF SYSTEM-WIDE DATA", lines[i-1]) && break
        end
    end
    s = 1                                         # index into RAW34_SECTIONS
    while i <= length(lines)
        line = lines[i]
        if startswith(lstrip(line), "@!") || s > length(RAW34_SECTIONS)
            push!(out, line); i += 1; continue
        end
        f, c = split_record(line)
        sec = RAW34_SECTIONS[s]
        if !isempty(f) && (f[1] == "0" || uppercase(f[1]) == "Q")    # end of a section
            push!(out, line); i += 1
            m = match(r"BEGIN (.+?) DATA", uppercase(c))
            nxt = m === nothing ? nothing : findfirst(==(strip(m.captures[1])), RAW34_SECTIONS)
            s = nxt === nothing ? s + 1 : nxt
            continue
        end
        if haskey(RAW34_FIELDS, sec)
            push!(out, join(reorder(f, RAW34_FIELDS[sec], RAW35_FIELDS[sec]), ", ")); i += 1
        elseif sec == "TRANSFORMER"
            three = length(f) >= 3 && parse_number(f[3]) !== nothing && parse_number(f[3]) != 0
            push!(out, lines[i]); push!(out, lines[i+1])              # line 1 and the impedances
            nw = three ? 3 : 1                                      # full winding lines
            for w in 1:nw
                fw, _ = split_record(lines[i+1+w])
                push!(out, join(winding34_to_35(fw), ", "))
            end
            three || push!(out, lines[i+3])                         # WINDV2, NOMV2
            i += three ? 5 : 4
        else
            push!(out, line); i += 1
        end
    end
    return out
end
