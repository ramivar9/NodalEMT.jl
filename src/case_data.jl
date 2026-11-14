# ==========================================================================
#  case_data.jl -- input data of a simulation case
#
#  A case is one JSON file with two parts:
#
#    "pfd"  PowerFlowData -- the solved power flow: buses, lines,
#           transformers, loads, shunts, synchronous generators and IBRs with
#           their dispatch.
#    "dyd"  DynamicData   -- the parameters of every dynamic model: GENROU
#           machine, SEXS exciter, TGOV1/HYGOV/GAST governors, IEEEST
#           stabilizer, REGCA/REECB/REPCA inverter models, and the PLL and
#           voltage filter of every bus.
#
#  A case is built from PSS/E files by load_psse (psse.jl), which can also
#  save it in this JSON format (`json` keyword) to be read again by load_case
#  without solving the power flow again.
#  Everything is stored as structs of plain vectors, so that the parameters
#  of each model can be read at a glance.  Conventions:
#    * all indices are 1-based;
#    * "generator i" is the i-th entry of pfd.gen_bus; GENROU machine i and
#      SEXS exciter i belong to generator i (ParaEMT's convention);
#    * governors and stabilizers are attached to generators by (bus, id).
# ==========================================================================

using JSON

# --------------------------------------------------------------------------
# Small JSON helpers
# --------------------------------------------------------------------------

jvec(d, key) = Float64.(d[key])                      # numeric vector
jivec(d, key) = Int.(round.(Float64.(d[key])))       # integer vector
jsvec(d, key) = String.(d[key])                      # string vector

"Complex vector, stored in the case file as `{\"re\": [...], \"im\": [...]}` (or `[]` when empty)."
function jcvec(d, key)
    v = d[key]
    v isa AbstractDict || (isempty(v) ? (return ComplexF64[]) : error("$key: expected {\"re\": ..., \"im\": ...}"))
    return ComplexF64.(Float64.(v["re"]), Float64.(v["im"]))
end

# --------------------------------------------------------------------------
# Power-flow data
# --------------------------------------------------------------------------

"""
    PowerFlowData

Solved power flow of the system: powers in MW and Mvar, voltages in pu,
angles in rad, impedances in pu on `basemva`.
"""
struct PowerFlowData
    basemva::Float64
    ws::Float64                    # synchronous angular frequency 2*pi*60 [rad/s]
    # buses
    bus_num::Vector{Int}
    bus_name::Vector{String}
    bus_type::Vector{Int}
    bus_Vm::Vector{Float64}
    bus_Va::Vector{Float64}
    bus_basekV::Vector{Float64}
    # synchronous generators
    gen_bus::Vector{Int}
    gen_id::Vector{String}
    gen_MW::Vector{Float64}
    gen_Mvar::Vector{Float64}
    gen_MVA_base::Vector{Float64}
    # inverter-based resources
    ibr_bus::Vector{Int}
    ibr_id::Vector{String}
    ibr_MW::Vector{Float64}
    ibr_Mvar::Vector{Float64}
    ibr_MVA_base::Vector{Float64}
    # lines (series R + jX and total charging B)
    line_from::Vector{Int}
    line_to::Vector{Int}
    line_id::Vector{String}
    line_RX::Vector{ComplexF64}
    line_chg::Vector{Float64}
    # transformers (series R + jX)
    xfmr_from::Vector{Int}
    xfmr_to::Vector{Int}
    xfmr_id::Vector{String}
    xfmr_RX::Vector{ComplexF64}
    xfmr_k::Vector{Float64}
    # loads
    load_bus::Vector{Int}
    load_id::Vector{String}
    load_MW::Vector{Float64}
    load_Mvar::Vector{Float64}
    # fixed shunts (G + jB in MW/Mvar at 1 pu voltage)
    shnt_bus::Vector{Int}
    shnt_id::Vector{String}
    shnt_gb::Vector{ComplexF64}
    # switched shunts (Mvar at the actual bus voltage)
    shnt_sw_bus::Vector{Int}
    shnt_sw_gb::Vector{ComplexF64}
end

function PowerFlowData(d::AbstractDict)
    return PowerFlowData(
        Float64(d["basemva"]), Float64(d["ws"]),
        jivec(d, "bus_num"), jsvec(d, "bus_name"), jivec(d, "bus_type"),
        jvec(d, "bus_Vm"), jvec(d, "bus_Va"), jvec(d, "bus_basekV"),
        jivec(d, "gen_bus"), jsvec(d, "gen_id"), jvec(d, "gen_MW"), jvec(d, "gen_Mvar"), jvec(d, "gen_MVA_base"),
        jivec(d, "ibr_bus"), jsvec(d, "ibr_id"), jvec(d, "ibr_MW"), jvec(d, "ibr_Mvar"), jvec(d, "ibr_MVA_base"),
        jivec(d, "line_from"), jivec(d, "line_to"), jsvec(d, "line_id"), jcvec(d, "line_RX"), jvec(d, "line_chg"),
        jivec(d, "xfmr_from"), jivec(d, "xfmr_to"), jsvec(d, "xfmr_id"), jcvec(d, "xfmr_RX"), jvec(d, "xfmr_k"),
        jivec(d, "load_bus"), jsvec(d, "load_id"), jvec(d, "load_MW"), jvec(d, "load_Mvar"),
        jivec(d, "shnt_bus"), jsvec(d, "shnt_id"), jcvec(d, "shnt_gb"),
        jivec(d, "shnt_sw_bus"), jcvec(d, "shnt_sw_gb"),
    )
end

nbus(pfd::PowerFlowData) = length(pfd.bus_num)
ngen(pfd::PowerFlowData) = length(pfd.gen_bus)
nibr(pfd::PowerFlowData) = length(pfd.ibr_bus)

"Position of bus number `busnum` in the bus list."
function bus_index(pfd::PowerFlowData, busnum::Integer)
    k = findfirst(==(busnum), pfd.bus_num)
    k === nothing && error("bus $busnum not found in the power-flow data")
    return k
end

"Generator with bus number `busnum` and PSS/E id `id` (used to attach governors and stabilizers)."
function generator_index(pfd::PowerFlowData, busnum::Integer, id::AbstractString)
    for i in 1:ngen(pfd)
        pfd.gen_bus[i] == busnum && strip(pfd.gen_id[i]) == strip(id) && return i
    end
    error("generator bus=$busnum id='$id' not found in the power-flow data")
end

# --------------------------------------------------------------------------
# Dynamic-model parameters (one struct per model; one vector entry per unit)
# --------------------------------------------------------------------------

"Governor model types."
@enum GovernorType GAST = 0 HYGOV = 1 TGOV1 = 2

"GENROU round-rotor synchronous machine (one entry per generator, machine base)."
struct GenrouParams
    bus::Vector{Int}
    id::Vector{String}
    Td0p::Vector{Float64};  Td0pp::Vector{Float64}     # open-circuit time constants [s]
    Tq0p::Vector{Float64};  Tq0pp::Vector{Float64}
    H::Vector{Float64};     D::Vector{Float64}         # inertia [s] and damping
    Xd::Vector{Float64};    Xq::Vector{Float64}        # reactances [pu]
    Xdp::Vector{Float64};   Xqp::Vector{Float64}
    Xdpp::Vector{Float64};  Xl::Vector{Float64}
    Ra::Vector{Float64}                                # stator resistance [pu]
    X0::Vector{Float64}                                # zero-sequence reactance [pu]
end

"SEXS simplified excitation system (one entry per generator)."
struct SexsParams
    bus::Vector{Int}
    id::Vector{String}
    TA::Vector{Float64}     # lead time constant (= TB * TA/TB of the data sheet)
    TB::Vector{Float64}     # lag time constant
    K::Vector{Float64}      # gain
    TE::Vector{Float64}     # exciter time constant
    Emin::Vector{Float64}
    Emax::Vector{Float64}
end

"TGOV1 steam turbine-governor; `gen[k]` is the generator that unit k drives."
struct Tgov1Params
    gen::Vector{Int}
    R::Vector{Float64};  T1::Vector{Float64}
    Vmax::Vector{Float64}; Vmin::Vector{Float64}
    T2::Vector{Float64}; T3::Vector{Float64}; Dt::Vector{Float64}
end

"HYGOV hydro turbine-governor."
struct HygovParams
    gen::Vector{Int}
    R::Vector{Float64};  r::Vector{Float64}
    Tr::Vector{Float64}; Tf::Vector{Float64}; Tg::Vector{Float64}
    VELM::Vector{Float64}; GMAX::Vector{Float64}; GMIN::Vector{Float64}
    TW::Vector{Float64}; At::Vector{Float64}; Dturb::Vector{Float64}; qNL::Vector{Float64}
end

"GAST gas turbine-governor."
struct GastParams
    gen::Vector{Int}
    R::Vector{Float64}; T1::Vector{Float64}; T2::Vector{Float64}; T3::Vector{Float64}
    LdLmt::Vector{Float64}; KT::Vector{Float64}
    VMAX::Vector{Float64}; VMIN::Vector{Float64}; Dturb::Vector{Float64}
end

"IEEEST power system stabilizer with speed input."
struct IeeestParams
    gen::Vector{Int}
    A1::Vector{Float64}; A2::Vector{Float64}; A3::Vector{Float64}
    A4::Vector{Float64}; A5::Vector{Float64}; A6::Vector{Float64}
    T1::Vector{Float64}; T2::Vector{Float64}; T3::Vector{Float64}
    T4::Vector{Float64}; T5::Vector{Float64}; T6::Vector{Float64}
    KS::Vector{Float64}; LSMAX::Vector{Float64}; LSMIN::Vector{Float64}
    VCU::Vector{Float64}
end

"REGCA renewable generator/converter model (one entry per IBR)."
struct RegcaParams
    bus::Vector{Int}
    id::Vector{String}
    Tg::Vector{Float64}; Rrpwr::Vector{Float64}
    Brkpt::Vector{Float64}; Zerox::Vector{Float64}
    Volim::Vector{Float64}; Lvpnt1::Vector{Float64}; Lvpnt0::Vector{Float64}
    Tfltr::Vector{Float64}; Khv::Vector{Float64}
    Iqrmax::Vector{Float64}; Iqrmin::Vector{Float64}
    MVAbase::Vector{Float64}   # MVA base of the dynamic model
end

"REECB renewable electrical control model."
struct ReecbParams
    PFFLAG::Vector{Int}; VFLAG::Vector{Int}; QFLAG::Vector{Int}; PQFLAG::Vector{Int}
    Vdip::Vector{Float64}; Vup::Vector{Float64}; Trv::Vector{Float64}
    dbd1::Vector{Float64}; dbd2::Vector{Float64}; Kqv::Vector{Float64}
    Iqhl::Vector{Float64}; Iqll::Vector{Float64}; Vref0::Vector{Float64}
    Tp::Vector{Float64}; Qmax::Vector{Float64}; Qmin::Vector{Float64}
    Vmax::Vector{Float64}; Vmin::Vector{Float64}
    Kqp::Vector{Float64}; Kqi::Vector{Float64}; Kvp::Vector{Float64}; Kvi::Vector{Float64}
    Tiq::Vector{Float64}; dPmax::Vector{Float64}; dPmin::Vector{Float64}
    Pmax::Vector{Float64}; Pmin::Vector{Float64}; Imax::Vector{Float64}; Tpord::Vector{Float64}
end

"REPCA renewable plant controller model."
struct RepcaParams
    remote_bus::Vector{Int}; branch_From_bus::Vector{Int}; branch_To_bus::Vector{Int}
    VCFlag::Vector{Int}; RefFlag::Vector{Int}; FFlag::Vector{Int}
    Tfltr::Vector{Float64}; Kp::Vector{Float64}; Ki::Vector{Float64}
    Tft::Vector{Float64}; Tfv::Vector{Float64}; Vfrz::Vector{Float64}
    Rc::Vector{Float64}; Xc::Vector{Float64}; Kc::Vector{Float64}
    emax::Vector{Float64}; emin::Vector{Float64}; dbd1::Vector{Float64}; dbd2::Vector{Float64}
    Qmax::Vector{Float64}; Qmin::Vector{Float64}; Kpg::Vector{Float64}; Kig::Vector{Float64}
    Tp::Vector{Float64}; fdbd1::Vector{Float64}; fdbd2::Vector{Float64}
    femax::Vector{Float64}; femin::Vector{Float64}; Pmax::Vector{Float64}; Pmin::Vector{Float64}
    Tg::Vector{Float64}; Ddn::Vector{Float64}; Dup::Vector{Float64}
end

"PLL and voltage-magnitude filter of every bus."
struct MeasurementParams
    pll_ke::Vector{Float64}   # PLL gain
    pll_te::Vector{Float64}   # PLL integrator time constant
    vm_te::Vector{Float64}    # voltage-magnitude filter time constant
end

"""
    DynamicData

All dynamic-model parameters of the case.  `gov_type[i]` is the governor
model of generator `i` and `gov_instance[i]` its position in that model's
parameter vectors; `pss_instance[i]` is 0 if generator `i` has no stabilizer.
"""
struct DynamicData
    genrou::GenrouParams
    sexs::SexsParams
    tgov1::Tgov1Params
    hygov::HygovParams
    gast::GastParams
    ieeest::IeeestParams
    regca::RegcaParams
    reecb::ReecbParams
    repca::RepcaParams
    meas::MeasurementParams
    gov_type::Vector{GovernorType}
    gov_instance::Vector{Int}
    pss_instance::Vector{Int}
end

function DynamicData(d::AbstractDict, pfd::PowerFlowData)
    n = ngen(pfd)
    p(k) = jvec(d, k)
    pi(k) = jivec(d, k)
    ps(k) = jsvec(d, k)

    genrou = GenrouParams(pi("gen_genrou_bus"), ps("gen_genrou_id"),
        p("gen_genrou_Td0p"), p("gen_genrou_Td0pp"), p("gen_genrou_Tq0p"), p("gen_genrou_Tq0pp"),
        p("gen_H"), p("gen_D"), p("gen_genrou_Xd"), p("gen_genrou_Xq"), p("gen_genrou_Xdp"),
        p("gen_genrou_Xqp"), p("gen_genrou_Xdpp"), p("gen_genrou_Xl"), p("gen_Ra"), p("gen_X0"))
    length(genrou.bus) == n || error("expected one GENROU model per generator ($n), got $(length(genrou.bus))")
    all(genrou.bus .== pfd.gen_bus) || @warn "GENROU order differs from the power-flow generator order"

    TA_o_TB = p("exc_sexs_TA_o_TB"); TB = p("exc_sexs_TB")
    sexs = SexsParams(pi("exc_sexs_bus"), ps("exc_sexs_id"), TB .* TA_o_TB, TB, p("exc_sexs_K"),
        p("exc_sexs_TE"), p("exc_sexs_Emin"), p("exc_sexs_Emax"))
    length(sexs.bus) == n || error("expected one SEXS exciter per generator")

    # governors and stabilizers are listed by (bus, id)
    gmap(busk, idk) = [generator_index(pfd, b, i) for (b, i) in zip(jivec(d, busk), jsvec(d, idk))]
    tgov1 = Tgov1Params(gmap("gov_tgov1_bus", "gov_tgov1_id"), p("gov_tgov1_R"), p("gov_tgov1_T1"),
        p("gov_tgov1_Vmax"), p("gov_tgov1_Vmin"), p("gov_tgov1_T2"), p("gov_tgov1_T3"), p("gov_tgov1_Dt"))
    hygov = HygovParams(gmap("gov_hygov_bus", "gov_hygov_id"), p("gov_hygov_R"), p("gov_hygov_r"),
        p("gov_hygov_Tr"), p("gov_hygov_Tf"), p("gov_hygov_Tg"), p("gov_hygov_VELM"), p("gov_hygov_GMAX"),
        p("gov_hygov_GMIN"), p("gov_hygov_TW"), p("gov_hygov_At"), p("gov_hygov_Dturb"), p("gov_hygov_qNL"))
    gast = GastParams(gmap("gov_gast_bus", "gov_gast_id"), p("gov_gast_R"), p("gov_gast_T1"), p("gov_gast_T2"),
        p("gov_gast_T3"), p("gov_gast_LdLmt"), p("gov_gast_KT"), p("gov_gast_VMAX"), p("gov_gast_VMIN"),
        p("gov_gast_Dturb"))
    ieeest = IeeestParams(gmap("pss_ieeest_bus", "pss_ieeest_id"), p("pss_ieeest_A1"), p("pss_ieeest_A2"),
        p("pss_ieeest_A3"), p("pss_ieeest_A4"), p("pss_ieeest_A5"), p("pss_ieeest_A6"), p("pss_ieeest_T1"),
        p("pss_ieeest_T2"), p("pss_ieeest_T3"), p("pss_ieeest_T4"), p("pss_ieeest_T5"), p("pss_ieeest_T6"),
        p("pss_ieeest_KS"), p("pss_ieeest_LSMAX"), p("pss_ieeest_LSMIN"), p("pss_ieeest_VCU"))

    # every generator has exactly one governor; a stabilizer is optional
    gov_type = fill(TGOV1, n)
    gov_instance = zeros(Int, n)
    assigned = falses(n)
    for (T, prm) in ((TGOV1, tgov1), (HYGOV, hygov), (GAST, gast)), (k, g) in enumerate(prm.gen)
        assigned[g] && error("generator $g has more than one governor")
        gov_type[g] = T; gov_instance[g] = k; assigned[g] = true
    end
    all(assigned) || error("every generator must have exactly one governor (TGOV1, HYGOV or GAST)")
    pss_instance = zeros(Int, n)
    for (k, g) in enumerate(ieeest.gen)
        pss_instance[g] = k
    end

    regca = RegcaParams(pi("ibr_regca_bus"), ps("ibr_regca_id"), p("ibr_regca_Tg"), p("ibr_regca_Rrpwr"),
        p("ibr_regca_Brkpt"), p("ibr_regca_Zerox"), p("ibr_regca_Volim"), p("ibr_regca_Lvpnt1"),
        p("ibr_regca_Lvpnt0"), p("ibr_regca_Tfltr"), p("ibr_regca_Khv"), p("ibr_regca_Iqrmax"),
        p("ibr_regca_Iqrmin"), p("ibr_MVAbase"))
    length(regca.bus) == nibr(pfd) || error("expected one REGCA model per IBR")
    reecb = ReecbParams(pi("ibr_reecb_PFFLAG"), pi("ibr_reecb_VFLAG"), pi("ibr_reecb_QFLAG"), pi("ibr_reecb_PQFLAG"),
        p("ibr_reecb_Vdip"), p("ibr_reecb_Vup"), p("ibr_reecb_Trv"), p("ibr_reecb_dbd1"), p("ibr_reecb_dbd2"),
        p("ibr_reecb_Kqv"), p("ibr_reecb_Iqhl"), p("ibr_reecb_Iqll"), p("ibr_reecb_Vref0"), p("ibr_reecb_Tp"),
        p("ibr_reecb_Qmax"), p("ibr_reecb_Qmin"), p("ibr_reecb_Vmax"), p("ibr_reecb_Vmin"), p("ibr_reecb_Kqp"),
        p("ibr_reecb_Kqi"), p("ibr_reecb_Kvp"), p("ibr_reecb_Kvi"), p("ibr_reecb_Tiq"), p("ibr_reecb_dPmax"),
        p("ibr_reecb_dPmin"), p("ibr_reecb_Pmax"), p("ibr_reecb_Pmin"), p("ibr_reecb_Imax"), p("ibr_reecb_Tpord"))
    repca = RepcaParams(pi("ibr_repca_remote_bus"), pi("ibr_repca_branch_From_bus"), pi("ibr_repca_branch_To_bus"),
        pi("ibr_repca_VCFlag"), pi("ibr_repca_RefFlag"), pi("ibr_repca_FFlag"), p("ibr_repca_Tfltr"),
        p("ibr_repca_Kp"), p("ibr_repca_Ki"), p("ibr_repca_Tft"), p("ibr_repca_Tfv"), p("ibr_repca_Vfrz"),
        p("ibr_repca_Rc"), p("ibr_repca_Xc"), p("ibr_repca_Kc"), p("ibr_repca_emax"), p("ibr_repca_emin"),
        p("ibr_repca_dbd1"), p("ibr_repca_dbd2"), p("ibr_repca_Qmax"), p("ibr_repca_Qmin"), p("ibr_repca_Kpg"),
        p("ibr_repca_Kig"), p("ibr_repca_Tp"), p("ibr_repca_fdbd1"), p("ibr_repca_fdbd2"), p("ibr_repca_femax"),
        p("ibr_repca_femin"), p("ibr_repca_Pmax"), p("ibr_repca_Pmin"), p("ibr_repca_Tg"), p("ibr_repca_Ddn"),
        p("ibr_repca_Dup"))

    meas = MeasurementParams(p("pll_ke"), p("pll_te"), p("vm_te"))
    length(meas.pll_ke) == nbus(pfd) || error("expected one PLL entry per bus")

    return DynamicData(genrou, sexs, tgov1, hygov, gast, ieeest, regca, reecb, repca, meas,
                       gov_type, gov_instance, pss_instance)
end

# --------------------------------------------------------------------------
# The case
# --------------------------------------------------------------------------

"Everything that describes the system: `pfd` (power flow) and `dyd` (dynamic parameters)."
struct CaseData
    name::String
    pfd::PowerFlowData
    dyd::DynamicData
end

"""
    load_case(path) -> CaseData

Reads a case JSON file, as written by `load_psse(raw, dyr; json=path)`.
"""
function load_case(path::AbstractString)
    raw = JSON.parsefile(path)
    pfd = PowerFlowData(raw["pfd"])
    dyd = DynamicData(raw["dyd"], pfd)
    name = haskey(raw, "meta") && haskey(raw["meta"], "name") ? String(raw["meta"]["name"]) : basename(path)
    return CaseData(name, pfd, dyd)
end

function Base.show(io::IO, c::CaseData)
    p = c.pfd
    print(io, "CaseData(\"", c.name, "\": ", nbus(p), " buses, ", ngen(p), " generators, ", nibr(p),
          " IBRs, ", length(p.line_from), " lines, ", length(p.xfmr_from), " transformers, ",
          length(p.load_bus), " loads)")
end
