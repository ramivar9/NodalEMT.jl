# Shared helpers of the tests.

using NodalEMT
using PowerSystems, PowerFlows

const PSSEDIR = normpath(joinpath(@__DIR__, "..", "data", "psse"))

const LOADED = Dict{Symbol,CaseData}()

"""
    test_case(:two_area) / test_case(:wecc)

A fresh copy of a bundled system, read from its .raw and .dyr files (with the
power flow solved) once per test run.  `vm_te = 0.1` is ParaEMT's voltage
filter of the two-area system.
"""
function test_case(name::Symbol)
    c = get!(LOADED, name) do
        name === :two_area ? load_psse(joinpath(PSSEDIR, "twoarea.raw"), joinpath(PSSEDIR, "twoarea.dyr"); vm_te=0.1) :
        name === :wecc ? load_psse(joinpath(PSSEDIR, "WECC240.raw"), joinpath(PSSEDIR, "WECC240.dyr")) :
        error("unknown test case $name")
    end
    return deepcopy(c)
end

maxabs(a, b) = isempty(a) ? 0.0 : maximum(abs.(a .- b))

"All generator states as one vector (per machine the GENROU states, then the controllers)."
function generator_state_vector(sim)
    x = sim.x; m = x.mac; dyd = sim.case.dyd
    v = Float64[]
    for i in eachindex(m.delta)
        append!(v, (m.delta[i], m.omega[i], m.id[i], m.iq[i], m.ifd[i], m.i1d[i], m.i1q[i], m.i2q[i],
                    m.ed[i], m.eq[i], m.psid[i], m.psiq[i], m.psifd[i], m.psi1q[i], m.psi1d[i], m.psi2q[i],
                    m.te[i], m.qe[i]))
    end
    for i in eachindex(m.delta)
        append!(v, (x.exc.v1[i], x.exc.EFD[i]))
    end
    for k in eachindex(dyd.tgov1.gen)
        append!(v, (x.tgov1.p1[k], x.tgov1.p2[k], x.tgov1.pm[k]))
    end
    for k in eachindex(dyd.hygov.gen)
        append!(v, (x.hygov.xe[k], x.hygov.xc[k], x.hygov.xg[k], x.hygov.xq[k], x.hygov.pm[k]))
    end
    for k in eachindex(dyd.gast.gen)
        append!(v, (x.gast.p1[k], x.gast.p2[k], x.gast.p3[k], x.gast.pm[k]))
    end
    for k in eachindex(dyd.ieeest.gen)
        s = x.pss
        append!(v, (s.y1[k], s.y2[k], s.y3[k], s.y4[k], s.y5[k], s.y6[k], s.y7[k], s.x1[k], s.x2[k], s.vs[k]))
    end
    return v
end

"Bus measurement states in ParaEMT's order (`x_bus_pv_1`)."
function bus_state_vector(sim)
    b = sim.x_bus
    v = Float64[]
    for k in eachindex(b.ze)
        append!(v, (b.ze[k], b.de[k], b.we[k], b.vt[k], b.vtm[k], b.dvtm[k]))
    end
    return v
end

"""
Mask that excludes GENROU's `te` and `qe`: output-only quantities whose values
in ParaEMT (Numba) differ from a direct evaluation at the 1e-6 level and that
never feed back into the dynamics.
"""
function dynamics_mask(ngen::Int, len::Int)
    m = trues(len)
    for g in 1:ngen
        m[(g - 1) * 18 + 17] = false
        m[(g - 1) * 18 + 18] = false
    end
    return m
end

"Voltage magnitude of bus index `k` in column `j` of the results."
vmag(res, k, nb, j) = sqrt(2 / 3 * (res.V[k, j]^2 + res.V[k+nb, j]^2 + res.V[k+2nb, j]^2))
