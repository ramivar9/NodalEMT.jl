# ==========================================================================
#  results.jl -- recorded results and CSV output
#
#  Results are stored every `record_every` steps in preallocated matrices,
#  one column per recorded instant.
# ==========================================================================

using Printf

"Recorded trajectories; column j of every matrix belongs to time `t[j]`."
mutable struct Results
    record_every::Int
    t::Vector{Float64}
    V::Matrix{Float64}          # node voltages (3 nbus x nt) [pu, instantaneous]
    vtm::Matrix{Float64}        # measured voltage magnitude of every bus [pu]
    freq::Matrix{Float64}       # PLL frequency of every bus [pu]
    omega::Matrix{Float64}      # generator speed [rad/s]
    delta::Matrix{Float64}      # generator rotor angle [rad]
    pe::Matrix{Float64}         # generator electrical power [pu, machine base]
    pm::Matrix{Float64}         # generator mechanical power
    EFD::Matrix{Float64}        # generator field voltage
    ibr_P::Matrix{Float64}      # IBR active power [pu, IBR base]
    ibr_Q::Matrix{Float64}      # IBR reactive power
    n::Int                      # number of columns filled so far
end

function Results(N::Int, nbus::Int, ngen::Int, nibr::Int, nt::Int, record_every::Int)
    return Results(record_every, zeros(nt), zeros(N, nt), zeros(nbus, nt), zeros(nbus, nt),
                   zeros(ngen, nt), zeros(ngen, nt), zeros(ngen, nt), zeros(ngen, nt), zeros(ngen, nt),
                   zeros(nibr, nt), zeros(nibr, nt), 0)
end

"Trims the preallocated matrices to the recorded instants."
function finalize!(r::Results)
    n = r.n
    r.t = r.t[1:n]
    for f in (:V, :vtm, :freq, :omega, :delta, :pe, :pm, :EFD, :ibr_P, :ibr_Q)
        setfield!(r, f, getfield(r, f)[:, 1:n])
    end
    return r
end

"""
    save_results(r, sim, dir)

Writes the results as CSV files in `dir`: `voltages.csv` (all node
voltages), `buses.csv` (voltage magnitude and frequency), `generators.csv`
and, if the system has IBRs, `ibr.csv`.
"""
function save_results(r::Results, sim, dir::AbstractString)
    mkpath(dir)
    pfd = sim.case.pfd
    nb = nbus(pfd)
    write_csv(joinpath(dir, "voltages.csv"), r.t,
              ["bus$(pfd.bus_num[k])_v$(ph)" for ph in ("a", "b", "c") for k in 1:nb], r.V)
    write_csv(joinpath(dir, "buses.csv"), r.t,
              vcat(["bus$(pfd.bus_num[k])_vtm" for k in 1:nb], ["bus$(pfd.bus_num[k])_freq" for k in 1:nb]),
              vcat(r.vtm, r.freq))
    gname = ["gen$(pfd.gen_bus[i])_$(strip(pfd.gen_id[i]))" for i in 1:ngen(pfd)]
    write_csv(joinpath(dir, "generators.csv"), r.t,
              vcat(gname .* "_omega", gname .* "_delta", gname .* "_pe", gname .* "_pm", gname .* "_EFD"),
              vcat(r.omega, r.delta, r.pe, r.pm, r.EFD))
    if nibr(pfd) > 0
        iname = ["ibr$(pfd.ibr_bus[i])_$(strip(pfd.ibr_id[i]))" for i in 1:nibr(pfd)]
        write_csv(joinpath(dir, "ibr.csv"), r.t, vcat(iname .* "_P", iname .* "_Q"), vcat(r.ibr_P, r.ibr_Q))
    end
    return dir
end

function write_csv(path::AbstractString, t::Vector{Float64}, names::Vector{String}, data::Matrix{Float64})
    open(path, "w") do io
        println(io, "time,", join(names, ","))
        for j in eachindex(t)
            @printf(io, "%.7e", t[j])
            for i in axes(data, 1)
                @printf(io, ",%.10e", data[i, j])
            end
            println(io)
        end
    end
    return path
end
