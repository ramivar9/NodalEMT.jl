# NodalEMT.jl

[![CI](https://github.com/ramivar9/NodalEMT.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ramivar9/NodalEMT.jl/actions/workflows/CI.yml)

NodalEMT.jl is an open-source Julia simulator of electromagnetic transients (EMT) for batches of
fault, outage and control scenarios of one power system on a multicore CPU. It implements the
network, machine, exciter, governor, stabilizer and inverter-based resource (IBR) models of NREL's
[ParaEMT](https://github.com/NatLabRockies/ParaEMT_public) in compiled Julia, and adds:

* **Faults and trips by compensation.** Bus faults and generator trips are applied with the
  Sherman–Morrison–Woodbury identity, so the network matrix `G` and its LU factorization never
  change during a simulation.
* **Scenario-parallel batches.** A batch of scenarios is split among worker threads. Each worker
  factorizes `G` once and runs complete scenarios with it; all other constants are stored once
  for the batch, and each scenario owns its states and compensation matrices.
* **π+Bergeron lines.** Every line whose travel time exceeds the time step is a lossy Bergeron
  model (EMTP's R/4, R/2, R/4 resistance split) with fractional-delay interpolation in packed
  ring buffers, which decouples the buses at its two ends in `G`; shorter lines stay π sections.
* **Exact periodic initialization** of the discrete Bergeron line from the power flow, so long
  lines start without a transient.
* **PSS/E input, end to end** (added after the paper). `load_psse` reads a `.raw` and a `.dyr`
  file, solves the AC power flow with [PowerFlows.jl](https://github.com/NREL-Sienna/PowerFlows.jl)
  and initializes the EMT simulation from it, so no external power-flow program is needed.

A time step allocates no memory. With the ParaEMT models, π lines and ParaEMT's own input data,
NodalEMT.jl reproduces ParaEMT to round-off (paper, Section V-B).

NodalEMT.jl is the software of the paper

> R. T. Alla and A. R. Ramapuram Matavalam, "NodalEMT.jl: An Open-Source Simulator for Fast
> Multi-Scenario EMT Studies on Multicore CPUs," submitted to the 2027 IEEE PES General Meeting.

The repository is `NodalEMT.jl`; the Julia package and module are `NodalEMT`.

## Installation

Julia 1.10 or later:

```julia
using Pkg
Pkg.add(url="https://github.com/ramivar9/NodalEMT.jl")
```

This also installs PowerSystems.jl and PowerFlows.jl, which read the PSS/E files and solve the
power flow. Installation precompiles the package and runs `load_psse` on both bundled systems,
which takes a few minutes once. After that, `using NodalEMT` takes about 4 s and the first
`load_psse` of a session about 3 s.

Or, from a clone of this repository:

```
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

## Quick start

One scenario: a three-phase fault at bus 8 of the Kundur two-area system, read from its PSS/E
files with the power flow solved by PowerFlows.jl.

```julia
using NodalEMT

psse = joinpath(pkgdir(NodalEMT), "data", "psse")
case = load_psse(joinpath(psse, "twoarea.raw"), joinpath(psse, "twoarea.dyr"); vm_te=0.1)
sim  = build_simulation(case)                        # 50 us step, Bergeron lines, updated models
schedule_bus_fault!(sim; bus=8, t_on=1.0, t_off=1.1, Zflt=1e-3)
schedule_generator_trip!(sim; t=3.0, gen=1)
res  = run!(sim, 5.0; record_every=10)               # results every 0.5 ms
save_results(res, sim, "results_2area")              # CSV files
```

(`vm_te = 0.1` s is ParaEMT's bus voltage filter for this system; the default 0.02 s is that of
the WECC system.)

A batch of scenarios on all threads (start Julia with `--threads=N`):

```julia
case  = load_psse(joinpath(psse, "WECC240.raw"), joinpath(psse, "WECC240.dyr"))
batch = build_batch(case, 100)                        # 100 scenarios, one factorization per worker
for s in 1:100
    schedule_bus_fault!(batch[s]; bus=case.pfd.bus_num[2s], t_on=0.2, t_off=0.26, Zflt=0.01)
end
results = run!(batch, 1.0; record_every=20)           # one Results per scenario
```

## Options

`build_simulation(case; kwargs...)` and `build_batch(case, S; kwargs...)` accept:

| Keyword | Default | Meaning |
|---|---|---|
| `ts` | `50e-6` | time step [s] |
| `variant` | `:updated` | model set: `:updated` (the paper's) or `:paraemt` (the ParaEMT models) |
| `bergeron_threshold` | `1.0` | lines with a travel time above this many steps are Bergeron lines; `Inf` keeps every line a π section |
| `network_events` | `:compensate` | `:compensate` keeps `G`; `:refactorize` adds each event to `G` and factorizes again (single simulations only) |

`build_batch` also accepts `nworkers` (default `Threads.nthreads()`) and `factorization`:
`:per_worker` (default), `:per_scenario` or `:shared`.

Disturbances: `schedule_bus_fault!(sim; bus, t_on, t_off, Zflt, phases)`,
`schedule_generator_trip!(sim; t, gen)` and
`schedule_reference_step!(sim; t, gen, kind=:governor or :exciter, dsp)`. Buses are given by
their number in the case, generators by their position in `case.pfd.gen_bus`.

### Model sets

The updated model set (`variant = :updated`) differs from the ParaEMT models in five points:
the numerical damping resistors are ten times lighter for lines and removed for loads and shunts;
the SEXS exciter includes the derivative of the stabilizer output in its lead term; and REECB,
REPCA and motoring machines are initialized in equilibrium. `variant = :paraemt` reproduces the
ParaEMT models, including ParaEMT's companion model of series capacitors.

## PSS/E input

`load_psse(raw, dyr; kwargs...)` returns the case that `build_simulation` and `build_batch` take:

1. PowerSystems.jl reads the `.raw` file (versions 33, 34 and 35; NodalEMT.jl rewrites a v34
   file in the v35 layout that PowerSystems.jl reads).
2. PowerFlows.jl solves the AC power flow.
3. NodalEMT.jl reads the `.dyr` file and builds the network of the EMT model from the solved
   system. It then checks that this network reproduces the power-flow balance.

| Keyword | Default | Meaning |
|---|---|---|
| `solve_power_flow` | `true` | `false` keeps the solution stored in the `.raw` file |
| `reactive_power_limits` | `false` | PV buses at a reactive limit become PQ buses |
| `solver` | `:newton_raphson` | or `:trust_region`, `:levenberg_marquardt` |
| `tol`, `maxiter` | `1e-9`, `50` | power-flow tolerance [pu] and iteration limit |
| `pll_ke`, `pll_te`, `vm_te` | `0.1`, `0.05`, `0.02` | bus PLL and voltage filter, not part of PSS/E data |
| `json` | `nothing` | a path: also save the case as JSON, readable by `load_case` without solving the power flow again |

Supported data:

* **Network:** lines with their charging and line-end shunts; two- and three-winding
  transformers with off-nominal ratios and magnetizing branches; fixed and switched shunts
  (capacitors, reactors and conductances); loads with any ZIP mix.
* **Dynamic models:** each unit is either a synchronous machine with GENROU, SEXS, one
  governor (TGOV1, HYGOV or GAST) and optionally IEEEST, or an IBR with REGCA1, REECB1 and
  REPCA1. Synchronous condensers (PG = PT = PB = 0) are machines too. As in ParaEMT, every
  machine needs a governor model, including condensers.
* **Rejected:** any other model attached to an in-service unit, phase-shifting transformers,
  HVDC lines and FACTS devices, each with a message naming them.

Where the conventions of PowerSystems.jl and PSS/E differ, `load_psse` follows PSS/E:

* A PV bus without an in-service unit is solved as PQ.
* Transformer vector groups (VECGRP) add no phase shift.
* Every unit keeps its scheduled P and Q. Only the change that the power flow makes at its
  bus (P at the slack bus, Q at PV buses) is shared among the bus's units by MVA base. For a
  `.raw` file that PSS/E has solved, this gives PSS/E's unit outputs.

Reading the WECC system takes about 3 s the first time in a session and 0.2 s after that. With
`json = path` the case is saved once; `load_case(path)` then reads it in about a second, without
solving the power flow.

PowerFlows.jl stores the admittance matrix in single precision. Its solution therefore balances
the exact network only to about 1e-7 of the largest branch admittance: 1e-5 pu for the
two-area system and 7e-4 pu for the WECC system. The five-decimal solution stored in the WECC
`.raw` file misses by 5e-2 pu.

PowerSystems.jl drops line-end shunt conductances (GI, GJ), so neither the power flow nor the
EMT model contains them.

## Examples

[`examples/README.md`](examples/README.md) lists the examples and gives short recipes for
reading cases, scheduling disturbances, choosing options, running batches and reading results.

```
julia --project=. examples/basic_simulation.jl
julia --project=. --threads=8 examples/contingency_screening.jl
julia --project=. examples/own_system.jl my.raw my.dyr
julia --project=. examples/two_area_fault.jl
julia --project=. --threads=16 examples/wecc_batch.jl 500
```

* `basic_simulation.jl` -- one simulation of the two-area system, from the PSS/E files to the
  results.
* `contingency_screening.jl` -- a fault at every 500-kV bus of the WECC system as a batch on all
  threads, ranked by post-fault voltage.
* `own_system.jl` -- a template for your own `.raw`/`.dyr` files.
* `two_area_fault.jl` -- the fault of the paper's PSCAD/EMTDC comparison (bus 8, 1 Ω, 100 ms),
  with π+Bergeron lines and with a single π section per line. With Bergeron lines the
  peak fault current is 8.19 kA (PSCAD: 8.06 kA); with π lines it is 12.69 kA.
* `wecc_batch.jl` -- the WECC contingency batch of the paper, timed on one worker and on all
  threads with each factorization policy. With `:shared`, Julia prints "waiting for UmfpackLU's
  lock": the substitutions of one shared factorization are serialized. In the paper, 500
  scenarios on 16 workers ran 12.2 times faster than on one worker with a factorization per
  worker, and 7.6 times faster with a shared factorization. At low thread counts the policies
  perform alike.

## Verification results

[`docs/verification.md`](docs/verification.md) shows, with figures:

* the power flow against PSS/E;
* the exact initialization of the Bergeron lines;
* the fault of the PSCAD/EMTDC comparison;
* compensation against refactorization;
* the undisturbed start from the power flow;
* a typical run;
* the speedup of scenario batches.

`docs/make_figures.jl` regenerates the figures and numbers.

## Tests

```
julia --project=. --threads=4 -e "using Pkg; Pkg.test()"
```

The tests run on the two bundled systems, read from their PSS/E files. They check:

* **Simulator:**
  * A Bergeron line starts in the exact steady state of the discrete model.
  * Compensation equals refactorization.
  * Every scenario of a batch equals an independent simulation bit for bit, for each
    factorization policy and worker count.
  * Each change of the updated model set acts where it should.
  * A time step allocates no memory.
* **PSS/E input:**
  * The `.dyr` parameter order agrees with the independent PSS/E mapping of PowerSystems.jl.
  * For both systems, the power flow of PowerFlows.jl reproduces the solution PSS/E stored in
    the `.raw` file to 2e-5 pu, and an undisturbed EMT simulation stays at the operating point.
  * A case with off-nominal taps, magnetizing branches, reactors, a three-winding transformer
    and a synchronous condenser gives an EMT network with the same admittance matrix as the
    power flow, and an undisturbed simulation stays at the operating point.
  * A case saved as JSON reads back identically.

## Source files

| File | Content |
|---|---|
| `src/case_data.jl` | power flow and dynamic parameters of a case; reading a saved JSON case |
| `src/bergeron.jl` | lossy Bergeron model: ring buffers, fractional delay, exact initialization |
| `src/network.jl` | companion models, the network matrix `G`, history currents |
| `src/compensation.jl` | faults and generator trips without changing `G` |
| `src/machine.jl` | GENROU synchronous machine as a Norton source |
| `src/controls.jl` | swing equation, SEXS, TGOV1, HYGOV, GAST, IEEEST |
| `src/bus_measurement.jl` | PLL and voltage filter of every bus |
| `src/ibr.jl` | IBRs: REGCA, REECB, REPCA |
| `src/results.jl` | recording and CSV output |
| `src/simulation.jl` | initialization, disturbances, one time step, time loop |
| `src/batch.jl` | scenario batches on worker threads |
| `src/psse.jl` | the `.dyr` reader, RAW v34 conversion, case assembly |
| `src/power_flow.jl` | `load_psse`: the `.raw` file and its power flow through PowerSystems.jl and PowerFlows.jl |
| `src/precompile.jl` | compiles reading, power flow and simulation when the package is installed |

The header of each file explains its model and refers to the equations of the paper.

## Test systems

`data/psse/` holds the two systems of the paper as PSS/E files:

* `twoarea.raw` and `twoarea.dyr` -- the Kundur two-area system (11 buses, 4 generators).
* `WECC240.raw` and `WECC240.dyr` -- the reduced WECC 240-bus system (243 buses, 111
  generators, 91 IBRs).

The `.raw` files are those of the public release of ParaEMT (commit d79d735): `twoarea.raw` is
RAW v33, `WECC240.raw` RAW v34. The `.dyr` files contain ParaEMT's dynamic data in PSS/E
format.

The paper's results were computed from ParaEMT's exported data, whose power-flow solution and
impedances PSS/E had rounded. Read from these files with the power flow re-solved, the systems
give the same fault simulations within 2e-4 of the peak voltage.

## Citation

If you use NodalEMT.jl, please cite the paper above (see `CITATION.cff`).

## License and acknowledgment

NodalEMT.jl is released under the BSD 3-Clause License (`LICENSE`). It reimplements the models of
ParaEMT, developed at the National Laboratory of the Rockies (formerly NREL) and distributed under
the BSD 3-Clause License; its notice and the provenance of the bundled data are in
`THIRD_PARTY_NOTICES.md`. This project is not endorsed by the National Laboratory of the Rockies.
