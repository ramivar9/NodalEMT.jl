# Examples

Run each example from the repository root, in the package's own environment:

```
julia --project=. -e 'using Pkg; Pkg.instantiate()'                  # once, in a clone
julia --project=. examples/basic_simulation.jl
```

`Pkg.instantiate()` installs and precompiles the package once, which takes a few minutes. After
that, `using NodalEMT` and reading a system each take a few seconds. Results are written to
`examples/results_*/`, which git ignores.

| Example | What it shows | Command |
|---|---|---|
| [`basic_simulation.jl`](basic_simulation.jl) | One simulation of the two-area system: read the PSS/E files, apply a fault and a governor step, run, read and save the results. | `julia --project=. examples/basic_simulation.jl` |
| [`contingency_screening.jl`](contingency_screening.jl) | A batch on all threads: a fault at every 500-kV bus of the WECC system, with a unit trip where there is one, ranked by post-fault voltage. | `julia --project=. --threads=8 examples/contingency_screening.jl` |
| [`own_system.jl`](own_system.jl) | A template for your own `.raw`/`.dyr` files: power flow with a solver fallback, checks, saving the case as JSON, a steady-state check and a fault. | `julia --project=. examples/own_system.jl my.raw my.dyr` |
| [`two_area_fault.jl`](two_area_fault.jl) | The fault of the paper's PSCAD/EMTDC comparison, with Bergeron lines and with single π sections (peak current 8.19 kA vs. 12.69 kA; PSCAD: 8.06 kA). | `julia --project=. examples/two_area_fault.jl` |
| [`wecc_batch.jl`](wecc_batch.jl) | The contingency batch of the paper, timed on one worker and on all threads with each factorization policy. | `julia --project=. --threads=16 examples/wecc_batch.jl 500` |

`basic_simulation.jl` produces the run of Fig. 5 in [`docs/verification.md`](../docs/verification.md).

## Recipes

The snippets assume `using NodalEMT` and
`psse = joinpath(pkgdir(NodalEMT), "data", "psse")`.

**Read a case once, reuse it without solving the power flow again**

```julia
case = load_psse("my.raw", "my.dyr"; json="my_case.json")   # power flow solved, case saved
case = load_case("my_case.json")                            # later: about a second, no power flow
```

`load_psse` keywords: `solve_power_flow=false` keeps the solution stored in the `.raw` file,
`reactive_power_limits=true` turns PV buses at a reactive limit into PQ buses,
`solver=:trust_region` or `:levenberg_marquardt` if Newton–Raphson does not converge, and
`pll_ke`, `pll_te`, `vm_te` set the bus PLL and voltage filter, which PSS/E files do not contain.

**Disturbances** (times in s; buses by number; generators by their position in `case.pfd.gen_bus`)

```julia
sim = build_simulation(case)
schedule_bus_fault!(sim; bus=8, t_on=1.0, t_off=1.1, Zflt=0.01)               # three-phase, Zflt in pu
schedule_bus_fault!(sim; bus=7, t_on=2.0, t_off=2.05, Zflt=0.05, phases=[1])  # phase a to ground
schedule_generator_trip!(sim; t=3.0, gen=2)
schedule_reference_step!(sim; t=4.0, gen=1, kind=:governor, dsp=0.002)        # or kind=:exciter
```

Faults and trips are applied by compensation: `G` and its factorization never change.

**Options**

```julia
sim = build_simulation(case; ts=20e-6)                 # time step [s], default 50e-6
sim = build_simulation(case; bergeron_threshold=Inf)   # every line a π section
sim = build_simulation(case; variant=:paraemt)         # ParaEMT's models instead of the updated set
```

**Run, or step yourself**

```julia
res = run!(sim, 5.0; record_every=20)     # to t = 5 s, storing every 20 steps
while current_time(sim) < 6.0             # or advance step by step,
    V = step!(sim)                        # V: node voltages of this step
end                                       # (schedule more disturbances at any time)
```

**Results** -- column `j` of every matrix belongs to `res.t[j]`:

| Field | Content |
|---|---|
| `res.V` | instantaneous node voltages [pu]: phase a of all buses, then b, then c |
| `res.vtm`, `res.freq` | measured voltage magnitude [pu] and PLL frequency [pu] of every bus |
| `res.omega`, `res.delta` | generator speed [rad/s] and rotor angle [rad] |
| `res.pe`, `res.pm`, `res.EFD` | electrical and mechanical power [pu, machine base], field voltage |
| `res.ibr_P`, `res.ibr_Q` | IBR active and reactive power [pu, IBR base] |

`bus_index(case.pfd, 8)` gives the row of bus 8 in `vtm` and `freq`, and its phase-a row in `V`;
the phase-b and phase-c rows are that plus `nbus(case.pfd)` and twice that.
`save_results(res, sim, dir)` writes all of it as CSV.

**A batch of scenarios on all threads** (start Julia with `--threads=N`)

```julia
case  = load_psse(joinpath(psse, "WECC240.raw"), joinpath(psse, "WECC240.dyr"))
batch = build_batch(case, 100)                      # one LU factorization per worker (default)
# or: build_batch(case, 100; factorization=:per_scenario or :shared, nworkers=8)
for s in 1:100
    schedule_bus_fault!(batch[s]; bus=case.pfd.bus_num[s], t_on=0.1, t_off=0.16, Zflt=0.01)
end
results = run!(batch, 1.0; record_every=40)         # one Results per scenario
```
