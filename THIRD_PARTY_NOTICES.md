# Third-party notices

## ParaEMT

NodalEMT.jl reimplements in Julia the network, machine, exciter, governor, stabilizer and IBR models
of ParaEMT, developed at the National Laboratory of the Rockies (formerly the National Renewable
Energy Laboratory, NREL) and released at https://github.com/NatLabRockies/ParaEMT_public.

* `data/psse/twoarea.raw` and `data/psse/WECC240.raw` are the PSS/E power-flow files of the
  public ParaEMT release (commit d79d735, `models/2area_psse/` and `models/240bus_psse/`).
* `data/psse/twoarea.dyr` and `data/psse/WECC240.dyr` contain the dynamic data of ParaEMT's
  spreadsheets `twoarea.xls` and `240bus.xls` from the same release, written in the PSS/E
  .dyr format by the authors of NodalEMT.jl.

## PowerSystems.jl and PowerFlows.jl

NodalEMT.jl depends on PowerSystems.jl and PowerFlows.jl (NREL Sienna, BSD 3-Clause License),
which `load_psse` uses to read `.raw` files and solve the power flow. Julia's package manager
installs them from their own repositories; they are not distributed with NodalEMT.jl.

ParaEMT is distributed under the following license:

```
Copyright (c) Year, Alliance for Sustainable Energy, LLC

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```
