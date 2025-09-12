# Numerical experiments

This directory contains all source code required to reproduce the numerical experiments presented in the paper. It is developed for Julia v1.11.3.

To reproduce the numerical experiments, clone this repository and start Julia with the project set to the local directory:

```bash
git clone https://github.com/andrewwinters5000/2025_nonlinear_bndy_flux.git
cd 2025_nonlinear_bndy_flux
julia --project=.
```

The mesh used in the channel flow computations is provided in the file `mesh_swe_channel_wiggles.inp`.
One can reproduce this mesh using HOHQMesh.jl and the script `interactive_swe_channel.jl`.

All the two-dimensional simulation were visualized in ParaView where the HDF5 output files were converted
using `Trixi2Vtk.jl` with a command like
```julia
trixi2vtk("out/solution_*.h5", output_directory="out", nvisnodes=8)
```

## Burgers'

The provided elixir for Trixi.jl uses the stable version. It can be run with
```julia
include("code/elixir_burgers_mms.jl");
```

To obtain the weakly unstable behavior of the entropy conservative flux execute
```julia
trixi_include("code/elixir_burgers_mms.jl", boundary_conditions = (x_neg = boundary_condition_inflow_ec,
                                                                   x_pos = boundary_condition_outflow_ec));
```

To obtain the local Lax-Friedrichs result use
```julia
trixi_include("code/elixir_burgers_mms.jl", boundary_conditions = (x_neg = boundary_condition_inflow_llf,
                                                                   x_pos = boundary_condition_outflow_llf));
```

One can use `plot(sol)` to obtain the figures in the manuscript.

## Subcritical manufactured solution

For the subcritical manufactured solution the nonlinearly stable fluxes are the default.
To generate this result use
```julia
trixi_include("code/elixir_shallowwater_subcritical_mms.jl", tspan = (0.0, 4.1));
```

The local Lax-Friedrichs result that crashes around t ≈ 4.1 uses
```julia
trixi_include("code/elixir_shallowwater_subcritical_mms.jl",
              surface_flux = (flux_lax_friedrichs, flux_nonconservative_wintermeyer_etal),
              boundary_condition = Dict(:Right  => boundary_condition_slip_wall,
                                        :Left   => boundary_condition_slip_wall,
                                        :Bottom => boundary_condition_subcritical_inflow_llf,
                                        :Top    => boundary_condition_subcritical_outflow_llf),
                                        tspan = (0.0, 4.1));
```

The HLL result runs through with internal dissipation via the command
```julia
trixi_include("code/elixir_shallowwater_subcritical_mms.jl",
              surface_flux = (flux_hll, flux_nonconservative_wintermeyer_etal),
              boundary_condition = Dict(:Right  => boundary_condition_slip_wall,
                                        :Left   => boundary_condition_slip_wall,
                                        :Bottom => boundary_condition_subcritical_inflow_hll,
                                        :Top    => boundary_condition_subcritical_outflow_hll),
                                        tspan = (0.0, 4.1));
```

## Supercritical manufactured solution

For the subcritical manufactured solution the nonlinearly stable fluxes are the default.
To generate this result use
```julia
trixi_include("code/elixir_shallowwater_supercritical_mms.jl", tspan = (0.0, 6.0));
```

The local Lax-Friedrichs result comes from
```julia
trixi_include("code/elixir_shallowwater_supercritical_mms.jl",
              boundary_condition = Dict(:Right  => boundary_condition_slip_wall,
                                        :Left   => boundary_condition_slip_wall,
                                        :Bottom => boundary_condition_supercritical_inflow_llf,
                                        :Top    => boundary_condition_supercritical_outflow_llf),
                                        tspan = (0.0, 6.0));
```

The HLL result comes from
```julia
trixi_include("code/elixir_shallowwater_supercritical_mms.jl",
              boundary_condition = Dict(:Right  => boundary_condition_slip_wall,
                                        :Left   => boundary_condition_slip_wall,
                                        :Bottom => boundary_condition_supercritical_inflow_hll,
                                        :Top    => boundary_condition_supercritical_outflow_hll),
                                        tspan = (0.0, 6.0));
```

## Geostrophic adjustment

The result using the Riemann invariant boundary conditions from the linear analysis
```julia
trixi_include("code/elixir_shallowwater_geostrophic_adjustment_riemann_bcs.jl");
```

The result from the nonlinearly stable outflow boundary result
```julia
trixi_include("code/elixir_shallowwater_geostrophic_adjustment_nonlinear_bcs.jl");
```

To obtain the result that using local Lax-Friedrichs together
with Riemann invariant boundary state construction that crashes at t ≈ 53 use
```julia
trixi_include("code/elixir_shallowwater_geostrophic_adjustment_riemann_bcs.jl",
              surface_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal));
```
