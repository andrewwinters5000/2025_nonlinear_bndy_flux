using OrdinaryDiffEqLowStorageRK
using Trixi
using Plots

###############################################################################
# Construct the equations and initial condition

equations = InviscidBurgersEquation1D()

# manufactured solution and its corresponding source term
function initial_condition_mms(x, t, equations::InviscidBurgersEquation1D)
    2.0 + sin(pi*(x[1] - t) - 0.7) |> SVector
end

function source_term_mms(u, x, t, equations::InviscidBurgersEquation1D)
    return SVector( pi*cos(pi*(x[1] - t) - 0.7) * (-1 + 2.0 + sin(pi*(x[1] - t) - 0.7)) )
end


function inflow(x, t, equations::InviscidBurgersEquation1D)
    return initial_condition_mms(coordinate_min, t, equations)
end

###############################################################################
# Various inflow boundary routines for a nonperiodic problem

function boundary_condition_inflow_ec(u_inner, orientation, normal_direction, x, t,
                                      surface_flux_function, equations::InviscidBurgersEquation1D)
    # Get the external solution data
    u_ext = inflow(x, t, equations)

    # entropy conservative flux (dissipation-free)
    return flux_ec(u_ext, u_inner, orientation, equations)
end

function boundary_condition_inflow_llf(u_inner, orientation, normal_direction, x, t,
                                       surface_flux_function, equations::InviscidBurgersEquation1D)
    # Get the external solution data
    u_ext = inflow(x, t, equations)

    # local Lax-Friedrichs flux
    return flux_lax_friedrichs(u_ext, u_inner, orientation, equations)
end

function boundary_condition_inflow_stable(u_inner, orientation, normal_direction, x, t,
                                          surface_flux_function, equations::InviscidBurgersEquation1D)
    # Get the external solution data
    u_ext = inflow(x, t, equations)

    # provably nonlinearly stable inflow flux
    return SVector((2 * u_ext[1] * sqrt(abs(u_ext[1]) * abs(u_inner[1])) - u_inner[1]^2 / 2) / 3)
end

###############################################################################
# Various outflow boundary routines for a nonperiodic problem

function boundary_condition_outflow_ec(u_inner, orientation, normal_direction, x, t,
                                       surface_flux_function, equations::InviscidBurgersEquation1D)
    # Get the external solution data
    u_ext = initial_condition_mms(coordinate_max, t, equations)

    # entropy conservative flux (dissipation-free) that is not a Riemann solver
    return flux_ec(u_inner, u_ext, orientation, equations)
end

function boundary_condition_outflow_llf(u_inner, orientation, normal_direction, x, t,
                                        surface_flux_function, equations::InviscidBurgersEquation1D)
    # Get the external solution data
    u_ext = initial_condition_mms(coordinate_max, t, equations)

    return flux_lax_friedrichs(u_inner, u_ext, orientation, equations)
end

function boundary_condition_outflow(u_inner, orientation, normal_direction, x, t,
                                    surface_flux_function, equations::InviscidBurgersEquation1D)
    # Calculate the boundary flux entirely from the internal solution state
    return Trixi.flux(u_inner, orientation, equations)
end

initial_condition = initial_condition_mms

# Use the provably stable inflow boundary flux
boundary_conditions = (x_neg = boundary_condition_inflow_stable,
                       x_pos = boundary_condition_outflow)

###############################################################################
# Build the DGSEM approximation space that uses the entropy conserving flux for Burgers'
volume_flux = flux_ec

solver = DGSEM(polydeg=7, surface_flux=flux_ec,
               volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Create the 1D mesh with a specified number of elements
coordinate_min = (-1.0,)
coordinate_max = (1.0,)
cells_per_dimension = (5,)
mesh = StructuredMesh(cells_per_dimension, coordinate_min, coordinate_max, periodicity=false)

###############################################################################
# Wrap everything together into a semidiscretization (basically the RHS operator)
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions=boundary_conditions,
                                    source_terms=source_term_mms)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 120.0) # exhibits weak instability when `flux_ec` is used at inflow
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100000
analysis_callback = AnalysisCallback(semi, interval=analysis_interval,
                                     extra_analysis_errors=(:l2_error_primitive,
                                                            :linf_error_primitive))

alive_callback = AliveCallback(analysis_interval=analysis_interval)

stepsize_callback = StepsizeCallback(cfl=0.75)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            ode_default_options()..., callback = callbacks,
            adaptive = false, dt = 1.0);
