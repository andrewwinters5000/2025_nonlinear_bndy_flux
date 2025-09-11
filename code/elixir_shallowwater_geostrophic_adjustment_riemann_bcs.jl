
using OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# semidiscretization of the shallow water equations with Coriolis source term

equations = ShallowWaterEquations2D(gravity = 1.0)

# Initial condition for a geostrophic adjustment test with elliptical mass imbalance.
function initial_condition_geostrophic_adjustment(x, t, equations::ShallowWaterEquations2D)
    # Parameters of the initial perturbation
    a_0 = 0.5    # amplitude
    r_e = 0.1    # edge width
    r_i = 1.0    # radius
    lambda = 2.5 # aspect ratio

    h = 1 + a_0 / 2 * (1 - tanh((sqrt((sqrt(lambda) * x[1])^2 + (x[2] / sqrt(lambda))^2) - r_i) / r_e))
    h_v1 = 0.0
    h_v2 = 0.0
    b = 0.0

    return SVector(h, h_v1, h_v2, b)
end

initial_condition = initial_condition_geostrophic_adjustment

"""
    source_terms_coriolis(u, x, t, equations::ShallowWaterEquationsWetDry2D)

This source term can used to account for the Coriolis forces in the [`ShallowWaterEquationsWetDry2D`](@ref).
The Coriolis parameter `f` is computed with a beta-plane approximation (f ≈ f0 + beta * x[2]), where
`f` is approximated by the first two terms of its Taylor expansion around x[2] = 0.
"""
@inline function source_terms_coriolis(u, x, t, equations::ShallowWaterEquations2D)
    # Compute the Coriolis parameter with a beta-plane approximation
    beta = 0.0
    f0 = 1.0

    f = f0 + beta * x[2]
    return SVector(0, f * u[3], -f * u[2], 0)
end

###############################################################################
# extrapolating Riemann invariant strategy for the subcritical outflow boundary
@inline function boundary_condition_outflow_riemann_invariant(u_inner, normal_direction,
                                                              direction, x, t,
                                                              surface_flux_functions,
                                                              equations::ShallowWaterEquations2D)
    _, nonconservative_flux_function = surface_flux_functions

    h_inner = Trixi.waterheight(u_inner, equations)
    v1_inner, v2_inner = velocity(u_inner, equations)
    h_outer = 1.0

    # get the appropriate normal components from the `direction`
    if direction == 1 # x_negaitve
        v1_outer = v1_inner + 2 * sqrt(equations.gravity) * (sqrt(h_outer) - sqrt(h_inner))
        u_boundary = SVector(h_outer, h_outer*v1_outer, u_inner[3], u_inner[4])
    elseif direction == 2 # x_positive
        v1_outer = v1_inner - 2 * sqrt(equations.gravity) * (sqrt(h_outer) - sqrt(h_inner))
        u_boundary = SVector(h_outer, h_outer*v1_outer, u_inner[3], u_inner[4])
    elseif direction == 3 # y_negative
        v2_outer = v2_inner + 2 * sqrt(equations.gravity) * (sqrt(h_outer) - sqrt(h_inner))
        u_boundary = SVector(h_outer, u_inner[2], h_outer*v2_outer, u_inner[4])
    elseif direction == 4 # y_positive
        v2_outer = v2_inner - 2 * sqrt(equations.gravity) * (sqrt(h_outer) - sqrt(h_inner))
        u_boundary = SVector(h_outer, u_inner[2], h_outer*v2_outer, u_inner[4])
    end

    # Calculate boundary flux

    # Note, flip sign of normal to make it outward pointing, then flip the sign of the normal flux back
    # to be inward pointing on the -x and -y sides due to the orientation convention used by StructuredMesh
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux = flux_lax_friedrichs(u_inner, u_boundary, normal_direction, equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux = -flux_lax_friedrichs(u_inner, u_boundary, -normal_direction, equations)
    end

    # For this test case, the nonconservative flux will be a zero vector.
    noncons_flux = nonconservative_flux_function(u_inner, u_boundary, normal_direction,
                                                 equations)

    return flux, noncons_flux
end

# Use the Riemann invariant subsonic outflow boundary (must use StructuredMesh)
boundary_conditions = (x_neg = boundary_condition_outflow_riemann_invariant,
                       x_pos = boundary_condition_outflow_riemann_invariant,
                       y_neg = boundary_condition_outflow_riemann_invariant,
                       y_pos = boundary_condition_outflow_riemann_invariant)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (flux_lax_friedrichs, flux_nonconservative_wintermeyer_etal)
polydeg = 8
solver = DGSEM(polydeg = polydeg,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Create the StructuredMesh
coordinates_min = (-10.0, -10.0)
coordinates_max = (10.0, 10.0)

trees_per_dimension = (16, 16)
mesh = StructuredMesh(trees_per_dimension, coordinates_min, coordinates_max, periodicity = false)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_conditions,
                                    source_terms = source_terms_coriolis)

###############################################################################
# ODE solver

tspan = (0.0, 100.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_integrals = (lake_at_rest_error, energy_total),
                                     save_analysis = true)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 0.25,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 0.5)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            ode_default_options()..., callback = callbacks,
            adaptive = false, dt = 1.0);
