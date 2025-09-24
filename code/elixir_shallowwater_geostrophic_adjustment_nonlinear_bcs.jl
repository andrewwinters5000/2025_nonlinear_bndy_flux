
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

# Novel boundary flux for subcritical outflow
@inline function boundary_condition_subcritical_outflow_novel(u_inner, normal_direction::AbstractVector,
                                                              x, t,
                                                              surface_flux_functions,
                                                              equations::ShallowWaterEquations2D)
    _, nonconservative_flux_function = surface_flux_functions
    # Use the `normal_vector` to match how derived
    s_hat = Trixi.norm(normal_direction)
    # Normalize the vector without using `normalize` since we need to multiply by the `s_hat` later
    normal = normal_direction / s_hat

    # Pull the gravity constant and compute some convenience constants
    g = equations.gravity
    s2g = 1 / (2 * g)
    s4g = 1 / (4 * g)

    # Scaling variable
    alfa = -1 + sqrt(3)

    # Get the boundary data
    u_boundary = initial_condition_geostrophic_adjustment(x, t, equations)
    h_ext, u_ext, v_ext, _ = cons2prim(u_boundary, equations)

    # Compute wave speed and normal velocity
    c_ext = sqrt(g * h_ext)
    u_n_ext = u_ext * normal[1] + v_ext * normal[2]

    # Extract the internal data
    h, u, v, _ = cons2prim(u_inner, equations)

    c = sqrt(g * h)
    u_n = u * normal[1] + v * normal[2]

    # Geometric mean of the water height and the eigenvalue
    h_geom = sqrt(h * h_ext)
    lambda1_geom = sqrt((c - u_n) * (c_ext - u_n_ext))

    # Compute numerical flux components
    f1 = (0.5 * alfa * h * u_n + (1 - alfa) * h * c + alfa * c * u_n^2 * s2g
          - lambda1_geom * alfa * c_ext * (alfa * c_ext - u_n_ext) * s2g)

    f2 = ((0.25 * alfa + 0.5) * h * u * u_n + 0.5 * (1 - alfa) * h * c * u + alfa * c * u * u_n^2 * s4g
          + 0.5 * (1 - alfa) * g * h^2 * normal[1] + 0.5 * h * u_n * ((1 + alfa) * c - u_n) * normal[1]
          - lambda1_geom * c_ext * (alfa * c_ext - u_n_ext) * (alfa * u - 2 * c * normal[1]) * s4g)

    f3 = ((0.25 * alfa + 0.5) * h * v * u_n + 0.5 * (1 - alfa) * h * c * v + alfa * c * v * u_n^2 * s4g
          + 0.5 * (1 - alfa) * g * h^2 * normal[2] + 0.5 * h * u_n * ((1 + alfa) * c - u_n) * normal[2]
          - lambda1_geom * c_ext * (alfa * c_ext - u_n_ext) * (alfa * v - 2 * c * normal[2]) * s4g)

    # For this test case, the nonconservative flux will be a zero vector.
    noncons_flux = nonconservative_flux_function(u_inner, u_boundary, normal_direction,
                                                 equations)

    # Return the novel version of the flux
    return SVector(f1, f2, f3, zero(eltype(u_inner))) * s_hat, noncons_flux
end

# Use the novel boundary flux everywhere
boundary_conditions = Dict(:x_neg => boundary_condition_subcritical_outflow_novel,
                           :x_pos => boundary_condition_subcritical_outflow_novel,
                           :y_neg => boundary_condition_subcritical_outflow_novel,
                           :y_pos => boundary_condition_subcritical_outflow_novel)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
polydeg = 8
solver = DGSEM(polydeg = polydeg,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Create the StructuredMesh
coordinates_min = (-10.0, -10.0)
coordinates_max = (10.0, 10.0)

trees_per_dimension = (32, 32)
mesh = P4estMesh(trees_per_dimension,
                 polydeg = 1, initial_refinement_level = 0,
                 coordinates_min = coordinates_min, coordinates_max = coordinates_max,
                 periodicity = false)

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

save_solution = SaveSolutionCallback(dt = 1.0,
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
