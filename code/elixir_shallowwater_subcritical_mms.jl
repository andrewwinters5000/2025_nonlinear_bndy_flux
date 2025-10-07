using OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# semidiscretization of the shallow water equations

equations = ShallowWaterEquations2D(gravity = 9.81)

###############################################################################
# routines for the MMS

@inline function initial_condition_moving_blob(x, t, equations::ShallowWaterEquations2D)
    # Parameters for the manufactured solution

    # Corresponds to a Froude number ≈ 0.25
    h0 = 32.0

    A0 = 1.0    # amplitude of perturbation to water height
    sigma = 8   # sharpness of Gaussian blob

    # pi/4 Rotated blob center where the Cartesian center was [-2.0; 1.5]
    x0 = -2.0 * sqrt(2.0)
    y0 = -1.0 / sqrt(2.0)

    # Calculate primitive variables
    # Note, the scaling by 1/(2g) makes controlling the Froude number easier

    # Constant velocities (original) # for subcritical
    v1 = 1.0 / sqrt(2.0)
    v2 = 1.0 / sqrt(2.0)

    # For moving blob testcase
    H  = (h0 + A0 * exp(-sigma * ((x[1] - v1 * t - x0)^2 + (x[2] - v2 * t - y0)^2))) / (2.0 * equations.gravity)
    b  = zero(eltype(x))

    return prim2cons(SVector(H, v1, v2, b), equations)
end

initial_condition = initial_condition_moving_blob

function source_term_moving_blob(u, x, t, equations::ShallowWaterEquations2D)
    # Parameters for the manufactured solution
    # Note, must match those from `initial_condition_moving_blob`

    # Corresponds to a Froude number ≈ 0.25
    h0 = 32.0

    A0 = 1.0    # amplitude of perturbation to water height
    sigma = 8   # sharpness of Gaussian blob

    # pi/4 Rotated blob center where the Cartesian center was [-2.0; 1.5]
    x0 = -2.0 * sqrt(2.0)
    y0 = -1.0 / sqrt(2.0)

    # Constant velocities
    v1 = 1.0 / sqrt(2.0)
    v2 = 1.0 / sqrt(2.0)

    # Convenience variable for the common term across H and its derivatives
    funk = -sigma * ((x[1] - v1 * t - x0)^2 + (x[2] - v2 * t - y0)^2)

    # Water height and its derivatives
    # Note, scaling by 1/(2g) makes controlling the Froude number easier
    H   = (h0 + A0 * exp(funk)) / (2.0 * equations.gravity)
    H_x = (-2.0 * A0 * sigma * (x[1] - v1 * t - x0) * exp(funk)) / (2.0 * equations.gravity)
    H_y = (-2.0 * A0 * sigma * (x[2] - v2 * t - y0) * exp(funk)) / (2.0 * equations.gravity)

    # MMS source terms assuming that velocity field are constants
    s1 = zero(eltype(u))
    s2 = equations.gravity * H * H_x
    s3 = equations.gravity * H * H_y

    return SVector(s1, s2, s3, zero(eltype(u)))
end

###############################################################################
# "standard" boundary condition routines where the Riemann solver, given the
# external and internal solution states, should extract the appropriate information
# OBS! Not provably stable for the nonlinear problem.

# Subcritical outflow boundary using LLF
function boundary_condition_subcritical_outflow_llf(u_inner, normal_direction::AbstractVector,
                                                    x, t, surface_flux_functions,
                                                    equations::ShallowWaterEquations2D)
    _, nonconservative_flux_function = surface_flux_functions

    # Get the external solution state
    u_boundary = initial_condition_moving_blob(x, t, equations)

    # local Lax-Friedrichs solver to impose the outflow boundary state
    flux = flux_lax_friedrichs(u_inner, u_boundary, normal_direction, equations)

    noncons_flux = nonconservative_flux_function(u_inner, u_ext, normal_direction,
                                                 equations)

    return flux, noncons_flux
end

# Subcritical outflow boundary using HLL
function boundary_condition_subcritical_outflow_hll(u_inner, normal_direction::AbstractVector,
                                                    x, t, surface_flux_functions,
                                                    equations::ShallowWaterEquations2D)
    _, nonconservative_flux_function = surface_flux_functions

    # Get the external solution state
    u_boundary = initial_condition_moving_blob(x, t, equations)

    # Harten-Lax-van Leer (HLL) solver to impose the outflow boundary state
    flux = flux_hll(u_inner, u_boundary, normal_direction, equations)

    noncons_flux = nonconservative_flux_function(u_inner, u_ext, normal_direction,
                                                 equations)

    return flux, noncons_flux
end

# Subcritical inflow boundary using LLF
function boundary_condition_subcritical_inflow_llf(u_inner, normal_direction::AbstractVector,
                                                   x, t, surface_flux_functions,
                                                   equations::ShallowWaterEquations2D)
    _, nonconservative_flux_function = surface_flux_functions

    # Get the external solution state
    u_boundary = initial_condition_moving_blob(x, t, equations)

    # local Lax-Friedrichs solver to impose the inflow boundary state
    flux = flux_lax_friedrichs(u_inner, u_boundary, normal_direction, equations)

    noncons_flux = nonconservative_flux_function(u_inner, u_ext, normal_direction,
                                                 equations)

    return flux, noncons_flux
end

# Subcritical inflow boundary using HLL
function boundary_condition_subcritical_inflow_hll(u_inner, normal_direction::AbstractVector,
                                                   x, t, surface_flux_functions,
                                                   equations::ShallowWaterEquations2D)
    _, nonconservative_flux_function = surface_flux_functions

    # Get the external solution state
    u_boundary = initial_condition_moving_blob(x, t, equations)

    # Harten-Lax-van Leer (HLL) solver to impose the inflow boundary state
    flux = flux_hll(u_inner, u_boundary, normal_direction, equations)

    noncons_flux = nonconservative_flux_function(u_inner, u_ext, normal_direction,
                                                 equations)

    return flux, noncons_flux
end

###############################################################################
# Novel nonlinearly stable subcritical outflow and inflow boundary fluxes
#   - 2 BCs at inflow and 1 BC at outflow when abs(un) < c

# Nonlinearly stable subcritical outflow boundary condition
@inline function boundary_condition_subcritical_outflow_stable(u_inner, normal_direction::AbstractVector,
                                                               x, t, surface_flux_functions,
                                                               equations::ShallowWaterEquations2D)
    _, nonconservative_flux_function = surface_flux_functions

    # Normalize the vector without using `normalize` since we need to multiply by the `s_hat` later
    s_hat = Trixi.norm(normal_direction)
    normal = normal_direction / s_hat

    # Pull the gravity constant and compute some convenience constants
    g = equations.gravity
    s2g = 1 / (2 * g)
    s4g = 1 / (4 * g)

    # Scaling variable
    alfa = -1 + sqrt(3)

    # Get the boundary data
    u_boundary = initial_condition_moving_blob(x, t, equations)
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

    # Compute the nonconservative piece; zero for this test case because bottom is 0
    noncons_flux = nonconservative_flux_function(u_inner, u_boundary, normal_direction,
                                                 equations)

    # Return the novel version of the flux
    return SVector(f1, f2, f3, zero(eltype(u_inner))) * s_hat, noncons_flux
end

# Nonlinearly stable subcritical inflow boundary condition
@inline function boundary_condition_subcritical_inflow_stable(u_inner, normal_direction::AbstractVector,
                                                              x, t, surface_flux_functions,
                                                              equations::ShallowWaterEquations2D)
    _, nonconservative_flux_function = surface_flux_functions

    # Normalize the vector without using `normalize` since we need to multiply by the `s_hat` later
    s_hat = Trixi.norm(normal_direction)
    normal = normal_direction / s_hat

    # Pull the gravity constant and compute some convenience constants
    g = equations.gravity
    s2g = 1 / (2 * g)
    s4g = 1 / (4 * g)

    # Scaling variable
    alfa = -1 + sqrt(3)

    # Get the boundary data
    u_boundary = initial_condition_moving_blob(x, t, equations)
    h_ext, u_ext, v_ext, _ = cons2prim(u_boundary, equations)

    # Compute wave speed as well as normal and tangential velocities
    c_ext = sqrt(g * h_ext)
    u_n_ext = u_ext * normal[1] + v_ext * normal[2]
    u_t_ext = -u_ext * normal[2] + v_ext * normal[1]

    # Extract the internal data
    h, u, v, _ = cons2prim(u_inner, equations)

    c = sqrt(g * h)
    u_n = u * normal[1] + v * normal[2]

    # Geometric mean of the water height and the eigenvalues
    h_geom = sqrt(h * h_ext)
    lambda1_geom = sqrt((abs(u_n) + c) * (abs(u_n_ext) + c_ext))
    lambda2_geom = sqrt(abs(u_n) * abs(u_n_ext))

    # Compute numerical flux components
    f1 = (0.5 * alfa * h * u_n + (1 - alfa) * h * c + alfa * c * u_n^2 * s2g
          - lambda1_geom * alfa * c_ext * (alfa * c_ext - u_n_ext) * s2g)

    f2 = ((0.25 * alfa - 0.5) * h * u * u_n + 0.5 * (1 - alfa) * h * c * u + alfa * c * u * u_n^2 * s4g
          + 0.5 * (1 - alfa) * g * h^2 * normal[1] + 0.5 * h * u_n * ((1 + alfa) * c + u_n) * normal[1]
          - lambda1_geom * c_ext * (alfa * c_ext - u_n_ext) * (alfa * u - 2 * c * normal[1]) * s4g
          + lambda2_geom * h_geom * normal[2] * u_t_ext)

    f3 = ((0.25 * alfa - 0.5) * h * v * u_n + 0.5 * (1 - alfa) * h * c * v + alfa * c * v * u_n^2 * s4g
          + 0.5 * (1 - alfa) * g * h^2 * normal[2] + 0.5 * h * u_n * ((1 + alfa) * c + u_n) * normal[2]
          - lambda1_geom * c_ext * (alfa * c_ext - u_n_ext) * (alfa * v - 2 * c * normal[2]) * s4g
          - lambda2_geom * h_geom * normal[1] * u_t_ext)

    # Compute the nonconservative piece; zero for this test case because bottom is 0
    noncons_flux = nonconservative_flux_function(u_inner, u_boundary, normal_direction,
                                                 equations)

    # Return the novel version of the flux
    return SVector(f1, f2, f3, zero(eltype(u_inner))) * s_hat, noncons_flux
end

###############################################################################

# Setup the boundary conditions
boundary_condition = Dict(:Right  => boundary_condition_slip_wall,
                          :Left   => boundary_condition_slip_wall,
                          :Bottom => boundary_condition_subcritical_inflow_stable,
                          :Top    => boundary_condition_subcritical_outflow_stable)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
polydeg = 5
solver = DGSEM(polydeg = polydeg,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Get the unstructured quad mesh from a file
mesh_file = joinpath(@__DIR__, "mesh_swe_channel_wiggles.inp")
mesh = P4estMesh{2}(mesh_file)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition,
                                    source_terms = source_term_moving_blob)

###############################################################################
# ODE solvers, callbacks, etc.

# Gaussian pulse should enter, propagate across the channel length, and exit
tspan = (0.0, 11.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = true)

alive_callback = AliveCallback(analysis_interval = 10000)

save_solution = SaveSolutionCallback(dt = 0.1,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl=0.9)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation  with CFL based time stepping

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            ode_default_options()..., callback = callbacks,
            adaptive = false, dt = 1.0);
