
using HOHQMesh

# Create a new HOHQMesh model project. The project name
# "mesh_swe_channel_wiggles" will be the name of the mesh file
# saved in the directory "out".
channel_mesh = newProject("mesh_swe_channel_wiggles", "out")

# Reset polynomial order and output formats
setPolynomialOrder!(channel_mesh, 3)
setMeshFileFormat!(channel_mesh, "ABAQUS") # for P4estMesh{2}
setPlotFileFormat!(channel_mesh, "sem")

# A background grid is required for the mesh generation.
addBackgroundGrid!(channel_mesh, [0.8655, 0.8655, 0.0])

# Create all the outer boundary "curves".
# Note, the curve names are those that will be present in the mesh file.

# Convenience variables to define the bounds of the domain
r = 0.0
R = 3.0

# Cosine and sine values of rotation angle pi/4
c = 1.0 / sqrt(2.0)
s = 1.0 / sqrt(2.0)

# Right side
right = newEndPointsLineCurve("Right", [r , r, 0.0],
                                       [2.0*R*c - s*r, 2.0*R*s + c*r, 0.0])

# Top side
xEqn = "x(t) = -1.0/sqrt(2)*(3*t-6) - 1.0/sqrt(2)*(-0.3*cos(4*pi*t)*sin(2*pi*t))"
yEqn = "y(t) = 1.0/sqrt(2)*(3*t+6) - 1.0/sqrt(2)*(-0.25*sin(4*pi*t)*cos(2*pi*t))"
zEqn = "z(t) = 0.0"
top = newParametricEquationCurve("Top", xEqn, yEqn, zEqn)


# Left side
left = newEndPointsLineCurve("Left", [2.0*R*c - s*R, 2.0*R*s + c*R, 0.0],
                                     [r*c - s*R, r*s + c*R, 0.0])

# rotated parametric curve to force the inflow boundary to be "wobbly".
# OBS! we use 1-t in order to preserve mesh handedness

# Bottom side
xEqn = "x(t) = -1.0/sqrt(2)*(3*(1-t)) - 1.0/sqrt(2)*(-0.265*cos(4*pi*(1-t))*sin(2*pi*(1-t)))"
yEqn = "y(t) = 1.0/sqrt(2)*(3*(1-t)) - 1.0/sqrt(2)*(-0.25*cos(4*pi*(1-t))*sin(2*pi*(1-t)))"
zEqn = "z(t) = 0.0"

bottom = newParametricEquationCurve("Bottom", xEqn, yEqn, zEqn)

# Add outer boundary curves in counter-clockwise order.
addCurveToOuterBoundary!(channel_mesh, bottom)
addCurveToOuterBoundary!(channel_mesh, right)
addCurveToOuterBoundary!(channel_mesh, top)
addCurveToOuterBoundary!(channel_mesh, left)

# Generate the mesh. Saves the mesh file to the directory "out".
generate_mesh(channel_mesh)
