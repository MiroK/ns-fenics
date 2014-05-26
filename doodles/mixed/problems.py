'Problems used to test mixed solvers.'

from dolfin import Constant, Mesh, Expression, MeshFunction, SubDomain,\
    Rectangle, Circle, Point, CellFunction, cells, refine,\
    FacetFunction, near, plot
from math import pi, cos, sin, sqrt
import os

# Problem is specified by mesh, domains where Dirichlet bcs should be
# prescribed, Dirichlet values (this includes the inflow profile which
# is a time-dependent expression), the Reynolds number and forcing

# Directory where meshes are stored
mesh_dir = 'meshes'

# Mesh directory aware loading of mesh
mesh_path = lambda mesh: os.path.join(mesh_dir, mesh)

# Build testing mesh for cylinder flow
# geometric parameters
x_min, x_max = 0, 2.2
y_min, y_max = 0, 0.41          # note that cylinder is a bit off center
c_x, c_y, r = 0.2, 0.2, 0.05


# Convenience function for cylinder domain building in FEnicS
def refine_cylinder(mesh):
    'Refine mesh by cutting cells around the cylinder.'
    h = mesh.hmin()
    center = Point(c_x, c_y)
    cell_f = CellFunction('bool', mesh, False)
    for cell in cells(mesh):
        if cell.midpoint().distance(center) < r + h:
            cell_f[cell] = True
    mesh = refine(mesh, cell_f)

    return mesh

# Create the mesh
# First define the domain
rect = Rectangle(x_min, y_min, x_max, y_max)
circ = Circle(c_x, c_y, r)
domain = rect - circ

# Mesh the domain
mesh = Mesh(domain, 45)

# Refine mesh n-times
n = 1
for i in range(n):
    mesh = refine_cylinder(mesh)


class InflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], x_min)


class OutflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], x_max)


class NoslipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        dx = x[0] - c_x
        dy = x[1] - c_y
        dr = sqrt(dx**2 + dy**2)
        return on_boundary and (near(x[1]*(y_max - x[1]), 0) or dr < r + 1E-3)

f_f = FacetFunction('size_t', mesh, 0)
InflowBoundary().mark(f_f, 13)
NoslipBoundary().mark(f_f, 12)
OutflowBoundary().mark(f_f, 15)


# Width of the L domain
def compute_width_L_cylinder(unit):
    return 4*unit + 4*unit*cos(pi/4)


# Width of the O domain
def compute_width_O_cylinder(unit):
    return 4*unit


class InflowProfileConstant(Expression):
    'Constant direction but pulsating magnitude.'
    def __init__(self, U_max, t, unit, s, K, compute_width):
        '''Unit and s characterize domain, t, K - time parameters,
        U_max is the maximal magnitude of inflow velocity.'''
        self.U_max, self.t, self.K = U_max, t, K
        # Channel width
        self.w = compute_width(unit)
        # Midpoint on the inflow
        M = [4*unit*cos(5*pi/4) + s*sin(5*pi/4) - 0.5*self.w*cos(5*pi/4),
             4*unit*sin(5*pi/4) - s*cos(5*pi/4) - 0.5*self.w*sin(5*pi/4)]
        self.M = M

    def eval(self, values, x):
        U_max, w, t, M, K = self.U_max, self.w, self.t, self.M, self.K
        d = sqrt((x[0] - M[0])**2 + (x[1] - M[1])**2)
        mag = U_max*(w/2-d)**2*sin(pi*t/K)/(w/2)**2/sqrt(2)
        values[0] = mag
        values[1] = -mag

    def value_shape(self):
        return (2, )


class InflowProfilePeriodic(Expression):
    'Profile with pulsating magnitide and changing direction.'
    def __init__(self, U_max, t, unit, s, K, compute_width):
        '''Unit and s characterize domain, t, K - time parameters,
        U_max is the maximal magnitude of inflow velocity.'''
        self.U_max, self.t, self.K = U_max, t, K
        # Channel width
        self.w = compute_width(unit)
        # Midpoint on the inflow
        M = [4*unit*cos(5*pi/4) + s*sin(5*pi/4) - 0.5*self.w*cos(5*pi/4),
             4*unit*sin(5*pi/4) - s*cos(5*pi/4) - 0.5*self.w*sin(5*pi/4)]
        self.M = M

    def eval(self, values, x):
        U_max, w, t, M, K = self.U_max, self.w, self.t, self.M, self.K
        d = sqrt((x[0] - M[0])**2 + (x[1] - M[1])**2)
        mag = U_max*(w/2-d)**2*sin(pi*t/K)/(w/2)**2/sqrt(2)
        values[0] = mag*abs(sin(2*pi*t/K))
        values[1] = -mag*abs(cos(2*pi*t/K))

    def value_shape(self):
        return (2, )

# ------------------------------------------------------------------------------


class CylinderFlow(object):
    'Flow past a cylinder'
    name = 'cylinder'
    # Forcing
    f = Constant((0., 0., 0.))
    # versus Constant((0, 0)) gives 11:41 vs 10:56, TODO worth it?

    # Mesh and function marking facets
    mesh = mesh
    # mesh = Mesh(mesh_path('cylinder.xdmf'))
    f_f = f_f
    # f_f =MeshFunction('size_t', mesh, mesh_path('cylinder_facet_region.xdmf'))

    # Inflow and Noslip domains to be used for BC construction
    inflow = [f_f, 13]
    noslip = [f_f, 12]

    # Duration
    T = 8

    Re = Constant(1000)
    U_max = 1.5
    u_in = Expression(('4*Um*(x[1]*(ymax-x[1]))*sin(pi*t/K)/(ymax*ymax)',
                       '0.0'),
                      Um=U_max, ymax=0.41, t=0, K=T)


class LCylinderFlowConstant(object):
    'Flow past a cylinder in the bend of L(V) shaped domain. Constant force.'
    name = 'l-cylinder-constant'
    # Forcing
    f = Constant((0., 0., 0.))
    # Versus Constant((0, 0)) gives 11:41 vs 10:56, TODO worth it?

    # Mesh and function marking facets
    mesh = Mesh(mesh_path('l-cylinder.xdmf'))
    f_f = MeshFunction('size_t', mesh, mesh_path('l-cylinder_facet_region.xdmf'))

    # Inflow and Noslip domains to be used for BC construction
    inflow = [f_f, 17]
    noslip = [f_f, 16]

    # Duration
    T = 10

    # Inflow profile
    U_max = 3.5
    unit = 0.1  # Basic scaling unit for mesh definition
    s = 3       # Arm length of the channel, ALWAYS SYNC THIS VALUE WITH .geo

    u_in = InflowProfileConstant(U_max=U_max, t=0, unit=unit, s=s, K=T,
                                 compute_width=compute_width_L_cylinder)

    # Reynolds number
    Re = Constant(100)


class LCylinderFlowPeriodic(LCylinderFlowConstant):
    'Flow past a cylinder in V turn. Periodic forcing'
    name = 'l-cylinder-periodic'
    u_in = InflowProfilePeriodic(U_max=LCylinderFlowConstant.U_max,
                                 t=0,
                                 unit=LCylinderFlowConstant.unit,
                                 s=LCylinderFlowConstant.s,
                                 K=LCylinderFlowConstant.T,
                                 compute_width=compute_width_L_cylinder)


class OCylinderFlowConstant(object):
    'Flow past a cylinder in the bend of O shaped turn. Constant force.'
    name = 'o-cylinder-constant'
    # Forcing
    f = Constant((0., 0., 0.))
    # Versus Constant((0, 0)) gives 11:41 vs 10:56, TODO worth it?

    # Mesh and function marking facets
    mesh = Mesh(mesh_path('o-cylinder.xdmf'))
    f_f = MeshFunction('size_t', mesh, mesh_path('o-cylinder_facet_region.xdmf'))

    # Inflow and Noslip domains to be used for BC construction
    inflow = [f_f, 18]
    noslip = [f_f, 17]

    # Duration
    T = 10

    # Inflow profile
    U_max = 3.5
    unit = 0.1  # Basic scaling unit for mesh definition
    s = 3       # Arm length of the channel, ALWAYS SYNC THIS VALUE WITH .geo

    u_in = InflowProfileConstant(U_max=U_max, t=0, unit=unit, s=s, K=T,
                                 compute_width=compute_width_O_cylinder)

    # Reynolds number
    Re = Constant(100)


class OCylinderFlowPeriodic(OCylinderFlowConstant):
    'Flow past a cylinder in the bend of O shaped turn. Periodic force.'
    name = 'o-cylinder-periodic'
    u_in = InflowProfilePeriodic(U_max=LCylinderFlowConstant.U_max,
                                 t=0,
                                 unit=LCylinderFlowConstant.unit,
                                 s=LCylinderFlowConstant.s,
                                 K=LCylinderFlowConstant.T,
                                 compute_width=compute_width_O_cylinder)

all_problems = [CylinderFlow,
                LCylinderFlowConstant, LCylinderFlowPeriodic,
                OCylinderFlowConstant, OCylinderFlowPeriodic]

if __name__ == '__main__':
    for problem in all_problems:
        plot(problem.mesh)
        plot(problem.f_f, interactive=True)
