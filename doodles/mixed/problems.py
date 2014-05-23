'Problems used to test mixed solvers.'

from dolfin import Constant, Mesh, Expression, MeshFunction, SubDomain,\
    Rectangle, Circle, plot, Point, CellFunction, cells, refine,\
    FacetFunction, near, sqrt
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

# Define domains
class InflowBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return on_boundary and near(x[0], x_min)

class NoslipBoundary(SubDomain):
  def inside(self, x, on_boundary):
    dx = x[0] - c_x
    dy = x[1] - c_y
    dr = sqrt(dx**2 + dy**2)
    return on_boundary and (near(x[1]*(y_max - x[1]), 0) or dr < r + 1E-3)

f_f = FacetFunction('size_t', mesh, 0)
InflowBoundary().mark(f_f, 13)
NoslipBoundary().mark(f_f, 12)


class CylinderFlow(object):
    'Flow past a cylinder'
    name = 'cylinder'
    # Forcing
    f = Constant((0., 0.))

    # Mesh and function marking facets
    mesh = mesh  # Mesh(mesh_path('cylinder.xml'))
    f_f = f_f    # MeshFunction('size_t', mesh, mesh_path('cylinder_facet_region.xml'))

    # Inflow and Noslip domains to be used for BC construction
    inflow = [f_f, 13]
    noslip = [f_f, 12]

    Re = Constant(1000)
    U_max = 1.5
    u_in = Expression(('4*Um*(x[1]*(ymax-x[1]))*sin(pi*t/8.0)/(ymax*ymax)',
                       '0.0'),
                      Um=U_max, ymax=0.41, t=0)


class LCylinderFlow(object):
    'Flow past a cylinder in the bend of L(V) shaped domain'
    name = 'l-cylinder'
    # Forcing
    f = Constant((0., -1.))

    # Mesh and function marking facets
    mesh = Mesh(mesh_path('l-cylinder.xml'))
    f_f = MeshFunction('size_t', mesh, mesh_path('l-cylinder_facet_region.xml'))

    # Inflow and Noslip domains to be used for BC construction
    inflow = [f_f, 17]
    noslip = [f_f, 16]

    # Inflow profile
    U_max = 3.5
    unit = 0.1  # Basic scaling unit for mesh definition
    s = 3       # Arm length of the channel, ALWAYS SYNC THIS VALUE WITH .geo

    class InflowProfileLCylinder(Expression):
        def __init__(self, U_max, t, unit, s):
            Expression.__init__(self)
            self.U_max, self.t = U_max, t
            # Channel width
            w = 4*unit + 4*unit*cos(pi/4)  # Channel width
            self.width = w
            # Midpoint on the inflow
            M = [4*unit*cos(5*pi/4) + s*sin(5*pi/4) - 0.5*w*cos(5*pi/4),
                 4*unit*sin(5*pi/4) - s*cos(5*pi/4) - 0.5*w*sin(5*pi/4)]
            self.M = M

        def eval(self, values, x):
            U_max, width, t, M = self.U_max, self.width, self.t, self.M
            d = sqrt((x[0] - M[0])**2 + (x[1] - M[1])**2)
            values[0] = U_max*(width/2-d)**2*sin(pi*t/8.0)/(width/2)**2
            values[1] = -U_max*(width/2-d)**2*sin(pi*t/8.0)/(width/2)**2

        def value_shape(self):
            return (2, )

    u_in = InflowProfileLCylinder(U_max=U_max, t=0, unit=unit, s=s)

    # Reynolds number
    Re = Constant(10)

all_problems = [CylinderFlow, LCylinderFlow]
