__author__ = 'Anders Logg <logg@simula.no>'
__date__ = '2008-03-19'
__copyright__ = 'Copyright (C) 2008-2010 ' + __author__
__license__  = 'GNU GPL version 3 or any later version'

# Modified by Miroslav Kuchta, 2014

from problembase import *

# Constants related to the geometry
x_min, x_max = 0, 2.2            # The domain is rect [x_min, x_max] x 
y_min, y_max = 0, 0.41           # [y_min, y_max] with cylinder at of radius r
c_x, c_y, r = 0.2, 0.2, 0.05     # centered at c_x, c_y

# Problem definition
class Problem(ProblemBase):
  'Flow past a cylinder.' 
  def __init__(self, options):
    ProblemBase.__init__(self, options)

    # Mesh is created by refining the region around cylinder 3 times
    # the initial resolution is given by 25 + 10*refinement_level
    refinement_level = options['refinement_level']
    resolution = 25 + 10*refinement_level

    rect = Rectangle(x_min, y_min, x_max, y_max)
    circ = Circle(c_x, c_y, r)
    domain = rect - circ
    mesh = Mesh(domain, resolution)

    for i in range(3):
      mesh = refine_cylinder(mesh, c_x, c_y, r)
    self.mesh = mesh

    # Create right-hand side function
    self.f =  Constant((0, 0))

    # Set viscosity (Re = 1000)
    self.nu = 1.0 / 1000.0

    # Characteristic velocity in the domain (used to determinde timestep)
    self.U = 3.5

    # Set end time
    self.T  = 8.0

  def initial_conditions(self, V, Q):
    u0 = Constant((0, 0))
    p0 = Constant(0)
    return u0, p0

  def boundary_conditions(self, V, Q, t):
    # Create inflow boundary condition
    self.g0 = Expression(('4*Um*(x[1]*(ymax-x[1]))*sin(pi*t/8.0)/(ymax*ymax)',\
                          '0.0'),
                         Um=1.5, ymax=y_max, t=t)
    bcu_inflow = DirichletBC(V, self.g0, InflowBoundary())

    # Create no-slip boundary condition
    bcu_outflow = DirichletBC(V, Constant((0., 0.)), NoslipBoundary())

    # Create outflow boundary condition for pressure
    bcp_outflow = DirichletBC(Q, Constant(0.), OutflowBoundary())

    # Collect boundary conditions
    bcs_u = [bcu_inflow, bcu_outflow]
    bcs_p = [bcp_outflow]

    return bcs_u, bcs_p, None

  def update(self, t, u, p, f):
    'Update the expression for inlet velocity.'
    self.g0.t = t

  def functional(self, t, u, p):
    'Pressure difference between center points at front and back of cylinder.'
    if t < self.T:
      return 0.0
    else:
      front = Point(c_x - r - DOLFIN_EPS, c_y)
      back = Point(c_x + r + DOLFIN_EPS, c_y)
      p_f = self.parallel_eval(p, front)
      p_b = self.parallel_eval(p, back)
      return p_f - p_b

  def reference(self, t):
    if t < self.T:
      return 0.0
    else:
      return -0.111444953719

  def __str__(self):
    return 'Cylinder'

#---------------------------------------------------------------------

def refine_cylinder(mesh, c_x, c_y, r):
  'Refine cells around the cylinder.'
  h = mesh.hmin()
  center = Point(c_x, c_y)
  cell_f = CellFunction('bool', mesh, False)
  for cell in cells(mesh):
    if cell.midpoint().distance(center) < r + h:
      cell_f[cell] = True
  mesh = refine(mesh, cell_f)

  return mesh

#---------------------------------------------------------------------

class InflowBoundary(SubDomain):
  'Left wall.'
  def inside(self, x, on_boundary):
    return on_boundary and near(x[0], x_min)

#---------------------------------------------------------------------

class OutflowBoundary(SubDomain):
  'Right wall.'
  def inside(self, x, on_boundary):
    return on_boundary and near(x[0], x_max)

#---------------------------------------------------------------------

class NoslipBoundary(SubDomain):
  'Top, bottom wall and cylinder surface.'
  def inside(self, x, on_boundary):
    dx = x[0] - c_x
    dy = x[1] - c_y
    dr = sqrt(dx**2 + dy**2)
    return on_boundary and (near(x[1]*(y_max - x[1]), 0) or dr < r + 1E-3)
