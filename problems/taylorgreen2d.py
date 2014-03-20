__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2008-09-15"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Harish Narayanan, 2009.
# Modified by Anders Logg, 2010.
# Modified by Miroslav Kuchta, 2014

from math import pi
from problembase import *

# Problem definition
class Problem(ProblemBase):
  "2D periodic problem with known analytical solution."
  def __init__(self, options):
    ProblemBase.__init__(self, options)

    # Create mesh
    N = options["N"]
    self.mesh = RectangleMesh(-1, -1, 1, 1, N, N)
   
    # Set the constrained domain
    self.constrained_domain = PeriodicDomain()

    # The body force term
    self.f = Constant((0, 0))

    # Set viscosity
    self.nu = 1.0 / 100.0

    # Set the final time
    self.T = 0.5

    # Exact solutions, strings for expression
    self.Ux = ('-(cos(pi*(x[0]))*sin(pi*(x[1]))) * exp(-2.0*nu*pi*pi*t)',
              ' (cos(pi*(x[1]))*sin(pi*(x[0]))) * exp(-2.0*nu*pi*pi*t)')
    self.Px = '-0.25*(cos(2*pi*(x[0])) + cos(2*pi*(x[1]))) * exp(-4*nu*pi*pi*t)'

  def initial_conditions(self, V, Q):
    # Use analytical solutions at t = 0 as initial values
    self.exact_u = Expression(self.Ux, nu=self.nu, t=0., degree=3)
    self.exact_p = Expression(self.Px, nu=self.nu, t=0., degree=3)

    # Construct the initial conditions by setting t = 0 in the exact
    self.exact_u.t = 0.
    self.exact_p.t = 0.

    return self.exact_u, self.exact_p

  def boundary_conditions(self, v, q, t):
    # There are no bcs.
    return [], [], None 

  def update(self, t, u, p):
    self.exact_u.t = t
    self.exact_p.t = t
    pass

  def functional(self, t, u, p):
    if t < self.T:
        return 0
    else:
        return 0.5*norm(u)**2

  def reference(self, t):
    if t < self.T:
        return 0
    else:
        self.exact_u.t = t
        return 0.5*norm(self.exact_u,  mesh=self.mesh)**2

  def __str__(self):
    return "Taylor-Green 2D"

#------------------------------------------------------------------------------

class PeriodicDomain(SubDomain):
  # Using Mikael's code, avoid bottom right and top left corners
  # In general target is left and bottom wall
  def inside(self, x, on_boundary):
    is_lr = near(x[0], 1) and near(x[1], -1) and on_boundary
    is_ul = near(x[0], -1) and near(x[1], 1) and on_boundary

    return (near(x[0], -1) or near(x[1], -1)) and not (is_lr or is_ul)

  # Map to target
  def map(self, x, y):
    if near(x[0], 1) and near(x[1], 1):
      y[0] = x[0] - 2.
      y[1] = x[1] - 2.
    elif near(x[0], 1):
      y[0] = x[0] - 2.
      y[1] = x[1]
    else:
      y[0] = x[0]
      y[1] = x[1] - 2.
