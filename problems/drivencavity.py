__author__ = 'Anders Logg <logg@simula.no>'
__date__ = '2008-03-19'
__copyright__ = 'Copyright (C) 2008-2010 ' + __author__
__license__  = 'GNU GPL version 3 or any later version'

# Modified by Miroslav Kuchta 2014

from problembase import *
from numpy import cos, pi

# Problem definition
class Problem(ProblemBase):
  '2D lid-driven cavity test problem with known reference value.'
  def __init__(self, options):
    ProblemBase.__init__(self, options)

    # Create mesh and skew towards walls  as in Oasis
    N = options['N']
    self.mesh = UnitSquare(N, N)
    x = self.mesh.coordinates()
    x[:] = (x - 0.5) * 2
    x[:] = 0.5*(cos(pi*(x-1.) / 2.) + 1.) 
    
    # Create right-hand side function
    self.f = Constant((0, 0))

    # Set viscosity (Re = 1000)
    self.nu = 1.0 / 1000.0
    self.U = 1.0

    # Set end-time
    self.T = 2.5

  def initial_conditions(self, V, Q):
    u0 = Constant((0, 0))
    p0 = Constant(0)

    return u0, p0

  def boundary_conditions(self, V, Q, t):
    no_slip_value = Constant((0., 0.))
    driven_value = Constant((1., 0.))
    bc_no_slip = DirichletBC(V, no_slip_value, no_slip_domain)
    bc_driven = DirichletBC(V, driven_value, driven_domain)

    bcs_u = [bc_no_slip, bc_driven]
    bcs_p = []
    # the last return value is meant for surface pressure value with some
    # solvers 
    return bcs_u, bcs_p, None

  def functional(self, t, u, p):
    # Only check final time
    if t < self.T:
      return 0
    else:
    # Compute stream function and report minimum
      psi = stream_function(u)
      psi_min = psi.vector().min()

      print 'Stream function has minimal value' , psi_min

      return psi_min

  def reference(self, t):
    # Only check final time
    if t < self.T:
      return 0.0

    return -0.061076605

  def __str__(self):
      return 'Driven cavity'

#--------------------------------------------------------------------

def no_slip_domain(x, on_boundary):
  'The walls of the driven cavity problem.'
  return on_boundary and (near(x[1], 0) or near(x[0]*(1 - x[0]), 0))

#--------------------------------------------------------------------

def driven_domain(x, on_boundary):
  'The lid of the driven cavity problem.'
  return on_boundary and near(x[1], 1)

#---------------------------------------------------------------------

def stream_function(u):
  'Compute stream function of given 2-d velocity vector.'
  V = u.function_space().sub(0).collapse()

  if V.mesh().topology().dim() != 2:
    raise ValueError('Only stream function in 2D can be computed.')

  psi = TrialFunction(V)
  phi = TestFunction(V)

  a = inner(grad(psi), grad(phi))*dx
  L = inner(u[1].dx(0) - u[0].dx(1), phi)*dx
  bc = DirichletBC(V, Constant(0.), DomainBoundary())

  A, b = assemble_system(a, L, bc)
  psi = Function(V)
  solve(A, psi.vector(), b)

  return psi
