__author__ = 'Anders Logg <logg@simula.no>'
__date__ = '2008-03-19'
__copyright__ = 'Copyright (C) 2008-2010 ' + __author__
__license__  = 'GNU GPL version 3 or any later version'

# Modified by Kent-Andre Mardal, 2008.
# Modified by Miroslav Kuchta, 2014

from dolfin import *
from math import *

class ProblemBase:
  'Base class for all problems.'
  def __init__(self, options):
    # Store options
    self.options = options

    # Parameters must be defined by subclass
    self.mesh = None     # mesh
    self.f = None        # force
    self.bcs_u = []      # boundary conditions on velocity
    self.bcs_p = []      # boundary conditions on pressure
    self.nu = None       # kinematic viscosity
    self.t = 0           # start time
    self.T = None        # end time
    self.dt = None       # time step
    self.u0 = None       # initial condition for velocity
    self.p0 = None       # initial condition for pressure
    self.u = None        # output velocity
    self.p = None        # output pressure
    self.U = 1.0         # maximal velocity, for time step computation
    self.output_location = ''

  def update_problem(self, t, u, p, f):
    'Update problem at time t.'
    # Update state
    self.t = t
    self.u = u
    self.p = p

    # Call problem-specific update
    self.update(t, u, p, f)

  def update(self, t, u, p, f):
    'Problem-speficic update at time t.'
    pass

  def functional(self, t, u, p):
    'Return value of functional of interest.'
    return 0.0

  def reference(self, t):
    'Return reference value for functional.'
    return None

  def tolerance(self, problem):
    'Return tolerance (used as local convergence criterion).'
    if str(problem) == 'Channel':
        return 1e-11
    elif str(problem) == 'Cylinder':
        return 1e-7
    else:
        return 1e-6
