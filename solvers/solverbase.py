__author__ = 'Anders Logg <logg@simula.no>'
__date__ = '2008-04-03'
__copyright__ = 'Copyright (C) 2008-2010 ' + __author__
__license__  = 'GNU GPL version 3 or any later version'

# Modified by Miroslav Kuchta, 2014

from dolfin import *

from math import *
from numpy import linspace
from time import time
from commands import getoutput
import os

# Common solver parameters
maxiter = default_maxiter = 200
tolerance = default_tolerance = 1e-4

class SolverBase:
  'Base class for all solvers.'
  def __init__(self, options):
    # Store options
    self.options = options

    # Reset some solver variables
    self._time = None
    self._cputime = 0.0
    self._timestep = 0

    # Reset files for storing solution
    self._ufile = None
    self._pfile = None

    # Reset storage for functional values and errors
    self._t = []
    self._M = []
    self._m = []
    self._e = []

  def get_timestep(self, problem):
    'Return time step and number of time steps for problem.'
    T = problem.T
    U = problem.U
    nu = problem.nu
    h  = problem.mesh.hmin()

    # Set the time step from problem if available
    if not problem.dt is None:
      dt = problem.dt
      print dt
      n  = int(T / dt)
      print 'Using problem.dt'
    else:
      dt =  0.2*(h / U)
      n  = int(T / dt + 1.0)
      dt = T / n
      print 'Computing time step according to stability criteria'

    # Compute range
    t_range = linspace(0, T, n+1)[1:]

    # Report time step
    print ' '
    print 'Number of timesteps:' , len(t_range)
    print 'Size of timestep:' , dt
    print ' '

    return dt, t_range[0], t_range

  def getMyMemoryUsage(self):
    mypid = os.getpid()
    mymemory = getoutput('ps -o rss %s' % mypid).split()[1]
    return mymemory

  def start_timing(self):
    '''Start timing, will be paused automatically during update
    and stopped when the end-time is reached.'''
    self._time = time()

  def solve(self, problem, dt, plot_solution=True):
    'Solve problem.'
    raise NotImplementedError

  def prefix(self, problem):
    '''Return file prefix for output files:
    results_dir/problem.output_locations/problem/solver.'''
    root = self.options['results_dir']
    p = problem.__module__.split('.')[-1].lower()
    s = self.__module__.split('.')[-1].lower()
    return os.path.join(root, problem.output_location, p, s)

  def update(self, problem, t, u, p):
    'Update problem at time t.'
    # Add to accumulated CPU time
    timestep_cputime = time() - self._time
    self._cputime += timestep_cputime

    # Compute divergence
    if self.options['compute_divergence']:
      check_divergence(u, p.function_space())

    problem.update_problem(t, u, p)

    # Evaluate functional and error
    m = problem.reference(t)
    M = problem.functional(t, u, p)
    if m is None:
      e = None
      print 'M = %g (missing reference value)' % M
    else:
      e = abs(M - m)
      print 'M = %g (reference %g), error = %g (maximum %g)' % (M, m, e, max([e] + self._e))

    # Store values
    self._t.append(t)
    self._M.append(M)
    self._m.append(m)
    self._e.append(e)

    # If there is going to be any saving make the directories
    if any(key for key in self.options.keys() if 'save' in key.lower()):
      results_dir = self.prefix(problem) 
      if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    
    # Save solution
    if self.options['save_solution']:
      # Save velocity and pressure
      frequency = self.options['save_frequency']
      refinement = self.options['refinement_level']
      if (self._timestep - 1) % frequency == 0:
        # Create files for saving
        if self._ufile is None:
          u_fname = os.path.join(results_dir, 'refinement_level_' + str(refinement) + '_u.pvd')
          self._ufile = File(u_fname)
        if self._pfile is None:
          p_fname = os.path.join(results_dir, 'refinement_level_' + str(refinement) + '_p.pvd')
          self._pfile = File(p_fname)
        self._ufile << u
        self._pfile << p

    # Save solution at t = T
    if self.options['save_solution_at_t=T']:
      if t >= problem.T:
        refinement = self.options['refinement_level']
        # Create files for saving, only happens once
        u_fname = os.path.join(results_dir, 'refinement_level_' + str(refinement) + '_at_end' +'_u.pvd')
        ufile = File(u_fname)
        p_fname = os.path.join(results_dir, 'refinement_level_' + str(refinement) + '_at_end' + '_p.pvd')
        pfile = File(p_fname)
        ufile << u
        pfile << p

    # Plot solution
    if self.options['plot_solution']:
      # Plot velocity and pressure
      plot(u, title='Velocity', rescale=True)
      plot(p, title='Pressure', rescale=True)

    # Check memory usage
    if self.options['check_mem_usage']:
      if (self._timestep - 1) % self.options['check_frequency'] == 0:
        print 'Memory usage is:' , self.getMyMemoryUsage()

    # Print progress
    print ''
    s = 'Time step %d finished in %g seconds, %g%% done (t = %g, T = %g).' \
        % (self._timestep, timestep_cputime, 100.0*(t / problem.T), t, problem.T)
    print s + '\n' + len(s)*'-'

    # Increase time step and record current time
    self._timestep += 1
    self._time = time()

  def eval(self):
    '''Return last functional value and maximum error in functional
    value on [0, T].'''
    # Plot values
    if self.options['plot_functional']:
      import matplotlib.pyplot as plt
      fig = plt.figure()
      plt.plot(self._t, self._M)
      plt.xlabel(r'$t$')
      plt.ylabel(r'Functional')
      plt.grid(True)
      plt.show()

    # Return value
    if self._e[0] is None:
      return self._M[-1], None
    else:
      return self._M[-1], max([0.0] + self._e)

  def cputime(self):
    'Return accumulated CPU time.'
    return self._cputime

#------------------------------------------------------------------------------

def sigma(u, p, nu):
  'Return stress tensor.'
  return 2*nu*sym(grad(u)) - p*Identity(u.cell().d)

#------------------------------------------------------------------------------

def has_converged(r, itr, method, maxiter=default_maxiter, tolerance=default_tolerance):
  'Check if solution has converged.'
  print 'Residual = ', r
  if r < tolerance:
    print '%s iteration converged in %d iteration(s).' % (method, itr + 1)
    return True
  elif itr == maxiter - 1:
    raise RuntimeError, '%s iteration did not converge.' % method
  return False

#------------------------------------------------------------------------------

def check_divergence(u, Q):
  'Check divergence of velocity.'
  # Compute L2 norm of divergence
  print '||div u||_L2 =', norm(u, 'Hdiv0')

  # Compute projection of div u into Q_0
  pdivu = project(div(u), Q)
  zero = Constant(Q.mesh(), 0.0)
  bc = DirichletBC(Q, zero, DomainBoundary())
  bc.apply(pdivu.vector())

  # Compute 'weak' L2 norm of divergence
  print '||div u||_w  =', sqrt(abs(assemble(pdivu*div(u)*dx, mesh=Q.mesh())))
