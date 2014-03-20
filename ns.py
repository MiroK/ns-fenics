#!/usr/bin/env python

__author__ = 'Anders Logg <logg@simula.no>'
__date__ = '2008-04-11'
__copyright__ = 'Copyright (C) 2008-2010 ' + __author__
__license__  = 'GNU GPL version 3 or any later version'

# Modified by Miroslav Kuchta, 2014

import sys, time, os
from dolfin import set_log_active, parameters, list_timings
from problems import Problem, problems
from solvers import Solver, solvers

# List of mesh sizes
mesh_sizes = [8, 11, 16, 23, 32, 45, 64]

# Default options
OPTIONS = {'refinement_level' :      3, # which mesh to use, 3d has internal 
           'viscosity_time_step' :   True, # t-step criterion uses viscosity
           'dt_division' :           2, # refine the timestep
           'save_solution' :         True,
           'save_frequency' :        1,
           'check_mem_usage' :       False,
           'check_frequency' :       10,
           'save_solution_at_t=T' :  True,
           'plot_solution' :         True,
           'plot_functional' :       True,
           'compute_divergence' :    False,
  	       'debug' :                 False,
  	       'max_steps':              None,
           'results_dir' :           './results',
           'krylov_solver_params' :  {'absolute_tolerance': 1e-25,
                                      'relative_tolerance': 1e-12,
                                      'monitor_convergence': False},
           'ffc_compiler_params' :   {'optimize': True, 
                                      'eliminate_zeros': True, 
                                      'precompute_basis_const': True, 
                                      'precompute_ip_const': True}
          }

def save_results(problem, solver, num_dofs, cputime, wct, functional, error):
  'Save results to file.'
  # Print summary
  print ''
  print 'Problem    |', problem
  print 'Solver     |', solver
  print 'Unknowns   |', num_dofs
  print 'CPU time   |', cputime
  print 'WCT time   |', wct
  print 'Overhead   |', wct - cputime
  print 'Functional |', functional
  print 'Error      |', error

  # Print DOLFIN summary
  set_log_active(True)
  list_timings()
  
  # Append to file
  results_dir = problem.options['results_dir']
  filename = os.path.join(results_dir, 'results.log')

  # create the dir for results if needed
  if not os.path.exists(os.path.dirname(filename)):
    os.makedirs(os.path.dirname(filename))

  with open(filename, 'a') as f:
    f.write('%s, %s, %s, %d, %.15g, %.15g, %.15g, %s\n' %
           (time.asctime(), problem, solver, num_dofs, cputime, wct, functional, str(error)))

#------------------------------------------------------------------------------

def usage():
  'Print usage'
  print '''\
Usage: ns problem solver

Available problems:

%s

Available solvers:

%s
''' % ('\n'.join('  ' + p for p in problems),
     '\n'.join('  ' + s for s in solvers))

#------------------------------------------------------------------------------

def main(args):
  'Parse command-line arguments and run solver'

  # Check arguments
  if not len(args) >= 2:
    usage()
    return 2

  # Get problem and solver
  problem_name, solver_name = args[:2]

  # Get options
  options = OPTIONS.copy()
  for arg in args[2:]:
    try:
      key, value = arg.split('=')
      try:
        options[key] = eval(value)
      except:
        options[key] = str(value)
    except:
      print 'Warning: Unhandled command-line argument', arg
  
  # Set cpp_compiler options #TODO What exact does this do?
  parameters["form_compiler"]["cpp_optimize"] = True

  # Set debug level
  set_log_active(options['debug'])

  # Set refinement level
  options['N'] = mesh_sizes[options['refinement_level']]

  # Create problem and solver
  problem = Problem(problem_name, options)
  solver = Solver(solver_name, options)
  print 'Problem: ' + str(problem)
  print 'Solver:  ' + str(solver)

  # Solve problem with solver
  wct = time.time()
  u, p = solver.solve(problem)

  # Compute elapsed time
  wct = time.time() - wct

  # Compute number of degrees of freedom
  num_dofs = u.vector().size() + p.vector().size()

  # Get functional value and error
  functional, error = solver.eval()

  # Save results
  save_results(problem, solver, num_dofs, solver.cputime(), wct, functional, error)

  return 0

#------------------------------------------------------------------------------

if __name__ == '__main__':
  sys.exit(main(sys.argv[1:]))
