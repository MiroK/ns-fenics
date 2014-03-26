#!/usr/bin/env python

__author__ = 'Anders Logg <logg@simula.no>'
__date__ = '2008-04-11'
__copyright__ = 'Copyright (C) 2008-2010 ' + __author__
__license__  = 'GNU GPL version 3 or any later version'

# Modified by Miroslav Kuchta, 2014

import sys, time, os, copy
from dolfin import set_log_active, parameters, list_timings, MPI
from problems import Problem, problems
from solvers import Solver, solvers
from bench import print_color

# List of mesh sizes
mesh_sizes = [8, 11, 16, 23, 32, 45, 64]

# Default options
OPTIONS = {'refinement_level' :      [0, 1, 2],\
           # which mesh to use, 3d has internal, either int - use given mesh
           # or list - use sequence of meshes, suitable for convergence rate
           # studies
           'viscosity_time_step' :   False,       # t-step criterion uses viscosity
           'dt_division' :           0,\
           # refine the timestep, either int - use given time step
           # or list - use sequence of refined time steps, suitable for
           # convergence rate studies
           'verbose' :               False, # print progress info
           'save_solution' :         True,
           'save_frequency' :        100,
           'check_mem_usage' :       False,
           'check_frequency' :       10,
           'save_solution_at_t=T' :  True,
           'plot_solution' :         False,
           'plot_functional' :       False,
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

def save_results(problem, solver, num_dofs, mesh_size, time_step, functional, error):
  'Save results to file.'
  # Print summary
  if MPI.process_number() == 0 :
    print ''
    print 'Problem    |', problem
    print 'Solver     |', solver
    print 'Unknowns   |', num_dofs
    print 'Mesh size  |', mesh_size
    print 'Time step  |', time_step
    print 'Functional |', functional
    print 'Error      |', error

    # Print DOLFIN summary
    set_log_active(True)
    list_timings()
    
    # Append to file, let each dx, dt have its own log file
    results_dir = problem.options['results_dir']
    dx = problem.options['refinement_level']
    dt = problem.options['dt_division']

    name = '%s_%s_results_dx%d_dt%d.log' % (str(problem), str(solver), dx, dt)
    filename = os.path.join(results_dir, name)

    # create the dir for results if needed
    if not os.path.exists(os.path.dirname(filename)):
      os.makedirs(os.path.dirname(filename))

    with open(filename, 'a') as f:
      f.write('%s, %s, %s, %d, %.15g, %.15g, %.15g, %s\n' %
             (time.asctime(), problem, solver, num_dofs, mesh_size, time_step, functional, str(error)))

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
  
  # Parse the spatial and time stepping options.
  # If either refinement_level or dt_division are lists feed the items of
  # list one by one
  type_dx = type(options['refinement_level']) 
  type_dt = type(options['dt_division'])

  if type_dx is list:
    if type_dt is int:
      # Do a dx convergence study
      if MPI.process_number() == 0:
        print_color('Dx convergence study, dt=%d' % options['dt_division'], 'cyan')
      refinement_levels = copy.copy(options['refinement_level'])
      
      # loop over dx
      for refinement_level in refinement_levels:
        if MPI.process_number() == 0:
          print_color('\tRefinement level %d' % refinement_level, 'yellow')
        options['refinement_level'] = refinement_level
        # solve with fixed dt and changing dx
        solve(solver_name, problem_name, options)
    
    else:
      raise ValueError('Use fixed dt_division with changing refinement_level!')
  
  elif type_dx is int:
    if type_dt is list:
      # Do a dt convergence study
      if MPI.process_number() == 0:
        print_color('Dt convergence study, dx=%d' % options['refinement_level'],'cyan')
      dt_divisions = copy.copy(options['dt_division'])

      # loop over dt
      for dt_division in dt_divisions:
        if MPI.process_number() == 0:
          print_color('\tDt division %d' % dt_division, 'yellow')  
        options['dt_division'] = dt_division
        #solve with fixed dx and changing dt
        solve(solver_name, problem_name, options)
    
    elif type(options['dt_division']) is int:
      if MPI.process_number() == 0:
        print_color('Fixed dx and dt', 'cyan')
      # solve with fixed dx and fixed dt
      solve(solver_name, problem_name, options)
    
    else:
      raise ValueError("options['dt_division'] must be int of list!")
  else:
      raise ValueError("options['refinement_level'] must be int of list!")

  return 0

#------------------------------------------------------------------------------

def solve(solver_name, problem_name, options):
  'Solve the problem by solver with options.'

  # Set cpp_compiler options 
  parameters["form_compiler"]["cpp_optimize"] = True

  # Set debug level
  set_log_active(options['debug'])

  # Set refinement level
  options['N'] = mesh_sizes[options['refinement_level']]

  # Create problem and solver
  problem = Problem(problem_name, options)
  solver = Solver(solver_name, options)

  time_step = solver.get_timestep(problem)[0]

  if MPI.process_number() == 0 and options['verbose']:
    print 'Problem: ' + str(problem)
    print 'Solver:  ' + str(solver)

  # Solve problem with solver
  wct = time.time()
  u, p = solver.solve(problem)

  # Compute elapsed time
  wct = time.time() - wct

  # Compute number of degrees of freedom
  num_dofs = u.vector().size() + p.vector().size()

  # Get the mesh size
  mesh_size = u.function_space().mesh().hmin()

  # Get functional value and error
  functional, error = solver.eval()

  # Save results
  cpu_time = solver.cputime()
  save_results(problem, solver, num_dofs, mesh_size, time_step, functional, error)

  return 0

#------------------------------------------------------------------------------

if __name__ == '__main__':
  sys.exit(main(sys.argv[1:]))
