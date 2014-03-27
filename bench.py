#!/usr/bin/env python

__author__ = 'Anders Logg <logg@simula.no>'
__date__ = '2009-10-12'
__copyright__ = 'Copyright (C) 2009-2010 ' + __author__
__license__  = 'GNU GPL version 3 or any later version'

# Modified by Miroslav Kuchta, 2014

import sys, commands
from problems import problems as bench_problems
from solvers import solvers as bench_solvers

def usage():
  'Print usage'
  print '''\
Usage: bench          (run all problems/solvers)
       bench problem  (run all solvers for given problem)
       bench solver   (run all problems for given solver)
'''

#------------------------------------------------------------------------------

def main(args):
  'Parse command-line arguments and run solver.'
  # Extract arguments (note: don't change the order here!)
  parameters = [arg for arg in args if '=' in arg] # exatract ns params
  args = [arg for arg in args if not '=' in arg]#no= -> can be solver or problem

  # Decide list of problems and solvers based on command-line arguments
  if len(args) == 0:
    problems = bench_problems
    solvers = bench_solvers
  if len(args) == 1 and args[0] in bench_problems:
    problems = args
    solvers = bench_solvers
  elif len(args) == 1 and args[0] in bench_solvers:
    problems = bench_problems
    solvers = args
  elif len(args) > 0:
    usage()
    exit(1)

  # Search for num processes in parameters
  num_processes = 1        # one process is used if not more given
  for param in parameters:
    if 'num_processes' in param:
      num_processes = int(param.split('=')[1])
      parameters.remove(param)
      break
  print_color('Using %d processes.' % num_processes, 'yellow')

  # Search for results_dir in parameters
  use_default_results_dir = True
  for param in parameters:
    if 'results_dir' in param:
      use_default_results_dir = False # result dir is fed to ns as cmdl arg
      break

  # If not found take it from ns
  if use_default_results_dir:
    from ns import OPTIONS as options
    try:
      options['results_dir']          # result dir will be found from options
    except KeyError:
      print 'Set resuts_dir in ns.OPTIONS or provide is command-line argument!'
      exit()

  # Iterate over problems and solvers
  for (i, problem) in enumerate(problems):
    print_color('  ' + problem, 'green')
    for (j, solver) in enumerate(solvers):
      print_color('    ' + solver, 'blue')

      # Solve problem with solver
      output = \
      commands.getoutput('mpirun -np %d ./ns.py %s %s %s' % (num_processes,
                                                             problem,
                                                             solver,
                                                             ' '.join(parameters)))
      #try:
      lines = output.split('\n')
      kws = ['Error', 'Time step', 'Mesh size']
      results = { kw : [] for kw in kws}
      for line in lines:
        for kw in kws:
          if kw in line:
            results[kw].append(line.split('|')[1].strip())

      n = len(results[kws[0]])
      for i in range(n):
        output = []
        for kw in kws:
          output.append('%s = %8s' % (kw, results[kw][i]))
        message = ' | '.join(output)
        print_color('       %s' % message, 'cyan')

      #except:
      #  print_color('\t\t\t FAILED', 'red')

  return 0

#------------------------------------------------------------------------------

def print_color(string, color):
  'Print string in given color.'
  if color == 'blue':
    print '\033[1;37;34m%s\033[0m' % string
  elif color == 'green':
    print '\033[1;37;32m%s\033[0m' % string
  elif color == 'red':
    print '\033[1;37;31m%s\033[0m' % string
  elif color == 'yellow':
    print '\033[1;33;33m%s\033[0m' % string
  elif color == 'cyan':
    print '\033[0;36;36m%s\033[0m' % string
  elif color == 'pink':
    print '\033[1;31;31m%s\033[0m' % string
  else:
    print string

#------------------------------------------------------------------------------

if __name__ == '__main__':
  sys.exit(main(sys.argv[1:]))
