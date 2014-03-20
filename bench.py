#!/usr/bin/env python

__author__ = 'Anders Logg <logg@simula.no>'
__date__ = '2009-10-12'
__copyright__ = 'Copyright (C) 2009-2010 ' + __author__
__license__  = 'GNU GPL version 3 or any later version'

# Modified by Miroslav Kuchta, 2014

import sys, commands, time
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

  # Search for results_dir in parameters
  use_default_results_dir = True
  for param in parameters:
    if 'results_dir' in param:
      resuts_dir = param.split('=')[1]
      use_default_results_dir = False
      break
  
  # If not found take it from ns
  if use_default_results_dir:
    from ns import OPTIONS as options
    try:
      results_dir = options['results_dir']
    except KeyError:
      print 'Set resuts_dir in ns.OPTIONS or provide is command-line argument!'
      exit()
      
  # Collect results in tables
  cputimes = {}
  errors = {}

  # Iterate over problems and solvers
  for (i, problem) in enumerate(problems):
    print_color(problem, 'green')
    for (j, solver) in enumerate(solvers):
      print_color('  ' + solver, 'blue')

      # Solve problem with solver
      output = commands.getoutput('./ns.py %s %s %s' % (problem,
                                                        solver,
                                                        ' '.join(parameters)))
      
      # Collect results
      try:
        cputime = float(output.split('CPU time   |')[1].split('\n')[0])
        wct = float(output.split('WCT time   |')[1].split('\n')[0])
        error = float(output.split('Error      |')[1].split('\n')[0])
        cputimes[(i, j)] = (problem, solver, cputime)
        errors[(i, j)] = (problem, solver, error)
        print '  Finished in %g seconds' % wct
      except:
        cputimes[(i, j)] = (problem, solver, '***')
        errors[(i, j)] = (problem, solver, '***')
        print_color('  *** Failed!', 'red')

      log_file = '/'.join([results_dir, 'bench.log'])
      with open('results/bench.log', 'a') as f:
        f.write('\n' + time.asctime() + '\n')

  # Report results
  print_table(cputimes, 'CPU times')
  print_table(errors, 'Errors')

  return 0

#------------------------------------------------------------------------------

def print_table(values, title):
  'Print nicely formatted table.'

  m = max([key[0] for key in values]) + 2
  n = max([key[1] for key in values]) + 2

  table = []
  for i in range(m):
    table.append(['' for j in range(n)])

  for i in range(m - 1):
    table[i + 1][0] = str(values[(i, 0)][0]).split(' ')[0]

  for j in range(n - 1):
    table[0][j + 1] = str(values[(0, j)][1]).split(' ')[0]

  for i in range(m - 1):
    for j in range(n - 1):
      value = values[(i, j)][2]
      if isinstance(value, float):
        value = '%.5g' % value
      table[i + 1][j + 1] = value

  table[0][0] = title

  column_sizes = [max([len(table[i][j]) for i in range(m)]) for j in range(n)]
  row_size = sum(column_sizes) + 3*(len(column_sizes) - 1) + 2

  print ''
  for i in range(m):
    print ' ' + '-'*row_size
    print '|',
    for j in range(n):
      print table[i][j] + ' '*(column_sizes[j] - len(table[i][j])),
      print '|',
    print ''
  print ' ' + '-'*row_size
  print ''

#------------------------------------------------------------------------------

def print_color(string, color):
  'Print string in given color.'
  if color == 'blue':
    print '\033[1;37;34m%s\033[0m' % string
  elif color == 'green':
    print '\033[1;37;32m%s\033[0m' % string
  elif color == 'red':
    print '\033[1;37;31m%s\033[0m' % string
  else:
    print string

#------------------------------------------------------------------------------

if __name__ == '__main__':
  sys.exit(main(sys.argv[1:]))
