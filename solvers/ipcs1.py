__author__ = 'Miroslav Kuchta <mirok@math.uio.no>'
__date__ = '2008-03-03'
__copyright__ = 'Copyright (C) 2008-2010 ' + __author__
__license__  = 'GNU GPL version 3 or any later version'

# Modified by Miroslav Kuchta, 2014

from solverbase import *

class Solver(SolverBase):
  '''Incremental pressure correction with nonlinearity treated by Adam-Bashword
  + Crank-Nicolson.'''
  def __init__(self, options):
    SolverBase.__init__(self, options)

  def solve(self, problem):
    # Get problem parameters
    mesh = problem.mesh
    dt, t, t_range = self.get_timestep(problem)

    # Define function spaces, see if the problem has constrained domains
    if hasattr(problem, 'constrained_domain'):
      cd = problem.constrained_domain
      message = 'Problem with constrained domain'
      print '\033[1;35;35m%s\033[0m' % message
    else:
      cd = None

    V = VectorFunctionSpace(mesh, 'CG', 2, constrained_domain=cd)
    Q = FunctionSpace(mesh, 'CG', 1, constrained_domain=cd)

    # Get initial and boundary conditions
    u0, p0 = problem.initial_conditions(V, Q)
    bcs_u, bcs_p, _ = problem.boundary_conditions(V, Q, t)

    # Define test and trial functions
    u = TrialFunction(V)
    v = TestFunction(V)
    p = TrialFunction(Q)
    q = TestFunction(Q)
    
    # Define functions
    nu = Constant(problem.nu)
    k  = Constant(dt)
    f  = problem.f
    n = FacetNormal(mesh)

    u0 = interpolate(u0, V)  # velocity at previous time step
    u1 = interpolate(u0, V)  # velocity two time steps back
    u_ = Function(V)         # current velocity
    
    p0 = interpolate(p0, Q)  # previous pressure
    p_ = Function(Q)         # current pressure
    
    # Now that u0, p0 are functions, make sure that they comply with boundary
    # conditions.
    bcs = {'u' : bcs_u, 'p' : bcs_p}
    ics = {'u' : u0, 'p' : p0}
    self.apply_bcs_to_ics(bcs, ics)
    
    # Tentative velocity, solve to u1
    U = 0.5*(u + u0)
    U_ = 1.5*u0 - 0.5*u1
    
    nonlinearity = inner(dot(grad(U), U_), v)*dx 

    F0 = (1./k)*inner(u - u0, v)*dx + nonlinearity\
         + nu*inner(grad(U), grad(v))*dx + inner(grad(p0), v)*dx\
         - inner(f, v)*dx     # solve to u_
    a0, L0 = system(F0)

    # Projection
    F1 = inner(grad(p - p0), grad(q))*dx + (1./k)*q*div(u_)*dx
    a1, L1 = system(F1)                                      # solve to p_ 

    # Finalize
    F2 = (1./k)*inner(u - u_, v)*dx + inner(grad(p_ - p0), v)*dx 
    a2, L2 = system(F2)                                      # solve to p_

    # Assemble matrices
    A0 = self.assemble(a0)
    A1 = self.assemble(a1)
    A2 = self.assemble(a2)

    # Create solvers; solver02 for tentative and finalize
    #                 solver1 for projection
    solver02 = KrylovSolver('gmres', 'hypre_euclid')

    solver1 = KrylovSolver('cg', 'hypre_amg')
    # Get the nullspace if there are no pressure boundary conditions
    foo = Function(Q)     # auxiliary vector for setting pressure nullspace
    if not bcs_p:
      null_vec = Vector(foo.vector())
      Q.dofmap().set(null_vec, 1.0)
      null_vec *= 1.0/null_vec.norm('l2')
      null_space = VectorSpaceBasis([null_vec])
      solver1.set_nullspace(null_space)

    # apply global options for Krylov solvers
    options = self.options
    if 'krylov_solver_params' in options:
      for solver in [solver02, solver1]:
        self.apply_krylov_solver_options(solver, options['krylov_solver_params'])

    # Time loop
    self.start_timing()
    for t in t_range:
      # Update boundary conditions
      bcs_u, bcs_p, _ = problem.boundary_conditions(V, Q, t)

      b = assemble(L0)
      [bc.apply(A0, b) for bc in bcs_u]
      solver02.solve(A0, u_.vector(), b)

      b = assemble(L1)
      [bc.apply(A1, b) for bc in bcs_p]
      if not bcs_p:
        null_space.orthogonalize(b);
      solver1.solve(A1, p_.vector(), b)

      b = assemble(L2)
      [bc.apply(A2, b) for bc in bcs_u]
      solver02.solve(A2, u_.vector(), b)

      # Update
      self.update(problem, t, u_, p_, f)
      u1.assign(u0)
      u0.assign(u_)
      p0.assign(p_)

    return u_, p_

  def __str__(self):
    return 'IPCS1'
