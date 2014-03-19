__author__ = "Anders Logg <logg@simula.no> and Kent-Andre Mardal <kent-and@simula.no>"
__date__ = "2008-03-03"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

class Solver(SolverBase):
  "Original pressure-correction scheme by Chorin and Temam."
  def __init__(self, options):
    SolverBase.__init__(self, options)

  def solve(self, problem):
    # Get problem parameters
    mesh = problem.mesh
    dt, t, t_range = self.get_timestep(problem)

    # Define function spaces
    V = VectorFunctionSpace(mesh, "CG", 2)
    Q = FunctionSpace(mesh, "CG", 1)
    DG = FunctionSpace(mesh, "DG", 0)

    # Get initial and boundary conditions
    u0, p0 = problem.initial_conditions(V, Q)
    bcs_u, bcs_p, _ = problem.boundary_conditions(V, Q, t)

    # Define test and trial functions
    v = TestFunction(V)
    q = TestFunction(Q)
    u = TrialFunction(V)
    p = TrialFunction(Q)

    # Define functions
    nu = Constant(problem.nu)
    k  = Constant(dt)
    f  = problem.f
    u0 = interpolate(u0, V)  
    p0 = interpolate(p0, Q)
    us = Function(V)
    u1 = Function(V)
    p1 = Function(Q)

    # Tentative velocity
    F0 = (1./k)*inner(u - u0, v)*dx + inner(dot(grad(u0), u0), v)*dx\
         + nu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
    a0, L0 = system(F0)

    # Projection
    F1 = inner(grad(p), grad(q))*dx + (1./k)*q*div(us)*dx
    a1, L1 = system(F1)

    # Finalize
    F2 = (1./k)*inner(u - us, v)*dx + inner(grad(p1), v)*dx 
    a2, L2 = system(F2)

    # Assemble matrices
    A0 = assemble(a0)
    A1 = assemble(a1)
    A2 = assemble(a2)

    # Create solvers; solver02 for tentative and finalize
    #                 solver1 for projection
    solver02 = KrylovSolver('gmres', 'ilu')

    solver1 = KrylovSolver('cg', 'petsc_amg')
    # Get the nullspace if there are no pressure boundary conditions
    if not bcs_p:
      null_vec = Vector(p0.vector())
      Q.dofmap().set(null_vec, 1.0)
      null_vec *= 1.0/null_vec.norm("l2")
      null_space = VectorSpaceBasis([null_vec])
      solver1.set_nullspace(null_space)

    # apply global options for Krylov solvers
    options = self.options
    for solver in [solver02, solver1]:
      solver.parameters['absolute_tolerance'] = \
                options['krylov_solver_absolute_tolerance']
      solver.parameters['relative_tolerance'] = \
                options['krylov_solver_relative_tolerance']
      solver.parameters['monitor_convergence'] = \
                options['krylov_solver_monitor_convergence']

    # Time loop
    self.start_timing()
    for t in t_range:
      # Update boundary conditions
      bcs_u, bcs_p, _ = problem.boundary_conditions(V, Q, t)

      b = assemble(L0)
      [bc.apply(A0, b) for bc in bcs_u]
      solver02.solve(A0, us.vector(), b)

      b = assemble(L1)
      [bc.apply(A1, b) for bc in bcs_p]
      if not bcs_p:
        null_space.orthogonalize(b);
      solver1.solve(A1, p1.vector(), b)

      b = assemble(L2)
      [bc.apply(A2, b) for bc in bcs_u]
      solver02.solve(A2, u1.vector(), b)

      # Update
      self.update(problem, t, u1, p1)
      u0.assign(u1)

    return u1, p1

  def __str__(self):
    return "Chorin"
