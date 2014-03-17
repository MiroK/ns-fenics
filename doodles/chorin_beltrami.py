'''
Solve the beltrami problem with Chorin method. We consider [-1, 1]^3 domain
with velocity boundary conditions everywhere. Initial conditions for pressure
and velocity are given by exact solution. For other parameters such as kinematic
viscosity, time step, simulation time see below. For testing, we are interested
in the relative error of velocity in L2 norm.
'''

from dolfin import *

set_log_level(WARNING)

U_exact = \
  ('-((a*(pow(E,a*x[2])*cos(a*x[0] + d*x[1]) + pow(E,a*x[0])*sin(a*x[1] + d*x[2])))/pow(E,pow(d,2)*t*etabyrho))',
   '-((a*(pow(E,a*x[0])*cos(a*x[1] + d*x[2]) + pow(E,a*x[1])*sin(d*x[0] + a*x[2])))/pow(E,pow(d,2)*t*etabyrho))',
   '-((a*(pow(E,a*x[1])*cos(d*x[0] + a*x[2]) + pow(E,a*x[2])*sin(a*x[0] + d*x[1])))/pow(E,pow(d,2)*t*etabyrho))')

P_exact = \
  ('-(rho/2.0)*(pow(a,2)*(pow(E,2*a*x[0]) + pow(E,2*a*x[1]) + pow(E,2*a*x[2]) + 2*pow(E,a*(x[1] + x[2]))*cos(d*x[0] + a*x[2])*sin(a*x[0] + d*x[1]) + 2*pow(E,a*(x[0] + x[1]))*cos(a*x[1] + d*x[2])*sin(d*x[0] + a*x[2]) + 2*pow(E,a*(x[0] + x[2]))*cos(a*x[0] + d*x[1])*sin(a*x[1] + d*x[2])))/(pow(E,pow(d,2)*t*etabyrho))')

u_params = {'a': pi/4.0, 'd': pi/2.0, 'E': e,             'etabyrho': 1.0, 't': 0.0}
p_params = {'a': pi/4.0, 'd': pi/2.0, 'E': e, 'rho': 1.0, 'etabyrho': 1.0, 't': 0.0}

#---------------------------------------------------------------------

def chorin_beltrami(n_cells):
  '''Parameter n_cells determins mesh size.'''

  mesh = BoxMesh(-1, -1, -1, 1, 1, 1, n_cells, n_cells, n_cells)
  h = mesh.hmin()
 
  #! problem specific
  f = Constant((0., 0., 0))      # force
  nu = Constant(1.)              # kinematic viscosity
  dt = 0.2*h/1.                  # time step CFL with 1 = max. velocity
  k = Constant(dt)               # time step 
  T = 0.5                        # total simulation time

  #! solver specific
  V = VectorFunctionSpace(mesh, 'CG', 2)
  Q = FunctionSpace(mesh, 'CG', 1)

  u = TrialFunction(V)
  v = TestFunction(V)
  p = TrialFunction(Q)
  q = TestFunction(Q)

  u_exact = Expression(U_exact, degree=3, **u_params)
  p_exact = Expression(P_exact, degree=3, **p_params)
  
  u0 = interpolate(u_exact, V)  
  p0 = interpolate(p_exact, Q)
  us = Function(V)
  u1 = Function(V)
  p1 = Function(Q)

  # tentative velocity
  F0 = (1./k)*inner(u - u0, v)*dx + inner(dot(grad(u0), u0), v)*dx\
       + nu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
  a0, L0 = system(F0)

  # projection
  F1 = inner(grad(p), grad(q))*dx + (1./k)*q*div(us)*dx
  a1, L1 = system(F1)

  # finalize
  F2 = (1./k)*inner(u - us, v)*dx + inner(grad(p1), v)*dx 
  a2, L2 = system(F2)

  A0 = assemble(a0)
  A1 = assemble(a1)
  A2 = assemble(a2)
 
  solver02 = KrylovSolver('gmres', 'ilu')

  solver1 = KrylovSolver('cg', 'petsc_amg')
  null_vec = Vector(p0.vector())
  Q.dofmap().set(null_vec, 1.0)
  null_vec *= 1.0/null_vec.norm("l2")
  null_space = VectorSpaceBasis([null_vec])
  solver1.set_nullspace(null_space)

  t = 0
  while t < T:
    t += dt
    
    # boundary conditions are timedependent
    u_exact.t = t
    bcs = DirichletBC(V, u_exact, DomainBoundary())
    bcs = [bcs]
    
    b = assemble(L0)
    [bc.apply(A0, b) for bc in bcs]
    solver02.solve(A0, us.vector(), b)

    b = assemble(L1)
    null_space.orthogonalize(b);
    solver1.solve(A1, p1.vector(), b)

    b = assemble(L2)
    [bc.apply(A2, b) for bc in bcs]
    solver02.solve(A2, u1.vector(), b)

    u0.assign(u1)
    print t

  e_norm = sqrt(assemble(inner(u_exact - u0, u_exact - u0)*dx, mesh=mesh))
  u_norm = sqrt(assemble(inner(u_exact, u_exact)*dx, mesh=mesh))

  return u0, interpolate(u_exact, V), error

if __name__ == '__main__':
  u, u_exact, error = chorin_beltrami(8)

  print error
  
  plot(u)
  plot(u_exact)
  interactive()
