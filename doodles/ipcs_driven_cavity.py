'''
Solve the driven cavity problem with IPCS method. We consider [0, 1]^2 domain
with no slip boundary conditions everywhere apart from y = 1 where the velocity
is (1, 0). Initial conditions for pressure and velocity are 0. For other pareme-
ters such as kinematic viscosity, time step, simulation time see below. For
testing, we are interested in the stream function at each time instance. Time
discretization is only first order.
'''

from dolfin import *

set_log_level(WARNING)

def stream_function(u):
  '''Compute stream function of given 2-d velocity vector.'''
  V = u.function_space().sub(0).collapse()

  if V.mesh().topology().dim() != 2:
    raise ValueError("Only stream function in 2D can be computed.")

  psi = TrialFunction(V)
  phi = TestFunction(V)

  a = inner(grad(psi), grad(phi))*dx
  L = inner(u[1].dx(0) - u[0].dx(1), phi)*dx
  bc = DirichletBC(V, Constant(0.), DomainBoundary())

  A, b = assemble_system(a, L, bc)
  psi = Function(V)
  solve(A, psi.vector(), b)

  return psi

#---------------------------------------------------------------------

def ipcs_driven_cavity(n_cells):
  '''Parameter n_cells determins mesh size.'''

  mesh = UnitSquareMesh(n_cells, n_cells)
  h = mesh.hmin()
  n = FacetNormal(mesh)
 
  #! problem specific
  f = Constant((0., 0.))         # force
  nu = Constant(1./1000)         # kinematic viscosity
  dt = 0.2*h/1.                  # time step CFL with 1 = max. velocity
  k = Constant(dt)               # time step 
  T = 2.5                        # total simulation time
  u0 = Constant((0., 0.))        # initial velocity
  p0 = Constant(0.)              # initial pressure

  no_slip_value = Constant((0., 0.)) 
  
  def no_slip_domain(x, on_boundary):
    return on_boundary and (near(x[1], 0) or near(x[0]*(1 - x[0]), 0))

  driven_value = Constant((1., 0.))

  def driven_domain(x, on_boundary):
    return on_boundary and near(x[1], 1)

  #! solver specific
  V = VectorFunctionSpace(mesh, 'CG', 2)
  Q = FunctionSpace(mesh, 'CG', 1)

  u = TrialFunction(V)
  v = TestFunction(V)
  p = TrialFunction(Q)
  q = TestFunction(Q)

  u0 = interpolate(u0, V)  
  p0 = interpolate(p0, Q)
  u1 = Function(V)
  p1 = Function(Q)

  # tentative velocity
  U = 0.5*(u + u0)
  F0 = (1./k)*inner(u - u0, v)*dx + inner(dot(grad(u0), u0), v)*dx\
       + inner(grad(p0), v)*dx + nu*inner(grad(U), grad(v))*dx \
       - inner(f, v)*dx
  a0, L0 = system(F0)

  # projection
  F1 = inner(grad(p - p0), grad(q))*dx + (1./k)*q*div(u1)*dx
  a1, L1 = system(F1)

  # finalize
  F2 = (1./k)*inner(u - u1, v)*dx + inner(grad(p1 - p0), v)*dx 
  a2, L2 = system(F2)

  # boundary conditions
  bc_no_slip = DirichletBC(V, no_slip_value, no_slip_domain)
  bc_driven = DirichletBC(V, driven_value, driven_domain)
  bcs = [bc_no_slip, bc_driven]

  A0 = assemble(a0)
  A1 = assemble(a1)
  A2 = assemble(a2)
 
  solver02 = KrylovSolver('gmres', 'ilu')

  solver1 = KrylovSolver('cg', 'petsc_amg')
  null_vec = Vector(p0.vector())
  Q.dofmap().set(null_vec, 1.0)
  null_vec *= 1.0/null_vec.norm('l2')
  null_space = VectorSpaceBasis([null_vec])
  solver1.set_nullspace(null_space)

  t = 0
  while t < T:
    t += dt
    
    b = assemble(L0)
    [bc.apply(A0, b) for bc in bcs]
    solver02.solve(A0, u1.vector(), b)

    b = assemble(L1)
    null_space.orthogonalize(b);
    solver1.solve(A1, p1.vector(), b)

    b = assemble(L2)
    [bc.apply(A2, b) for bc in bcs]
    solver02.solve(A2, u1.vector(), b)

    u0.assign(u1)
    p0.assign(p1)
    print t

  phi0 = stream_function(u0)
  return u0, phi0

if __name__ == '__main__':
  u, phi = ipcs_driven_cavity(32)

  print phi.vector().min()
  
  plot(u)
  plot(phi)
  interactive()
