'''
Solve the driven cavity problem with GRPC method. We consider [0, 1]^2 domain
with no slip boundary conditions everywhere apart from y = 1 where the velocity
is (1, 0). Initial conditions for pressure and velocity are 0. For other pareme-
ters such as kinematic viscosity, time step, simulation time see below. For
testing, we are interested in the stream function at each time instance.
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

#------------------------------------------------------------------------------

def sigma(u, p, nu):
  return 2*nu*sym(grad(u)) - p*Identity(u.cell().d)

#------------------------------------------------------------------------------

def grpc_driven_cavity(n_cells):
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
  u1 = interpolate(u0, V)
  p1 = interpolate(p0, Q)

  solver0 = KrylovSolver('gmres', 'ilu')        # solver for M+dt*N inverse
  solver1 = KrylovSolver('cg', 'petsc_amg')     # solver for Schur inverse
  null_vec = Vector(p0.vector())
  Q.dofmap().set(null_vec, 1.0)
  null_vec *= 1.0/null_vec.norm("l2")
  null_space = VectorSpaceBasis([null_vec])
  solver1.set_nullspace(null_space)

  maxiter = 100 # inner loop max iteration
  tau_y1 = 2.0  # params in pressure update
  tau_y2 = 2.0

  # Boundary conditions
  bc_no_slip = DirichletBC(V, no_slip_value, no_slip_domain)
  bc_driven = DirichletBC(V, driven_value, driven_domain)
  bcs = [bc_no_slip, bc_driven]

  # Velocity and pressure residuals
  U = 0.5*(u0 + u1)
  P = p1
 
  # residual of first equation of system
  Ru = inner(v, u1 - u0)*dx + k*inner(v, (dot(grad(U), U)))*dx \
     + k*inner(sym(grad(v)), sigma(U, P, nu))*dx \
     - k*nu*inner(v, dot(n, grad(U)))*ds - k*inner(v, f)*dx
  
  # residual of the second eq.
  Rp = k*q*div(U)*dx

  # Assemble preconditioners
  ax  = inner(v, u)*dx + 0.5*k*inner(v, (dot(grad(u), u0)))*dx \
      + 0.5*k*2*nu*inner(sym(grad(v)), sym(grad(u)))*dx \
      - 0.5*k*nu*inner(v, dot(n, grad(u)))*ds

  ay1 = k**2*(inner(grad(q), grad(p)))*dx
  ay2 = k**2*((1.0/(nu*k)) * q*p)*dx

  Kx  = assemble(ax)
  Ky1 = assemble(ay1)
  Ky2 = assemble(ay2)
  [bc.apply(Kx) for bc in bcs] # there are no bcs for pressure

  # Get solution vectors
  x = u1.vector()
  y = p1.vector()
  delta_x = Vector(x.size())
  delta_y = Vector(y.size())

  # Time loop
  t = 0
  while t < T:
    t += dt
      
    # GRPC iteration
    for iter in range(maxiter):
      # Velocity update
      rx = assemble(Ru)
      [bc.apply(rx, x) for bc in bcs]
      delta_x.zero()
      solver0.solve(Kx, delta_x, rx)  #  dx = A-1*U
      x.axpy(-1.0, delta_x)           # U = U - dx

      # Pressure update
      ry = assemble(Rp)
      
      delta_y.zero()
      solver1.solve(Ky1, delta_y, ry)
      y.axpy(-tau_y1, delta_y)

      delta_y.zero()
      solve(Ky2, delta_y, ry, 'cg', 'jacobi')
      y.axpy(-tau_y2, delta_y)

      # Check for convergence
      r = sqrt(norm(rx)**2 + norm(ry)**2)
      if r < 1E-6:
        break

    print t
    u0.assign(u1)

  phi0 = stream_function(u0)
  return u0, phi0

#------------------------------------------------------------------------------

if __name__ == '__main__':
  u, phi = grpc_driven_cavity(64)

  print phi.vector().min()
  
  plot(u)
  plot(phi)
  interactive()
