'''
Solve the Taylor-Green vortex problem with Chorin method. We consider
[-1, 1]^2 domain periodic boundary conditions. Initial conditions for
pressure and velocity are given by the exact solution. For other parame-
ters such as kinematic viscosity, time step, simulation time see below. For
testing, we are the kinetic energy.
'''

from dolfin import *

U_exact = ('-(cos(pi*(x[0]))*sin(pi*(x[1]))) * exp(-2.0*nu*pi*pi*t)',
           ' (cos(pi*(x[1]))*sin(pi*(x[0]))) * exp(-2.0*nu*pi*pi*t)')

P_exact = '-0.25*(cos(2*pi*(x[0])) + cos(2*pi*(x[1]))) * exp(-4.0*nu*pi*pi*t)'

class PeriodicDomain(SubDomain):
  # Left boundary is "target domain" G
  def inside(self, x, on_boundary):
      return bool((near(x[0], -1) or near(x[1], -1)) and\
              (not ((near(x[0], -1) and near(x[1], 1)) or 
                    (near(x[0], 1) and near(x[1], -1)))) and on_boundary)

  def map(self, x, y):
    if near(x[0], 1) and near(x[1], 1):
        y[0] = x[0] - 2.
        y[1] = x[1] - 2.
    elif near(x[0], 1):
        y[0] = x[0] - 2.
        y[1] = x[1]
    else:   # near, 1)
        y[0] = x[0]
        y[1] = x[1] - 2.

constrained_domain = PeriodicDomain()

set_log_level(WARNING)

def chorin_taylor_green(n_cells):
  '''Parameter n_cells determins mesh size.'''

  mesh = RectangleMesh(-1, -1, 1, 1, n_cells, n_cells)
  h = mesh.hmin()
 
  #! problem specific
  f = Constant((0., 0.))         # force
  nu = Constant(1./100)         # kinematic viscosity
  dt = 0.2*h/1.                  # time step CFL with 1 = max. velocity
  k = Constant(dt)               # time step 
  T = 0.5                        # total simulation time

  #! solver specific
  V = VectorFunctionSpace(mesh, 'CG', 2, constrained_domain=constrained_domain)
  Q = FunctionSpace(mesh, 'CG', 1, constrained_domain=constrained_domain)

  u = TrialFunction(V)
  v = TestFunction(V)
  p = TrialFunction(Q)
  q = TestFunction(Q)

  u0 = Expression(U_exact, degree=3, t=0, nu=float(nu))
  p0 = Expression(P_exact, degree=3, t=0, nu=float(nu))

  u0 = interpolate(u0, V)  
  p0 = interpolate(p0, Q)
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
  null_vec *= 1.0/null_vec.norm('l2')
  null_space = VectorSpaceBasis([null_vec])
  solver1.set_nullspace(null_space)

  t = 0
  while t < T:
    t += dt
    
    b = assemble(L0)
    solver02.solve(A0, us.vector(), b)

    b = assemble(L1)
    null_space.orthogonalize(b);
    solver1.solve(A1, p1.vector(), b)

    b = assemble(L2)
    solver02.solve(A2, u1.vector(), b)

    u0.assign(u1)
    print t

  energy = 0.5*assemble(inner(u0, u0)*dx)

  return u0, energy

#----------------------------------------------------------------------

if __name__ == '__main__':
  u, energy = chorin_taylor_green(64)

  print energy
  
  plot(u)
  interactive()
