'''
Solve the flow driven by pressure gradient with Chorin method. We consider
[0, 1]^2 domain with no slip boundary conditions at y = 1, y = 0 and neumann
+ set pressure at x = 0, and x = 1. Initial conditions for pressure and velocity
are 0 and 1-x respectively. For other paremeters such as kinematic viscosity,
time step, simulation time see below. For testing, we are interested in the
value of velocity at (x=1, y=0.5)
'''

from dolfin import *

set_log_level(WARNING)

class InflowBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[0], 0) and on_boundary

class OutflowBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[0], 1) and on_boudnary

class NoslipBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[1]*(1 - x[1]), 0) and on_boundary

#---------------------------------------------------------------------

def chorin_channel(n_cells):
  '''Parameter n_cells determins mesh size.'''

  mesh = UnitSquareMesh(n_cells, n_cells)
  h = mesh.hmin()
 
  #! problem specific
  f = Constant((0., 0.))         # force
  nu = Constant(1./8.)           # kinematic viscosity
  dt = 0.2*h/1.                  # time step CFL with 1 = max. velocity
  k = Constant(dt)               # time step 
  T = 0.5                        # total simulation time
  u0 = Constant((0., 0.))        # initial velocity
  p0 = Expression('1 - x[0]')    # initial pressure

  #! solver specific
  V = VectorFunctionSpace(mesh, 'CG', 2)
  Q = FunctionSpace(mesh, 'CG', 1)

  u = TrialFunction(V)
  v = TestFunction(V)
  p = TrialFunction(Q)
  q = TestFunction(Q)

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

  # boundary conditions
  b_v = DirichletBC(V, Constant((0.0, 0.0)), NoslipBoundary())
  b_p0 = DirichletBC(Q, Constant(1.), InflowBoundary())
  b_p1 = DirichletBC(Q, Constant(0.),  OutflowBoundary())
  bcs_v = [b_v]
  bcs_p = [b_p0, b_p1]

  A0 = assemble(a0)
  A1 = assemble(a1)
  A2 = assemble(a2)
 
  solver02 = KrylovSolver('gmres', 'ilu')
  solver1 = KrylovSolver('cg', 'petsc_amg')

  t = 0
  while t < T:
    t += dt
    
    b = assemble(L0)
    [bc.apply(A0, b) for bc in bcs_v]
    solver02.solve(A0, us.vector(), b)

    b = assemble(L1)
    [bc.apply(A1, b) for bc in bcs_p]
    solver1.solve(A1, p1.vector(), b)

    b = assemble(L2)
    [bc.apply(A2, b) for bc in bcs_v]
    solver02.solve(A2, u1.vector(), b)

    u0.assign(u1)
    print t

  value = u0(1, 0.5)
  return u0, value

#----------------------------------------------------------------------

if __name__ == '__main__':
  u, value = chorin_channel(64)

  print value
  
  plot(u)
  interactive()
