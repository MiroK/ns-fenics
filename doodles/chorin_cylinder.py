'''
Solve the flow around the cylinder with Chorin method. We consider domain
[0, 2.2]x[0, 0.41] with cyliner or radius 0.05 centered at [0.2, 0.2].No slip
boundary conditions at y = 0, y = 0.41. Velocity profile is presscribed at
inlet while zero 0 is mainained at outlet. Initial conditions for pressure
and velocity are 0. For other paremeters such as kinematic viscosity,
time step, simulation time see below. For testing, we are interested in the
pressure difference between front and back of the cylinder.
'''

from dolfin import *

def refine_cylinder(mesh, c_x, c_y, r):
  h = mesh.hmin()
  center = Point(c_x, c_y)
  cell_f = CellFunction('bool', mesh, False)
  for cell in cells(mesh):
    if cell.midpoint().distance(center) < r + h:
      cell_f[cell] = True
  mesh = refine(mesh, cell_f)

  return mesh

#------------------------------------------------------------------------------

set_log_level(WARNING)

#---------------------------------------------------------------------

def chorin_cylinder(resolution):
  '''Parameter resolution determins mesh size.'''
  # geometric params and mesh creation
  x_min, x_max = 0, 2.2
  y_min, y_max = 0, 0.41
  c_x, c_y, r = 0.2, 0.2, 0.05

  rect = Rectangle(x_min, y_min, x_max, y_max)
  circ = Circle(c_x, c_y, r)
  domain = rect - circ

  mesh = Mesh(domain, resolution)

  for i in range(3):
    mesh = refine_cylinder(mesh, c_x, c_y, r)

  h = mesh.hmin()

  # define domains
  class InflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
      return on_boundary and near(x[0], x_min)

  class OutflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
      return on_boundary and near(x[0], x_max)

  class NoslipBoundary(SubDomain):
    def inside(self, x, on_boundary):
      dx = x[0] - c_x
      dy = x[1] - c_y
      dr = sqrt(dx**2 + dy**2)
      return on_boundary and (near(x[1]*(y_max - x[1]), 0) or dr < r + 1E-3)

  #! problem specific
  f = Constant((0., 0.))  # force
  nu = Constant(1./1000)  # kinematic viscosity
  dt = 0.2*h/3.5          # time step CFL with 3.5 = max. velocity
  k = Constant(dt)        # time step 
  T = 8.0                 # total simulation time
  u0 = Constant((0., 0.)) # initial velocity
  p0 = Constant(0.)       # initial pressure

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
  
  # Expression for inlet bc
  g0 = Expression(('4*Um*(x[1]*(ymax-x[1]))*sin(pi*t/8.0)/(ymax*ymax)', '0.0'),
                           Um=1.5, ymax=y_max, t=0)
  
  bc_noslip = DirichletBC(V, Constant((0., 0.)), NoslipBoundary())
  bc_p = DirichletBC(Q, Constant(0.), OutflowBoundary())

  A0 = assemble(a0)
  A1 = assemble(a1)
  A2 = assemble(a2)
 
  solver02 = KrylovSolver('gmres', 'ilu')
  solver1 = KrylovSolver('cg', 'petsc_amg')
  
  t = 0
  while t < T:
    t += dt
    # boundary condition update
    g0.t = t
    bc_inlet = DirichletBC(V, g0, InflowBoundary())

    bcs_v = [bc_noslip, bc_inlet]
    bcs_p = [bc_p]
    
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
  
  front = Point(c_x - r - DOLFIN_EPS, c_y)
  back = Point(c_x + r + DOLFIN_EPS, c_y)
  value = p1(front) - p1(back)
  return p1, value

#----------------------------------------------------------------------

if __name__ == '__main__':
  p, value = chorin_cylinder(50)

  print value
  
  plot(p)
  interactive()
