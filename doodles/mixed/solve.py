from elements import all_elements, make_mixed_function_space
from problems import all_problems
from dolfin import *
import os

set_log_level(WARNING)

def mixed_solve(problem, element, solver_name):
    'Mixed solver.'
    print 'Solving %s problem with %s element' % (problem.name, solver_name)

    #Make directory for results
    results_root = 'results'
    results_dir = os.path.join(results_root, solver_name, problem.name)
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    return None

    # Extract information from problem
    mesh = problem.mesh

    noslip_boundary = problem.noslip_boundary
    inflow_boundary = problem.inflow_boundary
    u_in = problem.u_in
    Re = problem.Re

    # Create function space with element
    M = make_mixed_function_space(mesh, element)

    # Noslip bc is not time-dependent and can be created now
    bc_noslip = DirichletBC(M.sub(0), Constant((0., 0.)), noslip_boundary)

    up = TrialFunction(M)
    u, p = split(up)

    vq = TestFunction(M)
    v, q = split(vq)

    # Current solution
    up0 = Function(M)
    u0, p0 = split(up0)

    # Previous solution
    up1 = Function(M)
    u1, p1 = split(up1)

    # Time step
    h = mesh.hmin()
    dt = Constant(0.005)
    k = dt**-1

    # Loads
    f0 = interpolate(f, V)
    f1 = interpolate(f, V)

    u_cn = 0.5*(u + u0)

    # Form for the first time step
    F0 = k*inner(u - u0, v)*dx + inner(dot(grad(u), u0), v)*dx +\
        Re**-1*inner(grad(u_cn), grad(v))*dx - inner(p, div(v))*dx -\
        inner(q, div(u))*dx - inner(f0, v)*dx
    a0, L0 = system(F0)

    # Form for other time steps
    u_ab = 1.5*u0 - 0.5*u1
    f_ab = 1.5*f0 - 0.5*f1
    F = k*inner(u - u0, v)*dx + inner(dot(grad(u_cn), u_ab), v)*dx +\
        Re**-1*inner(grad(u_cn), grad(v))*dx + inner(p, div(v))*dx +\
        inner(q, div(u))*dx - inner(f_ab, v)*dx
    a, L = system(F)

    t = 0
    T = 0.2
    step = 0

    while t < T:
        t += float(dt)
        step += 1
        print t

        u_in.t = t
        bc_inflow = DirichletBC(M.sub(0), u_in, inflow_boundary)
        bcs = [bc_noslip, bc_inflow]

        if step == 1:
            solve(a0 == L0, up0, bcs)
        else:
            solve(a == L, up0, bcs)
            up1.assign(up0)
            f0.assign(f1)

        if not step%10:
            pass

    # Check global mass conservation
    n = FacetNormal(mesh)
    global_div = assemble(dot(u0, n)*ds)

    # Check local divergence as Garth does it
    local_div_garth = sqrt(assemble(inner(div(u0), div(u0))*dx))

    # Check local divergence following Mikael
    v_degree = V.ufl_element().degree()
    if v_degree:
        Div_space = FunctionSpace(mesh, 'CG', v_degree-1)
    else:
        Div_space = FunctionSpace(mesh, 'DG', v_degree-1)

    div_u = project(div(u0), Div_space)
    local_div_mikael = norm(div_u, 'L2')

    print 'Global: %g, Mikael: %g, Garth: %g' %\
        (global_div, local_div_mikael, local_div_garth)

if __name__ == '__main__':
    names = all_elements.keys()
    element, solver_name = all_elements[names[0]], names[0]
    problem = all_problems[0]

    mixed_solve(problem, element, solver_name)
