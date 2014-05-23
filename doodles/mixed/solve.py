from elements import all_elements, make_function_spaces
from problems import all_problems
from dolfin import *
import os

set_log_level(WARNING)

def mixed_solve(problem, element):
    'Mixed solver.'
    print 'Solving %s problem with %s element' % (problem.name, element.name)

    # Make directory for results
    results_root = 'results'
    results_dir = os.path.join(results_root, element.name, problem.name)
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    # Extract information from problem
    mesh = problem.mesh
    noslip_boundary = problem.noslip
    inflow_boundary = problem.inflow
    u_in = problem.u_in
    U_max = problem.U_max
    Re = problem.Re
    f = problem.f
    T = problem.T

    # Create function space with element
    V, Q, M = make_function_spaces(mesh, element)

    # Create boundary conditions
    bc_noslip = DirichletBC(M.sub(0), Constant((0., 0.)), *noslip_boundary)
    bc_inflow = DirichletBC(M.sub(0), u_in, *inflow_boundary)
    bcs = [bc_noslip, bc_inflow]

    up = TrialFunction(M)
    u, p = split(up)

    vq = TestFunction(M)
    v, q = split(vq)

    # Solution at time level t - dt
    up0 = Function(M)
    u0, p0 = split(up0)

    #  Solution as time level t - 2dt
    up1 = Function(M)
    u1, p1 = split(up1)

    # Time step
    h = mesh.hmin()
    dt = Constant(0.25*h/U_max)
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
        Re**-1*inner(grad(u_cn), grad(v))*dx - inner(p, div(v))*dx -\
        inner(q, div(u))*dx - inner(f_ab, v)*dx
    a, L = system(F)

    # Solution at current level t
    uph = Function(M)
    uh, ph = uph.split()  # Current state components

    # Pickard loop variables
    up_ = Function(M)
    u_, p_ = up_.split()  # Previous state components

    t = 0
    step = 0

    u_out = XDMFFile(os.path.join(results_dir, 'u_%d.xdmf' % float(Re)))
    u_out.parameters['rewrite_function_mesh'] = False
    u_plot = Function(V)
    while t < T:
        t += float(dt)
        step += 1
        print '\nstep number =', step, ', time =', t

        u_in.t = t

        # Pickard loop
        iter = 0
        iter_max = 40
        tol = 1E-4
        converged = False
        e0 = None
        E = None
        while not converged and iter < iter_max:
            iter += 1
            print '\titer number =', iter

            # Remeber the old solution
            up_.assign(uph)

            # Get new solution with relaxation
            if step == 1:
                solve(a0 == L0, uph, bcs,
                      solver_parameters={'linear_solver': 'mumps'})
                uph.vector()[:] = 0.5*(uph.vector()[:] + up_.vector()[:])
                up0.assign(uph)
            else:
                solve(a == L, uph, bcs,
                      solver_parameters={'linear_solver': 'mumps'})
                uph.vector()[:] = 0.5*(uph.vector()[:] + up_.vector()[:])
                up1.assign(up0)
                up0.assign(uph)

            # Compute the error and decide convergence
            e = norm(uh, 'l2')
            if iter < 2:
                e0 = e
                converged = False
            else:
                E = abs(e - e0)/(e + e0)
                converged = E < tol
                e0 = e

            print '\t error =', E

        # TODO forces!!
        f0.assign(f1)

        if not step % 10:
            u_plot.assign(up0.split(True)[0])
            plot(u_plot, title='%s @ %g' % (element.name, t))
            u_out << u_plot, t

    # Check global mass conservation
    n = FacetNormal(mesh)
    global_div = assemble(dot(u0, n)*ds)

    # Check local divergence as Garth does it
    local_div_g = sqrt(assemble(inner(div(u0), div(u0))*dx))

    # Check local divergence following Mikael
    v_degree = V.ufl_element().degree()
    if v_degree-1:
        Div_space = FunctionSpace(mesh, 'CG', v_degree-1)
    # Capture projecting to degree 0
    else:
        Div_space = FunctionSpace(mesh, 'DG', v_degree-1)

    div_u = project(div(u0), Div_space)
    local_div_m = norm(div_u, 'L2')

    return '%s %s Global: %g, Mikael: %g, Garth: %g' %\
        (problem.name, element.name, global_div, local_div_m, local_div_g)

if __name__ == '__main__':

    for element in all_elements:
        for problem in all_problems:
            result = mixed_solve(problem, element)
            with open('data.txt', 'a') as f:
                f.write('%s\n' % result)
