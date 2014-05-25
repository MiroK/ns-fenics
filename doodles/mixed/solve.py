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

    # Time step
    h = mesh.hmin()
    dt = Constant(0.25*h/U_max)
    k = dt**-1

    # Loads
    f0 = interpolate(f, M)
    F0 = f0.vector()

    # Forms for assembling matrices on lhs
    w = inner(u, v)*dx                     # mass matrix form
    n0 = inner(dot(grad(u), u0), v)*dx     # nonlinearity
    s = Constant(0.5)*Re**-1*inner(grad(u), grad(v))*dx # ~stiffness matrix form
    bt = -inner(p, div(v))*dx              # divergence form
    b = -inner(q, div(u))*dx               # gradient form
    A = assemble(k*w + s + bt + b)         # part of lhs matrix which stays constant

    # Forms for assembling matrices/vector on rhs
    L = inner(Constant((0., 0.,)), v)*dx   # place holder
    W = assemble(w)                        # mass matrix
    K = assemble(k*w - s) # part of rhs matrix which stays constant
    b_ = assemble(L)                       # auxiliary vector
    b = Vector(b_)                         # rhs vector

    # Solution at current level t
    uph = Function(M)
    uh, ph = uph.split()  # Current state components
    UPH = uph.vector()

    # Pickard loop variables
    up_ = Function(M)
    u_, p_ = up_.split()  # Previous state components

    # Prepare output
    u_out = XDMFFile(os.path.join(results_dir, 'u_%d.xdmf' % float(Re)))
    u_out.parameters['rewrite_function_mesh'] = False
    u_plot = Function(V)

    # Create solver
    solver = LUSolver('mumps')
    solver.parameters['same_nonzero_pattern'] = True  # FIXME, advantage?

    t = 0
    step = 0
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
            UP0 = up0.vector()

            iter += 1
            print '\titer number =', iter

            # Remeber the old solution
            up_.assign(uph)
            UP_ = up_.vector()

            # Assemble the t-dep part of lhs matrix and
            # add to A to it yielding comple lhs matrix N0
            N0 = assemble(n0)
            N0.axpy(1, A, False)

            # Put together the rhs vector b
            W.mult(F0, b_)
            K.mult(UP0, b)
            b.axpy(1, b_)

            # Apply boundary conditions
            [bc.apply(N0, b) for bc in bcs]

            # Get new solution with relaxation
            solve(N0, UPH, b)
            UPH *= 0.5
            UPH.axpy(0.5, UP_)
            print uph.vector().norm('l2')

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

        # TODO handle transient forces!!

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

    print 'Problems:', [(i, problem) for i, problem in enumerate(all_problems)]
    print 'Elements', [(i, element) for i, element in enumerate(all_elements)]

    i_problem = int(raw_input('Select problem: '))
    i_element = int(raw_input('Select element: '))

    # import sys
    # assert len(sys.argv) == 3
    # i_problem = sys.argv[1]
    # i_element = sys.argv[2]

    if (-1 < i_problem < len(all_problems)) and\
            (-1 < i_element < len(all_elements)):
        problem = all_problems[i_problem]
        element = all_elements[i_element]
        result = mixed_solve(problem, element)
        with open('data.txt', 'a') as f:
            f.write('%s\n' % result)
    else:
        exit()
