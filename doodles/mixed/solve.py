from elements import all_elements, make_function_spaces
from problems import all_problems
import mpi4py.MPI as mpi
from dolfin import *
import os

comm = mpi.COMM_WORLD
rank = comm.Get_rank()


def p_print(arg):
    'Print only on process 0.'
    if rank == 0:
        print arg

set_log_level(WARNING)


def mixed_solve(problem, element):
    'Mixed solver.'
    p_print('Solving %s problem with %s element' % (problem.name, element.name))
    timer = Timer('Mixed solver')
    timer.start()

    # Make directory for results
    results_root = 'results/foo'
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
    UP0 = up0.vector()
    u0, p0 = split(up0)

    # Solution as time level t - 2dt
    up1 = Function(M)
    UP1 = up1.vector()
    u1, p1 = split(up1)

    # Time step
    h = mesh.hmin()
    h = MPI.min(h)
    dt = Constant(0.25*h/U_max)
    k = dt**-1

    # Loads
    f0 = interpolate(f, V)
    f1 = interpolate(f, V)

    # Forms for assembling lhs
    u_ab = 1.5*u0 - 0.5*u1

    w = inner(u, v)*dx               # mass matrix form
    n0 = inner(dot(grad(u), u0), v)*dx  # nonlinearity in 1st step
    n1 = Constant(0.5)*inner(dot(grad(u), u_ab), v)*dx  # nonlin. in 2nd step
    s = Constant(0.5)*Re**-1*inner(grad(u), grad(v))*dx  # ~stiffness m. form
    bt = -inner(p, div(v))*dx        # divergence form
    b = -inner(q, div(u))*dx         # gradient form
    A = assemble(k*w + s + bt + b)   # part of lhs matrix which stays constant

    # Forms for assembling rhs
    f_ab = 1.5*f0 - 0.5*f1

    L0 = inner(f0, v)*dx
    L1 = inner(f_ab, v)*dx
    K = assemble(k*w - s)

    # Aux form to get consisten b0, b1 vectors
    L = inner(Constant((0, 0)), v)*dx
    b0 = assemble(L)
    b1 = assemble(L)

    # Solution at current level t
    uph = Function(M)
    UPH = uph.vector()
    uh, ph = uph.split()  # Current state components

    # Pickard loop variables
    up_ = Function(M)
    UP_ = up_.vector()
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
        p_print('\nstep number = %d, time = %g' % (step, t))

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
            p_print('\titer number = %d' % iter)

            # Remeber the old solution
            UP_.zero()
            UP_.axpy(1, UPH)

            if step == 1:
                # Assemble the t-dep part of lhs matrix and
                # add to A to it yielding comple lhs matrix N0
                N0 = assemble(n0)
                N0.axpy(1, A, False)

                # Put together the rhs vector b
                assemble(L0, tensor=b0) # b0 =inner(f0, v)*dx
                K.mult(UP0, b1)         # b1 = k*inner(u0, v) - 0.5*inner(grad(u0), grad(v))*dx
                b0.axpy(1, b1)

                # Apply boundary conditions
                [bc.apply(N0, b0) for bc in bcs]

                # Get new solution with relaxation
                solve(N0, UPH, b0)
                UPH *= 0.5
                UPH.axpy(0.5, UP_)

                # Assign UPH to UP0
                UP0.zero()
                UP0.axpy(1, UPH)
            else:
                # Assemble the t-dep part of lhs matrix and
                # add to A to it yielding comple lhs matrix N0
                N1 = assemble(n1)
                N1.mult(UP0, b1)        # b1 = 0.5*inner(dot(grad(u0), u_ab), v)*dx

                N1.axpy(1, A, False)

                # Put together the rhs vector b
                assemble(L1, tensor=b0) # b0 = inner(f_ab, v)*dx
                b0.axpy(-1, b1)         # b0 = b0 - b1
                K.mult(UP0, b1)  # b1 = k*inner(u0, v) - 0.5*inner(grad(u0), grad(v))*dx
                b1.axpy(1, b0)

                # Apply boundary conditions
                [bc.apply(N1, b1) for bc in bcs]

                # Get new solution with relaxation
                solve(N1, UPH, b1)
                UPH *= 0.5
                UPH.axpy(0.5, UP_)

                # Assign UP0 to UP1, UPH to UP0
                UP1.zero()
                UP1.axpy(1, UP0)

                UP0.zero()
                UP0.axpy(1, UPH)

            # Compute the error and decide convergence
            e = norm(uh, 'l2')
            E = -1
            if iter < 2:
                e0 = e
                converged = False
            else:
                E = abs(e - e0)/(e + e0)
                converged = E < tol
                e0 = e

            p_print('\t error = %g' % E)

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

    # Collect timer info
    time = timer.stop()
    time = MPI.sum(time)

    return '%s %s Global: %g, Mikael: %g, Garth: %g, time: %g, proc.: %d' %\
        (problem.name, element.name,
         global_div, local_div_m, local_div_g,
         time, comm.Get_size())

if __name__ == '__main__':

    p_print('Problems:')
    p_print([(i, problem) for i, problem in enumerate(all_problems)])
    p_print('Elements')
    p_print([(i, element) for i, element in enumerate(all_elements)])

    if rank == 0:
        i_problem = int(raw_input('Select problem: '))
        i_element = int(raw_input('Select element: '))
    else:
        i_problem, i_element = None, None

    i_problem = comm.bcast(i_problem, root=0)
    i_element = comm.bcast(i_element, root=0)

    if (-1 < i_problem < len(all_problems)) and\
            (-1 < i_element < len(all_elements)):
        problem = all_problems[i_problem]
        element = all_elements[i_element]
        result = mixed_solve(problem, element)

        if rank == 0:
            with open('data.txt', 'a') as f:
                f.write('%s\n' % result)
    else:
        exit()
