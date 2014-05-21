'Problems used to test mixed solvers.'

from dolfin import Constant, Mesh, sqrt, near, Expression
import os

# Problem is specified by mesh, domains where Dirichlet bcs should be
# prescribed, Dirichlet values (this includes the inflow profile which
# is a time-dependent expression), the Reynolds number and forcing

# Directory where meshes are stored
mesh_dir = 'meshes'

# Mesh directory aware loading of mesh
mesh_path = lambda mesh: os.path.join(mesh_dir, mesh)


def on_circle(x, on_boundary, x0, y0, radius):
    'Return true for points on the surface of circle(x0, y0, radius).'
    dx = x[0] - x0
    dy = x[1] - y0
    dr = sqrt(dx**2 + dy**2)
    return on_boundary and dr < radius + 1E-3


class CylinderFlow(object):
    'Flow past a cylinder'
    name = 'cylinder'
    # Forcing
    f = Constant((0., 0.))

    # Flow past a cylinder
    mesh = Mesh(mesh_path('cylinder.xml'))

    x_min, x_max = 0, 2.2
    y_min, y_max = 0, 0.41
    circle = (0.2, 0.2, 0.05)

    @classmethod
    def inflow_boundary(cls, x, on_boundary):
        return on_boundary and near(x[0], cls.x_min)

    @classmethod
    def noslip_boundary(cls, x, on_boundary):
        return on_boundary and (near(x[1]*(cls.y_max - x[1]), 0) or
                                on_circle(x, on_boundary, *cls.circle))

    Re = Constant(1000)
    U_max = 3.5
    u_in = Expression(('4*Um*(x[1]*(ymax-x[1]))*sin(pi*t/8.0)/(ymax*ymax)',
                       '0.0'),
                      Um=U_max, ymax=y_max, t=0)

all_problems = [CylinderFlow]
