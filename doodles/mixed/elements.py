'Mixed elements for the test'

from dolfin import VectorFunctionSpace, FunctionSpace, MixedFunctionSpace


class TaylorHood(object):
    'Taylor-Hood element'
    name = 'taylor-hood'
    V = [['CG', 2]]
    Q = [['CG', 1]]


class Mini(object):
    'Mini element'
    name = 'mini'
    V = [['CG', 1], ['Bubble', 3]]
    Q = [['CG', 1]]


class CrouzeixRaviart(object):
    'Crouzeix-Raviart element'
    name = 'crouzeix-raviart'
    V = [['CG', 2], ['Bubble', 3]]
    Q = [['DG', 1]]


def make_function_spaces(mesh, element):
    'Create mixed function space on the mesh.'
    assert hasattr(element, 'V') and hasattr(element, 'Q')

    # Create velocity space
    n = len(element.V)
    assert n > 0
    V = VectorFunctionSpace(mesh, *element.V[0])
    for i in range(1, n):
        V += VectorFunctionSpace(mesh, *element.V[i])

    # Create pressure space
    assert len(element.Q) == 1
    Q = FunctionSpace(mesh, *element.Q[0])

    M = MixedFunctionSpace([V, Q])
    return V, Q, M

# Put all to list for import
all_elements = [TaylorHood, CrouzeixRaviart, Mini]
