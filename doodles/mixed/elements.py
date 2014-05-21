'Mixed elements for the test'

from dolfin import VectorFunctionSpace, FunctionSpace, MixedFunctionSpace


# Taylor-Hood element
th = {'V': [['CG', 2]], 'Q': [['CG', 1]]}

# Mini
mini = {'V': [['CG', 1], ['Bubble', 3]], 'Q': [['CG', 1]]}

# Crouzeix-Raviart
cr = {'V': [['CG', 2], ['Bubble', 3]], 'Q': [['DG', 1]]}


def make_function_spaces(mesh, element):
    'Create mixed function space on the mesh.'
    assert 'V' in element.keys() and 'Q' in element.keys()

    # Create velocity space
    n = len(element['Q'])
    assert n > 0
    V = VectorFunctionSpace(mesh, *element['V'][0])
    for i in range(1, n):
        V += VectorFunctionSpace(mesh, *element['V'][i])

    # Create pressure space
    assert len(element['Q']) == 1

    Q = FunctionSpace(mesh, *element['Q'][0])

    M = MixedFunctionSpace([V, Q])
    return V, Q, M 

# Put all to list for import
all_elements = {'taylor-hood': th,
                'crouzeix-raviart': cr,
                'mini': mini}
