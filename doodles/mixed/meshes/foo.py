from dolfin import Mesh, MeshFunction
import sys

if __name__  == '__main__':
    mesh = Mesh(sys.argv[1])
    f_f = MeshFunction('size_t', mesh, sys.argv[2])
