'''Convert xml to xdmf meshes.'''

from dolfin import Mesh, MeshFunction, File
import os

def is_mesh_file(file):
    name, ext = os.path.splitext(file)
    if 'facet_region' not in name or\



for f in os.listdir('.'):
    name, ext = os.path.splitext(f)
    if ext == '.xml':
        new_name = ''.join([name, '.xdmf'])

        # Convert facet function
        if 'facet' in name:
            print 'Converting {0} to {1}'.format(f, new_name)

            # Look for the mesh file to create mesh
            i = f.find('facet')
            mesh_file_name = ''.join([name[:(i-1)], '.xml'])
            mesh = Mesh(mesh_file_name)
            File(new_name) << MeshFunction('size_t', mesh, f)

            # Test
            foo = MeshFunction('size_t', mesh, new_name)
        # Convert mesh
        elif 'physical' not in name:
            print 'Converting {0} to {1}'.format(f, new_name)

            File(new_name) << mesh

            # Test
            bar = Mesh(new_name)
