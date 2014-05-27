'''Convert xml to xdmf meshes.'''

from dolfin import Mesh, MeshFunction, File
import sys
import os


def is_xml_mesh_file(file):
    'Check if file is .xml that could be mesh.'
    name, ext = os.path.splitext(file)
    if ext == '.xml':
        if not('facet_region' in name or 'physical_region' in name):
            return True
    return False


def get_facet_xml(mesh_xml):
    'Get a name of facet region file of given xml mehs file.'
    assert is_xml_mesh_file(mesh_xml), 'Not a .xml file for mesh!'

    name, xml_ext = os.path.splitext(mesh_xml)
    return ''.join([name, '_facet_region', xml_ext])


if __name__ == '__main__':
    if len(sys.argv) == 1:
        mesh_files = filter(is_xml_mesh_file, os.listdir('.'))
    else:
        mesh_files = sys.argv[1:]


    facet_files = map(get_facet_xml, mesh_files)

    for mesh_file, facet_file in zip(mesh_files, facet_files):
        print 'Converting', mesh_file, facet_file,
        # Convert from xml to xdmf
        mesh = Mesh(mesh_file)
        new_mesh_file = ''.join([os.path.splitext(mesh_file)[0], '.xdmf'])
        f_f = MeshFunction('size_t', mesh, facet_file)
        new_facet_file = ''.join([os.path.splitext(facet_file)[0], '.xdmf'])
        print 'to', new_mesh_file, new_facet_file
        print
        File(new_mesh_file) << mesh
        File(new_facet_file) << f_f

        # Test
        mesh = Mesh(new_mesh_file)
        f_f = MeshFunction('size_t', mesh, new_facet_file)
