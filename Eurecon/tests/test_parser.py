import os
import os.path as path
import sys
import unittest

import MDAnalysis as mda
import numpy as np
import open3d as o3d
from biopandas.pdb import PandasPdb

from eurecon import Conformation, Parser, Transform

base_dir = path.dirname(path.dirname(path.abspath(__file__)))
sys.path.append(base_dir)





class Tests(unittest.TestCase):

    def _get_parser(self):
        parser = Parser()
        parser.result_path = 'results/TEST'
        return parser

    def test_parse_transform_axes(self):
        axes_file_path = 'tests/tests_data/tesselation_axes/tesselation_vertices_layer_8.txt'
        rmsd = 2
        partition = 0.5

        test_obj = Parser().parse_transform(axes_file_path, rmsd, partition)
        result = np.loadtxt(axes_file_path)

        np.testing.assert_allclose(test_obj.axes.tolist(), result.tolist(), rtol=1e-8)
        self.assertEqual(test_obj.rmsd, rmsd)
        self.assertEqual(test_obj.partition, partition)

    def test_mol2_and_pdb(self):
        file_name = 'tests/tests_data/molecules/1b5e_1.mol2'

        data_object = mda.Universe(file_name)
        point_coords = data_object.select_atoms('protein')
        test_sample = Parser().parse_molecule(file_name)

        # Test data object
        print(str(test_sample))
        print(str(data_object))
        self.assertEqual(str(test_sample[0]), str(data_object))
        

        # Test coords
        points_coords = point_coords.atoms.positions.T
        self.assertEqual(points_coords.tolist(), test_sample[1].tolist())

        # Test mol length
        molecule_length = len(points_coords[0])
        self.assertEqual(molecule_length, test_sample[2])

    def test_data_object_pcd(self):
        file_name = 'tests/tests_data/molecules/Rookpts.pts'

        test_sample = Parser().parse_molecule(file_name)
        data_object = o3d.io.read_point_cloud(file_name)

        # Test data object
        self.assertEqual(str(test_sample[0]), str(data_object))

        # Test coords
        points_coords = np.transpose(np.asarray(data_object.points))
        self.assertEqual(points_coords.tolist(), test_sample[1].tolist())

        # Test mol length
        molecule_length = len(points_coords[0])
        self.assertEqual(molecule_length, test_sample[2])
    
    def test_data_object_mesh(self):
        file_name = 'tests/tests_data/molecules/Rookstl.stl'

        test_sample = Parser().parse_molecule(file_name)
        data_object = o3d.io.read_triangle_mesh(file_name)

        # Test object
        self.assertEqual(str(test_sample[0]), str(data_object))

        # Test coords
        pcdd = o3d.geometry.PointCloud()
        pcdd.points = data_object.vertices
        points_coords = np.transpose(np.asarray(data_object.vertices))
        self.assertEqual(points_coords.tolist(), test_sample[1].tolist())

        # Test mol length
        molecule_length = len(np.transpose(points_coords))
        self.assertEqual(molecule_length, test_sample[2])

    def test_create_path(self):
        data_file_name = '4d2i.pdb'
        path_name = 'TEST'
        counter = 1
        format = data_file_name.split('.')[-1]
        new_file_name = f'results/{path_name}/{counter}.{format}'
        sample = Parser()
        test_sample = sample.create_path(new_file_name)
        a = os.path.isdir('results/TEST')
        self.assertEqual(True, a)

    def test_write_mol2_pdb_conformation(self):
        file_name = 'tests/tests_data/molecules/1b5e_1.mol2'
        counter = 16
        data_file_name = '1b5e_1.mol2'
        sample = self._get_parser()

        data_object = mda.Universe(file_name)
        points_coords = data_object.atoms.positions
        molecule_length = len(points_coords[0])

        sample.write_conformation(data_object, points_coords, data_file_name, counter, molecule_length)

        a = os.path.isfile('results/TEST/16.mol2')
        self.assertEqual(True, a)

    def test_write_mesh_conformation(self):
        file_name = 'tests/tests_data/molecules/Rookstl.stl'
        counter = 16
        data_file_name = 'Rookstl.stl'
        sample = self._get_parser()

        data_object = o3d.io.read_triangle_mesh(file_name)
        np_vertices = np.array(data_object.vertices)
        molecule_length = len(np.transpose(np_vertices))

        sample.write_conformation(data_object, np_vertices, data_file_name, counter, molecule_length)
        a = os.path.isfile('results/TEST/16.stl')

        self.assertEqual(True, a)

    def test_write_point_conformation(self):
        file_name = 'tests/tests_data/molecules/Rookpts.pts'
        counter = 16
        data_file_name = 'Rookpts.pts'
        sample = self._get_parser()

        data_object = o3d.io.read_point_cloud(file_name)
        points_coords = data_object.points
        molecule_length = len(points_coords[0])

        sample.write_conformation(data_object, points_coords, data_file_name, counter, molecule_length)
        a = os.path.isfile('results/TEST/16.pts')

        self.assertEqual(True, a)

    def test_write_default(self):
        sample = self._get_parser()
        file_name = 'tests/tests_data/molecules/4d2i.pdb'
        counter = 16
        data_file_name = '4d2i.pdb'

        data_object = mda.Universe(file_name)
        data_object = data_object.select_atoms('protein')
        points_coords = data_object.atoms.positions.T
        molecule_length = len(points_coords[0])

        simple_object = Conformation(data_object, points_coords, molecule_length, file_name, None)
        test_object = sample.write_default(simple_object, points_coords)

        a = os.path.isfile('results/TEST/4d2i.pdb')
        self.assertEqual(True, a)

    def test_write_all_conformations(self):
        sample = self._get_parser()
        file_name = 'tests/tests_data/molecules/4d2i.pdb'
        data_file_name = '4d2i.pdb'

        data_object = mda.Universe(file_name)
        data_object1 = data_object.select_atoms('protein')
        points_coords = data_object1.atoms.positions.T
        molecule_length = len(points_coords[0])

        simple_object = Conformation(data_object1, points_coords, molecule_length, file_name, None)
        conf_list = [points_coords.T, points_coords.T]

        axes = np.loadtxt('tests/tests_data/tesselation_axes/tesselation_vertices_layer_8.txt')
        transform = Transform(axes, 2, 0.5)
        test_object = sample.write_all_conformations(conf_list, simple_object, transform)
        a = os.path.exists('tests/results/TEST/0.pdb')
        b = os.path.exists('tests/results/TEST/1.pdb')
        self.assertEqual(a, b)

      
      
if __name__=='__main__':
      unittest.main()
