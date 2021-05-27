import os
import os.path as path
import sys
import unittest

import numpy as np


from eurecon import Conformation
base_dir = path.dirname(path.dirname(path.abspath(__file__)))
sys.path.append(base_dir)



class Tests(unittest.TestCase):
    def test_calc_inertia_tensor(self):
        COM_moved = np.loadtxt('tests/tests_data/COM_moved_to_origin.txt', delimiter = ',')
        object_length = len(COM_moved)
        simple_object = Conformation([], COM_moved.T, object_length, 'whatever', None)
        result = np.loadtxt('tests/tests_data/inertia_tensor.txt', delimiter=',')
        np.testing.assert_allclose(simple_object.calc_inertia_tensor(), result, rtol=1e-8)

    def test_calc_mass_center(self):
        conf_not_moved = np.loadtxt('tests/tests_data/stl_input_array.txt', delimiter = ',')
        object_length = len(conf_not_moved.T)

        simple_object = Conformation([], conf_not_moved, object_length, 'whatever', None)
        result = np.loadtxt('tests/tests_data/COM.txt', delimiter = ',')
        self.assertEqual(simple_object.calc_mass_center(conf_not_moved).tolist(), result.tolist())

    def test_get_coords_moved_to_mass_center(self):
        conf_not_moved = np.loadtxt('tests/tests_data/stl_input_array.txt', delimiter = ',')
        object_length = len(conf_not_moved.T)

        simple_object = Conformation([], conf_not_moved, object_length, 'whatever', None)
        result = np.loadtxt('tests/tests_data/COM_moved_to_origin.txt', delimiter = ',')
        self.assertEqual(simple_object.get_coords_moved_to_mass_center(conf_not_moved).tolist(), result.tolist())

      
if __name__=='__main__':
      unittest.main()
