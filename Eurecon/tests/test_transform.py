import os.path as path
import sys
import unittest

import numpy as np

from eurecon import Transform

base_dir = path.dirname(path.dirname(path.abspath(__file__)))
sys.path.append(base_dir)


class Tests(unittest.TestCase):

    def test_axes_normalization(self):
        axes = np.loadtxt('tests/tests_data/tesselation_axes/tesselation_vertices_layer_8.txt')
        simple_object = Transform(axes, 2, 0.5)
        result = np.loadtxt('tests/tests_data/normalized_array.txt', delimiter=',')
        self.assertEqual(Transform.normalize_axes(self,axes).tolist(), result.tolist())
      
      
if __name__=='__main__':
      unittest.main()
