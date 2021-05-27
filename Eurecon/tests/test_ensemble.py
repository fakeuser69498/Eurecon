import os
import os.path as path
import sys
import unittest

import numpy as np

from eurecon import Ensemble, Eurecon, Parser

base_dir = path.dirname(path.dirname(path.abspath(__file__)))
sys.path.append(base_dir)






class TestEnsemble(unittest.TestCase):

    def gen_ens(
            self,
            data_path = 'tests/tests_data/molecules/4d2i.pdb',
            tess_path = 'tests/tests_data/tesselation_axes/tesselation_vertice.txt',
            rmsd = 2,
            partition = 0.5
    ):
        parser = Parser()
        transform = parser.parse_transform(
            tess_path,
            rmsd,
            partition
        )

        coords, base_conformation = parser.parse_base_conformation(
            data_path,
            None
        )

        ensemble = Ensemble(
            base_conformation,
            transform,
            parser
        )
        return ensemble

    def test_generate_ensemble(self):
        ensemble = self.gen_ens()

        generated_ensemble = np.loadtxt('tests/tests_data/generated_ensemble.txt', delimiter=' ')
        ensemble.generate_ensemble(False)
        np.testing.assert_allclose(ensemble.conformations[0], generated_ensemble, rtol=1e-3)

    TEST_ANGLE = 4.714494162841217223e-02

    def test_count_rotation_matrix(self):
        ensemble = self.gen_ens()

        res = ensemble.count_rotation_matrix(ensemble.transform.axes[0], self.TEST_ANGLE)
        generated_matrix = np.loadtxt('tests/tests_data/rotation_matrix.txt', delimiter=' ')
        np.testing.assert_allclose(res, generated_matrix, rtol=1e-4)

    def test_count_angle(self):
        ensemble = self.gen_ens()
        res, part = ensemble.rmsd_angle(ensemble.transform.axes[0])
        np.testing.assert_allclose(res, self.TEST_ANGLE, rtol=1e-3)

    def test_gen_conformation(self):
        ensemble = self.gen_ens()
        res = ensemble.generate_conformation(ensemble.transform.axes[0])
        generated_conf = np.loadtxt('tests/tests_data/generated_conformation.txt', delimiter=' ')
        np.testing.assert_allclose(res, generated_conf, rtol=1e-3)

    def test_trivial_rmsd(self):
        ensemble = self.gen_ens()
        conformation = ensemble.generate_conformation(ensemble.transform.axes[0]) - ensemble.base_conformation.center_of_mass
        res = ensemble.calc_trivial_rmsd(conformation)
        np.testing.assert_allclose(res, 2, rtol=1e-3)

if __name__=='__main__':
      unittest.main()
