#   Copyright 2019 1QBit
#   
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import unittest

from pyscf import gto

from openqemist.electronic_structure_solvers import FCISolver, CCSDSolver
from openqemist.problem_decomposition import DMETProblemDecomposition
from openqemist.problem_decomposition.electron_localization import iao_localization, meta_lowdin_localization

H10_RING = """
        H      0.970820393250   0.000000000000   0.000000000000
        H      0.785410196625   0.570633909777   0.000000000000
        H      0.300000000000   0.923305061153   0.000000000000
        H     -0.300000000000   0.923305061153   0.000000000000
        H     -0.785410196625   0.570633909777   0.000000000000
        H     -0.970820393250   0.000000000000   0.000000000000
        H     -0.785410196625  -0.570633909777   0.000000000000
        H     -0.300000000000  -0.923305061153   0.000000000000
        H      0.300000000000  -0.923305061153   0.000000000000
        H      0.785410196625  -0.570633909777   0.000000000000
        """

TETRA = """
        C   -0.580517    0.479710   -0.503842
        H   -1.269514    1.045226   -1.098983
        C   -0.462763   -0.528766    0.571790
        H   -1.009781   -1.152554    1.250252
        C    0.529158    0.565261    0.470571
        H    1.157329    1.232467    1.026343
        C    0.514143   -0.516046   -0.538653
        H    1.121845   -1.126087   -1.176808
        """


class DMETOnClassicalTest(unittest.TestCase):

    def test_hring_ml_fci_no_mf(self):
        mol = gto.Mole()
        mol.atom = H10_RING
        mol.basis = "3-21g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        solver = DMETProblemDecomposition()
        solver.electron_localization_method = meta_lowdin_localization
        solver.electronic_structure_solver = FCISolver()
        energy = solver.simulate(mol, [1,1,1,1,1,1,1,1,1,1])

        self.assertAlmostEqual(energy, -4.4989730244, places=6)

    def test_tetra_ml_fci_no_mf(self):
        mol = gto.Mole()
        mol.atom = TETRA
        mol.basis = "3-21g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        solver = DMETProblemDecomposition()
        solver.electron_localization_method = meta_lowdin_localization
        solver.electronic_structure_solver = CCSDSolver()
        energy = solver.simulate(mol, [1,1,1,1,1,1,1,1])

        self.assertAlmostEqual(energy, -153.1879732845, places=6)

if __name__ == "__main__":
    unittest.main()
