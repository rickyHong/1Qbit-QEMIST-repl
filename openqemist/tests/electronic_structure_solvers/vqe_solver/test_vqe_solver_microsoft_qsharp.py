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
from enum import Enum

from pyscf import gto, scf

from openqemist.electronic_structure_solvers import VQESolver, FCISolver
from openqemist.quantum_solvers.parametric_quantum_solver import ParametricQuantumSolver

H2 = """
   H 0.00 0.00 0.0
   H 0.00 0.00 0.74137727
   """


class VQESolverTest(unittest.TestCase):

    def test_h2_sto3g(self):
        from openqemist.quantum_solvers import MicrosoftQSharpParametricSolver

        mol = gto.Mole()
        mol.atom = H2
        mol.basis = "sto-3g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        solver = VQESolver()
        solver.hardware_backend_type = MicrosoftQSharpParametricSolver
        solver.ansatz_type = MicrosoftQSharpParametricSolver.Ansatze.UCCSD
        #solver.optimizer =

        energy = solver.simulate(mol)

        self.assertAlmostEqual(energy, -1.1372704178510415, places=3)

    @unittest.skip("After a sign change to the t2 amplitudes, may be quick enough to converge for unittests")
    def test_h2_321g(self):
        from openqemist.quantum_solvers import MicrosoftQSharpParametricSolver

        mol = gto.Mole()
        mol.atom = H2
        mol.basis = "3-21g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        solver = VQESolver()
        solver.hardware_backend_type = MicrosoftQSharpParametricSolver
        solver.ansatz_type = MicrosoftQSharpParametricSolver.Ansatze.UCCSD
        #solver.optimizer =

        energy = solver.simulate(mol)

        self.assertAlmostEqual(energy, -1.1478300615818977, places=3)


if __name__ == "__main__":
    unittest.main()
