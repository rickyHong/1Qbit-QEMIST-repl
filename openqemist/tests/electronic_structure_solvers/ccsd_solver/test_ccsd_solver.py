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

from pyscf import gto, scf

from openqemist.electronic_structure_solvers import CCSDSolver

H2 = """
   H 0.00 0.00 0.0
   H 0.00 0.00 0.74137727
   """

Be = """Be 0.0 0.0 0.0"""

class CCSDSolverTest(unittest.TestCase):

    def test_h2_no_mf(self):
        """ Test the CCSDSolver against result from reference implementation with
        mean field not calculated."""
        mol = gto.Mole()
        mol.atom = H2
        mol.basis = "3-21g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        solver = CCSDSolver()
        energy = solver.simulate(mol)

        self.assertAlmostEqual(energy, -1.1478300596229851, places=8)
        #TODO: test RDM here

    def test_h2_with_mf(self):
        """ Test the CCSDSolver against result from reference implementation with
        mean field calculated and passed in."""
        mol = gto.Mole()
        mol.atom = H2
        mol.basis = "3-21g"
        mol.charge = 0
        mol.spin = 0
        mol.build()


        mf = scf.RHF(mol)
        mf.verbose = 0
        mf.scf()

        solver = CCSDSolver()
        energy = solver.simulate(mol, mf)

        self.assertAlmostEqual(energy, -1.1478300596229851, places=8)
        #TODO: test RDM here

    def test_be_no_mf(self):
        """ Test the CCSDSolver against result from reference implementation with
        mean field not calculated."""
        mol = gto.Mole()
        mol.atom = Be
        mol.basis = "3-21g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        solver = CCSDSolver()
        energy = solver.simulate(mol)

        self.assertAlmostEqual(energy, -14.531416589890926, places=8)
        #TODO: test RDM here

    def test_be_with_mf(self):
        """ Test the CCSDSolver against result from reference implementation with
        mean field calculated and passed in."""
        mol = gto.Mole()
        mol.atom = Be
        mol.basis = "3-21g"
        mol.charge = 0
        mol.spin = 0
        mol.build()


        mf = scf.RHF(mol)
        mf.verbose = 0
        mf.scf()

        solver = CCSDSolver()
        energy = solver.simulate(mol, mf)

        self.assertAlmostEqual(energy, -14.531416589890926, places=8)
        #TODO: test RDM here

    def test_get_rdm_without_simulate(self):
        """Test that the runtime error is raised when user calls get RDM without
        first running a simulation."""
        solver = CCSDSolver()
        self.assertRaises(RuntimeError, solver.get_rdm);

if __name__ == "__main__":
    unittest.main()
