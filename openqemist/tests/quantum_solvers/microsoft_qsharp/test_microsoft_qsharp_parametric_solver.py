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
import numpy as np

from pyscf import gto, scf, mp

from openqemist.quantum_solvers import MicrosoftQSharpParametricSolver

def get_file_path_stub():
    """ Gets the path of the test files from anywhere in the test tree."

    The direcory structure should be $SOMETHING/QEMIST/qemist/tests/$SOMETHINGELSE
    so we trim after "tests", then add the path to the results files so we can
    run the tests from anywhere in the tree."""
    import os
    cwd = os.getcwd()
    tests_root = cwd[0:cwd.find("tests") + 5]
    return tests_root + "/quantum_solvers/microsoft_qsharp/data/"


H2 = """
        H       0.00        0.00        0.0
        H       0.00        0.00        0.74137727
"""

H4 = """
    H   0.7071067811865476   0.0                 0.0
    H   0.0                  0.7071067811865476  0.0
    H  -1.0071067811865476   0.0                 0.0
    H   0.0                 -1.0071067811865476  0.0
"""

# References for H2
amplitudes_H2 = [1.69971474e-05, 5.65855806e-02]
expected_1rdm_H2 = np.load(get_file_path_stub() + "h2_onerdm.npy")
expected_2rdm_H2 = np.load(get_file_path_stub() + "h2_twordm.npy")

# references for H4
amplitudes_H4 = [-3.00520142e-05, -3.41547577e-05,  7.61837556e-06 ,-2.24075399e-05,
                 1.12512690e-02,  3.42703399e-01,  3.44523818e-02,  1.46586868e-02,
                 7.69034155e-02,  7.99964875e-03, -1.81430817e-01, -1.06611015e-01,
                 1.12938142e-02, -3.75164050e-02]
expected_1rdm_H4 = np.load(get_file_path_stub() + "h4_onerdm.npy")
expected_2rdm_H4 = np.load(get_file_path_stub() + "h4_twordm.npy")

class MicrosoftQSharpParametricSolverTest(unittest.TestCase):

    def test_no_mf_H2(self):
        """Tests that all the values are set correctly in the constructor."""
        mol = gto.Mole()
        mol.atom = H2
        mol.basis = "sto-3g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        ansatz = MicrosoftQSharpParametricSolver.Ansatze.UCCSD
        solver = MicrosoftQSharpParametricSolver(ansatz, mol)

        # Test that the number of variational parameters is as expected
        self.assertEqual(solver.amplitude_dimension, 2)

        # Test "simulate"
        energy = solver.simulate(amplitudes_H2)
        self.assertAlmostEqual(energy, -1.13727, places=3)
        
        # Test "get_rdm"
        one_rdm, two_rdm = solver.get_rdm()
        
        for index, value_ref in np.ndenumerate(expected_1rdm_H2):
            self.assertAlmostEqual(one_rdm[index], value_ref, msg='one_rdm error at index ' + str(index), delta=1e-3)
        for index, value_ref in np.ndenumerate(expected_2rdm_H2):
            self.assertAlmostEqual(two_rdm[index], value_ref, msg='two_rdm error at index ' + str(index), delta=1e-3)

    @unittest.skip("A fix may be needed for this test to pass")
    def test_no_mf_H4(self):
        """Tests that all the values are set correctly in the constructor."""
        mol = gto.Mole()
        mol.atom = H4
        mol.basis = "sto-3g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        ansatz = MicrosoftQSharpParametricSolver.Ansatze.UCCSD
        solver = MicrosoftQSharpParametricSolver(ansatz, mol)

        # Test that the number of variational parameters is as expected
        self.assertEqual(solver.amplitude_dimension, 14)

        # Test "simulate"
        energy = solver.simulate(amplitudes_H4)
        self.assertAlmostEqual(energy, -1.9778376625571703, places=3)

        # Test "get_rdm"
        one_rdm, two_rdm = solver.get_rdm()
        
        for index, value_ref in np.ndenumerate(expected_1rdm_H4):
            self.assertAlmostEqual(one_rdm[index], value_ref, msg='one_rdm error at index ' + str(index), delta=1e-3)
        for index, value_ref in np.ndenumerate(expected_2rdm_H4):
            self.assertAlmostEqual(two_rdm[index], value_ref, msg='two_rdm error at index ' + str(index), delta=1e-3)
    
    def test_mf_H2(self):
        """Tests that all the values are set correctly in the constructor."""
        mol = gto.Mole()
        mol.atom = H2
        mol.basis = "sto-3g"
        mol.charge = 0
        mol.spin = 0
        mol.build()
        mean_field = scf.RHF(mol)
        mean_field.verbose = 0
        mean_field.scf()

        ansatz = MicrosoftQSharpParametricSolver.Ansatze.UCCSD
        solver = MicrosoftQSharpParametricSolver(ansatz, mol, mean_field)

        self.assertEqual(solver.amplitude_dimension, 2)

        energy = solver.simulate(amplitudes_H2)

        self.assertAlmostEqual(energy, -1.13727, places=3)

    def test_simulate_dimension_throw(self):
        """Tests that all the values are set correctly in the constructor."""
        mol = gto.Mole()
        mol.atom = H2
        mol.basis = "sto-3g"
        mol.charge = 0
        mol.spin = 0
        mol.build()
        mean_field = scf.RHF(mol)
        mean_field.verbose = 0
        mean_field.scf()

        ansatz = MicrosoftQSharpParametricSolver.Ansatze.UCCSD
        solver = MicrosoftQSharpParametricSolver(ansatz, mol, mean_field)

        # solver.amplitude_dimension = 2, this should throw.
        self.assertRaises(ValueError, solver.simulate, [0])

if __name__ == "__main__":
    unittest.main()
