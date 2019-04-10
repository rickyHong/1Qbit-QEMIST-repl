import unittest
from enum import Enum

from pyscf import gto, scf

from qemist.electronic_structure_solvers import VQESolver
from qemist.quantum_solvers.parametric_quantum_solver import ParametricQuantumSolver

H2 = """
   H 0.00 0.00 0.0
   H 0.00 0.00 0.74137727
   """

Be = """Be 0.0 0.0 0.0"""

class MockQuantumSolver(ParametricQuantumSolver):
    """ Class that mocks out the abstract ParametricQuantumSolver for testing.

    For the purpose of this test, all we need is that we return the same thing,
    causing the optimizer to converge. We assume that the optimizer is working
    correctly.
    """

    class Ansatze(Enum):
        UCCSD = 0

    def __init__(self, molecule, anstatz, mean_field=None):
        self.initial_wavefunction = None
        self.n_qubits = None
        self.amplitude_dimension = None

    def simulate(self, amplitudes):
        return 3

    def get_rdm(self):
        pass

class VQESolverTest(unittest.TestCase):

    def test_be_no_mf_default_optimizer(self):
        """ Test that the solver returns with the default optimizer from scipy.

        Since simlate is always returning the same number, this should force the
        optimizer to converge.
        """
        mol = gto.Mole()
        mol.atom = Be
        mol.basis = "3-21g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        solver = VQESolver()
        solver.hardware_backend_type = MockQuantumSolver
        solver.ansatz_type = MockQuantumSolver.Ansatze.UCCSD

        energy = solver.simulate(mol)

        # The converged energy value should be what the mock solver returned at
        # the simulation call.
        self.assertEqual(energy, 3)


if __name__ == "__main__":
    unittest.main()
