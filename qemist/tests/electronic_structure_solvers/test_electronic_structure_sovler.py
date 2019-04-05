import unittest

from qemist.electronic_structure_solvers.electronic_structure_solver import ElectronicStructureSolver

class MockElectronicStructureSolver(ElectronicStructureSolver):

    def simulate(self, molecule, mean_field=None):
        """ Implements the abstract function in the superclass. """
        pass

    def get_rdm(self):
        """ Implements the abstract function in the superclass. """
        pass

class ElectronicStructureSolverTest(unittest.TestCase):

    def test_init(self):
        solver = MockElectronicStructureSolver()
        self.assertTrue(solver.quantum_solver == None)

