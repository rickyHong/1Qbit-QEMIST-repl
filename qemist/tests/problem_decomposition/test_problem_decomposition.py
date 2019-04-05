import unittest

from qemist.problem_decomposition.problem_decomposition import ProblemDecomposition

class MockProblemDecomposition(ProblemDecomposition):

    def simulate(self, molecule, fragment_atoms, mean_field = None):
        """ Implements the abstract function in the superclass. """
        pass

class ProblemDecompositionTest(unittest.TestCase):

    def test_init(self):
        solver = MockProblemDecomposition()
        self.assertTrue(solver.electronic_structure_solver == None)
        self.assertTrue(solver.electron_localization_method == None)


