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

from openqemist.problem_decomposition.problem_decomposition import ProblemDecomposition

class MockProblemDecomposition(ProblemDecomposition):

    def simulate(self, molecule, fragment_atoms, mean_field = None):
        """ Implements the abstract function in the superclass. """
        pass

class ProblemDecompositionTest(unittest.TestCase):

    def test_init(self):
        solver = MockProblemDecomposition()
        self.assertTrue(solver.electronic_structure_solver == None)
        self.assertTrue(solver.electron_localization_method == None)


