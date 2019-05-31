import unittest

from pyscf import gto

from openqemist.electronic_structure_solvers import FCISolver, CCSDSolver
from openqemist.problem_decomposition import DMETProblemDecomposition
from openqemist.problem_decomposition.electron_localization import iao_localization, meta_lowdin_localization

H10_RING = """
H          1.6180339887          0.0000000000          0.0000000000
H          1.3090169944          0.9510565163          0.0000000000
H          0.5000000000          1.5388417686          0.0000000000
H         -0.5000000000          1.5388417686          0.0000000000
H         -1.3090169944          0.9510565163          0.0000000000
H         -1.6180339887          0.0000000000          0.0000000000
H         -1.3090169944         -0.9510565163          0.0000000000
H         -0.5000000000         -1.5388417686          0.0000000000
H          0.5000000000         -1.5388417686          0.0000000000
H          1.3090169944         -0.9510565163          0.0000000000
        """

H4_RING = """
H   0.7071067811865476   0.0                 0.0
H   0.0                  0.7071067811865476  0.0
H  -1.0071067811865476   0.0                 0.0
H   0.0                 -1.0071067811865476  0.0
        """

class DMETProblemDecompositionMicrosoftQSharpTest(unittest.TestCase):

    def test_h4ring_vqe_uccsd_microsoftqsharp(self):
        """
        DMET on H4 ring with fragment size one, using VQE-UCCSD backend
        from Microsoft.
        """

        mol = gto.Mole()
        mol.atom = H4_RING
        mol.basis = "minao"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        # Initialize VQE object with MicrosoftQSharp backend
        from openqemist.quantum_solvers import MicrosoftQSharpParametricSolver
        from openqemist.electronic_structure_solvers import VQESolver
        vqe = VQESolver()
        vqe.hardware_backend_type = MicrosoftQSharpParametricSolver
        vqe.ansatz_type = MicrosoftQSharpParametricSolver.Ansatze.UCCSD

        # Run DMET
        dmet = DMETProblemDecomposition()
        dmet.electron_localization_method = meta_lowdin_localization
        dmet.electronic_structure_solver = vqe
        energy_vqe = dmet.simulate(mol, [1,1,1,1])

        self.assertAlmostEqual(energy_vqe, -1.9916120594, places=3)

    def test_h10ring_vqe_uccsd_microsoftqsharp(self):
        """
           DMET on H10 ring with fragment size one, using VQE-UCCSD backend
           from Microsoft.
        """
        mol = gto.Mole()
        mol.atom = H10_RING
        mol.basis = "minao"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        # Initialize VQE object with MicrosoftQSharp backend
        from openqemist.quantum_solvers import MicrosoftQSharpParametricSolver
        from openqemist.electronic_structure_solvers import VQESolver
        vqe = VQESolver()
        vqe.hardware_backend_type = MicrosoftQSharpParametricSolver
        vqe.ansatz_type = MicrosoftQSharpParametricSolver.Ansatze.UCCSD

        # Run DMET
        dmet = DMETProblemDecomposition()
        dmet.electron_localization_method = meta_lowdin_localization
        dmet.electronic_structure_solver = vqe
        energy_vqe = dmet.simulate(mol, [1,1,1,1,1,1,1,1,1,1])

        self.assertAlmostEqual(energy_vqe, -5.367532792556453, places=3)

if __name__ == "__main__":
    unittest.main()
