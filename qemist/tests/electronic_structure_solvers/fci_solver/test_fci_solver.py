import unittest

from pyscf import gto, scf

from qemist.electronic_structure_solvers import FCISolver

H2 = """
    H 0.00 0.00 0.0
    H 0.00 0.00 0.74137727
    """

Be = """Be 0.0 0.0 0.0"""

class FCISolverTest(unittest.TestCase):

    def test_h2_no_mf(self):
        """ Test the FCISolver against result from reference implementation with
        mean field not calculated."""
        mol = gto.Mole()
        mol.atom = H2
        mol.basis = "3-21g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        solver = FCISolver()
        energy = solver.simulate(mol)

        self.assertAlmostEqual(energy, -1.1478300596229851, places=8)
        #TODO: test RDM here

    def test_h2_with_mf(self):
        """ Test the FCISolver against result from reference implementation with
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

        solver = FCISolver()
        energy = solver.simulate(mol, mf)

        self.assertAlmostEqual(energy, -1.1478300596229851, places=8)
        #TODO: test RDM here

    def test_be_no_mf(self):
        """ Test the FCISolver against result from reference implementation with
        mean field not calculated."""
        mol = gto.Mole()
        mol.atom = Be
        mol.basis = "3-21g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        solver = FCISolver()
        energy = solver.simulate(mol)

        self.assertAlmostEqual(energy, -14.531444379108095, places=8)
        #TODO: test RDM here

    def test_be_with_mf(self):
        """ Test the FCISolver against result from reference implementation with
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

        solver = FCISolver()
        energy = solver.simulate(mol, mf)

        self.assertAlmostEqual(energy, -14.531444379108095, places=8)
        #TODO: test RDM here

    def test_get_rdm_without_simulate(self):
        """Test that the runtime error is raised when user calls get RDM without
        first running a simulation."""
        solver = FCISolver()
        self.assertRaises(RuntimeError, solver.get_rdm);

if __name__ == "__main__":
    unittest.main()
