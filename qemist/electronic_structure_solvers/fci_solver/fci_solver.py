from pyscf import ao2mo, fci, scf
from functools import reduce
import numpy as np

from ..electronic_structure_solver import ElectronicStructureSolver

class FCISolver(ElectronicStructureSolver):

    def __init__(self):
        self.ci = None
        self.norb = None
        self.nelec = None
        self.cisolver = None

    def simulate(self, molecule, mean_field=None):
        # Calculate the mean field if the user has not already done it.
        if not mean_field:
            mean_field = scf.RHF(molecule)
            mean_field.verbose = 0
            mean_field.scf()

        # Check the convergence of the mean field
        if not mean_field.converged:
            warnings.warn("FCISolver simulating with mean field not converged.", RuntimeWarning)

        h1 = reduce(np.dot, (mean_field.mo_coeff.T, mean_field.get_hcore(),
            mean_field.mo_coeff))
        twoint = mean_field._eri
        self.norb = len(mean_field.mo_energy)
        eri = ao2mo.restore(8, twoint, self.norb)
        eri = ao2mo.incore.full(eri, mean_field.mo_coeff)
        eri = ao2mo.restore(1, eri, self.norb)
        self.cisolver = fci.direct_spin0.FCI(molecule)
        self.cisolver.verbose = 0
        self.nelec = molecule.nelectron
        energy, self.ci = self.cisolver.kernel(h1, eri, h1.shape[1], self.nelec, ecore=molecule.energy_nuc())

        return energy


    def get_rdm(self):
        if not self.norb or not self.nelec:
            raise RuntimeError("Cannot retrieve RDM because no simulation has been run.")

        fci_onerdm = self.cisolver.make_rdm1(self.ci, self.norb, self.nelec)
        fci_twordm = self.cisolver.make_rdm2(self.ci, self.norb, self.nelec)

        return fci_onerdm, fci_twordm


