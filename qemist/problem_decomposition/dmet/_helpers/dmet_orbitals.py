"""
Construct localized orbitals for DMET calculation
"""
from pyscf import scf, ao2mo
import numpy as np
from functools import reduce

class dmet_orbitals:
    """
    Localized the SCF orbitals
    Get the integrals and matrices for active and core regions in DMET
    """
    def __init__(self, mol, mf, active_space, localization_function):
        """
        Let us initialize the class
        Get energies and integrals for the entire system
        :param mol: Molecule object (PySCF)
        :param mf: Mean field object (PySCF)
        :param active_space: The active space in DMET calculation (All orbitals in the initial SCF calculation)
        :param localization_function: Function for electron localization.
        """
        # TODO: Is active space always claculated from the molecule?

        # Obtain the elements from the low-level SCF calculations
        self.mol_full = mol
        self.mf_full = mf
        self.low_scf_energy = mf.e_tot
        low_scf_dm = reduce(np.dot, (mf.mo_coeff, np.diag(mf.mo_occ), mf.mo_coeff.T))
        low_scf_twoint = scf.hf.get_veff(mf.mol, low_scf_dm, 0, 0, 1)
        self.low_scf_fock = mf.mol.intor('cint1e_kin_sph') + mf.mol.intor('cint1e_nuc_sph') + low_scf_twoint

        # Define the active space if possible
        self.dmet_active_orbitals = np.zeros([mf.mol.nao_nr()], dtype=int)
        self.dmet_active_orbitals[active_space] = 1
        self.number_active_orbitals = np.sum(self.dmet_active_orbitals)
        self.number_active_electrons = int(np.rint(mf.mol.nelectron - np.sum(mf.mo_occ[self.dmet_active_orbitals==0])))

        # Localize the orbitals (IAO)
        self.localized_mo = localization_function(mol, mf)

        # Define the core space if possible (Initial calculations treat the entire molecule ...)
        core_mo_dm = np.array(mf.mo_occ, copy=True)
        core_mo_dm[self.dmet_active_orbitals == 1] = 0
        core_ao_dm = reduce(np.dot, (mf.mo_coeff, np.diag(core_mo_dm), mf.mo_coeff.T))
        core_twoint = scf.hf.get_veff(mf.mol, core_ao_dm, 0, 0, 1)
        core_oneint = self.low_scf_fock - low_scf_twoint + core_twoint

        # Define the energies and matrix elements based on the localized orbitals
        self.core_constant_energy = mf.mol.energy_nuc() + np.einsum('ij,ij->', core_oneint - 0.5*core_twoint, core_ao_dm)
        self.active_oneint = reduce(np.dot, (self.localized_mo.T, core_oneint, self.localized_mo))
        self.active_fock = reduce(np.dot, (self.localized_mo.T, self.low_scf_fock, self.localized_mo))

    def dmet_fragment_hamiltonian(self, bath_orb, norb_high, onerdm_core):
        """
        Construct the Hamiltonian for a DMET fragment
        :param bath_orb: the bath orbitals
        :param norb_high: the number of orbitals in the fragment
        :param onerdm_core: the one particle RDM
        :return: one-electron integrals, fock matrix, two-electron integrals for a fragment
        """

        # one electron integrals
        frag_oneint = reduce(np.dot, (bath_orb[ : , : norb_high].T, self.active_oneint, bath_orb[ : , : norb_high]))

        # fock matrix
        density_matrix = reduce(np.dot, (self.localized_mo, onerdm_core, self.localized_mo.T))
        two_int = scf.hf.get_veff(self.mol_full, density_matrix, 0, 0, 1)
        new_fock = self.active_oneint + reduce(np.dot, ((self.localized_mo.T, two_int, self.localized_mo)))
        frag_fock = reduce(np.dot, (bath_orb[ : , : norb_high ].T, new_fock, bath_orb[ : , : norb_high]))

        # two electron integrals
        coefficients = np.dot(self.localized_mo, bath_orb[ : , : norb_high])
        frag_twoint = ao2mo.outcore.full_iofree(self.mol_full, coefficients, compact=False).reshape( \
                                                norb_high,  norb_high,  norb_high,  norb_high)

        return frag_oneint, frag_fock, frag_twoint

