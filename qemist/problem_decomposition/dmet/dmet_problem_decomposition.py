import warnings

import scipy
from pyscf import scf
from functools import reduce
import numpy as np

from qemist.electronic_structure_solvers import CCSDSolver
from ..problem_decomposition import ProblemDecomposition
from ..electron_localization import iao_localization

from . import _helpers as helpers

class DMETProblemDecomposition(ProblemDecomposition):

    def __init__(self):
        self.verbose = True
        self.electronic_structure_solver = CCSDSolver()
        self.electron_localization_method = iao_localization

    def simulate(self, molecule, fragment_atoms, mean_field = None):

        # Calculate the mean field if the user has not already done it.
        if not mean_field:
            mean_field = scf.RHF(molecule)
            mean_field.verbose = 0
            mean_field.scf()

        # Check the convergence of the mean field
        if not mean_field.converged:
            warnings.warn("DMET simulating with mean field not converged.", RuntimeWarning)

        # Construct orbital object
        orbitals = helpers._orbitals(molecule, mean_field,
                range(molecule.nao_nr()), self.electron_localization_method)

        # TODO: remove last argument, combiging fragments not supported
        orb_list, orb_list2, atom_list2 = helpers._fragment_constructor(molecule,
                fragment_atoms, 0)

        # Initialize the energy list and SCF procedure employing newton-raphson algorithm
        energy = []
        chemical_potential = 0.0
        chemical_potential = scipy.optimize.newton(self._oneshot_loop, chemical_potential, args = (orbitals, orb_list, orb_list2, energy))

        # Get the final energy value
        niter = len(energy)
        dmet_energy = energy[niter-1]

        if self.verbose:
            print(' \t*** DMET Cycle Done *** ')
            print(' \tDMET Energy ( a.u. ) = '+'{:17.10f}'.format(dmet_energy))
            print(' \tChemical Potential   = '+'{:17.10f}'.format(chemical_potential))

        return dmet_energy


    def _oneshot_loop(self, chemical_potential, orbitals, orb_list, orb_list2, energy_list):
        onerdm_low = helpers._low_rdm(orbitals.active_fock, orbitals.number_active_electrons)

        niter = len(energy_list)+1

        if self.verbose:
            print(" \tIteration = ", niter)
            print(' \t----------------')
            print(' ')

        number_of_electron = 0.0
        energy_temp = 0.0

        for i, norb in enumerate(orb_list):
            if self.verbose:
                print("\t\tFragment Number : # ", i+1)
                print('\t\t------------------------')

            t_list=[]
            t_list.append(norb)
            temp_list = orb_list2[i]

            # Construct bath orbitals
            bath_orb, e_occupied = helpers._fragment_bath(orbitals.mol_full, t_list, temp_list, onerdm_low)

            # Obtain one particle rdm for a fragment
            norb_high, nelec_high, onerdm_high = helpers._fragment_rdm(t_list, bath_orb, e_occupied, orbitals.number_active_electrons)

            # Obtain one particle rdm for a fragment
            one_ele, fock, two_ele = orbitals.dmet_fragment_hamiltonian(bath_orb, norb_high, onerdm_high)

            # Construct guess orbitals for fragment SCF calculations
            guess_orbitals = helpers._fragment_guess(t_list, bath_orb, chemical_potential, norb_high, nelec_high, orbitals.active_fock)

            # Carry out SCF calculation for a fragment
            mf_fragment, fock_frag_copy, mol_frag = helpers._fragment_scf(t_list, two_ele, fock, nelec_high, norb_high, guess_orbitals, chemical_potential)

            energy = self.electronic_structure_solver.simulate(mol_frag, mf_fragment)
            cc_onerdm, cc_twordm = self.electronic_structure_solver.get_rdm()

            # Compute the fargment energy
            fragment_energy, total_energy_rdm, one_rdm = self._compute_energy(mf_fragment, cc_onerdm, cc_twordm, fock_frag_copy, t_list, one_ele, two_ele, fock)

            # Sum up the energy
            energy_temp += fragment_energy

            # Sum up the number of electrons
            number_of_electron += np.trace(one_rdm[ : t_list[0], : t_list[0]])

            if self.verbose:
                print("\t\tFragment Energy                 = "+'{:17.10f}'.format(fragment_energy))
                print("\t\tNumber of Electrons in Fragment = "+'{:17.10f}'.format(np.trace(one_rdm)))
                print('')


        energy_temp += orbitals.core_constant_energy
        energy_list.append(energy_temp)

        return number_of_electron - orbitals.number_active_electrons

    def _compute_energy(self, mf_frag, onerdm, twordm, fock_frag_copy, t_list, oneint, twoint, fock):
        """
        Calculate fragment energy
        :param mf_frag: mean field object (for fragment) of pyscf
        :param cc_onerdm: one-particle reduced density matrix
        :param cc_twordm: two-particle reduced density matrix
        :param fock_frag_copy: Fock matrix with the chemical potential subtracted
        :param t_list: Number of fragment & bath orbitals
        :param oneint: One-electron integrals for the fragment
        :param twoint: Two-electron integrals for the fragment
        :param fock: The fock matrix for the fragment
        :return: Fragment energy and reduced density matrix from CCSD calculation
        """
        # Execute CCSD calculation
        norb = t_list[0]

        # Calculate the one- and two- RDM for DMET energy calculation (Transform to AO basis)
        one_rdm = reduce(np.dot, (mf_frag.mo_coeff, onerdm, mf_frag.mo_coeff.T))

        twordm = np.einsum('pi,ijkl->pjkl', mf_frag.mo_coeff, twordm)
        twordm = np.einsum('qj,pjkl->pqkl', mf_frag.mo_coeff, twordm)
        twordm = np.einsum('rk,pqkl->pqrl', mf_frag.mo_coeff, twordm)
        twordm = np.einsum('sl,pqrl->pqrs', mf_frag.mo_coeff, twordm)

        # Calculate the total energy based on RDMs
        total_energy_rdm = np.einsum('ij,ij->', fock_frag_copy, one_rdm ) + 0.5 * np.einsum('ijkl,ijkl->', twoint, twordm )

        # Calculate fragment expectation value
        fragment_energy_one_rdm = 0.25  * np.einsum('ij,ij->',     one_rdm[ : norb, :      ], fock[: norb, :     ] + oneint[ : norb, :     ]) \
                               + 0.25  * np.einsum('ij,ij->',     one_rdm[ :     , : norb ], fock[:     , : norb] + oneint[ :     , : norb])

        fragment_energy_twordm = 0.125 * np.einsum('ijkl,ijkl->', twordm[ : norb, :     , :     , :      ], twoint[ : norb, :     , :     , :     ]) \
                               + 0.125 * np.einsum('ijkl,ijkl->', twordm[ :     , : norb, :     , :      ], twoint[ :     , : norb, :     , :     ]) \
                               + 0.125 * np.einsum('ijkl,ijkl->', twordm[ :     , :     , : norb, :      ], twoint[ :     , :     , : norb, :     ]) \
                               + 0.125 * np.einsum('ijkl,ijkl->', twordm[ :     , :     , :     , : norb ], twoint[ :     , :     , :     , : norb])


        fragment_energy = fragment_energy_one_rdm + fragment_energy_twordm
        return fragment_energy, total_energy_rdm, one_rdm

