"""
Carry out SCF calculation for DMET calculation
"""
from pyscf import gto, scf, ao2mo
from functools import reduce
import numpy as np
import scipy

def dmet_fragment_scf( t_list, two_ele, fock, number_electrons, number_orbitals, guess_orbitals, chemical_potential ):
    """
    Carry out SCF calculation for the fragment
    :param t_list: Number of fragment & bath orbitals
    :param two_ele: Two-electron integrals for the fragment
    :param fock: The fock matrix for the fragment
    :param number_electrons: Number of electrons for the fragment
    :param number_orbitals: Number of orbitals for the fragment
    :param guess_orbitals: The guess orbitals for SCF calculation
    :param chemical_potential: The chemical potential of the previous iteration
    :return: mf_frag (mean field object (for fragment) of pyscf, the fock matrix subtracting the chemical potential, and mol_frag(molecule object (for fragment) of pyscf)
    """
    # Deep copy the fock matrix
    fock_frag_copy = fock.copy()

    # Subtract the chemical potential to make the number of electrons consistent
    if (chemical_potential != 0.0):
        for orb in range(t_list[0]):
            fock_frag_copy[orb, orb] -= chemical_potential

    # Determine the molecular space (set molecule object of pyscf)
    mol_frag = gto.Mole()
    mol_frag.build(verbose=0)
    mol_frag.atom.append(('C', (0, 0, 0)))
    mol_frag.nelectron = number_electrons
    mol_frag.incore_anyway = True

    # Perform SCF calculation (set mean field object of pyscf)
    mf_frag = scf.RHF(mol_frag)
    mf_frag.get_hcore = lambda *args: fock_frag_copy
    mf_frag.get_ovlp = lambda *args: np.eye(number_orbitals)
    mf_frag._eri = ao2mo.restore(8, two_ele, number_orbitals)
    mf_frag.scf(guess_orbitals)

    # Calculate the density matrix for the fragment
    dm_frag = reduce(np.dot, (mf_frag.mo_coeff, np.diag(mf_frag.mo_occ), mf_frag.mo_coeff.T))

    # Use newton-raphson algorithm if the above SCF calculation is not converged
    # print("SCF Converged ? = ", mf_frag.converged)
    if ( mf_frag.converged == False ):
        mf_frag.get_hcore = lambda *args: fock_frag_copy
        mf_frag.get_ovlp = lambda *args: np.eye(number_orbitals)
        mf_frag._eri = ao2mo.restore(8, two_ele, number_orbitals)
        mf_frag = scf.RHF(mol_frag).newton()
        energy = mf_frag.kernel(dm0 = dm_frag)
        dm_frag = reduce(np.dot, (mf_frag.mo_coeff, np.diag(mf_frag.mo_occ), mf_frag.mo_coeff.T))
    
    # print("SCF Energy = ", mf_frag.e_tot)

    # Calculate the density matrix
    number_occupied = int(number_electrons/2)
    fock_fragment = fock_frag_copy + np.einsum('ijkl,ij->kl', two_ele, dm_frag ) - 0.5 * np.einsum('ijkl,ik->jl', two_ele, dm_frag )
    eigenvalues, eigenvectors = scipy.linalg.eigh(fock_fragment)
    inew_index = eigenvalues.argsort()
    eigenvalues = eigenvalues[inew_index]
    eigenvectors = eigenvectors[ :, inew_index]
    dm_frag2 = np.dot(eigenvectors[:, :number_occupied], eigenvectors[:, :number_occupied].T) * 2

    # Check the difference between the RDMs obtained in two ways
    # print(" Density difference = ", scipy.linalg.norm(dm_frag - dm_frag2))
    
    return mf_frag, fock_frag_copy, mol_frag

