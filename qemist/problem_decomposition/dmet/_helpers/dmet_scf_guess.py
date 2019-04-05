"""
A code for guess orbital construction for fragment SCF calculation in DMET calculation
"""
import scipy
import numpy as np
from functools import reduce

def dmet_fragment_guess(t_list, bath_orb, chemical_potential, norb_high, number_active_electron, active_fock):
    """
    Construct the guess orbitals for fragment SCF calculation
    :param t_list: Number of fragment & bath orbitals
    :param bath_orb: The bath orbitals
    :param chemical_potential: The chemical potential of the previous iteration
    :param norb_high: The number of orbitals in the fragment
    :param number_active_electron: The number of active electrons in the fragment
    :param active_fock: The fock matrix from the low-level calcualtion
    :return: The guess orbitals
    """

    # Construct the fock matrix of the fragment (subtract the chemical potential for consistency)
    fock_fragment = reduce(np.dot, (bath_orb[ : , : norb_high].T, active_fock, bath_orb[ : , : norb_high]))
    norb = t_list[0]
    if(chemical_potential != 0):
        for i in range(norb):
            fock_fragment[i, i] -= chemical_potential

    # Diagonalize the fock matrix and get the eigenvectors
    eigenvalues, eigenvectors = scipy.linalg.eigh(fock_fragment)
    eigenvectors = eigenvectors[ : , eigenvalues.argsort()]

    # Extract the eigenvectors of the the occupied orbitals) as the guess orbitals
    frag_guess = np.dot(eigenvectors[ :, : int(number_active_electron/2)], eigenvectors[ :, : int(number_active_electron/2)].T) * 2
    
    return frag_guess

