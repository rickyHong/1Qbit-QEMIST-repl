"""
A code for constructing one-particle RDM for DMET calculation
"""
import numpy as np
from functools import reduce

def dmet_low_rdm(active_fock, number_active_electrons):
    """
    Calculate one-particle RDM from low-level calculation results
    :param active_fock: fock from low level calculation
    :param number_active_electrons: number of electrons in the entire system
    :return: one-particle RDM
    """

    # Extract the occupied part of the one-particle RDM
    num_occ = number_active_electrons / 2
    e, c = np.linalg.eigh(active_fock)
    new_index = e.argsort()
    e = e[new_index]
    c = c[ : , new_index]
    onerdm = np.dot(c[ : , : int(num_occ)], c[ : , : int(num_occ)].T) * 2

    return onerdm

def dmet_fragment_rdm(t_list, bath_orb, core_occupied, number_active_electrons):
    """
    construct the one-particle RDM for the core orbitals
    :param t_list: Number of fragment & bath orbitals
    :param bath_orb: Bath orbitals
    :param core_occupied: Core occupied part of the MO coefficients
    :param number_active_electrons: Number of active electrons in the entire system
    :return: number of orbitals, number of electrons in the fragments, and the core RDM
    """

    # Obtain number of active orbitals
    number_orbitals = t_list[0] + t_list[1]

    # Round the values above or below threshold
    for i, core in enumerate(core_occupied):
        if (core < 0.01):
            core_occupied[i] = 0.0
        elif (core > 1.99):
            core_occupied[i] = 2.0

    # Define the number of electrons in the fragment
    number_ele_temp = np.sum(core_occupied)
    number_electrons = int(round(number_active_electrons - number_ele_temp))

    # Obtain the one particle RDM for the fragment (core)
    core_occupied_onerdm = reduce(np.dot, (bath_orb, np.diag(core_occupied), bath_orb.T))
    
    return number_orbitals, number_electrons, core_occupied_onerdm

