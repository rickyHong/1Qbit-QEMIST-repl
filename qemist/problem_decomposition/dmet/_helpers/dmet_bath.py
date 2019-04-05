"""
A code for bath orbital construction for DMET calculation
"""

import numpy as np


def dmet_fragment_bath(mol, t_list, temp_list, onerdm_low):
    """
    Construct the bath orbitals for DMET fragment calculation
    :param mol: The molecule object from PySCF
    :param t_list: Number of fragment & bath orbitals
    :param temp_list: List of the minimum and maximum number for active orbitals
    :param onerdm_low: The one-particle RDM from the low-level calculation
    :return: bath orbitals, and orbital energies
    """

    # Extract the one-particle RDM for the active space
    onerdm_embedded = dmet_onerdm_embed(mol, temp_list, onerdm_low)

    # Diagonalize it
    e, c = np.linalg.eigh(onerdm_embedded)

    # Sort the eigenvectors with the eigenvalues
    e_sorted, c_sorted = dmet_bath_orb_sort(t_list, e, c)

    # Add the core contribution
    bath_orb, e_core = dmet_add_to_bath_orb(mol, t_list, temp_list, e_sorted, c_sorted)

    return bath_orb, e_core

def dmet_onerdm_embed(mol, temp_list, onerdm_before):
    """
    Extract the one particle RDM of the active space
    :param mol: The molecule object from PySCF
    :param temp_list: List of the minimum and maximum number for active orbitals
    :param onerdm_before: The one-particle RDM from the low-level calculation
    :return: The extracted one particle RDM
    """

    # Get the number of orbitals
    norbital_total = mol.nao_nr()

    # Reshape the RDM
    onerdm_matrix = np.reshape(onerdm_before, (norbital_total, norbital_total))

    if temp_list[0] == 0:
        # If it is the first fragment, just determine the maximum for extraction
        onerdm_temp  = onerdm_matrix[ : , temp_list[1]: ]
        onerdm_temp3 = onerdm_temp[temp_list[1]: , : ]
    else:
        # Determine the minimum and maximum orbitals for extraction
        onerdm_temp  = onerdm_matrix[ : , : temp_list[0]]
        onerdm_temp2 = onerdm_matrix[ : , temp_list[1]: ]
        onerdm_temp3 = np.hstack((onerdm_temp, onerdm_temp2))
        onerdm_temp  = onerdm_temp3[ : temp_list[0], : ]
        onerdm_temp2 = onerdm_temp3[temp_list[1]: , : ]
        onerdm_temp3 = np.vstack((onerdm_temp, onerdm_temp2))

    return onerdm_temp3        

def dmet_bath_orb_sort(t_list, e_before, c_before):
    """
    Sort the bath orbitals with the eigenvalues (orbital energies)
    :param t_list: Number of fragment & bath orbitals
    :param e_before: The orbitals before sort
    :param c_before: The orbitals energies before sort
    :return: The sorted eigenvalues and eigenvectors
    """

    # Sort the orbital energies (Occupation of 1.0 should come first...)
    new_index = np.maximum(-e_before, e_before - 2.0).argsort()

    # Throw away some orbitals above threshold
    thresh_orb = np.sum(-np.maximum(-e_before, e_before - 2.0)[new_index] > 1e-13)

    # Determine the number of bath orbitals
    norb = min(np.sum(thresh_orb), t_list[0])
    t_list.append(norb)

    # Sort the bath orbitals with its energies
    e_new = e_before[new_index]
    c_new = c_before[ : , new_index]
    
    return e_new, c_new      

def dmet_add_to_bath_orb( mol, t_list, temp_list, e_before, c_before ):
    """
    Add the frozen core part to the bath orbitals
    :param mol: The molecule object from PySCF
    :param t_list: Number of fragment & bath orbitals
    :param temp_list: List of the minimum and maximum number for active orbitals
    :param e_before: energy before addition of frozen core
    :param c_before: orbitals before addition of frozen core
    :return: the final bath orbitals and orbital energies
    """

    # Copy the bath orbitals and energies be fore adding the core
    add_e = - e_before[t_list[1]: ]
    add_c = c_before[ : , t_list[1]: ]
    new_index = add_e.argsort()

    # Sort the orbitals based on its energies
    c_before[ : , t_list[1]: ] = add_c[ : , new_index]
    add_e = - add_e[new_index]

    # The orbital energies with core part
    norbital_total = mol.nao_nr()
    e_occupied_core_orbitals = np.hstack((np.zeros([t_list[0] + t_list[1]]), add_e))

    # Add the core part in the orbitals
    for orb in range(0, t_list[0]):
        c_before = np.insert(c_before, orb, 0.0, axis=1) 
    i_temp = 0
    for orb_total in range( 0, norbital_total ):
        if ((orb_total >= temp_list[0]) and (orb_total < temp_list[1])):
            c_before = np.insert(c_before, orb_total, 0.0, axis=0) 
            c_before[orb_total, i_temp] = 1.0
            i_temp += 1
    
    return c_before, e_occupied_core_orbitals

