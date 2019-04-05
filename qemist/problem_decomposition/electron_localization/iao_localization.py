from pyscf import gto
from pyscf.lo import iao
from functools import reduce
from pyscf.lo import orth
import numpy as np
import scipy

def iao_localization(mol, mf):
    """
    Construct IAO (Localize the canonical orbitals obtained from
    low-level SCF calculations)for DMET Calculation!!!
    :param mol: The molecule object from PySCF
    :param mf: The mean field object from PySCF
    :return: The localized Intrinsic atomic orbitals
    """

    #   Construct IAO from occupied orbitals
    iao1 = _iao_occupied_orbitals(mol, mf)

    #   Construct IAO from complementary space
    iao2 = _iao_complementary_orbitals(mol, iao1)

    #   Gather two and Assign the IAOs to atoms, rearrange them
    iao_lo = _iao_atoms(mol, iao1, iao2)

    return iao_lo


def _iao_occupied_orbitals(mol, mf):
    """
    Get the IAO for occupied space
    :param mol: The molecule object from PySCF
    :param mf: The mean field object from PySCF
    :return: IAO for occupied space
    """

    #   Get MO coefficient of occupied MOs
    occupied_orbitals = mf.mo_coeff[:, mf.mo_occ > 0.5]

    #   Get mol data in minao basis
    min_mol = iao.reference_mol(mol)

    #   Calculate the Overlaps for total basis
    s1 = mol.intor_symmetric('int1e_ovlp')

    #   ... for minao basis
    s2 = min_mol.intor_symmetric('int1e_ovlp')

    #   ... between the two basis (and transpose)
    s12 = gto.mole.intor_cross('int1e_ovlp', mol, min_mol)
    s21 = s12.T

    #   Calculate P_12 = S_1^-1 * S_12 using Cholesky decomposition
    s1_sqrt = scipy.linalg.cho_factor(s1)
    s2_sqrt = scipy.linalg.cho_factor(s2)
    p12 = scipy.linalg.cho_solve(s1_sqrt, s12)

    #   C~ = second_half ( S_1^-1 * S_12 * first_half ( S_2^-1 * S_21 * C ) )
    c_tilde = scipy.linalg.cho_solve(s2_sqrt, np.dot(s21, occupied_orbitals))
    c_tilde = scipy.linalg.cho_solve(s1_sqrt, np.dot(s12, c_tilde))
    c_tilde = np.dot(c_tilde, orth.lowdin(reduce(np.dot, (c_tilde.T, s1, c_tilde))))

    #   Obtain C * C^T * S1 and C~ * C~^T * S1
    ccs1 = reduce(np.dot, (occupied_orbitals, occupied_orbitals.conj().T, s1))
    ctcts1 = reduce(np.dot, (c_tilde, c_tilde.conj().T, s1))

    #   Calculate A = ccs1 * ctcts1 * p12 + ( 1 - ccs1 ) * ( 1 - ctcts1 ) * p12
    iao_active = (p12 + reduce(np.dot, (ccs1, ctcts1, p12)) * 2 - np.dot(ccs1, p12) - np.dot(ctcts1, p12))

    #   Orthogonalize A
    iao_active = np.dot(iao_active, orth.lowdin(reduce(np.dot, (iao_active.T, s1, iao_active))))

    return iao_active


def _iao_complementary_orbitals(mol, iao_ref):
    """
    Get IAO for complementary space (Virtual Orbitals)
    :param mol: The molecule object from PySCF
    :param mf: The mean field object from PySCF
    :param iao_ref: The IAO for occupied
    :return:
    """

    #   Get the total number of AOs
    norbital_total = mol.nao_nr()

    #   Calculate the Overlaps for total basis
    s1 = mol.intor_symmetric('int1e_ovlp')

    #   Construct the complementary space AO
    number_iaos = iao_ref.shape[1]
    number_inactive = norbital_total - number_iaos
    iao_com_ref = _iao_complementary_space(iao_ref, s1, number_inactive)

    #   Get a list of active orbitals
    min_mol = iao.reference_mol(mol)
    norbital_active, active_list = _iao_count_active(mol, min_mol)

    #   Obtain the Overlap-like matrices
    s21 = s1[active_list, : ]
    s2 = s21[ : , active_list]
    s12 = s21.T

    #   Calculate P_12 = S_1^-1 * S_12 using Cholesky decomposition
    s1_sqrt = scipy.linalg.cho_factor(s1)
    s2_sqrt = scipy.linalg.cho_factor(s2)
    p12 = scipy.linalg.cho_solve(s1_sqrt, s12)

    #   C~ = orth ( second_half ( S_1^-1 * S_12 * first_half ( S_2^-1 * S_21 * C ) ) )
    c_tilde = scipy.linalg.cho_solve(s2_sqrt, np.dot(s21, iao_com_ref))
    c_tilde = scipy.linalg.cho_solve(s1_sqrt, np.dot(s12, c_tilde))
    c_tilde = np.dot(c_tilde, orth.lowdin(reduce(np.dot, (c_tilde.T, s1, c_tilde))))

    #   Obtain C * C^T * S1 and C~ * C~^T * S1
    ccs1 = reduce(np.dot, (iao_com_ref, iao_com_ref.conj().T, s1))
    ctcts1 = reduce(np.dot, (c_tilde, c_tilde.conj().T, s1))

    #   Calculate A = ccs1 * ctcts1 * p12 + ( 1 - ccs1 ) * ( 1 - ctcts1 ) * p12
    iao_comp = (p12 + reduce(np.dot, (ccs1, ctcts1, p12)) * 2 - np.dot(ccs1, p12) - np.dot(ctcts1, p12))
    iao_comp = np.dot(iao_comp, orth.lowdin(reduce(np.dot, (iao_comp.T, s1, iao_comp))))

    return iao_comp


def _iao_count_active(mol, min_mol):
    """
    Figure out the basis functions matching with MINAO
    :param mol: The molecule object for the basis set used
    :param min_mol: The molecule object for the minimal basis MINAO
    :return: Number of active orbitals, and the list of active orbitals
    """

    #   Initialize the list
    active_number_list = []

    #   Loop over all basis and see if there are labels matching with the MINAO ones
    itemp = 0
    for total_basis in mol.spheric_labels():
        count_yes = 0
        for min_basis in min_mol.spheric_labels():
            if min_basis == total_basis:
                count_yes = 1
        if count_yes == 0:
            active_number_list.append(itemp)

        itemp += 1

    #   Make the list a numpy array
    number_active = len(active_number_list)
    active_number_list = np.array(active_number_list)

    return number_active, active_number_list


def _iao_complementary_space(iao_ref, s, number_inactive):
    """
    Determine the complementary space
    :param iao_ref: IAO for occupied space
    :param s: the overlap matrix
    :param number_inactive: The number of inactive (core) orbitals
    :return: the inactive part of the MO coefficients in IAO
    """
    #   Construct the "density matrix" for active space
    density_active = np.dot(iao_ref, iao_ref.T)

    #   Get the MO Coefficient from the IAO density matrix
    a_mat = reduce(np.dot, (s, density_active, s))
    eigval, eigvec = scipy.linalg.eigh(a=a_mat, b=s)

    #   Extract inactive part of "MO Coefficient" and return it
    eigen_vectors = eigvec[:, : number_inactive]

    return eigen_vectors


#
def _iao_atoms(mol, iao1, iao2):
    """
    Assign IAO to atom centers and rearrange the IAOs.
    :param mol: The molecule object from PySCF
    :param mf: The mean field object from PySCF
    :param iao1: the IAO from occupied
    :param iao2: the IAO from virtual
    :return: the rearranged IAOs
    """

    # Calclate the integrals for assignment
    number_orbitals = mol.nao_nr()
    r_int1e = mol.intor('cint1e_r_sph', 3)
    iao_combine = np.hstack((iao1, iao2))

    # Calculate xyz for each orbital
    x = np.diag(reduce(np.dot,(iao_combine.T, r_int1e[0], iao_combine)))
    y = np.diag(reduce(np.dot,(iao_combine.T, r_int1e[1], iao_combine)))
    z = np.diag(reduce(np.dot,(iao_combine.T, r_int1e[2], iao_combine)))

    # Align the coordinates
    orbitals_temp = np.vstack((x, y, z))
    orbitals = orbitals_temp.T

    # Assign each orbital to atom center
    atom_list = _dmet_atom_list(mol, orbitals)

    # Prepare the orbital labels
    orb_list = _dmet_orb_list(mol, atom_list)

    # Rearrange the orbitals
    iao_combine = iao_combine[ : , orb_list]

    # Orthogonalize the orbitals
    s1 = mol.intor_symmetric('int1e_ovlp')
    iao_combine = np.dot(iao_combine, orth.lowdin(reduce(np.dot, (iao_combine.T, s1, iao_combine))))

    return iao_combine


def _dmet_atom_list(mol, orbitals):
    """
    Assign each orbital to an atom
    :param mol: The molecule object from PySCF
    :param orbitals: Coordinates for the orbitals centers
    :return: The list for atom assignment of the IAOs
    """

    # Initialize the list
    number_orbitals = mol.nao_nr()
    newlist = []

    # Calculate the distance from atom centers and determine the nearest
    for i in range(number_orbitals):
        i_temp = 0
        distance_temp = scipy.linalg.norm(orbitals[i, :] - mol.atom_coord(0))
        for j in range(1, mol.natm):
            distance = scipy.linalg.norm(orbitals[i, :] - mol.atom_coord(j))
            if (distance < distance_temp):
                distance_temp = distance
                i_temp = j
            else:
                pass
        newlist.append(i_temp)

    return newlist


def _dmet_orb_list(mol, atom_list):
    """
    Rearrange the orbital label
    :param mol: The molecule object from PySCF
    :param atom_list: The atom list of IAO assignment
    :return: The orbital list in new order
    """
    newlist = []
    for i in range(mol.natm):
        for j in range(mol.nao_nr()):
            if (atom_list[j] == i):
                newlist.append(j)

    return newlist
