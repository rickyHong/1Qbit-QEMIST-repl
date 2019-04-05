"""
Construct the orbital list for the fragments
from the atom list
"""
def dmet_fragment_constructor( mol, atom_list, number_fragment ):
    """
    Make a list of number of orbitals for each fragment
    while obtaining the list if we consider combining
    fragments
    :param mol: The molecule object from PySCF
    :param atom_list: The list describing number of atoms for each fragment
    :param number_fragment: The number of fragment for each DMET calculation
    :return: The lists for number of orbitals and the list defining the orbital space for each DMET calculation
    """

    # Make a new atom list based on how many fragments for DMET calculation
    if number_fragment == 0 :
        atom_list2 = atom_list
    else:
        # Calculate the number of DMET calculations
        number_new_fragment = int(len(atom_list)/(number_fragment+1)) # number of DMET calulation per loop
        atom_list2 = []

        # Define the number of atoms per DMET calculation
        for i in range( number_new_fragment ):
            num = 0
            for j in range( number_fragment + 1 ):
                k = (number_fragment+1)*i+j
                num += atom_list[ k ]
            atom_list2.append(num)

    # Initialize the list of the number of orbitals
    orb_list = []
    orb_list2 = []
    isum = 0
    isum2 = -1
    iorb = 0
    jorb = 0

    # Calculate the number of orbitals for each atom
    for i in atom_list2 : 
        itemp = 0
        isum2 += i
        for total_basis in mol.spheric_labels():
            item = total_basis.split()
            item0 = int(item[0])
            if ( ( item0 >= isum ) and ( item0 <= isum2 ) ):
               itemp+=1
        isum += i 
        jorb += itemp
        orb_list.append(itemp)
        orb_list2.append([iorb,jorb])
        iorb += itemp

    return orb_list, orb_list2, atom_list2
