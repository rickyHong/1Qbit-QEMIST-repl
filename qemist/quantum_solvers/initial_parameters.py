def mp2_intitial_amplitudes(molecule, mean_field):
    """ Computes and prepares the MP2 inital amplitudes.

    Computed the inital ampltudes with PySCF, then reorders the elements
    into the QEMIST convention.
    """
    from pyscf import mp
    import numpy as np
    import itertools
    mp2_fragment = mp.MP2(mean_field)
    mp2_fragment.verbose = 0
    mp2_correlation_energy, mp2_t2 = mp2_fragment.kernel()
    scf_energy = mean_field.e_tot
    mp2_total_energy = scf_energy + mp2_correlation_energy

    # Reordering the amplitudes into QEMIST convention
    n_spatial_orbitals = len(mean_field.mo_energy)
    n_occupied = int(np.ceil(molecule.nelectron / 2))
    n_virtual = n_spatial_orbitals - n_occupied

    singles = []
    doubles_1 = []
    doubles_2 = []

    icount = 0

    # Get singles and doubles amplitudes associated with one
    # spatial occupied-virtual pair
    for p, q in itertools.product(range(n_virtual), range(n_occupied)):
        virtual_spatial = p
        occupied_spatial = q

        # Get singles amplitude
        # Just get up amplitude, since down should be the same
        singles.append(2.0e-5)

        temp=-mp2_t2[q,q,p,p]/2.0
        if(abs(temp)<1.0e-15):
            doubles_1.append(0.0)
        else:
            doubles_1.append(-mp2_t2[q,q,p,p]/2.0)

    # Get doubles amplitudes associated with two spatial occupied-virtual pairs
    for (p, q), (r, s) in itertools.combinations(
            itertools.product(range(n_virtual), range(n_occupied)), 2):
        # Get indices of spatial orbitals
        virtual_spatial_1 = p
        occupied_spatial_1 = q
        virtual_spatial_2 = r
        occupied_spatial_2 = s

        # Get amplitude
        doubles_2.append(-mp2_t2[occupied_spatial_1,occupied_spatial_2,virtual_spatial_1,virtual_spatial_2])

    return singles + doubles_1 + doubles_2

