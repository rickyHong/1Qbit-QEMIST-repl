#   Copyright 2019 1QBit
#   
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import itertools

def complex_as_dict(re, im):
    import cmath
    x = complex(re, im)
    return {'Real': re, 'Imaginary': im, 'Magnitude': abs(x), 'Phase': cmath.phase(x)}

def alpha_spinorbital(MO):
    ''' Return the corresponding alpha spinorbital index given a molecular orbital (MO) index
    Args:
        MO(int):     molecular orbital index
    Returns:
        (2*MO)(int): alpha spin-orbital index
    '''
    return 2*MO

def beta_spinorbital(MO):
    ''' Return the corresponding beta spinorbital index given a molecular orbital (MO) index
    Args:
        MO(int):         molecular orbital index
    Returns:
        (2*MO + 1)(int): beta spin-orbital index
    '''
    return 2 * MO + 1

def count_amplitudes(n_spinorbitals,n_electrons):
    """Count the number of singles and doubles amplitudes for a given UCCSD-VQE run
    Args:
        n_spinorbitals(int): integer representing the number of spinorbitals (qubits) for
                             a given molecule and basis set
        n_electrons(int):    integer representing the number of electrons for a given molecule
    Returns:
        n_amplitudes(int):   integer representing the total number of amplitudes (MO basis)
    """
    # Compute the number of MOs and the number of occupied and virtual MOs
    n_MO = n_spinorbitals // 2
    n_occ = n_electrons // 2
    n_virt = n_MO - n_occ
    # Compute the number of singles and doubles amplitudes
    n_singles = n_occ*n_virt
    n_doubles = n_singles*(n_singles + 1) // 2
    return n_singles + n_doubles

def compute_cluster_operator(n_spinorbitals, n_electrons, amplitudes, multiply=False, operator = []):
    """Compute or update the cluster operator for a given UCCSD-VQE run


        n_spinorbitals(int): integer representing the number of spinorbitals (qubits) for
                             a given molecule and basis set
        n_electrons(int):    integer representing the number of electrons for a given molecule
        amplitudes(list):    list of the amplitudes, with the singles appearing first, followed
                             by the diagonal (i,i,a,a) doubles and then the off-diagonal (i,j,a,b)
                             doubles
        multiply(bool):      optional boolean to indicate whether we are performing an amplitude 
                             update (i.e. multiplying a new set of amplitudes by the corresponding 
                             operators) or not
        operator(list):      optional list of the contributions to the cluster operator 

    Returns:
        ref(tuple):          tuple of tuples representing the reference configuration
        t(list):             list of tuples representing the cluster operator
    """

    # Compute the number of MOs and the number of occupied and virtual MOs
    n_MO = n_spinorbitals // 2
    n_occ = n_electrons // 2
    n_virt = n_MO - n_occ

    # Extract the singles amplitudes and diagonal doubles amplitudes
    singles = amplitudes[:(n_occ*n_virt)]
    doubles_diag = amplitudes[(n_occ*n_virt):(2*n_occ*n_virt)]
    doubles_offdiag = amplitudes[(2*n_occ*n_virt):]

    # Spin indexing
    spin_index = [alpha_spinorbital, beta_spinorbital]

    # -------------------------------------------------------
    # Reference configuration
    # -------------------------------------------------------

    #Loop over occupied orbitals
    li = []
    j = 0
    for i in range(n_occ):
        # Define alpha and beta spinorbitals
        i_a = alpha_spinorbital(i)
        i_b = beta_spinorbital(i)

        if(multiply):
            li += [i_a,i_b]
        else:
            li += [('u',i_a), ('u', i_b)]
    # Define the reference state
    if(multiply):
        ref = ((1.0,0.0),li)
        #ref = (li,1)
    else:
        #ref = ((1.0, 0.0),(li,1))
        ref = (li,1)

    t = []

    # -------------------------------------------------------
    # Single excitations and diagonal double excitations
    # -------------------------------------------------------

    # Loop over occupied and virtual orbitals
    for i, (m,n) in enumerate(itertools.product(range(n_virt),range(n_occ))):
     
        # n labels virtual orbitals (offset of n_occ)
        m += n_occ

        # Loop over spin
        for spin in range(2):

            # Mapping of spatial orbitals to spin-orbitals
            ind_1 = spin_index[spin]
            ind_2 = spin_index[1-spin]

            # Define spin-orbital labels
            m_1 = ind_1(n)
            m_2 = ind_2(n)
            n_1 = ind_1(m)
            n_2 = ind_2(m)

            if(multiply):

                # Multiply the singles operators by the correct amplitudes
                t += [((singles[i], 0.0), operator[j][1])]
                #t += [(operator[j][1], complex_as_dict(singles[i],0.0))]
                j += 1
                
                # Multiply the diagonal doubles operators by the correct amplitudes
                if(m_1 != m_2 and n_1 != n_2):
                    t += [((doubles_diag[i], 0.0), operator[j][1])]
                    j += 1
                
            else:           

                # Generate the singles excitations in the proper format
                t += [(([('u',n_1),('d',m_1)],1), complex_as_dict(singles[i],0.0))]

                # Generate the diagonal doubles excitations in the proper format
                if(m_1 != m_2 and n_1 != n_2):
                    t += [(([('u',n_1),('u',n_2),('d',m_1),('d',m_2)],1), complex_as_dict(doubles_diag[i],0.0))]

    # Loop over unique off-diagonal doubles
    for i, ((m,u),(n,v)) in enumerate(itertools.combinations(
            itertools.product(range(n_virt),range(n_occ)),2)
            ):

        # m and n label virtual orbitals (offset of n_occ)
        m += n_occ
        n += n_occ

        # Loop over spin
        for (spin_1, spin_2) in itertools.product(range(2), repeat=2):

            # Mapping of spatial orbitals to spin-orbitals
            ind_1 = spin_index[spin_1]
            ind_2 = spin_index[spin_2]

            # Define spin-orbital labels
            m_1 = ind_1(u)
            m_2 = ind_2(v)
            n_1 = ind_1(m)
            n_2 = ind_2(n)

            if(multiply):

                # Multiply the off-diagonal doubles excitation operators by the correct amplitudes
                if(m_1 != m_2 and n_1 != n_2):
                    t += [((1.0*doubles_offdiag[i],0.0), operator[j][1])]
                    j += 1
            else:

                # Generate the off-diagonal doubles excitations in the proper format
               if(m_1 != m_2 and n_1 != n_2):
                   t += [(([('u',n_1),('u',n_2),('d',m_1),('d',m_2)],1), complex_as_dict(doubles_offdiag[i],0.0))]

    # If multiply is true, then append the reference configuration to t
    if (multiply): t += [ref]

    return ref, t
