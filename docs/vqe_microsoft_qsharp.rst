
The Variational Quantum Eigensolver (VQE)
=========================================

1. A brief overview of the VQE algorithm
----------------------------------------

The Variational Quantum Eigensolver (VQE)
[`Peruzzo_et_al.,_2014 <https://arxiv.org/abs/1304.3061>`__,
`McClean_et_al.,_2015 <https://arxiv.org/abs/1509.04279>`__] has been
introduced as a hybrid quantum–classical algorithm for simulating
quantum systems. Some examples of quantum simulation using VQE include
solving the molecular electronic Schrödinger equation and model systems
in condensed matter physics (e.g., Fermi– and Bose–Hubbard models). In
this notebook, we focus on VQE within the context of solving the
molecular electronic structure problem for the ground-state energy of a
molecular system. The second-quantized Hamiltonian of such a system,
within the Born-Oppenheimer approximation, assumes the following form:

.. raw:: latex

   \begin{equation}
   \hat{H} = h_{\text{nuc}} + \sum_{p,q} h^{p}_{q} \hat{a}^{\dagger}_p \hat{a}_q + \sum_{p,q,r,s} h^{pq}_{rs} \hat{a}^{\dagger}_p \hat{a}^{\dagger}_q \hat{a}_s \hat{a}_r\nonumber
   \end{equation}

Here, :math:`h_{\text{nuc}}` denotes the nuclear repulsion energy. The
coefficients :math:`h^{p}_{q}` and :math:`h^{pq}_{rs}` are obtained by
solving a mean-field problem. The Hamiltonian is then transformed into
the qubit basis (e.g., Jordan–Wigner, Bravyi–Kitaev). This means that it
is expressed entirely in terms of operators acting on qubits:

.. raw:: latex

   \begin{equation}
   \hat{H} = h_{\text{nuc}} + \sum_{\substack{p \\ \alpha}} h_{p}^{\alpha} \sigma_p^{\alpha} + \sum_{\substack{p,q \\ \alpha,\beta}} h_{pq}^{\alpha\beta}\sigma_p^{\alpha}\otimes\sigma_{q}^{\beta} + \sum_{\substack{p,q,r \\ \alpha,\beta,\gamma}}h_{pqr}^{\alpha\beta\gamma}\sigma_p^{\alpha}\otimes\sigma_{q}^{\beta}\otimes\sigma_r^{\gamma} + \ldots \nonumber
   \end{equation}

In this expression, the :math:`\sigma_p^\alpha` are Pauli matrices
(:math:`\alpha \in \{x,y,z\}`), acting on the :math:`p`-th qubit. We now
consider a trial wavefunction ansatz
:math:`\vert \Psi(\vec{\theta}) \rangle = U(\vec{\theta}) \vert 0 \rangle`
that depends on :math:`m` parameters defining
:math:`\vec{\theta}=(\theta_1, \theta_2, \ldots, \theta_m)`, which enter
a unitary operator that acts on the reference (i.e., mean-field) state
:math:`\vert 0 \rangle`. The variational principle dictates that we can
minimize the expectation value of the Hamiltonian,

.. raw:: latex

   \begin{equation}
   E = \min_{\vec{\theta}} \frac{\langle \Psi(\vec{\theta}) \vert \hat{H} \vert \Psi(\vec{\theta}) \rangle}{\langle \Psi(\vec{\theta}) \vert \Psi(\vec{\theta}) \rangle} \geq E_{\text{gs}}\nonumber
   \end{equation}

to determine the optimal set of variational parameters. The energy thus
computed will be an upper bound to the true ground-state energy
:math:`E_{\text{gs}}`. Once a suitable variational trial ansatz has been
chosen (e.g., a unitary coupled-cluster ansatz, a heuristic ansatz), we
must provide a suitable set of initial guess parameters. If our ansatz
is written in according to the second quantization picture, we must also
transform it into the qubit basis before proceeding. We must also apply
other approximations (e.g., Trotter–Suzuki) to render it amenable for
translation into a quantum circuit. The resulting qubit form of the
ansatz can then be translated into a quantum circuit and, thus, able to
be implemented on quantum hardware. Once the initial state has been
prepared using a quantum circuit, energy measurements are performed
using quantum hardware or an appropriate simulation tool. The energy
value obtained is the sum of the measurements of the expectation values
of each of the terms that contribute to the Hamiltonian (assuming the
wavefunction has been normalized to unity):

.. raw:: latex

   \begin{equation}
   E = \langle \Psi(\vec{\theta}) \vert \hat{H} \vert \Psi(\vec{\theta}) \rangle =\langle\hat{H}\rangle = h_{\text{nuc}} + \sum_{\substack{p \\ \alpha}} h_{p}^{\alpha} \langle\sigma_p^{\alpha}\rangle + \sum_{\substack{p,q \\ \alpha,\beta}} h_{pq}^{\alpha\beta}\langle\sigma_p^{\alpha}\otimes\sigma_{q}^{\beta}\rangle + \sum_{\substack{p,q,r \\ \alpha,\beta,\gamma}}h_{pqr}^{\alpha\beta\gamma}\langle\sigma_p^{\alpha}\otimes\sigma_{q}^{\beta}\otimes\sigma_r^{\gamma}\rangle + \ldots \nonumber
   \end{equation}

The computed energy is then input to a classical optimizer in order to
find a new set of variational parameters, which are then used to prepare
a new state (i.e., a quantum circuit) on the quantum hardware. The
process is repeated until convergence. The algorithm is illustrated
below.

.. image:: img/VQE_overview.png
    :align: center
    :width: 750pt

2. Computing the ground–state energy of H\ :math:`_{\text{2}}` with UCCSD-VQE
-----------------------------------------------------------------------------

The **Microsoft Quantum Development Kit (QDK)** provides a way to
simulate quantum circuits on classical hardware and quantum processors.
It uses the Microsoft Q# language, which was developed specifically to
handle hybrid quantum–classical workflows.

The `**Quantum Development Kit chemistry
library** <https://docs.microsoft.com/quantum/libraries/chemistry/>`__
provides key functionalities for tackling problems in quantum chemistry.
It is written in C#, and relies on Q# operations to implement various
quantum algorithms. This is an open source `GitHub
repository <https://github.com/Microsoft/QuantumLibraries>`__ that
accepts suggestions and contributions.

Although users are able to write and call their own code in Q# and C#,
this is not a requirement. This notebook uses Python exclusively. All
the functionality needed to execute the example that follows can be
accessed through the Quantum Development Kit Python interoperability
package for python, **qsharp**, available on pip. Further details about
this pip package are available at
https://docs.microsoft.com/quantum/install-guide/python .

This section shows how these functionalities can be used to compute the
ground state energy of H\ :math:`_{\mathrm{2}}` (the simplest molecule)
in a minimal basis set, using the unitary coupled-cluster ansatz with
single and double excitations (UCCSD) and compare the results with the
exact results obtained in this basis. A ball-stick model for
H\ :math:`_\text{2}` is shown below. The distance between the two
hydrogen atoms is called the bond length, and its value is set to
approximately 0.7414\ :math:`~`\ Å in this section.

.. image:: img/H2.png
    :align: center
    :width: 200pt


2.1 The Q# Python package
~~~~~~~~~~~~~~~~~~~~~~~~~

The cell below prepares the Q# environment and loads the useful
functionalities of the chemistry library through ``qsharp.chemistry``.
This notebook later details how each of these play a role in this
implementation of VQE.

.. code:: ipython3

    from qsharp.chemistry import load_broombridge, load_fermion_hamiltonian, load_input_state, encode

2.2 Input data
~~~~~~~~~~~~~~

Users need to provide quantities defining the target molecular system,
such as the following:

-  one- and two-electron integrals
-  nuclear repulsion energy

The use of VQE requires to specify extra input, such as the following:

-  the type of ansatz desired (UCCSD, for example)
-  the values for initial variational parameters
-  an initial state (a reference wavefunction, such as the Hartree–Fock
   wavefunction)

The `**Broombridge
format** <https://docs.microsoft.com/quantum/libraries/chemistry/schema/broombridge>`__,
created by Microsoft and PNNL, provides a way to store all the input
information in a human-readable .yaml file. Loading a pre-existing
Broombridge file containing the information of interest for the target
molecular system is the shortest way to get started with running VQE.
The ``qsharp.chemistry`` package allows Python packages to load and work
with the data stored in Broombridge files.

The following code snippet shows how to load existing data from a
Broombridge file (here for H\ :math:`_\text{2}` at a bond length of
0.7414), and explores the resulting data structure.

.. code:: ipython3

    # C# Chemistry library :: Loading molecular data (electronic integrals, etc.) from Broombridge                                                                                                    
    filename = 'data/hydrogen_0.2.yaml'                                                                                                                                                              
    broombridge_data =  load_broombridge(filename)

The data structure is easier to navigate when using a pretty-print
application or a proper IDE.

It is worth mentioning that users do not need a Broombridge file
describing the molecular system of interest in order to get started.
They could, for example, compute and provide their own data at runtime
using third-party libraries such as PySCF, and then be free to extract
and overwrite the information in the data structures produced by reading
any Broombridge file.

The instructions below show how users can read information stored in a
data structure (writing to the data structure is just as
straightforward).

.. code:: ipython3

    # Retrieve basis set and geometry used to generate the input data
    basis_set = broombridge_data.problem_description[0].basis_set
    geometry = broombridge_data.problem_description[0].geometry
    
    # Retrieve the nuclear repulsion and the one-electron integrals (Mulliken convention)
    nuclear_repulsion = broombridge_data.problem_description[0].coulomb_repulsion['Value']
    one_electron_integrals =  broombridge_data.problem_description[0].hamiltonian['OneElectronIntegrals']['Values']
    
    print("nuclear_repulsion = ", nuclear_repulsion)
    print("one_electron_integrals = ", one_electron_integrals)


.. parsed-literal::

    nuclear_repulsion =  0.713776188
    one_electron_integrals =  [([1, 1], -1.252477495), ([2, 2], -0.475934275)]


Note that users who are directly writing to the data structures should
be aware that the Python interop relies on JSON serialization, and
should use fundamental data types. They should make sure to pass lists
instead of NumPy arrays, or to cast their integer and floating point
values with the built-in **int** and **float** Python functions to avoid
JSON serialization errors at runtime.

2.3 Qubit Hamiltonian, UCCSD ansatz, and initial variational parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following section shows how to prepare the qubit Hamiltonian (also
referred to as the Pauli Hamiltonian) and access the information related
to one of the available ansatz for VQE: UCCSD.

Note that the underlying data structures may change in the future. The
code cells below encourage users to print their content by directly
accessing the available fields, exposed by the ``dir`` built-in Python
function.

The fermionic Hamiltonian can be built using the chemistry library, and
is returned to the Python context:

.. code:: ipython3

    ferm_hamiltonian = broombridge_data.problem_description[0].load_fermion_hamiltonian()
    print("ferm_hamiltonian ::", ferm_hamiltonian)
    print(dir(ferm_hamiltonian))


.. parsed-literal::

    ferm_hamiltonian :: <qsharp.chemistry.FermionHamiltonian object at 0x7f3e8c9a2940>
    ['__class__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'add_terms', 'system_indices', 'terms']


A Broombridge file can contain suggestions of initial states to use to
carry electronic computations of a molecule. In particular, they can be
used by the UCCSD ansatz to store information about the initial state
(i.e., the reference wavefunction) as well as initial values for the
variational parameters and the spin-orbital excitations to whic they
correspond.

Several initial states can be available and stored in a Broombridge file
as a result of classical computations from libraries such as NWChem, for
example. The user can specify which initial state to load with the
following code snippet:

.. code:: ipython3

    input_state = load_input_state(filename, "UCCSD |G>")
    print("input_state ::", input_state)
    print(dir(input_state))


.. parsed-literal::

    input_state :: <qsharp.chemistry.InputState object at 0x7f3e8c0ade80>
    ['Energy', 'MCFData', 'Method', 'SCFData', 'UCCData', '__class__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__']


Users can decide what excitations should be included in the ansatz and
how the values of variational parameters can be tied to specific
excitations, or enforce that a unique value should be tied to several
terms during the classical optimization later. The last entry in
``inputstate[Superposition]`` is the initial state, here showing a
Hartree-Fock state, with the two lower orbitals filled with one electron
each.

The chemistry library can now build the qubit Hamiltonian with a
transformation such as the Jordan–Wigner transformation.

.. code:: ipython3

    jw_hamiltonian = encode(ferm_hamiltonian, input_state)
    print("jw_hamiltonian :: \n", jw_hamiltonian)


.. parsed-literal::

    jw_hamiltonian :: 
     (4, ([([0], [0.17120128499999998]), ([1], [0.17120128499999998]), ([2], [-0.222796536]), ([3], [-0.222796536])], [([0, 1], [0.1686232915]), ([0, 2], [0.12054614575]), ([0, 3], [0.16586802525]), ([1, 2], [0.16586802525]), ([1, 3], [0.12054614575]), ([2, 3], [0.1743495025])], [], [([0, 1, 2, 3], [0.0, -0.0453218795, 0.0, 0.0453218795])]), (3, [((0.001, 0.0), [2, 0]), ((-0.001, 0.0), [3, 1]), ((-0.001, 0.0), [2, 3, 1, 0]), ((1.0, 0.0), [0, 1])]), -0.09883444600000002)


Please note that, currently, the underlying ``JordanWignerEncodingData``
data structure from the chemistry library is also used to store the
initial state for UCCSD as well as the variational parameters
representing the one- and two-body amplitudes (specified as the third
entry of the resulting ``jw_hamiltonian`` tuple object). In the future,
the objects may be kept separate and thus the ``inputState`` field may
not be required to compute the qubit Hamiltonian. Users can, however,
retrieve the values of the variational parameters directly from the data
structure, with a function such as the following:

.. code:: ipython3

    def get_var_params(jw_hamiltonian):
        """ Retrieve the values of variational parameters from the jw_hamiltonian object """
        _, _, input_state, _ = jw_hamiltonian
        _, var_params = input_state
        params = [param for ((param, _), _) in var_params]
        return params[:-1]
    
    var_params = get_var_params(jw_hamiltonian)
    print(var_params)


.. parsed-literal::

    [0.001, -0.001, -0.001]


2.4 Energy evaluation using the Q# quantum algorithms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The qsharp package can be used to directly call quantum algorithms
written in Q#. These can be user defined, or come from one of the
available Q# libraries.

The energy is computed as an expectation value
:math:`E(\theta) = \langle \Psi(\vec{\theta}) \vert \hat{H} \vert \Psi(\vec{\theta}) \rangle =\langle\hat{H}\rangle`,
which can be estimated by drawing many samples of the underlying
distribution (e.g., running the quantum circuit and measuring for each
sample). This approach is the one used on quantum hardware, and relies
on sampling to approach the expectation value, using the ``simulate``
function. The accuracy of the expectation value, and therefore the
result of the energy evaluation, directly correlates with the number of
samples used. The fast frequency estimator provided in the QDK allows
for the approximation of the result for a very large number of samples
without incurring longer runtimes.

.. code:: ipython3

    import qsharp
    
    # It is possible to create a Python object to represent a
    # Q# callable from the chemistry library
    estimate_energy = qsharp.QSharpCallable("Microsoft.Quantum.Chemistry.JordanWigner.VQE.EstimateEnergy", "")
    
    # The Q# operation can then be called through the simulate method
    # A large number of samples is selected for high accuracy
    energy = estimate_energy.simulate(jwHamiltonian=jw_hamiltonian, nSamples=1e18)
    
    print("Energy evaluated at {0} : {1} \n".format(var_params, energy))


.. parsed-literal::

    Energy evaluated at [0.001, -0.001, -0.001] : -1.1170458249057043 
    


2.5 Classical optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~

VQE is a quantum–classical hybrid algorithm that aims to compute
:math:`E = \min_{\vec{\theta}} \: \langle \Psi(\vec{\theta}) \vert \hat{H} \vert \Psi(\vec{\theta}) \rangle`.
This approach relies on solving an optimization problem, using a
classical optimizer to tune the values of the variational parameters
:math:`\{\theta_i\}_{i=1}^{m}`.

There are several Python libraries that provide implementations of
optimizers based on different heuristics, and SciPy is one that is
widely used. The optimizers in ``scipy.optimize`` have a common
interface that require users to provide the following:

-  A handle to a Python function to perform energy evaluations. It takes
   the variational parameters as its first input, leaving other
   parameters that are to be left out of the optimization process
   afterwards
-  Values for the initial parameters
-  Optional parameters used by our energy evaluation function, that
   should not be optimized
-  Optional parameters defining the behaviour and termination criteria
   for the chosen optimizer

The first item requires the user to provide a Python wrapper (here named
``energy_eval_wrapper``) with the expected signature, in order to call
the ``energy_evaluation`` operation available in the Quantum Development
Kit chemistry library. This wrapper requires the variational parameters
to be passed as a list or a NumPy array and, currently, an extra step is
needed to modify the data structure passed to the Q# context in order to
use the correct values (defined in ``set_var_params`` below).

.. code:: ipython3

    def set_var_params(var_params, jw_hamiltonian):
        """ Set variational parameters stored in the JW data-structure to the desired values"""
        # Unpack data structure
        a1, a2, input_state, a3 = jw_hamiltonian
        b1, amps = input_state
        # Unpack and overwrite variational parameters
        new_amps = [((var_params[i], 0.0), amps[i][1]) for i in range(len(var_params))]
        new_amps.append(amps[-1])
        # Re-pack the data structure
        input_state = (b1, new_amps)
        jw_hamiltonian = (a1, a2, input_state, a3)
        return jw_hamiltonian

.. code:: ipython3

    def energy_eval_wrapper(var_params, jw_hamiltonian, n_samples):
        """
            A wrapper whose signature is compatible with the use of scipy optimizers,
            calling the Q# energy_evalaution from the Microsoft Chemistry library
        """
        
        # NumPy arrays are currently not supported by the Python interops
        # This ensures that neither the user nor SciPy call the energy evaluation function with a NumPy array
        var_params = list(var_params)
    
        # Set the varational parameters to the right values in the jw_hamiltonian object
        jw_hamiltonian = set_var_params(var_params, jw_hamiltonian)
    
        # Estimate the energy
        energy = estimate_energy.simulate(jwHamiltonian=jw_hamiltonian, nSamples=1e18)
        
        print("Energy evaluated at {0} : {1} \n".format(var_params, energy))
        return energy

These two functions can then be used to run VQE. For simplicity, a
specific optimizer from the SciPy library is used, with given
hyperparameters such as tolerance or step size. Since accuracy of energy
evaluation is correlated to the number of samples drawn, it is important
to set it to a number large enough to guarantee that it is consistent
with the optimizers convergence criteria, and to ensure the correct
approximation of derivatives used by some optimizers. Setting a very
large number of samples would solve this issue.

.. code:: ipython3

    from scipy.optimize import minimize
    
    def VQE(initial_var_params, jw_hamiltonian, n_samples):
        """ Run VQE Optimization to find the optimal energy and the associated variational parameters """
        
        opt_result = minimize(energy_eval_wrapper,
                              initial_var_params,
                              args=(jw_hamiltonian, n_samples),
                              method="COBYLA",
                              tol=0.000001,
                              options={'disp': True, 'maxiter': 200,'rhobeg' : 0.05})
        
        return opt_result

.. code:: ipython3

    # Run VQE and print the results of the optimization process
    # A large number of samples is selected for higher accuracy
    opt_result = VQE(var_params, jw_hamiltonian, n_samples=1e18)
    print(opt_result)


.. parsed-literal::

    Energy evaluated at [0.001, -0.001, -0.001] : -1.1170458249805946 
    
    Energy evaluated at [0.051000000000000004, -0.001, -0.001] : -1.1150922913587564 
    
    Energy evaluated at [0.001, 0.049, -0.001] : -1.1151731333946984 
    
    Energy evaluated at [0.001, -0.001, 0.049] : -1.0951663452290825 
    
    Energy evaluated at [-0.003430545145800475, -0.005247198268165383, -0.05062188606433574] : -1.1309325802637815 
    
    Energy evaluated at [-0.010445845745391688, -0.011972188439066844, -0.09966839526426408] : -1.1367311459645775 
    
    Energy evaluated at [-0.04537630678071315, -0.04664899461950352, -0.10846509947574681] : -1.1331125319130892 
    
    Energy evaluated at [0.007178556191494, -0.029702738302247473, -0.09975816116182788] : -1.1363369902386204 
    
    Energy evaluated at [0.0066719657259913835, 0.011794797493765025, -0.14019136295048812] : -1.1359148381807342 
    
    Energy evaluated at [-0.0018869400097001523, -8.869547265090948e-05, -0.1199298791073761] : -1.1371913075695055 
    
    Energy evaluated at [0.010521493795355897, 0.019543610306031745, -0.1106776717859935] : -1.1367978572235846 
    
    Energy evaluated at [-0.015400536412154294, 0.003231884615372931, -0.1406990018657065] : -1.1358549240045668 
    
    Energy evaluated at [0.008149364771357121, -0.007271259306185439, -0.12191314992068644] : -1.1370698001316133 
    
    Energy evaluated at [-0.0036900235956029005, -0.002869720518257383, -0.11463107836777874] : -1.1372454767122933 
    
    Energy evaluated at [-0.0054437370217885125, -0.0048661142745314336, -0.11627562273894247] : -1.1372024903856963 
    
    Energy evaluated at [-0.004237819979936429, 0.0006665889202230337, -0.10950692627434655] : -1.137236549418017 
    
    Energy evaluated at [0.0015459922292730094, -0.005549589945731425, -0.11251803834607091] : -1.1372469700604788 
    
    Energy evaluated at [0.0032727534564713162, -0.003411541488510308, -0.11400554278600344] : -1.1372554240953492 
    
    Energy evaluated at [0.00323023473116795, -0.0024965224458383317, -0.11273970766263776] : -1.137260075113446 
    
    Energy evaluated at [0.004466434309022911, 0.00021904131835736055, -0.11181064951252091] : -1.137251797675844 
    
    Energy evaluated at [0.0044591327963416725, -0.003258520063994653, -0.11214761356087996] : -1.1372503916234913 
    
    Energy evaluated at [0.000280223928466486, -0.0024632307759734325, -0.11170919665897876] : -1.1372628703533096 
    
    Energy evaluated at [-1.8952911918186178e-05, -0.0009874144639177085, -0.11212616869790665] : -1.137268219995569 
    
    Energy evaluated at [-0.00024988152941975546, 0.0005129781387458419, -0.11175616834454018] : -1.1372674415612742 
    
    Energy evaluated at [-0.0007851360556661726, -0.0011227672846251054, -0.11205549644972532] : -1.1372669731581557 
    
    Energy evaluated at [0.00010327428680280016, -0.000692019058545139, -0.11365561576020648] : -1.1372694880750607 
    
    Energy evaluated at [0.0015218952725846444, -0.0007800424311977558, -0.1143045939066538] : -1.1372660496985976 
    
    Energy evaluated at [6.491968103281382e-05, 8.172189144339337e-05, -0.1135545938367093] : -1.1372700138895104 
    
    Energy evaluated at [-0.00032394453173831663, 6.687751905319809e-05, -0.1135885373082352] : -1.1372698900312619 
    
    Energy evaluated at [0.0003162121043161492, 0.000650294717502832, -0.114027804806539] : -1.1372684252598586 
    
    Energy evaluated at [6.401635903201671e-05, 0.00012684867376921196, -0.11277464876170319] : -1.1372702606142973 
    
    Energy evaluated at [0.00020192304224664857, 0.00046327088530851216, -0.11263185613958025] : -1.137269878230053 
    
    Energy evaluated at [0.00024399331007042308, 5.112414925927232e-05, -0.11277005897533424] : -1.1372702215426567 
    
    Energy evaluated at [-0.00013761953374534535, -0.0001873844480939262, -0.11265980799153413] : -1.137270098827592 
    
    Energy evaluated at [4.4394134442322125e-05, 0.00011490169258425586, -0.11296860548522411] : -1.1372703861649724 
    
    Energy evaluated at [2.6585387068326508e-05, 0.0001990321151588752, -0.11314396766236227] : -1.1372703703493303 
    
    Energy evaluated at [-4.8137261618236166e-05, -4.7159746994468276e-05, -0.11302623363113395] : -1.1372704072940851 
    
    Energy evaluated at [3.5720058904602486e-05, -8.982699713341738e-05, -0.11305239082372744] : -1.1372704078378484 
    
    Energy evaluated at [4.0280220670831035e-05, -0.0001084056463214515, -0.11300746617002783] : -1.13727039999582 
    
    Energy evaluated at [3.1096912977298147e-06, -0.00010465041286918662, -0.1131432400033304] : -1.1372703955107273 
    
    Energy evaluated at [8.020493777816505e-05, -2.95087939575659e-06, -0.11305561245103349] : -1.1372704091974448 
    
    Energy evaluated at [7.114371840639009e-05, 4.5405120530772794e-06, -0.11310300400248174] : -1.1372704076304005 
    
    Energy evaluated at [0.00010123568224228144, -1.3933152703004576e-05, -0.11306136952564164] : -1.1372704065487556 
    
    Energy evaluated at [4.6211439516522815e-05, 2.2048570222726553e-05, -0.11303104305165694] : -1.1372704102846574 
    
    Energy evaluated at [2.1928223809803936e-05, 3.551875273645271e-05, -0.11307120603721584] : -1.1372704122419846 
    
    Energy evaluated at [-8.074820719078531e-06, -2.9832792791550475e-06, -0.11306994245473628] : -1.1372704143439654 
    
    Energy evaluated at [-5.311916300356113e-05, -1.3016810251023318e-05, -0.1130858961776689] : -1.1372704106622367 
    
    Energy evaluated at [5.028837354391568e-06, -1.6777314190588513e-05, -0.11308524164405172] : -1.1372704131841944 
    
    Energy evaluated at [-3.5433477428735012e-06, -1.2590628666746656e-05, -0.1130479603930003] : -1.1372704135895173 
    
    Energy evaluated at [-2.0060919261995336e-05, -1.1931434713569711e-06, -0.1130714054883817] : -1.137270413855897 
    
    Energy evaluated at [-3.0239290622430287e-06, 8.097109608944007e-06, -0.11306909093403901] : -1.137270414081069 
    
    Energy evaluated at [-7.250959246873536e-06, -4.8526044800430895e-06, -0.1130641909524497] : -1.137270413991808 
    
    Energy evaluated at [-5.319460654529441e-06, -5.32021686171819e-06, -0.11307486176532157] : -1.1372704141290777 
    
    Energy evaluated at [-7.028080638375138e-06, -1.194519749966079e-07, -0.11307006915959633] : -1.1372704143441998 
    
    Energy evaluated at [-8.224046329610561e-06, 2.796524810104976e-07, -0.11307092862927882] : -1.137270414032172 
    
    Energy evaluated at [-4.943306254431374e-06, -7.861536878643887e-07, -0.11306794255424031] : -1.1372704143773436 
    
    Energy evaluated at [-5.9272894072404254e-06, -3.781668107744599e-07, -0.11306685002035084] : -1.1372704144161916 
    
    Energy evaluated at [-6.916674440001822e-06, 1.1830839211208511e-07, -0.11306400610484998] : -1.1372704140216614 
    
    Energy evaluated at [-6.323258813442331e-06, 8.980325296733387e-07, -0.11306758679026699] : -1.1372704142905068 
    
    Energy evaluated at [-7.021315898170545e-06, -1.2817650085998322e-06, -0.11306741120015319] : -1.137270414102867 
    
    Energy evaluated at [-5.2381128618438445e-06, -1.7869348471014356e-07, -0.1130661534242832] : -1.1372704140919954 
    
         fun: -1.1372704140919954
       maxcv: 0.0
     message: 'Optimization terminated successfully.'
        nfev: 61
      status: 1
     success: True
           x: array([-5.23811286e-06, -1.78693485e-07, -1.13066153e-01])


.. code:: ipython3

    # Print difference with exact FCI value known for this bond length
    fci_value = -1.1372704220924401
    print("Difference with exact FCI value :: ", abs(opt_result.fun - fci_value))


.. parsed-literal::

    Difference with exact FCI value ::  8.000444751132818e-09


3 Potential energy surface of H\ :math:`_\text{2}` with VQE, using the 1QBit OpenQEMIST package
-----------------------------------------------------------------------------------------------

The potential energy surface of this molecule can be obtained by
plotting the energy of the system as a function of the distance between
the hydrogen atoms.

This section shows how the 1QBit OpenQEMIST package allows users to run
VQE without relying on an input Broombridge file, or worrying about
modifying the data structures returned by the Python interop in the
previous section. Users can directly provide the geometry and basis set
of the target molecular system: OpenQEMIST computes the mean field and
electronic integrals using PySCF, generates the UCCSD one- and two-body
excitations, and provides good initial variational parameters using MP2
amplitudes.

OpenQEMIST provides several electronic structure solvers, such as VQE,
FCI, and CCSD. This package can be used to compute the
H\ :math:`_\text{2}` bond dissociation curve using VQE, with Microsoft
libraries, and compare it to the exact FCI values, computed on-the-fly.
Running the code cells in this section should yield a plot that closely
resembles the one below:

.. image:: img/h2_vqe.png
    :align: center
    :width: 600pt

.. code:: ipython3

    # Import the OpenQEMIST package from 1QBit and PySCF
    import openqemist
    import pyscf
    import numpy as np

.. code:: ipython3

    from pyscf import gto, scf
    from openqemist.electronic_structure_solvers import VQESolver, FCISolver
    from openqemist.quantum_solvers.parametric_quantum_solver import ParametricQuantumSolver
    from openqemist.quantum_solvers import MicrosoftQSharpParametricSolver
    
    # Iterate over different bond lengths
    bond_lengths = np.arange(0.4, 1.7, 0.1)
    energies_FCI, energies_VQE = [], []
    
    for bond_length in bond_lengths:
    
        # Create molecule object with PySCF
        H2 = [['H',[ 0, 0, 0]], ['H',[0,0, bond_length]]]
        mol = gto.Mole()
        mol.atom = H2
        mol.basis = "sto-3g"
        mol.charge = 0
        mol.spin = 0
        mol.build()
    
        # Compute FCI energy with PySCF, for reference
        solver = FCISolver()
        energy = solver.simulate(mol)
        energies_FCI += [energy]
        
        # Compute energy with VQE, instantiating a VQESolver object using the UCCSD ansatz
        solver = VQESolver()
        solver.hardware_backend_type = MicrosoftQSharpParametricSolver
        solver.ansatz_type = MicrosoftQSharpParametricSolver.Ansatze.UCCSD
        energy = solver.simulate(mol)
        energies_VQE += [energy]


.. parsed-literal::

    VQE : initial amplitudes
     [2e-05, 0.021503834487911277] 
    
    
    
    		Optimal UCCSD Singlet Energy: -0.9141502312031388
    		Optimal UCCSD Singlet Amplitudes: [-1.6596050e-06  2.9706155e-02]
    		Number of Function Evaluations :  21
    VQE : initial amplitudes
     [2e-05, 0.02513703925499215] 
    
    
    
    		Optimal UCCSD Singlet Energy: -1.0551590730735914
    		Optimal UCCSD Singlet Amplitudes: [4.72134759e-07 3.59485067e-02]
    		Number of Function Evaluations :  26
    VQE : initial amplitudes
     [2e-05, 0.02936700404572922] 
    
    
    
    		Optimal UCCSD Singlet Energy: -1.1162858146609462
    		Optimal UCCSD Singlet Amplitudes: [2.05541355e-05 4.35251734e-02]
    		Number of Function Evaluations :  34
    VQE : initial amplitudes
     [2e-05, 0.0341700987859853] 
    
    
    
    		Optimal UCCSD Singlet Energy: -1.1361890700497943
    		Optimal UCCSD Singlet Amplitudes: [1.17940598e-05 5.24413917e-02]
    		Number of Function Evaluations :  26
    VQE : initial amplitudes
     [2e-05, 0.039547719554227444] 
    
    
    
    		Optimal UCCSD Singlet Energy: -1.1341478045494089
    		Optimal UCCSD Singlet Amplitudes: [5.11612497e-06 6.27612372e-02]
    		Number of Function Evaluations :  32
    VQE : initial amplitudes
     [2e-05, 0.045541402725168836] 
    
    
    
    		Optimal UCCSD Singlet Energy: -1.1205606307437053
    		Optimal UCCSD Singlet Amplitudes: [3.11066004e-06 7.45854573e-02]
    		Number of Function Evaluations :  35
    VQE : initial amplitudes
     [2e-05, 0.05222992635439974] 
    
    
    
    		Optimal UCCSD Singlet Energy: -1.1011506003729807
    		Optimal UCCSD Singlet Amplitudes: [1.14960830e-05 8.80958506e-02]
    		Number of Function Evaluations :  27
    VQE : initial amplitudes
     [2e-05, 0.059706998024107984] 
    
    
    
    		Optimal UCCSD Singlet Energy: -1.0791918501736077
    		Optimal UCCSD Singlet Amplitudes: [7.01266147e-06 1.03432271e-01]
    		Number of Function Evaluations :  32
    VQE : initial amplitudes
     [2e-05, 0.0680584711024549] 
    
    
    
    		Optimal UCCSD Singlet Energy: -1.0567407729284206
    		Optimal UCCSD Singlet Amplitudes: [-2.16971436e-05  1.20663339e-01]
    		Number of Function Evaluations :  34
    VQE : initial amplitudes
     [2e-05, 0.07735157579471495] 
    
    
    
    		Optimal UCCSD Singlet Energy: -1.0351850444603288
    		Optimal UCCSD Singlet Amplitudes: [-1.96373416e-05  1.39659783e-01]
    		Number of Function Evaluations :  38
    VQE : initial amplitudes
     [2e-05, 0.0876349401540645] 
    
    
    
    		Optimal UCCSD Singlet Energy: -1.0154686955300325
    		Optimal UCCSD Singlet Amplitudes: [-1.03276798e-05  1.60160360e-01]
    		Number of Function Evaluations :  37
    VQE : initial amplitudes
     [2e-05, 0.0989425215637509] 
    
    
    
    		Optimal UCCSD Singlet Energy: -0.9981499394622826
    		Optimal UCCSD Singlet Amplitudes: [-3.64259658e-06  1.81669203e-01]
    		Number of Function Evaluations :  43
    VQE : initial amplitudes
     [2e-05, 0.11129638554843897] 
    
    
    
    		Optimal UCCSD Singlet Energy: -0.9834723123785654
    		Optimal UCCSD Singlet Amplitudes: [-1.93579249e-05  2.03597833e-01]
    		Number of Function Evaluations :  36


.. code:: ipython3

    import matplotlib.pyplot as plt
    %matplotlib inline
    
    plt.plot(bond_lengths, energies_FCI, color = 'black', label='Full CI')
    plt.plot(bond_lengths, energies_VQE, 'ro', label='UCCSD-VQE')
    plt.title("Potential energy surface of H2")
    plt.xlabel("Distance between hydrogen atoms (angstroms)")
    plt.ylabel("Energy (hartrees)")
    plt.legend()

