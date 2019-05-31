
OpenQEMIST-DMET with Microsoft Quantum Development Kit Libraries
================================================================

1. Introduction
---------------

One of the main objectives of quantum chemistry calculations in the area
of materials science is to solve the electronic structure problem,
:math:`H\Psi=E\Psi`, as accurately as possible, in order to accelerate
the materials design process. Here we solve the problem for the molecule
shown below.

.. image:: img/exact.png
    :align: center
    :width: 200pt

The computational cost for performing accurate calculations
of the electronic structure of molecules, however, is usually very
expensive. For example, the cost of performing the full CI calculation
scales exponentially on a classical computer as the size of the system
increases. Therefore, when we target large-sized molecules, those
relevant for industry problems, it becomes essential to employ an
appropriate strategy for reducing the computational cost, while
maintaining accuracy, when performing electronic structure calculations.
One such strategy is to decompose a molecular system into its
constituent fragments and its environment, for each fragment
independently, as appropriate for the problem. First, the environment
outside of a fragment is calculated using a less-accurate method than
will be used to calculate the electronic structure of a fragment. Then,
the electronic structure problem for a given fragment is solved to a
high degree of accuracy, which includes the quantum mechanical effects
of the environment. The quantum mechanical description is updated (i.e.,
solved iteratively as shown below) by incorporating the just-performed
highly accurate calculation. In the following schematic illustration, the molecule
shown above is decomposed into fragments. Each molecular fragment
CH\ :math:`_\mathrm{3}` and CH\ :math:`_\text{2}` are the fragments chosen for
the electronic structure calculation, with the rest of the molecular system being the
surrounding environment.

.. image:: img/iterations.png
    :align: center
    :width: 600pt

2. Density matrix embedding theory (DMET)
-----------------------------------------

2-A. Theory
~~~~~~~~~~~

One successful decomposition approach is the DMET
method\ :math:`^{1,2}`. The DMET method decomposes a molecule into
fragments, and each fragment is treated as an open quantum system that
is entangled with each of the other fragments, all taken together to be
that fragment’s surrounding environment (or “bath”). In this framework,
the electronic structure of a given fragment is obtained by solving the
following Hamiltonian, by using a highly accurate quantum chemistry
method, such as the full CI method or a coupled-cluster method.

.. math::  H_{I}=\sum^{\text{frag}+\text{bath}}_{rs}  \left[ h_{rs} + \sum_{mn} \left[ (rs|mn) - (rn|ms) \right] D^{\text{(mf)env}}_{mn} \right] a_{r}^{\dagger}a_{s} + \sum_{pqrs}^{\text{frag}+\text{bath}} (pq|rs) a_{p}^{\dagger}a_{r}^{\dagger}a_{s}a_{q} - \mu\sum_{r}^{\text{frag}} a_{r}^{\dagger}a_{r} 

The expression
:math:`\sum_{mn} \left[ (rs|mn) - (rn|ms) \right] D^{\text{(mf)env}}_{mn}`
describes the quantum mechanical effects of the environment on the
fragment, where :math:`D^{\text{(mf)env}}_{mn}` is the mean-field
electronic density obtained by solving the Hartree–Fock equation. The
quantum mechanical effects from the environment are incorporated through
the one-particle term of the Hamiltonian. The extra term
:math:`\mu\sum_{r}^{\text{frag}} a_{r}^{\dagger}a_{r}` ensures, through
the adjustment of :math:`\mu`, that the number of electrons in all of
the fragments, taken together, becomes equal to the total number of
electrons in the entire system.

2-B. Performance of DMET
~~~~~~~~~~~~~~~~~~~~~~~~

Here we provide examples of the performance results of DMET calculations
(using classical simulation), employing some organic hydrocarbons
(C\ :math:`_\mathrm{n}`\ H\ :math:`_\mathrm{n}`), below: tetrahedrane
(C\ :math:`_\mathrm{4}`\ H\ :math:`_\mathrm{4}`, top left), prismane
(C\ :math:`_\mathrm{6}`\ H\ :math:`_\mathrm{6}`, top centre), cubane
(C\ :math:`_\mathrm{8}`\ H\ :math:`_\mathrm{8}`, top right), cuneane
(C\ :math:`_\mathrm{8}`\ H\ :math:`_\mathrm{8}`, bottom left), pentaprismane
(C\ :math:`_\mathrm{10}`\ H\ :math:`_\mathrm{10}`, bottom centre), and diademane
(C\ :math:`_\mathrm{10}`\ H\ :math:`_\mathrm{10}`, bottom right). In these examples,
each fragment consists of only one atom, thereby largely reducing the
size of the electronic structure problem to be solved. Of the several
electronic structure solvers used in DMET calculation we select the CCSD
method (as it is the one most commonly used), the Meta-Löwdin as the
localization\ :math:`^{3}` scheme, and cc-pVDZ as the basis set. Visualizations
are done with py3Dmol\ :math:`^{4}`.

.. code:: ipython3

    import py3Dmol
    view = py3Dmol.view(width=600,height=600,viewergrid=(2,3))
    
    tetrahedrane = open('crd/tetrahedrane.xyz', 'r').read()
    prismane = open('crd/prismane.xyz', 'r').read()
    cubane = open('crd/cubane.xyz', 'r').read()
    cuneane = open('crd/cuneane.xyz', 'r').read()
    pentaprismane = open('crd/pentaprismane.xyz', 'r').read()
    diademane = open('crd/diademane.xyz', 'r').read()
    
    view.addModel(tetrahedrane,'xyz',viewer=(0,0))
    view.addModel(prismane,'xyz',viewer=(0,1))
    view.addModel(cubane,'xyz',viewer=(0,2))
    view.addModel(cuneane,'xyz',viewer=(1,0))
    view.addModel(pentaprismane,'xyz',viewer=(1,1))
    view.addModel(diademane,'xyz',viewer=(1,2))
    
    view.setStyle({'stick':{'colorscheme':'cyanCarbon'}})
    view.zoomTo()
    view.show()



.. raw:: html

    <div id="3dmolviewer_15591687007323346"  style="position: relative; width: 600px; height: 600px">
            <p id="3dmolwarning_15591687007323346" style="background-color:#ffcccc;color:black">You appear to be running in JupyterLab (or JavaScript failed to load for some other reason).  You need to install the 3dmol extension: <br>
            <tt>jupyter labextension install jupyterlab_3dmol</tt></p>
            </div>
    <script>
    
    var loadScriptAsync = function(uri){
      return new Promise((resolve, reject) => {
        var tag = document.createElement('script');
        tag.src = uri;
        tag.async = true;
        tag.onload = () => {
          resolve();
        };
      var firstScriptTag = document.getElementsByTagName('script')[0];
      firstScriptTag.parentNode.insertBefore(tag, firstScriptTag);
    });
    };
    
    if(typeof $3Dmolpromise === 'undefined') {
    $3Dmolpromise = null;
      $3Dmolpromise = loadScriptAsync('https://3dmol.csb.pitt.edu/build/3Dmol.js');
    }
    
    var viewer_15591687007323346 = null;
    var warn = document.getElementById("3dmolwarning_15591687007323346");
    if(warn) {
        warn.parentNode.removeChild(warn);
    }
    $3Dmolpromise.then(function() {
    var viewergrid_15591687007323346 = null;
    viewergrid_15591687007323346 = $3Dmol.createViewerGrid($("#3dmolviewer_15591687007323346"),{rows: 2, cols: 3, control_all: true},{backgroundColor:"white"});
    viewer_15591687007323346 = viewergrid_15591687007323346[0][0];
    	viewergrid_15591687007323346[0][0].addModel("8\n\nC               -0.580517    0.479710   -0.503842\nH               -1.269514    1.045226   -1.098983\nC               -0.462763   -0.528766    0.571790\nH               -1.009781   -1.152554    1.250252\nC                0.529158    0.565261    0.470571\nH                1.157329    1.232467    1.026343\nC                0.514143   -0.516046   -0.538653\nH                1.121845   -1.126087   -1.176808\n","xyz");
    	viewergrid_15591687007323346[0][1].addModel("12\n\nC       -0.779069    0.812384   -0.334875\nH       -1.513521    1.549153   -0.638932\nC       -0.779000   -0.696246   -0.536005\nH       -1.513519   -1.327946   -1.021958\nC       -0.779055   -0.116090    0.870880\nH       -1.513553   -0.221553    1.660871\nC        0.779033    0.812221   -0.335269\nH        1.513544    1.548876   -0.639457\nC        0.779051   -0.696501   -0.535677\nH        1.513486   -1.328469   -1.021412\nC        0.779039   -0.115700    0.870931\nH        1.513564   -0.220471    1.660987\n","xyz");
    	viewergrid_15591687007323346[0][2].addModel("16\n\nC        0.971771    0.891992    0.331675\nH        1.751319    1.607597    0.597981\nC       -0.171066    1.143910   -0.715642\nH       -0.308195    2.061588   -1.289876\nC        1.219679   -0.532602   -0.281012\nH        2.198034   -0.959893   -0.506708\nC       -0.076781    0.280861    1.328242\nH       -0.138252    0.506176    2.394154\nC        0.076780   -0.280863   -1.328242\nH        0.138275   -0.506175   -2.394154\nC       -1.219669    0.532598    0.281015\nH       -2.198031    0.959876    0.506714\nC        0.171050   -1.143913    0.715642\nH        0.308185   -2.061593    1.289873\nC       -0.971766   -0.891982   -0.331679\nH       -1.751323   -1.607581   -0.597978\n","xyz");
    	viewergrid_15591687007323346[1][0].addModel("16\n\nC        0.000002   -0.969623    0.784721\nC       -0.000002   -0.969594   -0.784761\nC        1.270774   -0.079806    0.761542\nC       -1.270806   -0.079860    0.761505\nC        1.270838   -0.079878   -0.761474\nC       -1.270806   -0.079824   -0.761511\nC        0.763631    1.126462   -0.000042\nC       -0.763632    1.126462    0.000014\nH        0.000005   -1.905675    1.341769\nH       -0.000007   -1.905624   -1.341842\nH        2.098533   -0.063036    1.463654\nH       -2.098600   -0.063159    1.463575\nH        2.098675   -0.063170   -1.463492\nH       -2.098609   -0.063048   -1.463570\nH        1.333415    2.048839   -0.000102\nH       -1.333415    2.048839    0.000036\n","xyz");
    	viewergrid_15591687007323346[1][1].addModel("20\n\nC        1.237614   -0.480729   -0.785182\nH        2.087357   -0.810965   -1.386440\nC        1.237658   -0.480550    0.785223\nH        2.087434   -0.810648    1.386508\nC        0.839579    1.028082   -0.785499\nH        1.416033    1.734163   -1.386929\nC        0.839623    1.028262    0.785219\nH        1.416109    1.734480    1.386457\nC       -0.074833   -1.325342   -0.785169\nH       -0.126099   -2.235639   -1.386290\nC       -0.074788   -1.325161    0.785474\nH       -0.126022   -2.235322    1.386805\nC       -0.718566    1.116073   -0.785360\nH       -1.211874    1.882779   -1.386576\nC       -0.718526    1.116252    0.785146\nH       -1.211800    1.883093    1.386217\nC       -1.283948   -0.338559   -0.785292\nH       -2.165341   -0.570962   -1.386760\nC       -1.283903   -0.338382    0.785441\nH       -2.165261   -0.570655    1.387008\n","xyz");
    	viewergrid_15591687007323346[1][2].addModel("20\n\nC        1.300981   -0.780133   -0.586795\nC        1.432197   -0.020041    0.713195\nC        1.322454    0.743435   -0.586557\nC       -0.017653   -1.516720   -0.586674\nC        0.000434    0.000018    1.336087\nC        0.025066    1.516763   -0.586282\nC       -0.733084   -1.230010    0.713600\nC       -1.326555   -0.736648   -0.586080\nC       -0.698638    1.250019    0.713309\nC       -1.305248    0.773436   -0.586318\nH        2.155120   -1.287031   -1.026652\nH        2.327511   -0.032430    1.327931\nH        2.190151    1.226966   -1.026190\nH       -0.033103   -2.509826   -1.026506\nH        0.000674   -0.000143    2.427208\nH        0.036626    2.509969   -1.026037\nH       -1.190202   -1.999741    1.328578\nH       -2.192656   -1.222918   -1.025867\nH       -1.136178    2.030929    1.328443\nH       -2.157662    1.283512   -1.025815\n","xyz");
    	viewergrid_15591687007323346[0][0].setStyle({"stick": {"colorscheme": "cyanCarbon"}});
    	viewergrid_15591687007323346[0][1].setStyle({"stick": {"colorscheme": "cyanCarbon"}});
    	viewergrid_15591687007323346[0][2].setStyle({"stick": {"colorscheme": "cyanCarbon"}});
    	viewergrid_15591687007323346[1][0].setStyle({"stick": {"colorscheme": "cyanCarbon"}});
    	viewergrid_15591687007323346[1][1].setStyle({"stick": {"colorscheme": "cyanCarbon"}});
    	viewergrid_15591687007323346[1][2].setStyle({"stick": {"colorscheme": "cyanCarbon"}});
    	viewergrid_15591687007323346[0][0].zoomTo();
    	viewergrid_15591687007323346[0][1].zoomTo();
    	viewergrid_15591687007323346[0][2].zoomTo();
    	viewergrid_15591687007323346[1][0].zoomTo();
    	viewergrid_15591687007323346[1][1].zoomTo();
    	viewergrid_15591687007323346[1][2].zoomTo();
    viewergrid_15591687007323346[1][2].render();
    viewergrid_15591687007323346[1][1].render();
    viewergrid_15591687007323346[1][0].render();
    viewergrid_15591687007323346[0][2].render();
    viewergrid_15591687007323346[0][1].render();
    viewergrid_15591687007323346[0][0].render();
    });
    </script>


2-B-a. Performance of DMET: Accuracy of calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This table shows the CCSD total energies (in a.u.), as well as the total
energy difference (in kcal/mol) of DMET, MP2, B3LYP (DFT), and HF with
respect to the reference CCSD value.

.. image:: img/Table_DMET_organic_compounds.png
    :align: center
    :width: 680pt

The total energy values of the DMET calculations agree with those
obtained from CCSD, with only a small error, even though the fragment
size in the DMET calculations is very small (i.e., there is only one
atom per fragment). The calculations require only about 5% of the
amplitudes (i.e., the parameters to be optimized) for tetrahedrane, and
only 0.1% for pentaprismane, compared to performing a CCSD calculation
of the full system. The number of terms of the Hamiltonian in DMET
calculations is 1.5% and 0.05% of the full system for tetrahedrane and
pentaprismane, respectively. A large amount of computational resources
will therefore be saved without affecting the accuracy of the
calculations.

2-B-b. Performance of DMET: Computational time
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The table below shows the computation time required for the full CCSD
and DMET calculations, and the computation time of the DMET calculation
per fragment (i.e., the DMET calculation time divided by the number of
fragments used to decompose the molecule). Although the present examples
are based on a serial implementation of DMET, the DMET calculation for
each fragment can be trivially parallelized. Therefore, the DMET calculation
time per fragment corresponds approximately to that of DMET executed in parallel.

.. image:: img/Table_DMET_time_organic_compounds.png
    :align: center
    :width: 400

As shown in the plot, the computation time of the parallellized DMET
calculation (blue) begins to demonstrate its advantage over the full
CCSD calculation (red) as the molecular size increases.

.. image:: img/Time_plot.png
    :align: center
    :width: 400

3. OpenQEMIST-DMET sample calculation (classical simulation): A ring of 10 hydrogen atoms
-----------------------------------------------------------------------------------------

3-A. Sample DMET script for OpenQEMIST
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here, we demonstrate how to perform DMET calculations using 1QBit’s
`OpenQEMIST (Quantum-Enabled Molecular ab Initio Simulation
Toolkit) <http://>`__ software package. Harnessing the power of emerging
quantum computing technologies combined with state-of-the-art classical
techniques, OpenQEMIST is able to deliver either state-of-the-art
classical solutions or, with the flip of a switch, map a computationally
challenging problem onto quantum computing architectures without the
need for any additional programming, neither by a user nor a developer.

We have selected a ring of 10 hydrogen atoms as a simple example of a
molecular system. The distance between adjacent hydrogen atoms has been
set to 1\ :math:`~`\ Å.

.. code:: ipython3

    H10='''
    H          1.6180339887          0.0000000000          0.0000000000
    H          1.3090169944          0.9510565163          0.0000000000
    H          0.5000000000          1.5388417686          0.0000000000
    H         -0.5000000000          1.5388417686          0.0000000000
    H         -1.3090169944          0.9510565163          0.0000000000
    H         -1.6180339887          0.0000000000          0.0000000000
    H         -1.3090169944         -0.9510565163          0.0000000000
    H         -0.5000000000         -1.5388417686          0.0000000000
    H          0.5000000000         -1.5388417686          0.0000000000
    H          1.3090169944         -0.9510565163          0.0000000000
    '''
    
    view = py3Dmol.view(width=400,height=400)
    view.addModel("10\n" + H10,'xyz',{'keepH': True})
    view.setStyle({'sphere':{}})
    view.setStyle({'model':0},{'sphere':{'colorscheme':'cyanCarbon','scale':'0.2'}})
    view.zoomTo()
    view.show()



.. raw:: html

    <div id="3dmolviewer_15591687007451837"  style="position: relative; width: 400px; height: 400px">
            <p id="3dmolwarning_15591687007451837" style="background-color:#ffcccc;color:black">You appear to be running in JupyterLab (or JavaScript failed to load for some other reason).  You need to install the 3dmol extension: <br>
            <tt>jupyter labextension install jupyterlab_3dmol</tt></p>
            </div>
    <script>
    
    var loadScriptAsync = function(uri){
      return new Promise((resolve, reject) => {
        var tag = document.createElement('script');
        tag.src = uri;
        tag.async = true;
        tag.onload = () => {
          resolve();
        };
      var firstScriptTag = document.getElementsByTagName('script')[0];
      firstScriptTag.parentNode.insertBefore(tag, firstScriptTag);
    });
    };
    
    if(typeof $3Dmolpromise === 'undefined') {
    $3Dmolpromise = null;
      $3Dmolpromise = loadScriptAsync('https://3dmol.csb.pitt.edu/build/3Dmol.js');
    }
    
    var viewer_15591687007451837 = null;
    var warn = document.getElementById("3dmolwarning_15591687007451837");
    if(warn) {
        warn.parentNode.removeChild(warn);
    }
    $3Dmolpromise.then(function() {
    viewer_15591687007451837 = $3Dmol.createViewer($("#3dmolviewer_15591687007451837"),{backgroundColor:"white"});
    	viewer_15591687007451837.addModel("10\n\nH          1.6180339887          0.0000000000          0.0000000000\nH          1.3090169944          0.9510565163          0.0000000000\nH          0.5000000000          1.5388417686          0.0000000000\nH         -0.5000000000          1.5388417686          0.0000000000\nH         -1.3090169944          0.9510565163          0.0000000000\nH         -1.6180339887          0.0000000000          0.0000000000\nH         -1.3090169944         -0.9510565163          0.0000000000\nH         -0.5000000000         -1.5388417686          0.0000000000\nH          0.5000000000         -1.5388417686          0.0000000000\nH          1.3090169944         -0.9510565163          0.0000000000\n","xyz",{"keepH": true});
    	viewer_15591687007451837.setStyle({"sphere": {}});
    	viewer_15591687007451837.setStyle({"model": 0},{"sphere": {"colorscheme": "cyanCarbon", "scale": "0.2"}});
    	viewer_15591687007451837.zoomTo();
    viewer_15591687007451837.render();
    });
    </script>


Here we give the steps of a sample DMET calculation script.

Import OpenQEMIST modules
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    import openqemist
    print(openqemist.__version__)


.. parsed-literal::

    0.0.1


Import DMET modules and localization schemes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    from openqemist.problem_decomposition import DMETProblemDecomposition
    from openqemist.problem_decomposition.electron_localization import iao_localization, meta_lowdin_localization

OpenQEMIST gives a user the ability to easily switch bewtween several
electronic structure solvers, regardless of whether it is a classical or
quantum solver. Here we present sample code using classical electronic
structure solvers. In this open source version of OpenQEMIST, the Full
CI and CCSD methods are currently available.

Import classical electronic structure solver modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    from openqemist.electronic_structure_solvers import FCISolver
    from openqemist.electronic_structure_solvers import CCSDSolver

In OpenQEMIST, the inputs to all items of type “object” in OpenQEMIST
are objects from the PySCF\ :math:`^{5}` program. First we create the
molecule object. Then we set up the OpenQEMIST objects.

Here we create a molecule object, and a problem decomposition object (a
OpenQEMIST object), which problem decomposition techniques in OpenQEMIST
require. The problem decomposition object holds an instance of an
electronic structure solver, in this case the classical CCSD solver. We
use the solver to perform the DMET simulation. An orbital localization
technique needs to be defined to execute the DMET simulation. We employ
the Meta-Löwdin localization method in this example DMET simulation.
The orbital were depicted with VMD\ :math:`^{6}`.

Build molecule object (using PySCF) for DMET calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    from pyscf import gto
    mol = gto.Mole() # Instantiate the molecule class in PySCF
    mol.atom = H10   # The coordinates of the atoms of the 10-hydrogen-atom ring are defined above
    mol.basis = "minao" # Use "minao" as the basis set
    mol.charge = 0 # Assign the charge of the molecule 
    mol.spin = 0 # Assign the spin of the molecule
    mol.build() # Build the molecule object




.. parsed-literal::

    <pyscf.gto.mole.Mole at 0x7faef51b7400>



Instantiate DMET class
^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    dmet_solver = DMETProblemDecomposition()

Instantiate CCSD class
^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    dmet_solver.electronic_structure_solver = CCSDSolver()

Select orbital localization technique (here Meta-Löwdin)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    dmet_solver.electron_localization_method = meta_lowdin_localization

We perform a DMET calculation with one atom per fragment, with the
localization of molecular orbitals being executed before entering the
DMET loop. The resulting orbitals localized on each fragment are
depicted here.

Perform DMET calculation
^^^^^^^^^^^^^^^^^^^^^^^^

The “simulate” function takes two arguments: 1. the molecule object; 2.
a list of the number of atoms each fragment has (i.e., in this case, 10
fragments, with one atom per fragment)

.. code:: ipython3

    energy = dmet_solver.simulate(mol, [1,1,1,1,1,1,1,1,1,1])
    
    print(energy)


.. parsed-literal::

    -5.3675327924452745


3-B. Results of DMET calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This plot shows the potential energy curves of the ring of hydrogen
atoms in atomic units (a.u.) for four methods.

After repeating the DMET calculations for the ring of hydrogen atoms by
symmetrically stretching the distance between them, we obtain discrete
sample points of the potential energy, which we plot alongside the
curves of the other methods. The energy has been plotted as the energy
per atom.

.. image:: img/H10_stretch.png
    :align: center
    :width: 500

The results obtained from the DMET-CCSD method (using problem
decomposition) are almost identical to those of the Full CI method
(without using problem decomposition). When we decompose the ring of
atoms into fragments, one of which includes only one hydrogen atom, the
DMET method creates a fragment orbital (left: the single orbital
distribution is shown in both pink and blue, with the colours depicting
the phases) and the bath orbital (right: the single orbital distribution
of the remaining nine hydrogen atoms is shown in both pink and blue,
with the colours depicting the phases).

.. image:: img/frag_and_bath.png
    :align: center
    :width: 450

Then, the DMET Hamiltonian will consist of only two electrons and two
(i.e., fragment and bath) orbitals. Therefore, the CCSD solver, which
treats single- and double-excitations, will provide results equivalent
to those of the Full CI solver. This is why the results of the DMET-CCSD
and Full CI methods almost coincide.

The DMET calculations require only about 0.6% of the amplitudes compared
to the CCSD calculation for the full system consisting of 10 hydrogen
atoms. The number of Hamiltonian terms in the DMET calculation in the
spacial orbital basis is only 20% (17% in the spin orbital basis) of the
full system. A large amount of computational resources are saved while
maintaining the same level of accuracy in the calculations.

4. DMET-VQE quantum simulation with `Microsoft Quantum Development Kit Libraries <https://docs.microsoft.com/en-us/quantum/?view=qsharp-preview>`__\ 
-----------------------------------------------------------------------------------------------------------------------------------------------------

4-A. Sample script of DMET calculation combining OpenQEMIST and MS QDK 1: A ring of 10 hydrogen atoms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here we describe how to perform DMET calculations using the variational
quantum eigensolver (VQE) as the electronic structure solver. We use the
UCCSD-VQE framework available in the Microsoft Quantum Development Kit
libraries.\ :math:`^{7}`

.. code:: ipython3

    from openqemist.quantum_solvers.parametric_quantum_solver import ParametricQuantumSolver
    from openqemist.quantum_solvers import MicrosoftQSharpParametricSolver
    from openqemist.electronic_structure_solvers import VQESolver
    
    vqe = VQESolver()
    vqe.hardware_backend_type = MicrosoftQSharpParametricSolver
    vqe.ansatz_type = MicrosoftQSharpParametricSolver.Ansatze.UCCSD
    
    dmet = DMETProblemDecomposition()
    dmet.electron_localization_method = meta_lowdin_localization
    dmet.electronic_structure_solver = vqe
    energy_vqe = dmet.simulate(mol, [1,1,1,1,1,1,1,1,1,1])

.. code:: ipython3

    print(energy_vqe)
    print(energy-energy_vqe)


.. parsed-literal::

    -5.368762305205017
    0.0012295127597425903


4-B. Results of DMET quantum simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The results of the DMET-CCSD calculation and those of the DMET-UCCSD-VQE
method almost coincide. It has been reported\ :math:`^{8}` that the
UCCSD method performs better than the CCSD method; however, as
discussed, both the CCSD and UCCSD solvers provide results equivalent to
those of the Full CI solver.

The DMET method decomposes the ring of 10 hydrogen atoms into 10
subproblems, each of which requires only four qubits to perform a
quantum simulation, whereas the CCSD calculation of the full system
requires 20 qubits, as shown.

.. image:: img/Table_DMET_qubits_H10.png
    :align: center
    :width: 430

4-C. Sample script of DMET calculation combining OpenQEMIST and MS QDK 2: Symmetric stretching of a ring of 10 hydrogen atoms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We repeat the DMET calculations for the ring of hydrogen atoms by
symmetrically stretching the distance between the atoms.

.. code:: ipython3

    energy_vqe_table = {}
    for x in range(1,22):
        HH = 0.5+((x-1)*0.1)
        H10 = open('crd/h10_'+str(x)+'.xyz', 'r').readlines()[1:]
        H10 = ''.join(H10)
        
        mol = gto.Mole()
        mol.atom = H10
        mol.basis = "minao"
        mol.charge = 0
        mol.spin = 0
        mol.build()
    
        dmet = DMETProblemDecomposition()
        dmet.electron_localization_method = meta_lowdin_localization
        dmet.electronic_structure_solver = vqe
        energy_vqe = dmet_solver.simulate(mol, [1,1,1,1,1,1,1,1,1,1])
        energy_vqe_table.update({str(HH):energy_vqe})


Results of DMET quantum simulation 2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The potential energy curve of the DMET-VQE and FCI methods almost
coincide, as shown.

It is also possible to estimate the quantum resources required based on
the classical DMET calculations shown above. For pentaprismane, shown in
an example classical simulation (see Section 2-B), DMET calculation
requires 56 qubits, whereas 380 qubits are necessary for CCSD
calculation of the full system. Thus, the DMET method can be a powerful
tool for greatly reducing the computational resources needed.

.. image:: img/H10_stretch_with_VQE.png
    :align: center
    :width: 500

5. References
-------------

1. Gerald Knizia and Garnet K.-L. Chan, “Density Matrix Embedding: A
   Simple Alternative to Dynamical Mean-Field Theory”, Phys. Rev. Lett.,
   109, 186404 (2012).
2. Sebastian Wouters, Carlos A. Jiménez-Hoyos, Qiming Sun, and Garnet
   K.-L. Chan, “A Practical Guide to Density Matrix Embedding Theory in
   Quantum Chemistry”, J. Chem. Theory Comput., 12, pp. 2706–2719
   (2016).
3. Qiming Sun and Garnet K.-L. Chan, “Exact and Optimal Quantum
   Mechanics/Molecular Mechanics Boundaries”, J. Chem. Theory Comp.,
   10, 3784--3790 (2014).
4. py3Dmol. https://github.com/3dmol/3Dmol.js/tree/master/py3Dmol
5. Qiming Sun, Timothy C. Berkelbach, Nick S. Blunt, George H. Booth,
   Sheng Guo, Zhendong Li, Junzi Liu, James D. McClain, Elvira R.
   Sayfutyarova, Sandeep Sharma, Sebastian Wouters, and Garnet Kin‐Lic
   Chan, “PySCF: the Python‐based simulations of chemistry framework”,
   Wiley Interdiscip. Rev. Comput. Mol. Sci., 8, e1340 (2017).
6. William Humphrey, Andrew Dalke, and Klaus Schulten, “VMD – Visual
   Molecular Dynamics”, J. Molec. Graphics, 14, pp. 33–38 (1996).
   http://www.ks.uiuc.edu/Research/vmd
7. Guang Hao Low, Nicholas P. Bauman, Christopher E. Granade, Bo Peng, Nathan
   Wiebe, Eric J. Bylaska, Dave Wecker, Sriram Krishnamoorthy, Martin Roetteler,
   Karol Kowalski, Matthias Troyer, Nathan A. Baker, “Q# and NWChem: Tools for
   Scalable Quantum Chemistry on Quantum Computers”, arXiv:1904.01131 (2019).
8. Michael Kühn, Sebastian Zanker, Peter Deglmann, Michael Marthaler,
   and Horst Weiß, “Accuracy and Resource Estimations for Quantum
   Chemistry on a Near-term Quantum Computer”, arXiv:1812.06814 (2018).

