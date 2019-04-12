import warnings
import itertools

import numpy as np
from pyscf import scf

from qemist.quantum_solvers.initial_parameters import mp2_intitial_amplitudes

from ..electronic_structure_solver import ElectronicStructureSolver

class VQESolver(ElectronicStructureSolver):
    """Estimates energy wih a variational quantum eigensolver algorithm.

    Uses the VQE algorithm to solve the electronic structure problem. By default
    an optimizer from scipy is used, but users can set any function whose first
    argument is the `simulate` function of the hardware backend and the second
    is the amplitudes to optimize over and returns a energy. See the
    implementation of `_default_optimizer` for a concrete example.
    Users should provide a hardware backend type that conforms to the interface
    of `qemist.quantum_solver.ParametricQuantumSolver` that the `VQESolver` will
    construct and use. Users should also provide an ansatze type that is
    supported by the backend.
    Users can also provide a function that takes a `pyscf.gto.Mole` as its first
    argument and `pyscf.scf.RHF` as is second and returns the inital amplitudes
    for the variational optimization. The user is responsible for ensuring that
    the dimension of the amplitudes vector is correct for the given molecule and
    andsatz choice.

    Attributes:
        hardware_backend_type (subclass of ParametricQuantumSolver): A type for
            the backend instance that is automatically constructed.
        ansatz_type (subclass of Enum): Type of ansatz that is supported by the
            backend.
        optimizer (function): Function that is called to optimize.
        initial_amplitude_function (function): Function that returns the initial
            amplitudes used for the optimization.
    """

    def __init__(self):
        self.hardware_backend_type = None
        self.hardware_backend = None
        self.ansatz_type = None
        self.optimizer = None
        self.inital_amplitude_function = mp2_intitial_amplitudes

    def simulate(self, molecule, mean_field=None):
        """Perform the simulation for the molecule.

        If the mean field is not provided it is automatically calculated.

        Args:
            molecule (pyscf.gto.Mole): The molecule to simulate.
            mean_field (pyscf.scf.RHF): The mean field of the molecule.
        """
        # Calculate the mean field if the user has not already done it.
        if not mean_field:
            mean_field = scf.RHF(molecule)
            mean_field.verbose = 0
            mean_field.scf()

        # Check the convergence of the mean field
        if not mean_field.converged:
            warnings.warn("VQESolver simulating with mean field not converged.",
                    RuntimeWarning)

        # Set up the instance of the hardware backend
        self.hardware_backend = self.hardware_backend_type(self.ansatz_type,
                molecule, mean_field)

        amplitudes = self.inital_amplitude_function(molecule, mean_field)

        # If the user didn't provide an optimizer, then we give them scipy's
        if not self.optimizer:
            self.optimizer = self._default_optimizer

        # Call the optimizer until it converges
        energy = self.optimizer(self.hardware_backend.simulate, amplitudes)

        return energy

    def get_rdm(self):
        """Returns the RDM from the hardware backend.

        Returns the reduced density matrices from the hardware backend. Does not
        catch and exceptions that the hardware backend raises if it is not in a
        state to return the RDM.

        Returns:
            tuple(numpy.array, numpy.array): The one- and two-element RDM
            matrices.

        Raises:
            RuntimeError: If no simulation has been run.
        """
        if not self.hardware_backend:
            raise RuntimeError("Cannot retrieve RDM because no simulation has been run.")

        return self.hardware_backend.get_rdm()

    def _default_optimizer(self, backend, amplitudes):
        """Function that can be set as a default optimizer.

        Funciton that is used by the class as a default optimizer when user does
        not provide one.
        """
        from scipy.optimize import minimize
        result = minimize(backend, amplitudes, method='SLSQP',
                options={'disp':False, 'maxiter': 15000})

        print("\n\t\tOptimal UCCSD Singlet Energy: {}".format(result.fun))
        print("\t\tOptimal UCCSD Singlet Amplitudes: {}".format(result.x))
        print("\t\tNumber of Iterations : ", result.nit)
        print("\t\tNumber of Function Evaluations : ", result.nfev)

        return result.fun
