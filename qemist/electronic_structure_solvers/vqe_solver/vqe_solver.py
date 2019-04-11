import warnings
import itertools

import numpy as np
from pyscf import scf

from qemist.quantum_solvers.initial_parameters import mp2_intitial_amplitudes

from ..electronic_structure_solver import ElectronicStructureSolver

class VQESolver(ElectronicStructureSolver):

    def __init__(self):
        self.hardware_backend_type = None
        self.ansatz_type = None
        self.optimizer = None
        self.inital_amplitude_function = mp2_intitial_amplitudes

    def simulate(self, molecule, mean_field=None):
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
        self.hardware_backend = self.hardware_backend_type(self.ansatz_type, molecule,
                mean_field)

        amplitudes = self.inital_amplitude_function(molecule, mean_field)

        # If the user didn't provide an optimizer, then we give them scipy's
        if not self.optimizer:
            self.optimizer = self._default_optimizer

        # Call the optimizer until it converges
        energy = self.optimizer(self.hardware_backend.simulate, amplitudes)

        return energy

    def get_rdm(self):
        return self.hardware_backend.get_rdm()

    def _default_optimizer(self, backend, amplitudes):
            from scipy.optimize import minimize
            result = minimize(backend, amplitudes, method='SLSQP',
                    options={'disp':False, 'maxiter': 15000})

            return result.fun
