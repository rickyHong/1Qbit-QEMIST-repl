import abc

from .. import quantum_solvers

class ElectronicStructureSolver(abc.ABC):

    def __init__(self):
        self.quantum_solver = None

    @abc.abstractmethod
    def simulate(self, molecule, mean_field=None):
        pass

    @abc.abstractmethod
    def get_rdm(self):
        pass
