import abc

class QuantumSolver(abc.ABC):

     def __init__(self):
         pass

     @abc.abstractmethod
     def simulate(self, molecule, mean_field=None):
         pass
