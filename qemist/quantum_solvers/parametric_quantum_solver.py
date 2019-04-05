import abc

class ParametricQuantumSolver(abc.ABC):
    """Performs energy estimations for a molecule with a parametric circuit.

    Performs energy estimations for a given molecule and a choice of ansatz
    circuit that is supported by the hardware. A concrete class will contain a
    nested Enum class that exposes the supported ansatze to the users. The state
    of the class is created in the costructor so that the `simulate` function
    can just accept an amplitude array and return an energy.

    Attributes:
        amplitude_dimension (int): The size of the amplitudes list that should
            be provided to the solver.
        initial_wavefunction (numpy.array): The initial wavefuntion for the
            simulations.
        n_qubits (int): The number of qubits in the circuit.
    """

    def __init__(self, molecule, anstatz):
        self.initial_wavefunction = None
        self.n_qubits = None
        self.amplitude_dimension = None
        pass

    @abc.abstractmethod
    def simulate(self, amplitudes):
        """Performs an energy simulation for the given amplitudes.

        Args:
            amplitudes (numpy.array): The amplitudes to use in the simulation.

        Returns:
            float: The energy from the simulation.

        Raises:
            ValueError: If amplitudes.shape doesn't equal
            (self.amplitude_dimension, 1).
        """
        pass

    @abc.abstractmethod
    def get_rdm(self):
        """ Returns the RDM from the simulation.

        This is intended to be called after the simulation loop of the calling
        code has converged to be passed on to problem decompositions.

        Returns:
            tuple(numpy.array, numpy.array): The one- and two-element RDM
            matrices.
        """
        pass
