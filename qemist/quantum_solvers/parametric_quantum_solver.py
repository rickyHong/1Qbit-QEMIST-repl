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

    def __init__(self, molecule, anstatz, mean_field=None):
        self.initial_wavefunction = None
        self.n_qubits = None
        self.amplitude_dimension = None

    @abc.abstractmethod
    def simulate(self, amplitudes):
        """Performs an energy simulation for the given amplitudes.

        Args:
            amplitudes (list): The ampltudes to use in the simulation.

        Returns:
            float: The energy from the simulation.

        Raises:
            ValueError: If len(ampltudes) doesn't equal amplitude_dimension.
        """
        pass

    @abc.abstractmethod
    def get_rdm(self):
        """ Returns the RDM from the simulation.

        This is intended to be called after the simulation loop of the calling
        code has converged to be passed on to problem decompositions.  In a
        concrete implementation, the `simulate` function would set the necessary
        internal state for the class so that this function can return the
        reduced density matrix.

        Returns:
            tuple(numpy.array, numpy.array): The one- and two-element RDM
            matrices.

        Raises:
            RuntimeError: If no simulation has been run.
        """
        pass
