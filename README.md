# QEMIST

This library aims at pushing the boundaries of advanced materials science by leveraging emerging quantum computing technologies together with state-of-the-art classical techniques in a hybrid framework.

It allows for the simulation of intractable electronic structure problems through problem decomposition and problem reduction techniques.

## Contents of the repository

Details the organization of this repository and the contents of each folder

- **cont_integration** :
Tools and script for continuous integration (versioning, automated testing and updating documentation)

- **documentation** :
Source code documentation and user documentation

- **examples** :
Examples and tutorials to learn how to use the different functionalities of the library

- **qemist** :
The source code of this quantum chemistry package

## Contents of QEMIST

QEMIST is organized into three layers: problem decomposition, electronic structure solvers, and hardware backends.
The problem decomposition layer is responsible for splitting the input molecule into smaller subproblems and submitting these to one (but conceivably more than one) eigenvalue solver, then processing these results into an output energy.
Some electronic structure solvers use classical methods, while other use wrappers over quantum hardware emulators from various hardware providers.
The quantum solver backend layer implements a common interface over libraries and emulators of quantum hardware.
