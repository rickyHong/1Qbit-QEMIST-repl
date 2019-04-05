# QEMIST
[![Build Status](https://travis-ci.com/1QB-Information-Technologies/openqemist.svg?token=zt4rNJ8MTUGcpVsToGyy&branch=master)](https://travis-ci.com/1QB-Information-Technologies/openqemist)

This library aims at pushing the boundaries of advanced materials science by leveraging emerging quantum computing technologies together with state-of-the-art classical techniques in a hybrid framework.

It allows for the simulation of intractable electronic structure problems through problem decomposition and problem reduction techniques.

## Installation

The package can be installed with pip or directly from source.
Add the python package to your ``PYTHONPATH``.
The dependencies for running the project are pyscf, numpy, and scipy.
For the most current list of dependencies, see the [Dockerfile](./docker_images/Dockerfile).

## Getting started

To get started install the package, then see the [Jupyter notebook](./examples/end_to_end.ipynb) for example usage.

## Contents of the repository

Details the organization of this repository and the contents of each folder

- **cont_integration** :
Tools and script for continuous integration (versioning, automated testing and updating documentation)

- **documentation** :
Source code documentation and user documentation

- **examples** :
Examples and tutorials to learn how to use the different functionalities of the library

- **qemist** :
The python package.

## Architechture of QEMIST

QEMIST is organized into three layers: problem decomposition, electronic structure solvers, and hardware backends.
The problem decomposition layer is responsible for splitting the input molecule into smaller subproblems and submitting these to one (but conceivably more than one) eigenvalue solver, then processing these results into an output energy.
Some electronic structure solvers use classical methods, while other use wrappers over quantum hardware emulators from various hardware providers.
The quantum solver backend layer implements a common interface over libraries and emulators of quantum hardware.

Currently DMET is the only implemented problem decomposition.
VQE, FCI, and CCSD are the implemented electronic structure solvers.

## Contributing
We welcome contributions to QEMIST! Please open an issue or submit a pull request on GitHub to start the process.

## Citing
If you use QEMIST in your research, please cite

Takeshi Yamazaki, Shunji Matsuura, Ali Narimani, Anushervon Saidmuradov, and Arman Zaribafiyan "Towards the Practical Application of Near-Term Quantum Computers in Quantum Chemistry Simulations: A Problem Decomposition Approach" Published on [arXiv](https://arxiv.org/abs/1806.01305) on Jun 4, 2018.
