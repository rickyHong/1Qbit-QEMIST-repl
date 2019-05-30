![logo](docs/img/1qbitlogo.png "1QBit is awesome!")
# OpenQEMIST
[![Build Status](https://travis-ci.com/1QB-Information-Technologies/openqemist.svg?token=zt4rNJ8MTUGcpVsToGyy&branch=master)](https://travis-ci.com/1QB-Information-Technologies/openqemist)

Harnessing the combined power of emerging quantum computing technologies and
state-of-the-art classical techniques, the Quantum-Enabled Molecular ab Initio
Simulation Toolkit, or QEMIST, is 1QBit’s innovative solution to a fundamental
and intractable problem in chemistry: ab initio simulation of molecules.

The accurate prediction of the electronic structure of a molecule is key to the
design of new materials, such as drug compounds and catalyst molecules, by
helping to anticipate a material’s properties before its synthesis in the lab.
However, obtaining this information using classical computers is computationally
intensive, and the resources required for an exact solution scale exponentially
with the size of the problem. Attempts to provide approximate approaches to this
problem on classical computers have been to date either limited to small-sized
systems or compromising on the accuracy of the simulation.

## Installation
To install the package with pip:

1. Follow the setup [instructions](https://docs.microsoft.com/en-us/quantum/install-guide/?view=qsharp-preview) for installing the .NET Core SDK and the Microsoft IQ# module.

2. Install OpenQEMIST with `pip install openqemist`.

To install OpenQEMIST from source, simply clone the GitHub repo and add the package
to your ``PYTHONPATH``. The dependencies for running the project are the Microsoft
.NET Core SDK, IQ#, and qsharp packages; pyscf; numpy; and scipy. The most current
list of dependencies, as well as dependencies for building the documentation can
be found in the [Dockerfile](./docker_images/Dockerfile).

## Getting started

To get started install the package, then see the [Jupyter notebooks](./examples/) for example usage.

## Contents of the repository

Details the organization of this repository and the contents of each folder

- **cont_integration** :
Tools and script for continuous integration (versioning, automated testing and updating documentation)

- **documentation** :
Source code documentation and user documentation

- **examples** :
Examples and tutorials to learn how to use the different functionalities of the library

- **openqemist** :
The python package.

## Architechture of OpenQEMIST

OpenQEMIST is organized into three layers: problem decomposition, electronic
structure solvers, and hardware backends. The problem decomposition layer is
responsible for splitting the input molecule into smaller subproblems and
submitting these to one (but conceivably more than one) eigenvalue solver, then
processing these results into an output energy. Some electronic structure
solvers use classical methods, while other use wrappers over quantum hardware
emulators from various hardware providers. The quantum solver backend layer
implements a common interface over libraries and emulators of quantum hardware.

Currently DMET is the only implemented problem decomposition.
VQE, FCI, and CCSD are the implemented electronic structure solvers.

## Contributing
We welcome contributions to OpenQEMIST! Please open an issue or submit a pull request on GitHub to start the process.

## Citing
If you use OpenQEMIST in your research, please cite

Takeshi Yamazaki, Shunji Matsuura, Ali Narimani, Anushervon Saidmuradov, and Arman Zaribafiyan "Towards the Practical Application of Near-Term Quantum Computers in Quantum Chemistry Simulations: A Problem Decomposition Approach" Published on [arXiv](https://arxiv.org/abs/1806.01305) on Jun 4, 2018.


Copyright 1QBit 2019. This software is released under the Apache Software License version 2.0.
