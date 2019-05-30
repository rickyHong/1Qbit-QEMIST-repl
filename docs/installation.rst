Installation Instructions
==================================

Pip
___
The simplest way to install OpenQEMIST is with pip.

1. Follow the setup instructions_ for installing the .NET Core SDK and the Microsoft IQ# module.

2. Install OpenQEMIST with ``pip install openqemist``.

.. _instructions: https://docs.microsoft.com/en-us/quantum/install-guide/?view=qsharp-preview

Installation from source
________________________
To install OpenQEMIST from source, simply clone the GitHub repo and add the package
to your ``PYTHONPATH``. The dependencies for running the project are the Microsoft
.NET Core SDK, IQ#, and qsharp packages; pyscf; numpy; and scipy. The most current
list of dependencies, as well as dependencies for building the documentation can
be found in the project's Dockerfile.
