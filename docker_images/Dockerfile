FROM fedora:28

RUN dnf -y update
RUN dnf -y install wget libgomp openblas-devel pandoc
RUN dnf clean all


# Openfermion and ProjectQ installation dependencies
RUN dnf -y install gcc redhat-rpm-config gcc-c++ python3-devel

# We do this because pyscf refuses to link to the system OpenBLAS in the pip
# install so we force load that version ahead of the one that it ships with.
ENV LD_PRELOAD=/usr/lib64/libopenblas.so

# Python modules for documentation, Jupyter notebook support, visualization
# and some computational packages
RUN pip3 install 'ipython==7.4.0' jupyter 'pyscf==1.6' 'openfermion==0.9.0' \
    openfermionprojectq numpy 'scipy==1.2' setuptools wheel sphinx py3Dmol \
    sphinx_rtd_theme nbsphinx


################### Microsoft QDK environment ###################

# Install QDK
WORKDIR /tmp/
RUN dnf clean all
RUN rpm --import https://packages.microsoft.com/keys/microsoft.asc
RUN wget -q https://packages.microsoft.com/config/fedora/27/prod.repo
RUN mv prod.repo /etc/yum.repos.d/microsoft-prod.repo
RUN chown root:root /etc/yum.repos.d/microsoft-prod.repo
RUN dnf install -y dotnet-sdk-2.2

# Install the IQSharp backend
RUN dotnet tool install -g Microsoft.Quantum.IQSharp --version 0.7.1905.3003-beta
RUN /root/.dotnet/tools/dotnet-iqsharp install --user

# Install python package
RUN pip3 install qsharp==0.7.1905.3003b1

########################### Finalize #############################

# Set the python path to find openqemist
ENV PYTHONPATH=$PYTHONPATH:/root/openqemist

WORKDIR /root/
