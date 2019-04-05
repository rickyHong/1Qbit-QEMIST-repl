FROM fedora:28

RUN dnf -y update
RUN dnf -y install wget ipython3 libgomp openblas-devel

# We do this because pyscf refuses to link to the system OpenBLAS in the pip
# install so we force load that version ahead of the one that it ships with.
ENV LD_PRELOAD=/usr/lib64/libopenblas.so

RUN pip3 install jupyter 'pyscf==1.6' 'openfermion==0.9.0' numpy scipy setuptools \
    wheel sphinx py3Dmol sphinx_rtd_theme
# This is needed because of a bug in jupyter.
# https://github.com/jupyter/jupyter_console/issues/158
RUN pip3 install 'prompt-toolkit==1.0.15'

# Install and configure the Q# python interopts
# Configure DNF repo for dotnet
WORKDIR /tmp/
RUN rpm --import https://packages.microsoft.com/keys/microsoft.asc
RUN wget -q https://packages.microsoft.com/config/fedora/27/prod.repo
RUN mv prod.repo /etc/yum.repos.d/microsoft-prod.repo
RUN chown root:root /etc/yum.repos.d/microsoft-prod.repo
RUN dnf install -y dotnet-sdk-2.2

# Install and configure the Q# beta stack
WORKDIR /root/ms_quantum
ADD ms_quantum.tar.gz /root/ms_quantum
RUN ls
RUN dotnet tool install --configfile NuGet.config -g Microsoft.Quantum.IQSharp
RUN pip3 install qsharp 
# We need the full path here because this isn't added to the path until the
# shell reloads.
RUN /root/.dotnet/tools/dotnet-iqsharp install

# Set the python path to find qemist
ENV PYTHONPATH=$PYTHONPATH:/root/QEMIST

WORKDIR /root/
