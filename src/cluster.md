# Setting up a Firedrake mini-cluster

This document briefly describes how to set up a Firedrake mini-cluster based on
the repository at https://github.com/CCMI-CDT/mini-cluster . The cluster is
based on the original repository at https://github.com/giovtorres/slurm-docker-cluster

## Tested platforms

This mini-cluster has been tested and runs on:

- Apple M3, 24GB RAM, running OS X Tahoe 26.3 with Docker Desktop 4.62.0
- Intel i7-12700, 32GB RAM, running Windows 11 Pro 25H2 with Docker Desktop 4.63.0

It has also been tested on an Intel N200 with 8GB RAM (Surface Go 4) and didn't
reliably run due to insufficient memory to start and actively use the suite of
containers. The 8GB RAM case seemed to be borderline - so 12GB and above 'may
be sufficient' but hasn't been tested.

## Installing Docker and git

You need a working local version of Docker and git. If you're confident setting
this up yourself, please do so.

If you're not confident about install Docker yourself, using Docker Desktop is
recommended, and available from: https://www.docker.com/products/docker-desktop/

Before continuing, make sure the Docker service is running - if you installed
Docker Desktop, you can just run the Docker Desktop instance. If you're using
Windows, and haven't used the Windows Subsystem for Linux (WSL), you may be
prompted to install it by running:

```
wsl --update
```

in a command prompt (terminal).

You can ensure that Docker's working as expected by running the following from
a command prompt:

```bash
docker run hello-world
```

Instructions for installing git are at: https://git-scm.com/install/

## Creating the cluster

To begin with, clone the mini-cluster repository:

```bash
git clone https://github.com/CCMI-CDT/mini-cluster
```

Change into the repository directory, and bring up the cluster:

```bash
cd mini-cluster/
docker compose up
```

This may take some time to download, particularly on slow connections,
but you should eventually see output indicating that the cluster is up
and running, and monitoring its own status.

You can either leave the terminal monitoring output, or press 'd'
to detach and return to the command prompt.

## Interacting with the cluster

To run a shell in the cluster, as a user, run:

```bash
docker exec -u user -it slurmctld bash
```

You should by default find yourself in /data , which is writeable
by you, persistent over cluster restarts, and available to all the
nodes in the cluster that your jobs will run on.

## Building Firedrake

To build firedrake in /data:

```bash
cd /data
python3 -m venv venv-firedrake
. venv-firedrake/bin/activate
export CC=mpicc CXX=mpicxx PETSC_DIR=/usr HDF5_MPI=ON HDF5_DIR=/usr
echo 'setuptools<81' > constraints.txt
export PIP_CONSTRAINT=constraints.txt
pip install --no-binary=h5py --no-binary=mpi4py vtk 'firedrake[check]'
```

There are some longer (on the scale of many minutes) pauses in output
expected during the pip install process - assuming your hardware is
roughly the same as or better than the tested hardware listed at the
beginning of this document, waiting is recommended. If you're unsure,
check your task manager to find out whether you're out of memory and
have containers in a stalled state.

## Verifying Firedrake

Check that firedrake imports without errors:

```bash
python -c "from firedrake import *"
```

Some warnings are expected here, but you shouldn't see errors.

## Verifying the cluster

To check that the cluster is working properly, create a small test:

```bash
cat > hello_mpi.py << EOF
from mpi4py import MPI
import socket

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
hostname = socket.gethostname()

print(f"Hello from rank {rank} of {size} on host {hostname}")
EOF
```

Then run the test:

```bash
srun -c 2 -n 2  python3 hello_mpi.py
```

You should see output of the form:

```
Hello from rank 0 of 1 on host c2
Hello from rank 0 of 1 on host c1
```

