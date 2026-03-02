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

## A primer on using cluster commands

This mini-cluster runs [Slurm](https://slurm.schedmd.com/), a workload manager
for HPC clusters. Slurm publish a [quickstart guide](https://slurm.schedmd.com/quickstart.html)
which is highly recommended for news users, but those wanting no more than the 
minimal operating commands for this cluster should be aware of the following:

### sinfo

View information about the cluster. Usage:

```
bash-5.1$ sinfo
PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST
normal*      up   infinite      2   idle c[1-2]
gpu          up   infinite      0    n/a 
bash5.1$
```

The output fields are:

* PARTITION : a collection of resources for running your work
* AVAIL: status of the partition
* TIMELIMIT: how long a single job of work is allowed to run for
* NODES: discrete physical or logical resources that run your work
* STATE: whether the resource is allocated to work or not
* NODELIST: the names of the resources in the partition

Note that this cluster has both normal and gpu parititions which jobs can be
submitted to, but by default the gpu partition is unconfigured.

Full documentation at: [https://slurm.schedmd.com/sinfo.html](https://slurm.schedmd.com/sinfo.html)

### squeue

View information about jobs. Usage:

```
bash-5.1$ squeue
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
                 4    normal    sleep     user  R       0:07      2 c[1-2]
bash-5.1$
```

Jobs, here, are units of work that you or other users have submitted to be run
by the cluster.

The output fields are:

* JOBID: an identifier for a particular job
* PARTITION: the collection of resources the job has been submitted to run on
* NAME: a human-readable name for the job
* USER: the submitter of the job
* ST: the job state - R indicates 'running', PD indicates 'pending'.
* TIME: how long the job has run for 
* NODES: how many discrete physical or logical resources are assigned to the job
* NODELIST: names of the resources assigned to the job

Full documentation at: [https://slurm.schedmd.com/squeue.html](https://slurm.schedmd.com/squeue.html)

### srun

Run an interactive job, where your terminal remains connected.

A minimal case, which returns the hostname of the node it runs on:

```
bash-5.1$ srun hostname
c1
bash-5.1$
```

A slightly more detailed case, running on two nodes (-n 2) and requesting two
cores on each node (-c 2):

```
bash-5.1$ srun -n 2 -c 2 hostname
c1
c2
bash-5.1$
```

Full documentation at: [https://slurm.schedmd.com/srun.html](https://slurm.schedmd.com/srun.html)

### sbatch

Run a non-interactive job, where the terminal returns immediately. Usage:

```
bash-5.1$ sbatch job.sh
Submitted batch job 11
bash-5.1$
```

Where 'job.sh' is a script that defines the job, such as:

```
#!/bin/bash

### Slurm specifications for the job
# Name of this job
#SBATCH --job-name=helmholtz
# Files to write stdout and stderr; %j is the job id
#SBATCH --output=helmholtz_%j.out
#SBATCH --error=helmholtz_%j.err
# Nodes requested for this job
#SBATCH --nodes=1
# CPUs requested per node
#SBATCH --ntasks-per-node=1       # Adjust to match your node's core count
# Time requested for this job
#SBATCH --time=00:05:00
# Partition the job will submit to
#SBATCH --partition=normal      # Change to your cluster's partition name

### Job commands

# Acctivate a venv
. /data/venv-firedrake/bin/activate

# Run the program
srun python3 helmholtz.py
```

sbatch should return immediately with the job ID of your submitted job.

The example above is the simplest resource allocation available - one core on
one node. Increase both parameters for parallel jobs, and be realistic with
your estimate of time. If your job exceeds the time you've specified it will
be terminated - this can be useful to catch stalled jobs, but frustrating if
your job runs slightly slower than expected and doesn't run to completion due
to being terminated by slurm.

Note in the above case - the batch script allows much more flexible and complex
work than the interactive 'srun', such as activating a venv.

Full documentation at: [https://slurm.schedmd.com/sbatch.html](https://slurm.schedmd.com/sbatch.html)

### scancel

Send signals to jobs, by default cancelling them. Usage:

```
bash-5.1$ scancel 13
bash-5.1$ 
```

In this case, the job cancelled has ID '13'. Use 'squeue' to find the ID of the
job you wish to cancel.

Full documentation at: [https://slurm.schedmd.com/scancel.html](https://slurm.schedmd.com/scancel.html)




