# Jug-ified Exposons Pipleine

Using the library [jug](https://github.com/luispedro/jug), I've created a parallel execution framework for computing exposons on large datasets (really the type of dataset they're most useful for).

## Installation

You should be able to install prerequisites (`mdtraj`, `enspara`, `jug`, etc.) with `requirements.txt`.

First, set up a virtual envronment:

```bash
virtualenv exposons
source exposons/bin/activate
```

```bash
pip install -r requirements.txt
```

If you're on MacOS, you'll need to clone and install `MDTraj` manually from `master`.


## Configuration

There's a configuration file, `proteins.json` in this repository. It's a dictionary mapping target name to configurations for that target:


```json
"protein_name": {
	"trajectories": [
		"/path/to/trajectory/trajectory00.xtc"
		"/path/to/trajectory/trajectory01.xtc",
		"/path/to/trajectory/trajectory02.xtc",
		"/path/to/trajectory/trajectory03.xtc",
	],
	"topology": "/path/to/trajectory/prot_masses.pdb",
	"cluster_radii": ["3.0"],
	"cluster_distance_metric": "euclidean",
	"lag_time": 150,
	"model_cluster_radius": "3.0",
	"kmedoids_updates": 10
}
```

This file is a dictionary of the form `target_name`: `{configurations}`. Each target needs the following keys set:

- `trajectories`: a list of paths to trajectories
- `topology`: the path to the topology (as a pdb, typically) that corresponds to `trajectories`
- `cluster_radii`: the cluster radii to choose
- `cluster_distance_metric`: the distance to use between SASA vectors during cluster
- `lag_time`: the lag time to build the msm
- `model_cluster_radius`: the radius for k-centers to cluster on
- `kmedoids_updates`: the number of kmedoids refinements to do.

## Checking your configuration

To check to see you've configured everything correctly, you can just run

```bash
jug status sasa-msms.py
```

That'll print something that looks like:

```
 Failed    Waiting      Ready   Complete     Active  Task name
----------------------------------------------------------------------------------------------
      0          1          0          0          0  lib.tasks.assemble_sasa_h5
      0          4          0          0          0  lib.tasks.condense_sparse_sidechain_sasas
      0          0          4          0          0  lib.tasks.map_sasa_sparse
..............................................................................................
      0          5          4          0          0  Total
```

Here, there are only four trajectories (you'll probably have more than that). There are more steps to the pipeline, but they'll appear later.

There are several other useful jug commands, which you can find at jug's documentation. In particular, I recommend `jug graph` to get a sense of how tasks relate to each other.

## Simple (All-at-Once) Execution

To run the pipeline serially, you can just run

```bash
jug execute sasa-msms.py
```

This will load each task (first `lib.tasks.map_sasa_sparse`), execute it, then write the result to disk, and start with the next task. This will continue until there are no more tasks to run.

The magic to jug is that you can run tasks in parallel though! So, for instance, with [GNU Parallel](https://www.gnu.org/software/parallel/) you could run:

```bash
parallel -N0 nice jug execute sasa-msms.py  ::: {1..4}
```

To run four processes, each running different tasks at a time.

## Complex (Step-by-Step) Execution

The reason to execute the tasks one-by-one is that you can apply the correct resouces

The tasks in the pipeline are:

- `lib.tasks.map_sasa_sparse` - compute the SASA of each atom in each trajectory and store the result as a sparse, compressed array (`npz` format). Parallelization uses `OMP_NUM_THREADS`, typically all availiable on linux, or 1 on macos (since clang doesn't support openmp).
- `lib.tasks.condense_sparse_sidechain_sasas` - condense the SASA of each atom into a single SASA for each residue, threads used: 1.
- `lib.tasks.assemble_sasa_h5` - combine each trajectory's residue SASA into a single `enspara` `RaggedArray`.
- `lib.tasks.cluster_features` - use k-medoids to cluster the set of 
- `lib.tasks.implied_timescales` - compute implied timescales plots. Parallelization uses `OMP_NUM_THREADS`, typically all availiable on linux, or 1 on macos (since clang doesn't support openmp).
- `lib.tasks.write_struct_ctrs` - write each center structure (from `cluster_features`) to a single trajectory on disk.

Depending on the task, on a high-performance computing clusters with slurm, I start with the following:

```bash
sbatch --array=0-9 --job-name sasa --exclusive --mem 8G --wrap 'jug execute sasa-msms.py --aggressive-unload --target lib.tasks.map_sasa_sparse'
```

This asks for 10 (`--array-0-9`) full (`--exclusive`) nodes each with 8GB (`--mem 8G`) or more or memory. Each will load one of the `lib.tasks.map_sasa_sparse` tasks, run it, and save the result, and then purge the result from memory (`--aggressive-unload`).

Then, because it's pretty low-compute, I'll do something like this for `lib.tasks.condense_sparse_sidechain_sasas` :

```bash
parallel -N0 nice jug execute sasa-msms.py  --target lib.tasks.condense_sparse_sidechain_sasas --aggressive-unload ::: {1..8}
```

There's only one `lib.tasks.assemble_sasa_h5` job per target, so that's easy:

```bash
jug execute sasa-msms.py --target lib.tasks.assemble_sasa_h5
```

Clustering is a bit more complicated. If you want to use more than one node (with MPI) you need MPI installed (e.g. MPICH) as well as `mpi4py`, which isn't installed by default in `requriements.txt`. If you don't need MPI, queue your clustering job with something like:

```bash
sbatch --job-name exposon-cluster --exclusive --wrap 'jug execute sasa-msms.py --target lib.tasks.cluster_features'
```

if you're using MPI, you can say something like:

```bash
sbatch --job-name exposon-cluster --exclusive --nodes=4 --cpus-per-task 24 --wrap 'jug execute sasa-msms.py --target lib.tasks.cluster_features'
```

to ask for 4 nodes (`--nodes=4`) with 24 cores each (`--cpus-per-task 24`). This will, of course, differ based on your queueing system
