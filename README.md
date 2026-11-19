# Polymer simulations of loop extrusion and enhancer-promoter tethering

This repository holds the simulation code used for the polymer simulations of the chrb manuscript
([citation]). It simulates a chromatin fibre as a bead-spring polymer in which loop-extruding cohesin and
sticky tether elements act together. Cohesin moves along the chain in a one-dimensional lattice model, and the
3D polymer is evolved by Langevin dynamics on the GPU between extrusion steps.

The code builds on the [polychrom](https://github.com/open2c/polychrom) polymer-simulation library of the
Open2C community and its loop-extrusion examples, and runs on the GPU molecular dynamics toolkit
[OpenMM](https://openmm.org). The simulations used polychrom 0.1.1 at commit `320f25e` (16 April 2025),
unmodified. The lab copy of that version, with its history, is [lab polychrom repository].

The branch `chrb-manuscript` (tag `chrb-manuscript-version`) is the code as run for the paper, plus the
example configs and environment files described below.

## Organisation

| Path | Contents |
|---|---|
| `simulation/main.py` | Entry point. Reads one config, validates it, and runs every replicate: loop extrusion, then molecular dynamics. Writes `config_resolved.json` beside the output. |
| `simulation/sim_config.py` | Every simulation setting, with its type, default and checks (pydantic model). Unknown keys are refused. |
| `simulation/config_loader.py` | Loads a YAML or JSON config, applies `--overrides`, and refuses to reuse an output folder made with different settings. |
| `simulation/run_sim_LE.py` | One-dimensional loop-extrusion model (numba): cohesin loading, two-sided extrusion, unloading, collisions between cohesins, and capture and hold at barrier sites. |
| `simulation/run_sim.py` | Runs the loop-extrusion step for a replicate, samples it at the frames passed to the 3D simulation, and optionally plots it. |
| `simulation/run_sim_MD.py` | Three-dimensional simulation with polychrom and OpenMM: chain bonds, bending stiffness, soft excluded volume, attraction between tether blocks, spherical confinement, and the cohesin bonds updated at each extrusion step. |
| `simulation/md_backend_slurm.py` | Optional: submits the molecular dynamics to a SLURM cluster instead of running locally. |
| `examples/chrb_manuscript/` | One config per condition simulated in the paper (below). |
| `examples/simulation/`, `debug/` | General examples (Rouse chain, crumpled polymer with and without extrusion) and a debugging config; not used in the paper. |
| `tests/` | Unit and integration tests (run with `pytest tests`). |
| `tutorial/` | Notebook and MATLAB script for loading simulation output; not used in the paper. |
| `environments/` | Conda environments (below). |

## Install

```
conda env create -f environments/environment.yml
conda activate chrb-sims
git clone <lab polychrom repository> polychrom
```
Put the polychrom source folder on the Python path, as in the runs:
- Windows (PowerShell): `$env:PYTHONPATH = "<path>\polychrom"`
- Linux / macOS: `export PYTHONPATH=<path>/polychrom`

`pip install -e <path>/polychrom` also works where a C++ compiler is available. It builds an extension that the
simulation does not use.

`environments/environment.yml` is a minimal conda-forge environment with the package versions the runs used
(Python 3.8, OpenMM 8.0.0, numpy 1.24.3, scipy 1.10.1, h5py 3.9.0, pydantic 2.10.6; numba 0.58.1, the nearest
conda-forge build to the 0.58.0 used). `environments/environment_as_run_windows.yml` is the full export of the
Windows environment the runs used (Windows Server 2019, CUDA 10.2, NVIDIA driver 497.09). It is a general-purpose
environment, and most of its packages are not used here.

## Run

```
python simulation/main.py --config examples/chrb_manuscript/config_WT.yaml
```
- Output goes to `out_dir/condition_name/<replicate>/`. `all_conformations.npy` holds the 3D coordinates in
  monomer diameters, shape (451, 757, 3): frame 0 is the end of equilibration, then one frame per
  loop-extrusion step. `LEFPositions_0.h5` holds the cohesin positions.
- Any setting can be changed on the command line, e.g.
  `--overrides num_replicates=10 start_idx_replicate=0 gpu_device='"1"'`.
- An interrupted run resumes from its last saved block when the same command is run again.
- One 450-step replicate, including equilibration, takes about 2 min on an NVIDIA T1000.

## Conditions simulated in the paper

Each config in `examples/chrb_manuscript/` is the base model with the change listed. Monomer positions are
0-based; the enhancer readout is monomer 260 and the promoter readout monomer 504. In the paper one
loop-extrusion step is 6 s, frame 0 is the start of nuclear cycle 14 and frame 450 is gastrulation.

| Config | Condition | Change from WT |
|---|---|---|
| `config_WT.yaml` | WT | none |
| `config_dNIPBL.yaml` | ΔNIPBL | no cohesin |
| `config_dCBS.yaml` | ΔCBS | CBS (monomer 261) capture probability 0 |
| `config_dTether.yaml` | ΔTether | enhancer tether deleted: no attraction, no barrier at 260 |
| `config_dTether_dNIPBL.yaml` | ΔTether + ΔNIPBL | no cohesin, no attraction |
| `config_wapl.yaml` | wapl²/+ | cohesin bound lifetime 200 steps |
| `config_wapl_dCBS.yaml` | wapl²/+ ΔCBS | both changes |
| `config_wapl_dTether.yaml` | wapl²/+ ΔTether | both changes |

Each example config loads to the same simulation settings as the config used for the paper's runs. Only the
bookkeeping fields (output folder, condition name, number of replicates, GPU, plotting) differ. The ΔTether +
ΔNIPBL runs used a config with other values in the barrier and sticky-block fields; with no cohesin and no
attraction these have no effect.

## Licence and credit

The repository is distributed under the MIT licence (`LICENSE`). `simulation/run_sim_MD.py` is adapted from the polychrom examples and is distributed under polychrom's MIT
licence, reproduced in its header. Code by Tee Udomlumleart and the Boettiger lab, Stanford University.
