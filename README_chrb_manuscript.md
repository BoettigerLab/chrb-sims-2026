# Simulations for the chrb manuscript

This branch (`chrb-manuscript`) holds the simulation code used for the polymer simulations of the chrb
manuscript, with example configs for every simulated condition and the environment the runs used.

The code builds on [polychrom](https://github.com/open2c/polychrom) (Open2C) and its loop-extrusion examples,
and runs on [OpenMM](https://openmm.org). The simulations used polychrom 0.1.1 at commit `320f25e`
(16 April 2025), unmodified; the lab copy of that version, with its history, is [lab polychrom repository].

## Contents

| Path | What it is |
|---|---|
| `simulation/main.py` | entry point: 1D loop extrusion, then 3D molecular dynamics, for one config |
| `examples/chrb_manuscript/config_<condition>.yaml` | settings of each simulated condition (below) |
| `environments/environment.yml` | minimal conda environment with the package versions of the runs |
| `environments/environment_as_run_windows.yml` | full export of the Windows environment the runs used |

## Install

```
conda env create -f environments/environment.yml
conda activate chrb-sims
git clone <lab polychrom repository> polychrom
```
Put the polychrom source folder on the Python path, as in the runs:
- Windows (PowerShell): `$env:PYTHONPATH = "<path>\polychrom"`
- Linux / macOS: `export PYTHONPATH=<path>/polychrom`

`pip install -e <path>/polychrom` also works where a C++ compiler is available; it builds an extension the
simulation does not use.

## Run

```
python simulation/main.py --config examples/chrb_manuscript/config_WT.yaml
```
- Output goes to `out_dir/condition_name/<replicate>/`. `all_conformations.npy` holds the 3D coordinates,
  shape (451, 757, 3): frame 0 is the end of equilibration, then one frame per loop-extrusion step.
- Any setting can be changed on the command line, e.g.
  `--overrides num_replicates=10 start_idx_replicate=0 gpu_device='"1"'`.
- One 450-step replicate takes about 2 min on an NVIDIA T1000.

## Conditions

All configs are the base model with the changes below. Monomer positions are 0-based; the enhancer readout is
monomer 260 and the promoter readout monomer 504.

| Config | Condition | Change from WT |
|---|---|---|
| `config_WT.yaml` | WT | none |
| `config_dNIPBL.yaml` | ΔNIPBL | `num_cohesin: 0` |
| `config_dCBS.yaml` | ΔCBS | CBS (261) capture probability 0 |
| `config_dTether.yaml` | ΔTether | enhancer tether deleted: no attraction, no barrier at 260 |
| `config_dTether_dNIPBL.yaml` | ΔTether + ΔNIPBL | no cohesin, no attraction |
| `config_wapl.yaml` | wapl²/+ | cohesin bound lifetime 200 LE steps |
| `config_wapl_dCBS.yaml` | wapl²/+ ΔCBS | both changes |
| `config_wapl_dTether.yaml` | wapl²/+ ΔTether | both changes |

Each example config loads to the same simulation settings as the config used for the paper's runs; only the
bookkeeping fields (output folder, condition name, number of replicates, GPU, plotting) differ. The ΔTether +
ΔNIPBL runs used a config with other values for the barrier and sticky-block fields; with no cohesin and no
attraction these have no effect.
