# ecearth-quests

Tools to set up, launch and monitor [EC-Earth4](https://ec-earth-4-docs.readthedocs.io/en/latest/) 
runs on HPC systems with SLURM. The repository collects the scripts used day to day to
generate and duplicate jobs, compile and update the model, build tuning
ensembles, follow the performance of running experiments, handle NEMO restarts,
and prepare NEMO grids and input fields.

Each script can be run directly from its folder.

## Installation

Clone the repository and create the conda environment:

```bash
conda env create -f environment.yaml
conda activate ecearth-quests
```

On HPC systems, load the required modules first with `. load_modules.sh`.

## Launching runs

The `launch/` folder contains the tools to prepare and submit EC-Earth4 runs:

| Script | Purpose |
|---|---|
| `generate_job.py` | Generate a job script from an experiment configuration |
| `duplicate_job.py` | Duplicate an existing job for a new experiment |
| `create_tuning_ensemble.py` | Build an ensemble of runs from a tuning parameter file |

Experiment configurations are in `configs/experiments/`. Start from
`config-example.yml` and keep your own `config.yml` local (it is ignored by git).
Tuning parameter files are in `configs/tuning/`.

## Other tools

| Folder | Content |
|---|---|
| `monitor/` | Run performance and HPC resource usage (SYPD, memory, TRES) |
| `restart/` | Rebuild, roll back and pack restart files |
| `scripts/` | various scripts for compile and update EC-Earth code |
| `tools/`   | NEMO and OIFS grids and fields (e.g. grid bounds, vertical interpolation) |

## Documentation

- [Tuning workflow](docs/tuning.md)
- [NEMO domain decomposition](docs/domain_decomposition.md)
- [NEMO tools](docs/nemo_tools.md)
- [Scaling tests](docs/scaling/README.md)

## License

See [LICENSE](LICENSE).