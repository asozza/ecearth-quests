# TUNING

Here we try to report the workflow behind the tuning.

## Forcings:
We added option in the `config.yaml` for selecting forcings (pre-industrial, historical, fixed year) and passing the information to `generate_job.py`.

## Single run with tuning
We added option for activating tuning for a single run, passing `tuning_params.yml` to `generate-job.py`. 
Remember also to activate `ecmean`.

## Ensemble runs
For running an ensemble of runs, we have a new python executable `tuning_ensemble.py`  for creating an ensemble of simulations with perturbed parameters.

## ECTUNER
We then run ECtuner (by Jost) to create sensivity and create tuned parameters. We tried with a set of 20 simulations (`lr**`) with pre-industrial forcing.