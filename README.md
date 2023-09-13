# Biodivine library for control of parametrised (partially specified) Boolean networks

A library to solve one-step, temporary and permanent source-target control of parametrised (partially-specified) Boolean networks.

The directory structure:

    .
    ├── auxiliary_scripts    # Scripts to do & process experiments
    ├── models               # Base experimental models
    ├── results              # Raw measured results from experiments
    ├── results_simple       # Raw measured results from experiments for simplified phenotype control procedure
    └── src                  # Library source code

### Auxiliary scripts

- `analyse_results.py` - A script showing quick statistics about the obtained experiment results
- `networks_sampler.py` - A script generating partially-specified samples of witness models
- `plot_results.ipynb` - A Jupyter notebook for visualization of the experiment results
- `run_groups.py` - A script for obtaining the experiment results, running the methods from library on the generated methods. Allows timeout specification.  

### Models

Base models for testing the library. Contains witness models from CellCollective platform and some parametrised version of the models.

### Results

The raw unprocessed outputs of experiments for both performance comparison and robustness metric.

### src

Source code of the library. Consists of following rust modules:

#### aeon module

Operations to perform base state-transition graph manipulations.

### control

Implementations of control algorithm on the perturbable graph.

### perturbation

Data structure representing state transition graph of Boolean network which is viable for perturbations.  

To run the basic experiments, execute `cargo run --release`
