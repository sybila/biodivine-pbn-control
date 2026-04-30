# Biodivine library for control of parametrised (partially specified) Boolean networks

A library to solve one-step, temporary and permanent source-target control of parametrised (partially-specified) Boolean networks.

The directory structure:

    .
    ├── Cargo.toml
    ├── Cargo.lock
    ├── README.md
    ├── auxiliary_scripts/    # Scripts to do & process experiments
    ├── models/               # Base experimental models
    ├── results/              # Raw measured results from experiments
    └── src/                  # Library source code
        ├── aeon/             # Simplified algorithms taken over from aeon-server
        ├── bin/              # Alternative console entry-points
        ├── control/          # Source-target control algorithms
        ├── perturbation/     # Perturbed state-transition graph implementation
        ├── phenotype_control # Phenotype control algorithms
        ├── lib.rs            # Library declaration
        └── main.rs           # Main console entry point

### Auxiliary scripts

- `analyse_results.py` - A script showing quick statistics about the obtained experiment results
- `networks_sampler.py` - A script generating partially-specified samples of witness models
- `plot_results.ipynb` - A Jupyter notebook for visualization of the experiment results
- `run_groups.py` - A script for obtaining the experiment results, running the methods from library on the generated methods. Allows timeout specification.  
- `run_phenotype.py` - A script to run phenotype-control experiments on parametrised Boolean-network models; supports timeouts and exports results for later analysis.

### Models

Base models for testing the library. Contains witness models from CellCollective platform and some parametrised version of the models.

### results

The raw unprocessed outputs of experiments for both performance comparison and robustness metric.

### src

Source code of the library. Consists of following rust modules:

#### aeon module

Operations to perform base state-transition graph manipulations.

### control module

Implementations of control algorithm on the perturbable graph.

### perturbation module

Data structure representing state transition graph of Boolean network which is viable for perturbations.  

### phenotype_control module

Implementations of phenotype control algorithm on the perturbable graph.


To run the basic experiments, execute `cargo run --release`