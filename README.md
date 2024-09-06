# ProkODE: A Computational Model to Unpakage Time-Series Protein Concentrations during Prokaryotic Development
Inputs:
- full genome sequence
- locations of annotated genes (indices within genome)
- args:
  - time frame to predict
  - list of proteins to show results

Outputs:
- graphs depicting protein concentrations across time frame (of listed proteins)

## Algorithm

\* first layer = file; second layer = function within file

1. /preprocessing/get_tfbs_decay.py
2. /diffeq_plot.py
  1. config_network_json() -> writes file that can be processed into individual ODEs (in plot_system())
  2. plot_system() -> computes entire ode system and plots graphs of desired protein concentrations


## Data Extraction Algorithms
### For Beta vals
