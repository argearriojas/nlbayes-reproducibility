# nlbayes-reproducibility
Code and instructions for reproducing the results for the NLBayes manuscript

## Linux setup

### System dependency
```bash
# Install GNU Scientific Library
sudo apt install -y libgsl-dev
```

### Python environment
```bash
# Prepare conda environment
conda create -n nlbayes-reproducibility python=3.11
conda activate nlbayes-reproducibility
pip install cython
pip install git+https://github.com/umbibio/nlbayes-python.git
pip install statsmodels scikit-learn seaborn ipykernel
```

### R environment
```bash
# Install system dependencies
sudo apt install -y \
    libcurl4-openssl-dev libxml2-dev libssl-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev

# Prepare R
R -q -e "install.packages(c('rjson', 'Rcpp', 'RcppProgress', 'devtools'))"
R -q -e "devtools::install_github('umbibio/nlbayes-r')"
```

## Experiments

### Fig 3. Simulations

For this experiment we use the command line tool provided in the python package. 
The `nlb-simulation` command automatically generates simulated networks and differential
expression, performs randomization as specified by named parameters and runs the inference algorithm to generate `csv` files with the results. More detailes can be found in the [`simulation.py`](https://github.com/umbibio/nlbayes-python/blob/main/nlbayes/commands/simulation.py) script.
```bash
#!/bin/bash

# run simulation with no randomizations
nlb-simulation --n_graphs 5 --n_replica 3 --net_seed  10 --evd_seed  20 --evd_rnd_p 0.00 --outdir ./simulations

# run simulations for data randomization experiments
nlb-simulation --n_graphs 5 --n_replica 3 --net_seed  30 --evd_seed  40 --evd_rnd_p 0.25 --outdir ./simulations
nlb-simulation --n_graphs 5 --n_replica 3 --net_seed  50 --evd_seed  60 --evd_rnd_p 0.50 --outdir ./simulations
nlb-simulation --n_graphs 5 --n_replica 3 --net_seed  70 --evd_seed  80 --evd_rnd_p 0.75 --outdir ./simulations
nlb-simulation --n_graphs 5 --n_replica 3 --net_seed  90 --evd_seed 100 --evd_rnd_p 1.00 --outdir ./simulations

# run simulations for network randomization experiments
nlb-simulation --n_graphs 5 --n_replica 3 --net_seed  30 --evd_seed  40 --net_rnd_p 0.25 --outdir ./simulations
nlb-simulation --n_graphs 5 --n_replica 3 --net_seed  50 --evd_seed  60 --net_rnd_p 0.50 --outdir ./simulations
nlb-simulation --n_graphs 5 --n_replica 3 --net_seed  70 --evd_seed  80 --net_rnd_p 0.75 --outdir ./simulations
nlb-simulation --n_graphs 5 --n_replica 3 --net_seed  90 --evd_seed 100 --net_rnd_p 1.00 --outdir ./simulations
```
Here we used only 3 replicas per experiment for brevity.

After generating the results files, we can process them and make Figure 3 with the methods `collect_simulations` and `make_figure_3`.
```python
#!/bin/env python

from helpers import collect_simulations
from figures import make_figure_3

results, metadata = collect_simulations('./simulations')
make_figure_3(results)
```

