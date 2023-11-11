# nlbayes-reproducibility
Code and instructions for reproducing the results for the NLBayes manuscript

## Linux setup

### How to install Python and R

To install Python we can use miniconda
```bash
#!/bin/bash

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
and follow the instructions.

To install R in a debian based Linux distribution like Ubuntu, we may use the following command:
```bash
#!/bin/bash

sudo apt install r-base
```
### System dependencies
The following commands will install several system libraries that are needed by the R and Python packages used:
```bash
#!/bin/bash

# Install R packages dependencies
sudo apt install -y \
    libcurl4-openssl-dev libxml2-dev libssl-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev

# Install git, build tools and the GNU Scientific Library
sudo apt install -y git build-essential libgsl-dev
```

### Python environment
It is recommended to create an independent Python environment for each project. Here we use the `conda` package manager.
```bash
#!/bin/bash

# Prepare conda environment
conda create -n nlbayes-reproducibility python=3.11
conda activate nlbayes-reproducibility

# Install NLBayes (python)
pip install cython
pip install git+https://github.com/umbibio/nlbayes-python.git

# Install additional packages used in analysis
pip install statsmodels scikit-learn seaborn ipykernel rpy2 decoupler upsetplot
```

### R environment
```bash
#!/bin/bash

# Prepare and install NLBayes
R -q -e "install.packages(c('rjson', 'Rcpp', 'RcppProgress', 'devtools'))"
R -q -e "devtools::install_github('umbibio/nlbayes-r')"

# Install packages used in analysis
R -q -e "install.packages(c('umap', 'BiocManager', 'ggplot2', 'ggrepel'))"
R -q -e "BiocManager::install(c('org.Hs.eg.db', 'viper', 'aracne.networks', 'GEOquery'))"
```

### Data download
```bash
#!/bin/bash

mkdir -p ./data

wget -O "data/three_tissue.rels.json" "https://umbibio.math.umb.edu/nlbayes/assets/data/networks/gtex_chip/homo_sapiens/tissue_independent/three_tissue.rels.json"
```

## Experiments

### Fig. 3. Simulations

For this experiment we use the command line tool provided in the python package. 
The `nlb-simulation` command automatically generates simulated networks and differential
expression, performs randomization as specified by named parameters and runs the inference algorithm to generate `csv` files with the results. More detailes can be found in the [`simulation.py`](https://github.com/umbibio/nlbayes-python/blob/main/nlbayes/commands/simulation.py) script.
```bash
#!/bin/bash

# run simulation with no randomizations
nlb-simulation --n_graphs 5 --n_replica 3 --net_seed  10 --evd_seed  20 --outdir ./simulations

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

from py_scripts.utils import collect_simulations
from figures import make_figure_3

results, metadata = collect_simulations('./simulations')
make_figure_3(results)
```
A new `PNG` file is saved in the `figures` subfolder.

### Fig. 4. NLBayes vs Viper

First we generate the results for the TF activity inference from the two methods.
The R method `compute.inference.comparison` generates three `CSV` files, one for each experiment (e2f3, c-myc, h-ras).
```R
#!/bin/env R

source('r_scripts/nlbayes_utils.R')
compute.inference.comparison()
```

Once the corresponding `CSV` files are available, we can make the figure.
```R
#!/bin/env R

source('r_scripts/figures.R')
make.figure.4()

```

To assess for agreement between both methods, we construct the corresponding 
confusion matrix and use it as a contingency matrix in a fisher exact test
with the 'greater' alternative hypothesis, i.e., the true odds ratio is greater
than one.
```R
#!/bin/env R

for (exp in c('e2f3', 'myc', 'ras')) {

    # read result
    filename <- paste0('data/oe_',exp,'_on_net_regulonbrca_nlbayes_and_viper.csv')
    df <- read.table(filename, sep = ',', header = TRUE)
    df$viper.pvalue[is.na(df$viper.pvalue)] <- 1

    # compute contingency table
    yy <- sum(df$posterior.p >= 0.2 & df$viper.pvalue <= 0.05)
    yn <- sum(df$posterior.p >= 0.2 & df$viper.pvalue > 0.05)
    ny <- sum(df$posterior.p < 0.2 & df$viper.pvalue <= 0.05)
    nn <- sum(df$posterior.p < 0.2 & df$viper.pvalue > 0.05)
    ctable <- matrix(c(yy, ny, yn, nn), nrow=2, ncol=2)

    # print test
    print(exp)
    print(ctable)
    print(fisher.test(ctable, alternative='greater'))
}
```

### Fig. 5. Comparison to other methods

To make figure 5, the `CSV` files from the previous step are needed.
```python
from figures import make_figure_5

make_figure_5()
```