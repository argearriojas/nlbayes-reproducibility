# NLBayes Reproducibility
This repository contains code and instructions for reproducing the results
of the paper titled: *"A Bayesian Noisy Logic Model for Inference of Transcription
Factor Activity from Single Cell and Bulk Transcriptomic Data"*
([Arriojas et. al. 2023](https://doi.org/10.1101/2023.05.03.539308)).

Please explore the notebooks [`figures.ipynb`](figures.ipynb)
and [`figures_with_outputs.ipynb`](figures_with_outputs.ipynb) to learn how to 
generate the corresponding files.

## Linux setup

### How to install Python and R

To install Python we can use miniconda
```bash
#!/bin/bash

# Download miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh -h
# usage: Miniconda3-latest-Linux-x86_64.sh [options]
# Installs Miniconda3 py311_23.9.0-0
# -b           run install in batch mode (without manual intervention),
#              it is expected the license terms (if any) are agreed upon

# Install miniconda
bash Miniconda3-latest-Linux-x86_64.sh -b

# Initialize the conda installation and make it available to use
$HOME/miniconda/bin/conda init
source $HOME/.bashrc
```

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
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    libudunits2-dev libgdal-dev

# Install git, build tools and the GNU Scientific Library
sudo apt install -y git build-essential libgsl-dev
```
### Data download and initialization
```bash
#!/bin/bash

Rscript download_data.R
Rscript init.R
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
# For R version <= 4.1, there may be broken dependencies related to ggpattern
#     you may consider installing the following specific package versions.
R -q -e "remotes::install_version('ggplot2', version='3.3.6', upgrade=FALSE)"
R -q -e "remotes::install_version('ggh4x', version='0.2.1', upgrade=FALSE)"
R -q -e "remotes::install_version('ggpattern', version='0.4.2', upgrade=FALSE)"
# Otherwise, you may do
R -q -e "install.packages(c('ggplot2', 'ggh4x', 'ggpattern'))"

R -q -e "install.packages(c('umap', 'BiocManager', 'ggrepel', 'mvtnorm', 'Seurat'))"
R -q -e "BiocManager::install(c('org.Hs.eg.db', 'viper', 'aracne.networks', 'GEOquery', 'glmGamPoi', 'clusterProfiler'), update=FALSE)"
```
