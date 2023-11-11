#!/bin/bash

mkdir -p ./data

wget -O "data/three_tissue.rels.json" "https://umbibio.math.umb.edu/nlbayes/assets/data/networks/gtex_chip/homo_sapiens/tissue_independent/three_tissue.rels.json"

wget -O "data/celllines_count_data.tar.gz" "https://poisson.math.umb.edu/~argenis/nlbayes/celllines_count_data.tar.gz"
tar -C ./data -zxvf data/celllines_count_data.tar.gz

wget -O "data/strandlab_count_data.tar.gz" "https://poisson.math.umb.edu/~argenis/nlbayes/strandlab_count_data.tar.gz"
tar -C ./data -zxvf data/strandlab_count_data.tar.gz
