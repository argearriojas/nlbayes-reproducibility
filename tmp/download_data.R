if (!dir.exists("data")) dir.create("data")

url <- "https://umbibio.math.umb.edu/nlbayes/assets/data/networks/gtex_chip/homo_sapiens/tissue_independent/three_tissue.rels.json"
destfile <- "data/three_tissue.rels.json"
if (!file.exists(destfile)) download.file(url, destfile, mode = "wb")


url <- "https://poisson.math.umb.edu/~argenis/nlbayes/celllines_count_data.tar.gz"
destfile <- "data/celllines_count_data.tar.gz"
if (!file.exists(destfile)) download.file(url, destfile, mode = "wb")
if (!dir.exists("data/celllines")) untar(destfile, exdir = "data")


url <- "https://poisson.math.umb.edu/~argenis/nlbayes/strandlab_count_data.tar.gz"
destfile <- "data/strandlab_count_data.tar.gz"
if (!file.exists(destfile)) download.file(url, destfile, mode = "wb")
if (!dir.exists("data/strandlab")) untar(destfile, exdir = "data")


if (!dir.exists("data/ancillary")) dir.create("data/ancillary")

url <- "https://poisson.math.umb.edu/~argenis/nlbayes/ancillary/d27.facs.rds"
destfile <- "data/ancillary/d27.facs.rds"
if (!file.exists(destfile)) download.file(url, destfile, mode = "wb")
