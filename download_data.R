suppressPackageStartupMessages({
  library(R.utils)
  library(RCurl)
})

src.urls <- c("https://github.com/argearriojas/nlbayes-reproducibility/releases/download/v1.0.0",
              "https://poisson.math.umb.edu/~argenis/nlbayes",
              "https://zenodo.org/records/10116664/files")

download <- function(destfile) {
  if (file.exists(destfile)) return(paste(destfile, "already exists. Skipping."))
  destdir <- dirname(destfile)
  if (!dir.exists(destdir)) dir.create(destdir, recursive = TRUE)

  filename <- basename(destfile)

  for (src.url in src.urls) {
    url <- file.path(src.url, filename)
    if (startsWith(url, "https://zenodo.org")) url <- paste0(url, "?download=1")

    if (!url.exists(url)) next
    download.file(url, destfile, mode = "wb")
    break
  }

  if (endsWith(destfile, ".tar.gz")) return(untar(destfile, exdir = "data"))
  if (endsWith(destfile, ".gz")) return(gunzip(destfile, remove = FALSE))
}


download("data/three_tissue.rels.json.gz")
download("data/celllines_count_data.tar.gz")
download("data/strandlab_count_data.tar.gz")
download("data/ancillary/d27.facs.rds")
