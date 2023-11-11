suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(gridExtra)
  library(reshape2)
  library(mvtnorm)
  library(stringr)
  library(dplyr)
  library(future)
})

plan("sequential")
options(future.globals.maxSize = 30 * 1024 * 1024^2)

# (BE) Basal Epithelia
markers.BE <- c("KRT14", "DST", "KRT15", "KRT5", "RGCC")
# (LE) Luminal Epithelia
markers.LE <- c("MSMB", "KLK3", "ACPP", "PLA2G2A", "KLK2")
# (OE1) Other Epithelia 1
markers.OE1 <- c("SCGB3A1", "LCN2", "PIGR", "WFDC2", "FCGBP")
# (OE2) Other Epithelia 2
markers.OE2 <- c("KRT13", "APOBEC3A", "CSTB", "LYPD3", "SERPINB1")
# (NE) Neuroendocrine
markers.NE <- c("CHGA", "GRP", "CALCA", "SCG2", "TPH1")
# (Fib) Fibroblast
markers.Fib <- c("APOD", "FBLN1", "PTGDS", "CFD", "DCN")
# (SM) Smooth muscle
markers.SM <- c("TPM2", "ACTA2", "RGS5", "MT1A", "MYH11")
# (Endo) Endothelia
markers.Endo <- c("IFI27", "ACKR1", "SELE", "HMOX1", "CLDN5")
# (Leu) Leukocyte
markers.Leu <- c("RGS1", "C1QA", "C1QB", "TYROBP", "C1QC")


markers <- c(
  markers.BE,
  markers.LE,
  markers.OE1,
  markers.OE2,
  markers.Fib,
  markers.SM,
  markers.Endo
)


markers.list <- list(
  "BE" = markers.BE,
  "LE" = markers.LE,
  "OE1" = markers.OE1,
  "OE2" = markers.OE2,
  "Fib" = markers.Fib,
  "SM" = markers.SM,
  "Endo" = markers.Endo
)

if (!file.exists("data/ancillary/features_intersection.rds")) {

  if (!dir.exists("data/ancillary")) dir.create("data/ancillary")

  # all cell lines share the same features
  counts <- Read10X(data.dir = "data/celllines/reanalyze_count_N1/outs/filtered_feature_bc_matrix")
  features <- rownames(counts)

  saveRDS(features, file = "data/ancillary/features_intersection_celllines.rds")

  # D17 and D17_FACS share the same features
  counts <- Read10X(data.dir = "data/strandlab/D17_FACS/")
  features <- intersect(features, rownames(counts))

  saveRDS(features, file = "data/ancillary/features_intersection.rds")

}
features.intersection <- readRDS("data/ancillary/features_intersection.rds")


substrRight <- function(x, n) substr(x, nchar(x) - n + 1, nchar(x))


celllines2.idents <- read.csv("data/celllines/FromLoupe_Celllines2_retouched.csv")
rownames(celllines2.idents) <- celllines2.idents[, 1]
celllines2.idents[, 1] <- NULL


load.data <- function(data.dir, project, load.fraction = 1., features = NULL, normalize.data = FALSE) {
  counts <- Read10X(data.dir = data.dir)
  if (!is.null(features)) counts <- counts[features, ]
  if (load.fraction < 1.) counts <- counts[, sample(seq_len(ncol(counts)), ncol(counts) * load.fraction)]
  if (normalize.data) {
    ## CLR normalization is too strong. It makes N1 and SFT merge to a single cluster
    ## In the merge with references, iHPF and pHPF also merge into a single cluster
    counts <- NormalizeData(counts, normalization.method = "CLR")
  }
  CreateSeuratObject(counts = counts, project = project)
}


pre.process <- function(obj, s.level = 1.5) {
  obj[["orig.ident"]] <- obj@project.name
  obj$Celllines2 <- celllines2.idents[Cells(obj), "Celllines2"]

  obj <- RenameCells(obj, add.cell.id = obj@project.name)
  obj[["tag"]] <- factor(substrRight(colnames(obj), 1))
  obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")

  obj <- CellCycleScoring(obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
  obj$CC.Difference <- obj$S.Score - obj$G2M.Score

  m1 <- mean(obj$percent.mt)
  s1 <- sd(obj$percent.mt)
  obj <- subset(obj, subset = percent.mt < m1+s.level*s1)

  m2 <- mean(obj$nCount_RNA)
  s2 <- sd(obj$nCount_RNA)
  m3 <- mean(obj$nFeature_RNA)
  s3 <- sd(obj$nFeature_RNA)

  obj <- subset(obj, subset = nCount_RNA < m2+s.level*s2 & nCount_RNA > m2-s.level*s2)
  obj <- subset(obj, subset = nFeature_RNA < m3+s.level*s3 & nFeature_RNA > m3-s.level*s3)

  obj <- prune.group.outliers(obj, ident = "orig.ident", dim1 = "nCount_RNA", dim2 = "nFeature_RNA", use.log = TRUE)

  obj
}


normalize.sct <- function(obj) {
  obj <- SCTransform(obj, do.scale = FALSE, do.center = FALSE, verbose = FALSE,
                     method = "glmGamPoi", conserve.memory = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, verbose = FALSE)
  obj
}


select.g1 <- function(obj) {
  m1 <- mean(obj$G2M.Score[obj$G2M.Score < 0])
  s1 <- sd(obj$G2M.Score[obj$G2M.Score < 0])
  m2 <- mean(obj$S.Score[obj$S.Score < 0])
  s2 <- sd(obj$S.Score[obj$S.Score < 0])
  t1 <- m1 + 1. * s1
  t2 <- m2 + 1. * s2

  obj <- subset(obj, subset = G2M.Score < t1 & S.Score < t2)
  obj
}


combine.merge <- function(objs, vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt", "CC.Difference"),
                          variable.features.n = 3000, npcs = 30) {
  for (i in seq_along(objs)) {
    DefaultAssay(objs[[i]]) <- "RNA"
    if ("SCT" %in% Assays(objs[[i]])) objs[[i]][["SCT"]] <- NULL
  }
  combined <- merge(objs[[1]], y = objs[-1], project = "merged_objs")
  combined <- SCTransform(combined, vars.to.regress = vars.to.regress, do.scale = FALSE, do.center = FALSE,
                          verbose = FALSE, method = "glmGamPoi", conserve.memory = FALSE,
                          variable.features.n = variable.features.n, assay = "RNA")
  combined <- RunPCA(combined, npcs = npcs, verbose = FALSE)
  combined <- RunUMAP(combined, reduction = "pca", dims = 1:npcs, verbose = FALSE)
  combined
}


prune.group.outliers <- function(obj, ident, dim1 = "UMAP_1", dim2 = "UMAP_2", level = 0.0001, use.log = FALSE) {
  df <- FetchData(obj, vars = c(ident, dim1, dim2))
  if (use.log) df[c(dim1, dim2)] <- log(df[c(dim1, dim2)])

  for (cid in unique(df[[ident]])) {
    m <- unlist(colMeans(df[df[[ident]] == cid, c(dim1, dim2)]), use.names = FALSE)
    s <- cov(df[df[[ident]] == cid, c(dim1, dim2)])
    df[df[[ident]] == cid, "in_cluster.density"] <- dmvnorm(df[df[[ident]] == cid, c(dim1, dim2)], m, s)
  }
  obj$in_cluster.density <- df[["in_cluster.density"]]
  obj <- subset(obj, subset = in_cluster.density > level)
  obj
}


percent.markers <- function(obj) {
  for (type in names(markers.list)) {
    markers.type <- markers.list[[type]]
    obj <- PercentageFeatureSet(obj, features = markers.type, col.name = paste0("percent.", type), assay = "RNA")
  }
  obj
}


assign.cell.types <- function(obj, ident) {
  # This method for cell type assignment is outdated. Please consider the medthod
  # described by Ianevski et. al. 2022 (https://doi.org/10.1038/s41467-022-28803-w)

  obj <- percent.markers(obj)
  cols <- paste0("percent.", names(markers.list))
  df <- FetchData(obj, vars = c(ident, cols))
  for (cid in unique(df[[ident]])) {
    cell.type <- NA
    pval <- 1.
    df.this <- df[df[[ident]] == cid, cols]
    df.othr <- df[df[[ident]] != cid, cols]
    for (type in names(markers.list)) {
      col <- paste0("percent.", type)
      x <- df.this[[col]]
      y <- df.othr[[col]]
      p <- t.test(x, y, alternative = "greater")$p.value
      if (p < 0.05 && p < pval) {
        cell.type <- type
        pval <- p
      }
    }
    df[df[[ident]] == cid, "cell.type"] <- cell.type
  }
  obj$cell.type <- df[["cell.type"]]
  obj
}


cluster.cells <- function(obj) {
  obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
  obj <- FindClusters(obj, verbose = FALSE)
  obj
}


find.markers.wrt <- function(obj, reference.ident, logfc.threshold = log2(1.5)) {
  df <- data.frame()
  for (ident in levels(Idents(obj))) {
    if (ident == reference.ident)
      next
    markers <- FindMarkers(obj, logfc.threshold = logfc.threshold, ident.1 = ident, ident.2 = reference.ident)
    markers$cluster <- ident
    markers$gene <- rownames(markers)
    df <- rbind(df, markers)
  }
  df
}


# Downloaded from www.housekeeping.unicamp.br/Housekeeping_GenesHuman.csv
hk.df1 <- read.table("data/celllines/Housekeeping_GenesHuman.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
# Read from https://www.genomics-online.com/resources/16/5049/housekeeping-genes/
hk.df2 <- read.table("data/celllines/hk_genes.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

remove.hk <- function(x) {
  x <- x[!startsWith(rownames(x), "MT-"), , drop=FALSE]
  x <- x[!startsWith(rownames(x), "RPL"), , drop=FALSE]
  x <- x[!startsWith(rownames(x), "RPS"), , drop=FALSE]
  x <- x[!(rownames(x) %in% hk.df1$Symbol), , drop=FALSE]
  x <- x[!(rownames(x) %in% hk.df2$Symbol), , drop=FALSE]
  x
}

filter.hk <- function(x) {
  x <- x[!startsWith(x, "MT-")]
  x <- x[!startsWith(x, "RPL")]
  x <- x[!startsWith(x, "RPS")]
  x <- x[!(x %in% hk.df1$Symbol)]
  x <- x[!(x %in% hk.df2$Symbol)]
  x
}

assign.identities <- function(obj, idents.in.order = c("orig.ident")) {
  df <- FetchData(obj, vars = idents.in.order)
  df$Identity <- NA
  for (ident in idents.in.order) {
    mask <- is.na(df$Identity) & !is.na(df[[ident]])
    if (any(mask)) {
      df[mask, "Identity"] <- df[mask, ident]
    }
  }
  obj$Identity <- df$Identity
  obj
}
