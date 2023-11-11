source("r_scripts/common_utils.R")

if (file.exists("data/ancillary/d27.facs.rds")) {
  stop("data/ancillary/d27.facs.rds already exists")
}

if (!dir.exists("data/ancillary")) {
  dir.create("data/ancillary")
}

obj <- load.data(data.dir = "data/strandlab/D27_FACS/", project = "d27.facs",
                 load.fraction = 0.3, features = features.intersection)

reference <- pre.process(obj, s.level = 1)
FeatureScatter(reference, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(reference, feature1 = "S.Score", feature2 = "G2M.Score")
FeatureScatter(reference, feature1 = "nCount_RNA", feature2 = "percent.mt")

reference <- normalize.sct(reference)
reference <- cluster.cells(reference)

reference <- assign.cell.types(reference, "seurat_clusters")
DimPlot(prune.group.outliers(reference, ident = "cell.type", level = 0.01), group.by = "cell.type")

reference <- prune.group.outliers(reference, ident = "cell.type", level = 0.01)

Idents(reference) <- "cell.type"
reference <- subset(reference, downsample = 500)
DimPlot(reference, group.by = "cell.type")

saveRDS(reference, file = "data/ancillary/d27.facs.rds")
