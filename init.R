suppressPackageStartupMessages({
  library(Seurat)
  library(org.Hs.eg.db)
})

if (!dir.exists("./data")) dir.create("./data")


if (!file.exists("data/annot.rds")) {
  annot.src <- select(org.Hs.eg.db, keys(org.Hs.eg.db, "ENTREZID"), "SYMBOL")
  annot <- as.list(annot.src$SYMBOL)
  names(annot) <- annot.src$ENTREZID
  saveRDS(annot, "data/annot.rds", compress = FALSE)

  write.table(data.frame(Gene.ID = names(annot),
                         Gene.symbol = unlist(annot, use.names = FALSE)),
              "data/annot.tsv", row.names = FALSE, sep = "\t")

}


if (!file.exists("data/tfactids.rds")) {
  tf.act.ids <- select(org.Hs.eg.db, c("GO:0003700"), "ENTREZID", keytype = "GOALL")$ENTREZID
  dna.bind.ids <- select(org.Hs.eg.db, c("GO:0003677"), "ENTREZID", keytype = "GOALL")$ENTREZID
  treg.act.ids <- select(org.Hs.eg.db, c("GO:0140110"), "ENTREZID", keytype = "GOALL")$ENTREZID
  regoft.ids <- select(org.Hs.eg.db, c("GO:0006355"), "ENTREZID", keytype = "GOALL")$ENTREZID

  tf.act.ids <- union(
    tf.act.ids,
    union(
      intersect(dna.bind.ids, treg.act.ids),
      intersect(dna.bind.ids, regoft.ids)
    )
  )
  saveRDS(tf.act.ids, "data/tfactids.rds", compress = FALSE)

}


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


if (!file.exists("data/ancillary/d27.facs.rds")) {
  source("r_scripts/common_utils.R")

  if (!dir.exists("data/ancillary")) dir.create("data/ancillary")

  obj <- load.data(data.dir = "data/strandlab/D27_FACS/", project = "d27.facs",
                   load.fraction = 0.3, features = features.intersection)

  reference <- pre.process(obj, s.level = 1)
  FeatureScatter(reference, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  FeatureScatter(reference, feature1 = "S.Score", feature2 = "G2M.Score")
  FeatureScatter(reference, feature1 = "nCount_RNA", feature2 = "percent.mt")

  reference <- normalize.sct(reference)
  reference <- cluster.cells(reference)

  reference <- assign.cell.types(reference, "seurat_clusters")
  reference <- prune.group.outliers(reference, ident = "cell.type", level = 0.01)

  Idents(reference) <- "cell.type"
  reference <- subset(reference, downsample = 500)

  saveRDS(reference, file = "data/ancillary/d27.facs.rds")

}
