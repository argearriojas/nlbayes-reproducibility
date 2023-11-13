suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
})


make.figure.4 <- function() {

  get.scatter <- function(experiment.name) {
    df <- read.csv(paste0("data/oe_", experiment.name, "_on_net_regulonbrca_nlbayes_and_viper.csv"))
    xmax <- max(- log10(df$viper.pvalue[!is.na(df$viper.pvalue)]))

    ggplot(df, aes(x = -log10(viper.pvalue),
                   y = posterior.p,
                   label = ifelse((posterior.p > 0.2 & viper.pvalue < 0.05) | rank(viper.pvalue) <= 3 | rank <= 3,
                                  symbol, ""))) +
      theme_bw(base_size = 14) +
      geom_vline(xintercept = -log10(c(0.05))) +
      geom_hline(yintercept = c(0.2)) +
      geom_point(color = ifelse((df$posterior.p >= 0.2 & df$viper.pvalue <= 0.05),
                                "red",
                                ifelse((df$posterior.p < 0.2 & df$viper.pvalue > 0.05),
                                       "gray",
                                       "black")),
        na.rm = TRUE
      ) +
      theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) +
      theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) +
      theme(strip.background = element_rect(colour = "black", fill = "white",
                                            size = 0.5, linetype = "solid")) +
      theme(strip.text = element_text(size = 14, face = "bold", angle = 0)) +

      coord_cartesian(xlim = c(NA, 0.5 + xmax), ylim = c(NA, 1.2)) +
      geom_text_repel(
        color = ifelse((df$posterior.p >= 0.2 & df$viper.pvalue <= 0.05),
                       "red",
                       ifelse((df$posterior.p < 0.2 & df$viper.pvalue > 0.05),
                              "gray",
                              "black")),
        force = 10,
        force_pull = 0.1,
        size = 2.5, fontface = "bold",
        box.padding = unit(0.6, "lines"),
        max.overlaps = 300,
        nudge_x = 0.25,
        nudge_y = 0.25,
        hjust = 0.25,
        segment.size = 0.1,
        na.rm = TRUE
      ) +
      labs(title = paste(paste0(toupper(experiment.name))), y = "TF posterior mean", x = "- log10( viper p-value )") +
      theme(
        plot.title = element_text(size = 14, face = "bold.italic", color = "red"),
        axis.title.x = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold")
      )
  }

  p1 <- get.scatter("e2f3")
  p2 <- get.scatter("myc")
  p3 <- get.scatter("ras")

  plot_title <- ggdraw() +
    draw_label("Results comparison between NLBayes and VIPER", fontface = "bold", x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))

  plot_row <- plot_grid(p1, p2, p3, labels = c("A)", "B)", "C)"), nrow = 1)
  f4 <- plot_grid(plot_title, plot_row, ncol = 1, rel_heights = c(0.1, 1))

  if (!dir.exists("figures")) {
    dir.create("figures")
  }

  ggsave2("figures/fig4.png", f4, height = 3, width = 10, bg = "white")

  f4
}


make.figure.6a <- function(obj.vect) {
  combined <- combine.merge(obj.vect,
                            vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt", "CC.Difference", "S.Score"))
  combined <- assign.identities(combined, idents.in.order = c("Celllines2", "orig.ident"))
  combined <- subset(combined, subset = Identity != "")
  combined$Phase <- sub("G2M", "G2/M", combined$Phase)
  f6a <- DimPlot(combined, reduction = "umap", group.by = "Phase", pt.size = 1, label = FALSE, label.size=8) +
    xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("") +
    theme(legend.position = c(0.75, 0.85), legend.text = element_text(size = 20),
          axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

  # skip save: ggsave2("figures/fig6a.png", plot = f6a, width = 12, height = 9, dpi = 150)

  f6a
}


make.figure.6b <- function(obj.vect, reference) {

  obj.list <- c(reference, obj.vect)
  obj.list <- sapply(obj.list, select.g1)
  obj.list <- sapply(obj.list, normalize.sct)

  features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
  obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)
  obj.list <- lapply(X = obj.list, FUN = RunPCA, features = features)

  anchors <- FindIntegrationAnchors(object.list = obj.list, reference = 1, k.anchor = 6, normalization.method = "SCT",
                                    anchor.features = features, dims = 1:30, reduction = "rpca")
  combined.integ <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)
  combined.integ <- RunPCA(combined.integ, verbose = FALSE)
  combined.integ <- RunUMAP(combined.integ, reduction = "pca", dims = 1:30)
  combined.integ <- assign.identities(combined.integ, idents.in.order = c("Celllines2", "cell.type", "orig.ident"))

  tmp <- subset(combined.integ, Identity != "Endo")
  identity.levels <- c("N1", "SFT1", "pHPF", "iHPF_A", "iHPF_B", "Fib", "SM", "BE", "LE", "OE1", "OE2")
  tmp[["Identity2"]] <- factor(tmp[[]]$Identity, levels = identity.levels)

  f6b <- DimPlot(tmp, reduction = "umap", group.by = "Identity2", pt.size = 1,
                 label = TRUE, label.size = 10, repel = TRUE) +
    xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("") +
    theme(axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30),
          legend.text = element_text(size = 24), legend.key.size = unit(26, "pt")) +
    guides(color = guide_legend(override.aes = list(size = 6)))

  # skip save: ggsave2("figures/fig6b.png", plot = f6b, width = 12, height = 9, dpi = 150)

  f6b
}


make.figure.6c <- function(celllines, top10.wrt.phpf) {
  top.absolute <- c()
  top10.list <- list()
  top10.deg.list <- list()
  for (ident in levels(Idents(celllines))) {
    count <- subset(celllines, idents = ident)[["SCT"]]@data
    total <- Matrix::rowSums(count) / ncol(count)
    x <- data.frame(total)
    x <- remove.hk(x)
    x <- head(x[order(x$total, decreasing = TRUE), , drop = FALSE], 20)
    top10.list[[ident]] <- rownames(head(x, 10))
    degs <- top10.wrt.phpf[top10.wrt.phpf$cluster == ident, ]$gene
    if (length(degs) > 0) top10.deg.list[[ident]] <- degs

    top.absolute <- c(top.absolute, rownames(x))
  }
  f6c <- DotPlot(celllines, features = unique(top.absolute), cols = c("white", "red3"), scale = FALSE) +
    theme(axis.text.x = element_text(size = 20, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 20, face = "bold", angle = 0, hjust = 1, vjust = 0.5))

  # skip save: ggsave2("figures/fig6c.png", plot = f6c, width = 12, height = 5, dpi = 150)

  f6c
}


make.figure.6d <- function(celllines, top10.wrt.phpf) {
  f6d <- DotPlot(celllines, features = unique(top10.wrt.phpf$gene), assay = "SCT", group.by = "Identity",
                 cols = "RdBu", scale = TRUE, col.min = -1., col.max = 1.) +
    theme(axis.text.x = element_text(size = 20, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 20, face = "bold", angle = 0, hjust = 1, vjust = 0.5))

  # skip save: ggsave2("figures/fig6d.png", plot = f6d, width = 12, height = 5, dpi = 150)

  f6d
}


make.figure.6e <- function(celllines) {
  source("r_scripts/nlbayes_utils.R")

  highly.expressed.genes <- c()
  for (cluster in c("iHPF_A", "iHPF_B", "N1", "SFT1")) {
    o <- subset(celllines, subset = Identity == cluster)
    x <- GetAssayData(o, slot = "data", assay = "SCT")
    y <- Matrix::rowSums(x) / ncol(x)
    y <- y[filter.hk(names(y))]
    highly.expressed.genes <- union(highly.expressed.genes, names(-sort(-y)[1:20]))
  }

  if (!file.exists("data/celllines_diff_expr_vs_pHPF.csv")) compute.sc.diff.expr()

  df <- read.table("data/celllines_diff_expr_vs_pHPF.csv", sep = ",", header = TRUE)
  df <- df[df$gene %in% highly.expressed.genes, ]
  df$gene <- factor(df$gene, levels = highly.expressed.genes)

  f6e <- ggplot(df, aes(x = gene, y = avg_log2FC)) +
    geom_bar(aes(fill = Identity), stat = "identity", position = position_dodge(0.6),
             color = ifelse(df$p_val_adj < 0.01, "black", "lightgray"), width = .1, alpha = .6) +
    geom_point(aes(size = -log10(p_val_adj), color=Identity), alpha = ifelse(df$p_val_adj < 0.01, 0.6, 0.),
               position = position_dodge(0.6)) + scale_size(range = c(1, 7)) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "red3") +
    geom_hline(yintercept = 0, color = "lightgray") +
    geom_vline(xintercept = 1:(length(highly.expressed.genes) - 1) + 0.5,
               linetype = "solid", color = "lightgray") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20),
          axis.text.y = element_text(face = "bold", size = 20),
          axis.title.y = element_text(face = "bold", size = 20),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18),
          legend.position = "bottom")

  # skip save: ggsave2("figures/fig6e.png", plot = f6e, width = 24, height = 5, dpi = 150)

  f6e
}


make.figure.6 <- function() {
  source("r_scripts/common_utils.R")

  data.dir <- "data/celllines/reanalyze_count_iHPF/outs/filtered_feature_bc_matrix"
  iHPF <- CreateSeuratObject(counts = Read10X(data.dir = data.dir), project = "iHPF")

  data.dir <- "data/celllines/reanalyze_count_pHPF/outs/filtered_feature_bc_matrix"
  pHPF <- CreateSeuratObject(counts = Read10X(data.dir = data.dir), project = "pHPF")

  data.dir <- "data/celllines/reanalyze_count_N1/outs/filtered_feature_bc_matrix"
  N1 <- CreateSeuratObject(counts = Read10X(data.dir = data.dir), project = "N1")

  data.dir <- "data/celllines/reanalyze_count_SFT1/outs/filtered_feature_bc_matrix"
  SFT1 <- CreateSeuratObject(counts = Read10X(data.dir = data.dir), project = "SFT1")


  iHPF <- pre.process(iHPF)
  pHPF <- pre.process(pHPF)
  N1 <- pre.process(subset(N1, downsample = 710))
  SFT1 <- pre.process(subset(SFT1, downsample = 850))

  obj.vect <- c(iHPF, pHPF, N1, SFT1)
  reference <- readRDS("data/ancillary/d27.facs.rds")

  f6a <- make.figure.6a(obj.vect)
  gc()
  f6b <- make.figure.6b(obj.vect, reference)
  gc()

  celllines <- obj.vect
  celllines <- sapply(celllines, select.g1)
  celllines <- combine.merge(celllines)
  celllines <- assign.identities(celllines, idents.in.order = c("Celllines2", "orig.ident"))
  Idents(celllines) <- "Identity"

  # keep track of the cells used in the analysis. This will be used for the RAW expression plots
  write.table(colnames(celllines), file = paste0("data/selected_barcodes.tsv"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  markers.wrt.phpf <- find.markers.wrt(celllines, reference.ident = "pHPF")
  top10.wrt.phpf <- markers.wrt.phpf %>% group_by(cluster) %>% top_n(10, avg_log2FC)

  f6c <- make.figure.6c(celllines, top10.wrt.phpf)
  f6d <- make.figure.6d(celllines, top10.wrt.phpf)
  f6e <- make.figure.6e(celllines)

  row1 <- plot_grid(f6a, f6b, ncol = 2, labels = c("A", "B"))
  row2 <- plot_grid(f6c, f6d, ncol = 2, labels = c("C", "D"))
  row3 <- f6e

  f6 <- plot_grid(row1, row2, row3, ncol = 1, rel_heights = c(4, 2, 2), labels = c("", "", "E"))

  ggsave2("figures/fig6.png", plot = f6, width = 24, height = 18, dpi = 150, bg = "white")

  f6
}


make.figure.7 <- function() {

  source("r_scripts/go_enrichment_utils.R")


  cellline.pHPF.markers <- read.csv("data/celllines_diff_expr_vs_pHPF.csv") %>%
    filter(gene %in% keys(org.Hs.eg.db, keytype = "SYMBOL"))

  logfc.level <- 1.
  top.n.terms <- 10
  enrichment.p.value.cutoff <- 0.01

  # the file name for the results of the enrichment analysis
  enrichment.results.rds.fn <- paste0("data/enrichment_results_logfc_", logfc.level, ".rds")

  ident.levels <- c("iHPF_A", "iHPF_B", "N1", "SFT1")
  ontology.levels <- c("BP", "MF", "CC")

  if (!file.exists(enrichment.results.rds.fn)) {
    enrichment.results <- data.frame()
    for (cellline.ident in ident.levels) {
      diff.expression <- get.diff.exp(cellline.pHPF.markers, cellline.ident, logfc.level)
      genes <- diff.expression$gene
      for (ontology in ontology.levels) {
        print(c(cellline.ident, ontology))
        result <- get.ego(genes, ontology)@result
        result$Cellline <- cellline.ident
        result$Ontology <- ontology
        enrichment.results <- rbind(enrichment.results, result)
      }
    }
    enrichment.results$Cellline <- factor(enrichment.results$Cellline, levels = ident.levels)
    enrichment.results$Ontology <- factor(enrichment.results$Ontology, levels = ontology.levels)
    levels(enrichment.results$Ontology) <- c("Biological Process", "Molecular Function", "Cellular Component")

    saveRDS(enrichment.results, enrichment.results.rds.fn)
  } else {
    enrichment.results <- readRDS(enrichment.results.rds.fn)
  }

  enrichment.results$Description <- str_replace(enrichment.results$Description, "regulation", "reg.")
  enrichment.results$Description <- str_replace(enrichment.results$Description, "signaling", "sig.")
  enrichment.results$Description <- str_replace(enrichment.results$Description, "serine", "ser")
  enrichment.results$Description <- str_replace(enrichment.results$Description, "threonine", "thr")
  enrichment.results$Description <- str_replace(enrichment.results$Description, "movement", "mov.")
  enrichment.results$Description <- str_replace(enrichment.results$Description, "rotational", "rot.")
  enrichment.results$Description <- enrichment.results$Description %>% substr(0, 60)

  filtered.enrichment.results <- (
    enrichment.results
    %>% filter(p.adjust < enrichment.p.value.cutoff)
    %>% group_by(Cellline, Ontology)
    %>% top_n(top.n.terms, -pvalue)
  )

  deg.counts.barplot <- get.deg.counts.barplot(cellline.pHPF.markers, logfc.level) +
    labs(title = paste0("Differential expression counts [ abs(log2FC) >", abs(logfc.level), " ]")) +
    theme(strip.text.y.left = element_text(angle = 0), plot.margin = unit(c(0.2, 0.2, 0., 0.4), "in"))

  enrichment.plot <- ggplot(filtered.enrichment.results,
                            aes(Cellline, Description, color = Ontology, size = -log(p.adjust))) +
    geom_point(alpha = 0.6) + facet_grid(Ontology ~ ., switch = "y", scales = "free", space = "free") +
    theme_bw() + ylab(NULL) + xlab(NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.3, size = 8, face = "bold"))

  if (logfc.level > 0)
    enrichment.plot <- enrichment.plot +
      labs(title = paste("GO terms enrichment on up-regulated genes")) +
      guides(color = "none")

  if (logfc.level < 0)
    enrichment.plot <- enrichment.plot +
      labs(title = paste("GO terms enrichment on down-regulated genes")) +
      guides(color = "none")

  enrichment.plot <- enrichment.plot +
    scale_y_discrete(position = "right") +
    force_panelsizes(cols=unit(0.8, "in")) +
    theme(legend.justification = c("left", "top"),
          legend.position = c(4, 0.2),
          plot.margin = unit(c(0.2, 0.5, 0.1, 0.4), "in"))

  f7 <- plot_grid(deg.counts.barplot, enrichment.plot,
                  labels = paste0(LETTERS[1:2], ")"),
                  ncol = 1, rel_heights = c(1, 4))

  ggsave2("figures/fig7.png", plot = f7, width = 5.5, height = 10, dpi = 150, bg = "white")
  f7
}

make.figure.8 <- function() {
  source("r_scripts/common_utils.R")
  suppressPackageStartupMessages({
    library(cowplot)
    library(stringr)
    library(ggpattern)
  })

  x <- AnnotationDbi::select(org.Hs.eg.db, keys(org.Hs.eg.db, "SYMBOL"), c("ENTREZID", "ENSEMBL"), "SYMBOL")
  entrez.from.ensembl <- x$ENTREZID
  names(entrez.from.ensembl) <- x$ENSEMBL

  selected.barcodes <- read.table("data/selected_barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)$V1
  selected.barcodes <- gsub("-.*", "", selected.barcodes)

  cell.identities <- celllines2.idents$Celllines2
  names(cell.identities) <- gsub("-.*", "", rownames(celllines2.idents))

  rna.expression.cutoffs <- list()
  df <- data.frame()
  show.n.tfs <- 0
  dfs <- list()
  for (cluster in c("iHPF_A", "iHPF_B", "N1", "SFT1")) {
    if (cluster %in% c("iHPF_A", "iHPF_B")) {
      dataset <- "iHPF"
    } else {
      dataset <- cluster
    }

    data.dir <- paste0("data/celllines/reanalyze_count_", dataset, "/outs/filtered_feature_bc_matrix")

    # read the raw count matrix
    counts <- Read10X(data.dir = data.dir, gene.column = 1)
    # strip lane information from the barcode
    colnames(counts) <- gsub("-.*", "", colnames(counts))
    # keep genes present in database
    counts <- counts[rownames(counts) %in% names(entrez.from.ensembl), ]
    # select cells from the cluster
    counts <- counts[, colnames(counts) %in% names(cell.identities[cell.identities == cluster])]
    # keep only the cells used in the analysis
    counts <- counts[, paste0(dataset, "_", colnames(counts)) %in% selected.barcodes]

    x <- rowMeans(counts)
    names(x) <- entrez.from.ensembl[names(x)]
    x <- sort(x, decreasing = TRUE)
    x <- x[!duplicated(names(x))]
    x <- x[order(names(x))]
    rna.expression <- x
    rna.expression.cutoff <- sort(x[x > 0])[as.integer(length(x[x > 0]) / 5)]

    dff <- read.csv(paste0("data/sc_", cluster, "_on_net_three_tissue_nlbayes.csv"), stringsAsFactors = FALSE)
    dff$id <- paste(dff$id)

    x <- rna.expression[dff$id]
    x[is.na(x)] <- 0
    dff$expression <- log10(x / min(x[x>0]))
    rna.expression.cutoff <- log10(rna.expression.cutoff / min(x[x>0]))

    dfs[[cluster]] <- dff
    rna.expression.cutoffs[[cluster]] <- rna.expression.cutoff

    dff$cluster <- cluster
    dff$rna.expression.cutoff <- rna.expression.cutoff
    df <- rbind(df, dff)
    show.n.tfs <- max(show.n.tfs, nrow(dff[dff$posterior.p > 0.2,]))

  }

  enrichment.cutoff <- -log10(0.05)
  df$psymbol <- str_pad(df$symbol, 12, "left")
  cluster_plots <- list()
  selected.tfs.union <- c()
  for (cluster in c("iHPF_A" ,"iHPF_B", "N1", "SFT1")) {
    rna.expression.cutoff <- rna.expression.cutoffs[[cluster]]
    d <- head(df[df$cluster == cluster,], show.n.tfs)
    if (cluster == "iHPF_B")
      d$psymbol <- str_pad(d$symbol, 14, "left")
    d$fsymbol <- factor(d$psymbol, levels = rev(d$psymbol))
    
    clear.x.axis <-  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    clear.y.axis <-  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    
    p0 <- ggplot(d) + 
      geom_blank(aes(y=fsymbol)) +
      theme_nothing() + 
      theme(
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_text(face="bold", size = 8, hjust = 1), axis.ticks.y = element_blank())
    p1 <- ggplot(d) + 
      geom_vline(xintercept = enrichment.cutoff) + 
      geom_col(aes(y=fsymbol, x= - log10(enrichment)), fill = ifelse(d$enrichment <= 0.05, "violetred", "lightgray")) + 
      theme_bw() +
      theme(
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_text(size = 9),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill="transparent", color=NA)) +
      clear.x.axis + 
      clear.y.axis
    p2 <- ggplot(d) + 
      geom_vline(xintercept = rna.expression.cutoff) + 
      geom_col(aes(y=fsymbol, x=expression), fill = ifelse(d$expression >= rna.expression.cutoff, "firebrick", "lightgray")) + 
      theme_bw() +
      theme(
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_text(size = 9),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill="transparent", color=NA)) +
      xlab("log10(RNA count)") + 
      clear.x.axis +
      clear.y.axis
    tfs.with.rna <- (d$expression >= rna.expression.cutoff)
    tfs.inferred <- (d$posterior.p >=0.2)
    selected.tfs <- tfs.with.rna & tfs.inferred
    selected.tfs.union <- union(selected.tfs.union, paste(d[selected.tfs, "symbol"]))
    fill <- ifelse(tfs.with.rna, ifelse(d$posterior.p >= 0.8, "forestgreen", ifelse(d$posterior.p >= 0.5, "steelblue", ifelse(d$posterior.p >= 0.2, "goldenrod", "lightgray"))), "lightgray")
    pattern <- ifelse(!tfs.inferred, "circle", "none")
    p3 <- ggplot(d) + 
      geom_vline(xintercept = c(0.2, 0.5, 0.8)) + 
      geom_col_pattern(aes(y=fsymbol, x=posterior.p), fill = fill, pattern = pattern, pattern_fill = "black", pattern_colour = "black", pattern_spacing = 0.05, pattern_density = 0.2) +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill="transparent", color=NA)) +
      xlab("Probability") + 
      clear.x.axis + 
      clear.y.axis + coord_cartesian(xlim=c(0., 1.))
    
    plot_title <- ggdraw() + 
      draw_label(paste("    "), fontface = "bold", x = 0, hjust = 0) +
      theme(plot.margin = margin(0, 0, 0, 7))
    
    plot_row <- plot_grid(p0, p3, p2, p1, nrow=1, align="h", rel_widths = c(1.1, 2, 2, 2))
    cluster_plots[[cluster]] <- plot_grid(plot_title, plot_row, ncol = 1, rel_heights = c(0.1, 1))
  }

  pg <- plot_grid(
    cluster_plots[["iHPF_A"]],
    cluster_plots[["iHPF_B"]],
    cluster_plots[["N1"]],
    cluster_plots[["SFT1"]],
    labels = c("     iHPF_A", "     iHPF_B", "        N1", "       SFT1")
  )

  d <- df[(df$symbol %in% selected.tfs.union)&(df$posterior.p >= 0.2),]
  f <- ifelse(
    d$expression < d$rna.expression.cutoff, "lightgray", ifelse(
    d$posterior.p >= 0.8, "forestgreen", ifelse(
    d$posterior.p >= 0.5, "steelblue", ifelse(
    d$posterior.p >= 0.2, "goldenrod","lightgray"))))
  p <- ggplot(d) +
    geom_point(aes(y=cluster, x=symbol, size=posterior.p), color=f) +
    scale_y_discrete(limits=rev) +
    labs(title = "Active TF inference") +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.y = element_text(face = "bold", size = 10, hjust = 1, vjust = 0.5), axis.title = element_blank(),
      axis.text.x = element_text(face = "bold", size = 10, angle = 90, hjust = 1, vjust = 0.5),
      title = element_text(face = "bold"),
      panel.background = element_rect(fill="transparent"),
      plot.background = element_rect(fill="transparent", color=NA)
    ) +
    coord_fixed(3/2)
  f8 <- plot_grid(pg, p, labels = c("", ""), ncol = 1, rel_heights = c(8, 3))

  ggsave2("figures/fig8.png", f8, height = 11, width = 9, dpi = 150, bg = "white")

  f8
}
