suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(reshape2)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(clusterProfiler)
  library(dplyr)
  library(ggh4x)
})


get.deg.counts.barplot <- function(markers.df, logfc.level) {

  df <- markers.df %>% filter(abs(avg_log2FC) > abs(logfc.level))
  df$`Diff. expression` <- ifelse(df$avg_log2FC > 0, "up regulated", "down regulated")
  nudge <- ((df %>% group_by(Identity, avg_log2FC > 0) %>% count)$n %>% max ) * 0.05
  ggplot(df, aes(y = `Diff. expression`, fill = `Diff. expression`)) +
    geom_bar() +
    geom_text(aes(label = after_stat(count)), stat = "count",nudge_x = nudge) +
    facet_grid(Identity~., switch = "y") +
    theme_bw() + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.5), "in")) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
}


get.diff.exp <- function(markers.df, cluster.name, logfc.level) {

  if (logfc.level > 0)
    diff.expression <- markers.df %>% filter(avg_log2FC > logfc.level) %>% filter(Identity == cluster.name)
  if (logfc.level < 0)
    diff.expression <- markers.df %>% filter(avg_log2FC < logfc.level) %>% filter(Identity == cluster.name)
  rownames(diff.expression) <- diff.expression$gene

  diff.expression$rank <- seq_len(nrow(diff.expression))
  diff.expression
}

get.down.diff.exp <- function(markers.df, cluster.name, logfc.level) {

  diff.expression <- markers.df %>% filter(avg_log2FC < -logfc.level) %>% filter(Identity == cluster.name)
  rownames(diff.expression) <- diff.expression$gene

  diff.expression$rank <- seq_len(nrow(diff.expression))
  diff.expression
}


get.ego <- function(genes, ontology, enrichment.p.value.cutoff = 0.01) {

  parse.eval <- function(r) base::eval(parse(text = r))

  universe <- readRDS("data/ancillary/features_intersection_celllines.rds")
  ego <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, universe = universe, keyType = "SYMBOL", ont = ontology,
                  pAdjustMethod = "BH", pvalueCutoff  = enrichment.p.value.cutoff, qvalueCutoff  = 0.05,
                  readable = FALSE)
  ego@result$`Fold enrichment` <- sapply(ego@result$GeneRatio, parse.eval) / sapply(ego@result$BgRatio, parse.eval)
  ego
}
