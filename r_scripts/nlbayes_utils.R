suppressPackageStartupMessages({
  library(dplyr)
  library(nlbayes)
  library(viper)
  library(aracne.networks)
  library(org.Hs.eg.db)
})

# > keytypes(org.Hs.eg.db)
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
# [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"
# [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"
# [16] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"
# [21] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"
# [26] "UNIPROT"


if (!file.exists("data/annot.rds")) {
  annot.src <- select(org.Hs.eg.db, keys(org.Hs.eg.db, "ENTREZID"), "SYMBOL")
  annot <- as.list(annot.src$SYMBOL)
  names(annot) <- annot.src$ENTREZID
  saveRDS(annot, "annot.rds", compress = FALSE)
}
annot <- readRDS("data/annot.rds")


if (!file.exists("data/annot.tsv"))
  write.table(data.frame(Gene.ID = names(annot),
                         Gene.symbol = unlist(annot, use.names = FALSE)),
              "data/annot.tsv", row.names = FALSE, sep = "\t")



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
tf.act.ids <- readRDS("data/tfactids.rds")


clean.tT <- function(tT, use.id = "Gene.ID") {

  # remove promiscuous probes
  tT <- tT[!grepl("///", tT[[use.id]], fixed = TRUE), ]

  # remove probes with no annotation
  tT <- tT[tT[[use.id]] != "", ]

  # remove duplicated probes. Keep those with greater variation
  tT <- tT[!duplicated(tT[[use.id]]), ]

  tT
}


rename.gset.features <- function(gset, tT, use.id = "Gene.ID") {
  # get the expression set subset with these selected probes
  gset <- gset[fData(gset)$ID %in% tT$ID, ]

  # rename features
  featureNames(gset) <- fData(gset)[[use.id]]

  gset
}


ents.rels.2.net <- function(ents, rels, remove.non.tfs = TRUE) {
  row.names(ents) <- ents$uid

  net <- list()
  for (srcuid in ents[ents$type == "Protein", "uid"]){
    srcid <- ents[srcuid, "id"]
    if (remove.non.tfs && !srcid %in% tf.act.ids) next
    if (nrow(rels[rels$srcuid == srcuid, ]) == 0) next
    trguids <- rels[rels$srcuid == srcuid, "trguid"]
    mors <- rels[rels$srcuid == srcuid, "type"]
    ids <- ents[trguids, "id"]
    names(mors) <- as.character(ids)
    net[[as.character(srcid)]] <- mors
  }
  net
}


net.2.ents.rels.evid <- function(net, evid) {
  srcs <- names(net)

  net2 <- net
  names(net2) <- NULL
  trgs <- unique(names(unlist(net2)))

  nsrc <- length(srcs)
  ntrg <- length(trgs)

  ents <- data.frame(
    uid = as.character(1:(nsrc + ntrg)),
    name = paste(annot[c(srcs, trgs)]),
    id = c(srcs, trgs),
    type = c(rep("Protein", nsrc), rep("mRNA", ntrg))
  )
  ents <- ents[ents$name != "NULL", ]

  src.ents <- ents[ents$type == "Protein", ]
  src.map <- src.ents$uid
  names(src.map) <- src.ents$id

  trg.ents <- ents[ents$type == "mRNA", ]
  trg.map <- trg.ents$uid
  names(trg.map) <- trg.ents$id

  rels <- data.frame()
  for (src in names(net)) {
    mor <- net[[src]]
    srcuid <- paste(src.map[[src]])
    if (srcuid == "NA") next
    df <- data.frame(uid = rep("", length(mor)),
                     srcuid = rep(srcuid, length(mor)),
                     trguid = paste(trg.map[names(mor)]),
                     type = mor)
    df <- df[df$trguid != "NA", ]

    rels <- rbind(rels, df)
  }
  rownames(rels) <- NULL
  rels$uid <- as.character(1:nrow(rels))

  val <- sign(evid$logFC)
  names(val) <- evid$Gene.ID
  trg.ents$val <- val[trg.ents$id]
  trg.ents[is.na(trg.ents$val), "val"] <- 0
  evid <- trg.ents[, c("uid", "val")]
  rownames(evid) <- NULL

  list(
    ents = ents,
    rels = rels,
    evid = evid
  )
}


read.network <- function(net.name) {
  ents <- read.table(paste0(net.name, "_ents.csv"), sep = ",", header = TRUE, as.is = TRUE)
  rels <- read.table(paste0(net.name, "_rels.csv"), sep = ",", header = TRUE, as.is = TRUE)

  ents.rels.2.net(ents, rels)
}


filter.regulon.1 <- function(regulon) {
  regulon <- regulon[intersect(tf.act.ids, names(regulon))]

  regulon
}


filter.regulon.2 <- function(regulon) {
  regulon <- filter.regulon.1(regulon)

  for (reg.id in names(regulon)) {
    # keep only edges with high confidence
    mask <- regulon[[reg.id]][["likelihood"]] > 0.5
    regulon[[reg.id]][["tfmode"]] <- regulon[[reg.id]][["tfmode"]][mask]
    regulon[[reg.id]][["likelihood"]] <- regulon[[reg.id]][["likelihood"]][mask]

    # discretize mode of regulation. Remove sign for edges with low correlation
    mask <- abs(regulon[[reg.id]][["tfmode"]]) < 0.5
    regulon[[reg.id]][["tfmode"]][mask] <- 0
    regulon[[reg.id]][["tfmode"]] <- sign(regulon[[reg.id]][["tfmode"]])
  }

  regulon
}


regulon2network <- function(regulon, keep.zeros = FALSE) {

  regulon <- filter.regulon.2(regulon)
  net <- list()
  for (src in names(regulon)) {
    x <- sign(regulon[[src]]$tfmode)
    if (!keep.zeros) x <- x[x != 0]
    net[[src]] <- x
  }

  net
}


regulon2long <- function(regulon, keep.zeros = FALSE) {

  dfs <- list()
  for (src in names(regulon)) {
    weight <- regulon[[src]]$tfmode
    target <- names(weight)
    n <- length(weight)
    dfs[[src]] <- data.frame(uid = paste(rep(src, n), target, sep = "-"),
                             source = rep(src, n),
                             target = target, weight = weight)
  }

  do.call(rbind, dfs)
}


net2long <- function(net, keep.zeros = FALSE) {

  long <- data.frame()
  for (src in names(net)) {
    weight <- net[[src]]
    target <- names(weight)
    n <- length(weight)
    df <- data.frame(uid = paste(rep(src, n), target, sep = "-"),
                     source = rep(src, n),
                     target = target, weight = weight)
    long <- rbind(long, df)
  }

  long
}


network2regulon <- function(net) {
  regulon <- list()
  for (src in names(net)) {
    regulon[[src]] <- list(tfmode = net[[src]],
                           likelihood = rep(1., length(net[[src]])))
  }
  regulon <- filter.regulon.1(regulon)
  class(regulon) <- "regulon"
  regulon
}


run.viper <- function(gset, regulon) {
  regulon <- filter.regulon.1(regulon)

  sig <- rowTtest(gset, "group", "pt", "wt")$statistic
  dnull <- ttestNull(gset, "group", "pt", "wt", per = 1000)
  mra <- msviper(sig, regulon, dnull, pleiotropy = FALSE)
  mra <- msviperAnnot(mra, annot)
  df <- data.frame(mra$es)
  df
}


run.nlbayes <- function(tT, net, uniform.t = FALSE, zn = 0., zy = 0., verbosity = 0, gr.level = 1.1,
                        logfc.trh = 1, max.degs = Inf, N = 20000, n.graphs = 5) {

  diff.exp <- head(tT[abs(tT$logFC) > logfc.trh & tT$adj.P.Val < 0.01, ], max.degs)
  evd <- sign(diff.exp$logFC)
  names(evd) <- diff.exp$Gene.ID

  inference.model <- ORNOR.inference(net, evidence = evd, uniform.t = uniform.t, N = N,
                                     s.leniency = 0.1, n.graphs = n.graphs, gr.level = gr.level,
                                     zn = zn, zy = zy, verbosity = verbosity, do.sample = TRUE,
                                     burnin = TRUE)
  inference.model <- postprocess.result(inference.model, annot)

  inference.model
}


compute.enrichment <- function(net, evd) {

  n.genes <- length(unique(unlist(sapply(net, names))))
  all.ydeg <- length(evd[evd != 0])
  all.ndeg <- n.genes - all.ydeg

  enrichment <- list()
  for (src in names(net)) {
    trg.ydeg <- sum(names(net[[src]]) %in% names(evd[evd != 0]))
    trg.ndeg <- length(net[[src]]) - trg.ydeg
    contingency <- matrix(c(trg.ydeg, trg.ndeg, all.ydeg - trg.ydeg, all.ndeg - trg.ndeg), ncol = 2)
    enrichment[[src]] <- fisher.test(contingency, alternative = "greater")$p.value
  }

  enrichment
}


explain.ornor <- function(net.t, posterior.p, evd, lvl = 0.5) {
  x <- posterior.p > lvl
  de.explained <- 0
  for (trg in names(evd)) {
    diff.exp <- evd[[trg]]
    if (diff.exp == 0) next
    mask <- x[names(net.t[[trg]])]
    y <- net.t[[trg]][mask]
    y <- y[y != 0]
    if (length(y) == 0) next
    ornor.pred <- min(y)
    de.explained <- de.explained + as.integer(ornor.pred == diff.exp)
  }

  de.explained
}


transpose.net <- function(net) {
  net.t <- list()
  for (src in names(net)) {
    trg.lst <- net[[src]]
    for (trg in names(trg.lst)) {
      val <- trg.lst[[trg]]
      src.lst <- net.t[[trg]]
      if (is.null(src.lst)) src.lst <- c()

      val <- c(val)
      names(val) <- c(src)
      net.t[[trg]] <- c(src.lst, val)
    }
  }

  net.t
}


get.chip.corrected.network <- function(net.name, remove.non.tfs = TRUE, keep.zeros = TRUE) {
  ents <- read.table("data/ChIPfilter.ents", header = TRUE)
  rels <- read.table(paste0("data/", net.name, ".rels"), header = TRUE)

  type.dict <- list("increase" = 1, "conflict" = 0, "decrease" = -1)
  rels$type <- unlist(type.dict[rels$type], use.names = FALSE)
  if (!keep.zeros) rels <- rels[rels$type != 0, ]

  net <- ents.rels.2.net(ents, rels, remove.non.tfs)

  net
}


get.oe.diff.expr <- function(experiment.name) {
  if (experiment.name %in% c("e2f3", "myc", "ras", "src", "bcat")) {
    suppressPackageStartupMessages({
      source(paste0("r_scripts/oe_", experiment.name, "_script.R"))
    })

    tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
    tT <- clean.tT(tT, use.id = "Gene.ID")
    gset <- rename.gset.features(gset, tT, use.id = "Gene.ID")

    filepath <- paste0("data/", experiment.name, "_experiment_tT.tsv")
    if (!file.exists(filepath))
      write.table(tT[, c("Gene.symbol", "Gene.ID", "t", "logFC", "P.Value", "adj.P.Val")],
                  filepath, sep = "\t", row.names = FALSE)

  } else if (experiment.name %in% c("tgfb", "cxcl12")) {
    tT <- read.table(paste0("cell_line_vehicle_", experiment.name, "_evidence_edgeR.txt"), header = TRUE)
    colnames(tT) <- c("Gene.ID", "adj.P.Val", "logFC")
    tT <- tT[order(tT$adj.P.Val), ]
    tT <- clean.tT(tT)

  } else if (experiment.name %in% c("iHPF_A", "iHPF_B", "N1", "SFT1")) {
    tT <- read.csv(paste0("data/de_markers_", experiment.name, "_vs_pHPF_tT.csv"), stringsAsFactors = FALSE)
    tT$Gene.ID <- paste(tT$Gene.ID)
    tT <- tT[order(tT$adj.P.Val), ]
    tT <- clean.tT(tT)
    gset <- NA

  } else {
    stop(paste("no definitions found for experiment", experiment.name))
  }


  list(tT = tT, gset = gset)
}


get.network <- function(network.name) {
  if (exists(network.name)) {
    regulon <- get(network.name)
    net <- regulon2network(regulon, keep.zeros = FALSE)

  } else if (file.exists(paste0("data/", network.name, ".rels.json"))) {
    net <- jsonlite::fromJSON(paste0("data/", network.name, ".rels.json"))
    net <- sapply(net, unlist)
    net <- sapply(net, function(x) x[x != 0])
    net <- net[names(net) %in% tf.act.ids]

  } else {
    stop(paste("network", network.name, "not available"))

  }

  net
}


get.regulon <- function(network.name) {
  if (exists(network.name)) {
    regulon <- get(network.name)

  } else if (file.exists(paste0("data/", network.name, ".rels.json"))) {
    net <- jsonlite::fromJSON(paste0("data/", network.name, ".rels.json"))
    net <- sapply(net, unlist)
    regulon <- network2regulon(net)

  } else {
    stop(paste("network", network.name, "not available"))

  }

  regulon
}


compute.inference.comparison <- function() {

  network.name <- "regulonbrca"
  net <- get.network(network.name)
  regulon <- get.regulon(network.name)

  for (experiment.name in c("myc", "e2f3", "ras")) {

    output.filename <- paste0("data/oe_", experiment.name, "_on_net_", network.name,
                              "_nlbayes_and_viper", ".csv")
    if (file.exists(output.filename)) next

    top.table <- get.oe.diff.expr(experiment.name)
    tT <- top.table[["tT"]]
    gset <- top.table[["gset"]]

    if (experiment.name == "ras") {
      logfc.trh <- 2.5
    } else {
      logfc.trh <- 1
    }

    inference.result <- run.nlbayes(tT, net, uniform.t = TRUE, verbosity = 2, logfc.trh = logfc.trh)
    nlbayes.df <- inference.result$result.info$tf.inference

    viper.df <- run.viper(gset, regulon)
    nlbayes.df$viper.pvalue <- viper.df[nlbayes.df$symbol, "p.value"]

    write.csv(nlbayes.df, output.filename, row.names = FALSE)
  }
}


compute.sc.diff.expr <- function() {
  source("r_scripts/common_utils.R")

  x <- AnnotationDbi::select(org.Hs.eg.db, keys(org.Hs.eg.db, "SYMBOL"), c("ENTREZID", "ENSEMBL"), "SYMBOL")
  entrez.from.symbol <- x$ENTREZID
  names(entrez.from.symbol) <- x$SYMBOL
  entrez.from.ensembl <- x$ENTREZID
  names(entrez.from.ensembl) <- x$ENSEMBL

  feature.table <- read.table("data/celllines/reanalyze_count_iHPF/outs/filtered_feature_bc_matrix/features.tsv.gz")
  ensembl.from.symbol <- feature.table$V1
  names(ensembl.from.symbol) <- feature.table$V2

  counts <- Read10X(data.dir = "data/celllines/reanalyze_count_iHPF/outs/filtered_feature_bc_matrix", gene.column = 2)
  iHPF <- CreateSeuratObject(counts = counts, project = "iHPF")

  counts <- Read10X(data.dir = "data/celllines/reanalyze_count_pHPF/outs/filtered_feature_bc_matrix", gene.column = 2)
  pHPF <- CreateSeuratObject(counts = counts, project = "pHPF")

  counts <- Read10X(data.dir = "data/celllines/reanalyze_count_N1/outs/filtered_feature_bc_matrix", gene.column = 2)
  N1 <- CreateSeuratObject(counts = counts, project = "N1")

  counts <- Read10X(data.dir = "data/celllines/reanalyze_count_SFT1/outs/filtered_feature_bc_matrix", gene.column = 2)
  SFT1 <- CreateSeuratObject(counts = counts, project = "SFT1")

  iHPF <- pre.process(iHPF)
  pHPF <- pre.process(pHPF)
  N1 <- pre.process(subset(N1, downsample = 710))
  SFT1 <- pre.process(subset(SFT1, downsample = 850))

  iHPF <- assign.identities(iHPF, idents.in.order = c("Celllines2", "orig.ident"))
  pHPF <- assign.identities(pHPF, idents.in.order = c("Celllines2", "orig.ident"))
  N1 <- assign.identities(N1, idents.in.order = c("Celllines2", "orig.ident"))
  SFT1 <- assign.identities(SFT1, idents.in.order = c("Celllines2", "orig.ident"))

  iHPF <- subset(iHPF, subset = Identity != "")
  pHPF <- subset(pHPF, subset = Identity != "")
  N1 <- subset(N1, subset = Identity != "")
  SFT1 <- subset(SFT1, subset = Identity != "")

  iHPF.pHPF <- combine.merge(c(iHPF, pHPF))
  Idents(iHPF.pHPF) <- "Identity"

  N1.pHPF <- combine.merge(c(N1, pHPF))
  Idents(N1.pHPF) <- "Identity"

  SFT1.pHPF <- combine.merge(c(SFT1, pHPF))
  Idents(SFT1.pHPF) <- "Identity"

  dataset.map <- list(iHPF_A = iHPF.pHPF, iHPF_B = iHPF.pHPF, N1 = N1.pHPF, SFT1 = SFT1.pHPF)

  celllines.diff.expr <- data.frame()
  for (ident.name in c("iHPF_A", "iHPF_B", "N1", "SFT1")) {

    dataset <- dataset.map[[ident.name]]

    tT <- FindMarkers(dataset, ident.1 = ident.name, ident.2 = "pHPF")
    tT$gene <- rownames(tT)

    tT$ensembl <- (ensembl.from.symbol[tT$gene])
    tT$entrez.e <- (entrez.from.ensembl[tT$ensembl])
    tT$entrez <- (entrez.from.symbol[tT$gene])
    tT <- tT[!(is.na(tT$entrez) & is.na(tT$entrez.e)), ]
    tmp <- tT[is.na(tT$entrez), ]
    tT[is.na(tT$entrez), "entrez"] <- tmp$entrez.e
    tT <- tT[!duplicated(tT$entrez), c("gene", "entrez", "p_val_adj", "avg_log2FC")]
    tT$Identity <- ident.name
    celllines.diff.expr <- rbind(celllines.diff.expr, tT)

    tT <- tT[, c("entrez", "p_val_adj", "avg_log2FC")]
    colnames(tT) <- c("Gene.ID", "adj.P.Val", "logFC")

    file.name <- paste0("data/de_markers_", ident.name, "_vs_pHPF_tT.csv")
    write.csv(tT, file.name, row.names = FALSE)
  }
  write.table(celllines.diff.expr, file = "data/celllines_diff_expr_vs_pHPF.csv",
              sep = ",", quote = FALSE, row.names = FALSE)
}


compute.sc.inference <- function() {

  network.name <- "three_tissue"
  net <- get.network(network.name)

  experiment.list <- c("iHPF_A", "iHPF_B", "N1", "SFT1")
  for (experiment.name in experiment.list) {

    output.filename <- paste0("data/sc_", experiment.name, "_on_net_", network.name,
                              "_nlbayes", ".csv")
    if (file.exists(output.filename)) next

    top.table <- get.oe.diff.expr(experiment.name)
    tT <- top.table[["tT"]]

    logfc.trh <- 1

    inference.result <- run.nlbayes(tT, net, uniform.t = TRUE, zn = 0., zy = 0., verbosity = 2, gr.level = 1.15,
                                    logfc.trh = logfc.trh, N = 1200)
    nlbayes.df <- inference.result$result.info$tf.inference

    write.csv(nlbayes.df, output.filename, row.names = FALSE)
  }
}
