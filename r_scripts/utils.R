suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(nlbayes)
  library(viper)
  library(aracne.networks)
})

# keytypes(org.Hs.eg.db)
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
# [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
# [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"         
# [16] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"        
# [21] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
# [26] "UNIPROT"


if(!exists('annot') && !file.exists('data/annot.rds')) {
  annot.src <- select(org.Hs.eg.db, keys(org.Hs.eg.db, 'ENTREZID'), 'SYMBOL')
  annot <- as.list(annot.src$SYMBOL)
  names(annot) <- annot.src$ENTREZID
  saveRDS(annot, "annot.rds", compress = FALSE)
} else if (file.exists('data/annot.rds')) {
  annot <- readRDS('data/annot.rds')
}

if(!file.exists('data/annot.tsv'))
  write.table(data.frame(Gene.ID=names(annot), Gene.symbol=unlist(annot, use.names=F)), 'data/annot.tsv', row.names = F, sep = '\t')



if(!exists('tf.act.ids') && !file.exists('data/tfactids.rds')) {
  tf.act.ids <- select(org.Hs.eg.db, c("GO:0003700"), "ENTREZID", keytype = 'GOALL')$ENTREZID
  dna.bind.ids <- select(org.Hs.eg.db, c("GO:0003677"), "ENTREZID", keytype = 'GOALL')$ENTREZID
  treg.act.ids <- select(org.Hs.eg.db, c("GO:0140110"), "ENTREZID", keytype = 'GOALL')$ENTREZID
  regoft.ids <- select(org.Hs.eg.db, c("GO:0006355"), "ENTREZID", keytype = 'GOALL')$ENTREZID
  
  tf.act.ids <- union(
    tf.act.ids,
    union(
      intersect(dna.bind.ids, treg.act.ids),
      intersect(dna.bind.ids, regoft.ids)
    )
  )
  saveRDS(tf.act.ids,'data/tfactids.rds', compress = FALSE)
} else if(file.exists('data/tfactids.rds')) {
  tf.act.ids <- readRDS('data/tfactids.rds')
}


clean.tT <- function(tT, use.id = "Gene.ID") {
  
  # remove promiscuous probes
  tT <- tT[!grepl('///', tT[[use.id]], fixed = TRUE),]
  
  # remove probes with no annotation
  tT <- tT[tT[[use.id]] != "",]
  
  # remove duplicated probes. Keep those with greater variation
  tT <- tT[!duplicated(tT[[use.id]]),]
  
  tT
}


rename.gset.features <- function(gset, tT, use.id = "Gene.ID") {
  # get the expression set subset with these selected probes
  gset <- gset[fData(gset)$ID %in% tT$ID,]

  # rename features
  featureNames(gset) <- fData(gset)[[use.id]]
  
  gset
}


ents.rels.2.net <- function(ents, rels, remove.non.tfs=TRUE) {
  row.names(ents) <- ents$uid

  net <- list()
  for(srcuid in ents[ents$type == 'Protein', 'uid']){
    srcid <- ents[srcuid, 'id']
    if (remove.non.tfs & !srcid %in% tf.act.ids) next
    if (nrow(rels[rels$srcuid == srcuid,]) == 0) next
    trguids <- rels[rels$srcuid == srcuid, 'trguid']
    mors <- rels[rels$srcuid == srcuid, 'type']
    ids <- ents[trguids, 'id']
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
    uid=as.character(1:(nsrc+ntrg)),
    name=paste(annot[c(srcs, trgs)]),
    id=c(srcs, trgs),
    type=c(rep('Protein', nsrc), rep('mRNA', ntrg))
  )
  ents <- ents[ents$name != 'NULL',]

  src.ents <- ents[ents$type == 'Protein',]
  src.map <- src.ents$uid
  names(src.map) <- src.ents$id
  
  trg.ents <- ents[ents$type == 'mRNA',]
  trg.map <- trg.ents$uid
  names(trg.map) <- trg.ents$id
  
  rels <- data.frame()
  for(src in names(net)) {
    mor <- net[[src]]
    srcuid <- paste(src.map[[src]])
    if (srcuid == "NA") next
    df <- data.frame(
      uid=rep("", length(mor)),
      srcuid=rep(srcuid, length(mor)),
      trguid=paste(trg.map[names(mor)]),
      type=mor)
    df <- df[df$trguid != "NA",]
    
    rels <- rbind(rels, df)
  }
  rownames(rels) <- NULL
  rels$uid <- as.character(1:nrow(rels))

  val <- sign(evid$logFC)
  names(val) <- evid$Gene.ID
  trg.ents$val <- val[trg.ents$id]
  trg.ents[is.na(trg.ents$val), 'val'] <- 0
  evid <- trg.ents[, c('uid', 'val')]
  rownames(evid) <- NULL
    
  list(
    ents=ents,
    rels=rels,
    evid=evid
  )
}


read.network <- function(net.name) {
  ents <- read.table(paste0(net.name, '_ents.csv'), sep=',', header = TRUE, as.is = TRUE)
  rels <- read.table(paste0(net.name, '_rels.csv'), sep=',', header = TRUE, as.is = TRUE)
  
  ents.rels.2.net(ents, rels)
}


filter.regulon.1 <- function(regulon) {
  regulon <- regulon[intersect(tf.act.ids, names(regulon))]

  regulon
}


filter.regulon.2 <- function(regulon) {
  regulon <- filter.regulon.1(regulon)

  for(reg.id in names(regulon)) {
    # keep only edges with high confidence
    mask <- regulon[[reg.id]][['likelihood']] > 0.5
    regulon[[reg.id]][['tfmode']] <- regulon[[reg.id]][['tfmode']][mask]
    regulon[[reg.id]][['likelihood']] <- regulon[[reg.id]][['likelihood']][mask]

    # discretize mode of regulation. Remove sign for edges with low correlation
    mask <- abs(regulon[[reg.id]][['tfmode']]) < 0.5
    regulon[[reg.id]][['tfmode']][mask] <- 0
    regulon[[reg.id]][['tfmode']] <- sign(regulon[[reg.id]][['tfmode']])
  }

  regulon
}


regulon2network <- function(regulon, keep.zeros=FALSE) {

  regulon <- filter.regulon.2(regulon)
  net <- list()
  for(src in names(regulon)) {
    x <- sign(regulon[[src]]$tfmode)
    if(!keep.zeros) x <- x[x!=0]
    net[[src]] <- x
  }

  net
}


regulon2long <- function(regulon, keep.zeros=FALSE) {
  
  long <- data.frame()
  for(src in names(regulon)) {
    weight <- regulon[[src]]$tfmode
    target <- names(weight)
    n <- length(weight)
    df <- data.frame(uid=paste(rep(src, n), target, sep='-'), source=rep(src, n), target=target, weight=weight)
    long <- rbind(long, df)
  }

  long
}


net2long <- function(net, keep.zeros=FALSE) {
  
  long <- data.frame()
  for(src in names(net)) {
    weight <- net[[src]]
    target <- names(weight)
    n <- length(weight)
    df <- data.frame(uid=paste(rep(src, n), target, sep='-'), source=rep(src, n), target=target, weight=weight)
    long <- rbind(long, df)
  }
  
  long
}


network2regulon <- function(net) {
  regulon <- list()
  for(src in names(net)) {
    regulon[[src]] <- list(
      tfmode=net[[src]],
      likelihood=rep(1., length(net[[src]]))
    )
  }
  regulon <- filter.regulon.1(regulon)
  class(regulon) <- 'regulon'
  regulon
}


run.viper <- function(gset, regulon) {
  regulon <- filter.regulon.1(regulon)

  sig <- rowTtest(gset, "group", "pt", "wt")$statistic
  dnull <- ttestNull(gset, "group", "pt", "wt", per=1000)
  mra <- msviper(sig, regulon, dnull, pleiotropy=FALSE)
  mra <- msviperAnnot(mra, annot)
  df <- data.frame(mra$es)
  df
}


run.nlbayes <- function(tT, net, uniform.t = FALSE, zn=0., zy=0., verbosity=0, gr.level = 1.10, logfc.trh = 1, max.degs=Inf, N = 20000) {

  diff.exp <- head(tT[abs(tT$logFC) > logfc.trh & tT$adj.P.Val < 0.01, ], max.degs)
  evd <- sign(diff.exp$logFC)
  names(evd) <- diff.exp$Gene.ID

  inference.model <- ORNOR.inference(net, evidence = evd, uniform.t = uniform.t, N = N, s.leniency = 0.1, n.graphs = 5, gr.level = gr.level, zn=zn, zy=zy, verbosity=verbosity, do.sample = TRUE, burnin = TRUE)
  inference.model <- postprocess.result(inference.model, annot)
  
  inference.model
}


compute.enrichment <- function(net, evd) {

  n.genes <- length(unique(unlist(sapply(net, names))))
  all.ydeg <- length(evd[evd!=0])
  all.ndeg <- n.genes - all.ydeg
  
  enrichment <- list()
  for(src in names(net)) {
    trg.ydeg <- sum(names(net[[src]]) %in% names(evd[evd != 0]))
    trg.ndeg <- length(net[[src]]) - trg.ydeg
    contingency <- matrix(c(trg.ydeg, trg.ndeg, all.ydeg - trg.ydeg, all.ndeg - trg.ndeg), ncol=2)
    enrichment[[src]] <- fisher.test(contingency, alternative = 'greater')$p.value
  }
  
  enrichment
}


explain.ornor <- function(net.t, posterior.p, evd, lvl=0.5) {
  x <- posterior.p > lvl
  de.explained = 0
  for(trg in names(evd)) {
    diff.exp <- evd[[trg]]
    if (diff.exp == 0) next
    mask <- x[names(net.t[[trg]])]
    y <- net.t[[trg]][mask]
    y <- y[y!=0]
    if (length(y) == 0) next
    ornor.pred <- min(y)
    de.explained = de.explained + as.integer(ornor.pred == diff.exp)
  }
  
  de.explained
}


transpose.net <- function(net) {
  net.t <- list()
  for(src in names(net)) {
    trg.lst <- net[[src]]
    for(trg in names(trg.lst)) {
      val = trg.lst[[trg]]
      src.lst <- net.t[[trg]]
      if(is.null(src.lst)) src.lst <- c()
      
      val <- c(val)
      names(val) <- c(src)
      net.t[[trg]] <- c(src.lst, val)
    }
  }
  
  net.t
}

get.chip.corrected.network <- function(net.name, remove.non.tfs=TRUE, keep.zeros=TRUE) {
  ents <- read.table('ChIPfilter.ents', header=TRUE)
  rels <- read.table(paste0(net.name, '.rels'), header=TRUE)

  type.dict <- list('increase'=1, 'conflict'=0, 'decrease'=-1)
  rels$type <- unlist(type.dict[rels$type], use.names = FALSE)
  if (!keep.zeros) rels <- rels[rels$type != 0,]

  net <- ents.rels.2.net(ents, rels, remove.non.tfs)

  net
}


get.oe.diff.expr <- function(experiment.name) {
  if(experiment.name %in% c('e2f3', 'myc', 'ras', 'src', 'bcat')) {
    suppressPackageStartupMessages({
        source(paste0('r_scripts/oe_',experiment.name,'_script.R'))
    })
    
    tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
    tT <- clean.tT(tT, use.id = "Gene.ID")
    gset <- rename.gset.features(gset, tT, use.id = "Gene.ID")
  } else if (experiment.name %in% c('tgfb', 'cxcl12')) {
    tT <- read.table(paste0('cell_line_vehicle_',experiment.name,'_evidence_edgeR.txt'), header = TRUE)
    colnames(tT) <- c('Gene.ID', 'adj.P.Val', 'logFC')
    tT <- tT[order(tT$adj.P.Val),]
    tT <- clean.tT(tT)
  } else if (experiment.name %in% c("iHPF_A", "iHPF_B", "N1", "SFT1")) {
    tT <- read.table(paste0('sc_markers_',experiment.name,'_vs_pHPF.csv'), sep = ',', header = TRUE)
    tT$Gene.ID <- paste(tT$Gene.ID)
    colnames(tT) <- c('Gene.ID', 'adj.P.Val', 'logFC')
    tT <- tT[order(tT$adj.P.Val),]
    tT <- clean.tT(tT)
    # df <- read.table('celllines_rna_mean_counts.csv', sep=',', header = TRUE, row.names = 1)
    # rna.expression <- df[[experiment.name]]
    # names(rna.expression) <- df[['Gene.ID']]
    # rna.expression.cutoff <- sort(rna.expression[rna.expression>0])[as.integer(length(rna.expression[rna.expression>0])/4)]
    # gene.is.expressed <- rna.expression >= rna.expression.cutoff
  } else {
    stop(paste("no definitions found for experiment", experiment.name))
  }

  filepath <- paste0('data/', experiment.name, '_experiment_tT.tsv')
  if(!file.exists(filepath))
    write.table(tT[,c('Gene.symbol', 'Gene.ID', 't', 'logFC', 'P.Value', 'adj.P.Val')], filepath, sep = '\t', row.names = F)

  list( tT=tT, gset=gset )
}


get.network <- function(network.name) {
  if(exists(network.name)) {
    regulon <- get(network.name)
    net <- regulon2network(regulon, keep.zeros = FALSE)
  } else if( file.exists(paste0('data/', network.name, '.rels.json'))) {
    net <- jsonlite::fromJSON(paste0('data/', network.name, '.rels.json'))
    net <- sapply(net, unlist)
  } else {
    stop(paste("network", network.name, "not ready for experiment", experiment.name))
  }

  net
}


get.regulon <- function(network.name) {
  if(exists(network.name)) {
    regulon <- get(network.name)
  } else if( file.exists(paste0('data/', network.name, '.rels.json'))) {
    net <- jsonlite::fromJSON(paste0('data/', network.name, '.rels.json'))
    net <- sapply(net, unlist)
    regulon <- network2regulon(net)
  } else {
    stop(paste("network", network.name, "not ready for experiment", experiment.name))
  }

  regulon
}


compute.inference.comparison <- function() {

  network.name <- 'regulonbrca'
  net <- get.network(network.name)
  regulon <- get.regulon(network.name)

  for (experiment.name %in% c('myc', 'e2f3', 'ras')) {
    top.table <- get.oe.diff.expr(experiment.name)
    tT <- top.table[['tT']]
    gset <- top.table[['gset']]

    if (experiment.name == 'ras') {
      logfc.trh = 2.5
    } else {
      logfc.trh = 1
    }

    inference.result <- run.nlbayes(tT, net, uniform.t = TRUE, verbosity=2, logfc.trh = logfc.trh)
    nlbayes.df <- inference.result$result.info$tf.inference

    viper.df <- run.viper(gset, regulon)
    nlbayes.df$viper.pvalue <- viper.df[nlbayes.df$symbol, 'p.value']

    output.filename <- paste0('data/oe_',experiment.name, '_on_net_', network.name,'_nlbayes_and_viper', '.csv')
    write.csv(nlbayes.df, output.filename, row.names = FALSE)
  }
}