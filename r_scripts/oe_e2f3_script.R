# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis with limma
suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(umap)
})

# load series and platform data from GEO

gset <- getGEO("GSE3151", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "1111111111XXXXXXXXXXXXXXXXXXXXXXXXXX000000000XXXXXXXXXX"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("pt","wt"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)


# tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
# tT <- clean.tT(tT, use.id = "Gene.ID")
# 
# tT <- tT[abs(tT$logFC)>1,]
# tT$adj.P.Val.2 <- p.adjust(tT$P.Value, 'fdr')
# 
# nrow(tT[tT$adj.P.Val<0.01,])
# nrow(tT[tT$adj.P.Val.2<0.01,])
# 
# plot(tT[tT$adj.P.Val<0.05, 'adj.P.Val'])
# plot(tT[tT$adj.P.Val.2<0.05, 'adj.P.Val.2'])
# 
