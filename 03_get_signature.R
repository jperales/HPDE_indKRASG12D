### TITLE : Extract a gene expression signature of oncogenic KRAS in HPDE ductal cells.
### AUTHOR : Perales-Paton - Javier (jperales@cnio.es)
### MOTIVATION : It is well described that KRAS_G12D is almost universal in PDAC.
### However, there are a lot 


library(Biobase)
library(limma)
library(ggplot2)
library(GSEABase)

### 1 Load the data ####
eset <- readRDS("./LVL3/eset.rds")    # Eset, already collapsed by gene symbol. RMA-quantile normalized.
indGene <- pData(eset)$indGene        # Vector of induced gene for each sample, in case of +DOX treatment
trt <- pData(eset)$trt                # +Dox OR -DOX treatment, for the induction of the expression of a given gene
groups <- pData(eset)$Group           # Biological group
MAT <- exprs(eset)                    # Gene expression matrix, collapsed by gene symbol. RMA-quantile normalized

# We check that the order of the groups are the name as the colnames in the expression matrix
stopifnot(all(colnames(MAT)==rownames(pData(eset))))

### 2 Make contrasts of interest ######
if(!dir.exists("./DiffExpr")) dir.create("./DiffExpr")
design <- model.matrix(~ 0+groups)
colnames(design) <- gsub("^groups","",colnames(design))
fit <- lmFit(MAT, design)
cont.matrix <- makeContrasts("KRASG12D"=KRASG12D.posDOX - KRASG12D.negDOX,
                             "KRASG12D_woKRASWT_woGFP"=(KRASG12D.posDOX - KRASG12D.negDOX) - (KRASWT.posDOX - KRASWT.negDOX) - (GFP.posDOX - GFP.negDOX),
                             "KRASWT"=KRASWT.posDOX - KRASWT.negDOX,
                             "KRASG12D_woGFP"=(KRASG12D.posDOX - KRASG12D.negDOX) - (GFP.posDOX - GFP.negDOX),
                             "KRASG12D_woKRASWT"=(KRASG12D.posDOX - KRASG12D.negDOX) - (KRASWT.posDOX - KRASWT.negDOX),
                             "KRASWT_woGFP"=(KRASWT.posDOX - KRASWT.negDOX) - (GFP.posDOX - GFP.negDOX),
                             "GFP"=(GFP.posDOX - GFP.negDOX),
                             levels = design)

fit2 <- contrasts.fit(fit, cont.matrix)

eBay <- eBayes(fit2)
topTab <- topTable(eBay,number = Inf)

# plot(KRASG12D ~ KRASG12D_woKRASWT,as.data.frame(eBay$coefficients))

saveRDS(eBay,"./DiffExpr/eBay.rds")
saveRDS(topTab,"./DiffExpr/topTab.rds")

# Just to see how different are the different contrasts
png("./DiffExpr/scatter_moderated-t_across_contrasts.png",width = 1200*3,height = 1200*3,res=280)
pairs(eBay$t,
      labels = gsub("_","\n",colnames(eBay$t)),
      cex.labels = 1.7,
      xlim=c(-40,40),ylim=c(-40,40))
dev.off()

# The DEGs for each contrast:
DEGs <- apply(decideTests(eBay),2, function(z) table(z)[c("-1","0","1")])
DEGs[is.na(DEGs)] <- 0
rownames(DEGs) <- c("Down","nonDiff","Up")
write.table(DEGs,"./DiffExpr/DEGs_across_contrasts.tsv",sep="\t",col.names = NA,row.names = TRUE,quote=FALSE)


# Finally, generate the gene set collection ######
if(!dir.exists("./Signatures")) dir.create("./Signatures")

GSC <- GeneSetCollection(list(
  # KRASG12D over-expression in HPDE_p53Rb1.loss_N150_UP
  GeneSet(setName="HPDE_KRASG12D.OE_p53Rb1.Loss_N150_UP",
          shortDescription="Up-regulated genes (n=150) - HPDE with loss of p53 and inactivation of the Rb pathway, and Over-Expression of KRAS G12D - source: GSE58055.",
          geneIds=names(sort(eBay$t[,"KRASG12D"],decreasing = TRUE)[1:150])),
  GeneSet(setName="HPDE_KRASG12D.OE_p53Rb1.Loss_N150_DN",
          shortDescription="Down-regulated genes (n=150) - HPDE with loss of p53 and inactivation of the Rb pathway, and Over-Expression of KRAS G12D - source: GSE58055.",
          geneIds=names(sort(eBay$t[,"KRASG12D"],decreasing = FALSE)[1:150])),
  # KRASG12D over-expression adjusted by the diff. expressed genes by KRASwt over-expression
  GeneSet(setName="HPDE_KRASG12D.OE_adjKRASwt.OE_p53Rb1.Loss_N150_UP",
          shortDescription="Up-regulated genes (n=150) - HPDE with loss of p53 and inactivation of the Rb pathway, and Over-Expression of KRAS G12D but substracting those differences due to Over-expression of KRAS wt - source: GSE58055.",
          geneIds=names(sort(eBay$t[,"KRASG12D_woKRASWT"],decreasing = TRUE)[1:150])),
  GeneSet(setName="HPDE_KRASG12D.OE_adjKRASwt.OE_p53Rb1.Loss_N150_DN",
          shortDescription="Down-regulated genes (n=150) - HPDE with loss of p53 and inactivation of the Rb pathway, and Over-Expression of KRAS G12D but substracting those differences due to Over-expression of KRAS wt - source: GSE58055.",
          geneIds=names(sort(eBay$t[,"KRASG12D_woKRASWT"],decreasing = FALSE)[1:150])),
  # KRASwt over-expression
  GeneSet(setName="HPDE_KRASwt.OE_p53Rb1.Loss_N150_UP",
          shortDescription="Up-regulated genes (n=150) - HPDE with loss of p53 and inactivation of the Rb pathway, and Over-Expression of KRASwt - source: GSE58055.",
          geneIds=names(sort(eBay$t[,"KRASWT"],decreasing = TRUE)[1:150])),
  GeneSet(setName="HPDE_KRASwt.OE_p53Rb1.Loss_N150_DN",
          shortDescription="Down-regulated genes (n=150) - HPDE with loss of p53 and inactivation of the Rb pathway, and Over-Expression of KRASwt - source: GSE58055.",
          geneIds=names(sort(eBay$t[,"KRASWT"],decreasing = FALSE)[1:150]))
))

toGmt(GSC,"./Signatures/HPDE_KRASG12D_p53Rb1PathwayLoss.gmt")
saveRDS(GSC,"./Signatures/HPDE_KRASG12D_p53Rb1PathwayLoss.GSC.rds")
