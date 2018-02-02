### TITLE : Get Gene Expression Profiles for the HPDE cell lines induced for the over-expression of KRAS_G12D, KRAS_wt, GFP
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es

### 0 Load libraries #########
library(Biobase)
library(GEOquery)
library(limma)
library(WGCNA)

### 1 read ###########
gset <- getGEO(filename = "./LVL1/GSE58055_series_matrix.txt.gz", GSEMatrix =TRUE, AnnotGPL=FALSE)

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

### 2 Define the biological groups ########
sample.groups <- setNames(as.character(pData(gset)$description.1),sampleNames(gset))
sample.groups <- gsub("HPDE-","",sample.groups)
sample.groups <- gsub("\\+","pos",sample.groups)
sample.groups <- gsub("\\-","neg",sample.groups)
sample.groups <- gsub("^Ind","",sample.groups)

groups1 <- factor(sapply(sample.groups,function(z) strsplit(z,split=" ")[[1]][1]),levels=c("KRASWT","KRASG12D","GFP"))
                    
groups2 <- factor(sapply(sample.groups,function(z) strsplit(z,split=" ")[[1]][2]),levels=c("negDOX","posDOX"))
groups3 <- factor(paste(as.character(groups1),as.character(groups2),sep="."))
groups3 <- relevel(groups3,ref="KRASG12D.negDOX")

# Set the groups
pData(gset)$indGene <- groups1
pData(gset)$trt <- groups2
pData(gset)$Group <- groups3


### 2 Save the LVL3 data (measures by probes)
if(!dir.exists("./LVL2")) dir.create("./LVL2");
saveRDS(gset,file = "./LVL2/eset.rds")

### 3 Collapse probe intensities by gene symbol using the MaxMean

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { 
  warning("The expression was log2-transformed.")
  ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex)
  }

stopifnot(all(rownames(fData(gset)$GENE_SYMBOL) == featureNames(gset)))
gene.symbols <- setNames(as.character(fData(gset)$GENE_SYMBOL),featureNames(gset))
gene.symbols[gene.symbols==""] <- "N/A"

cr_obj.MaxMean <- collapseRows(exprs(gset),
                               rowGroup = gene.symbols,
                               rowID = rownames(gset),
                               method="MaxMean")
datETcollapsed <- cr_obj.MaxMean$datETcollapsed
MAT <- datETcollapsed[which(rownames(datETcollapsed)!="N/A"),]


### 4 Create new eset #####
eset.sym <- ExpressionSet(assayData = MAT,
                          phenoData = new("AnnotatedDataFrame", data=pData(gset)),
                          annotation = annotation(gset))
### 5 Save the new eset into the lvl3
if(!dir.exists("./LVL3")) dir.create("./LVL3")

saveRDS(eset.sym,"./LVL3/eset.rds")
