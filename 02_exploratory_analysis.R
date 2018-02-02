### TITLE : Exploratory analysis of the dataset
### AUTHOR : Perales-Paton - Javier (jperales@cnio.es)

library(ggplot2)
library("GSVA")
library(GSEABase)


### 1 Load the data ####
eset <- readRDS("./LVL3/eset.rds")    # Eset, already collapsed by gene symbol. RMA-quantile normalized.
indGene <- pData(eset)$indGene        # Vector of induced gene for each sample, in case of +DOX treatment
trt <- pData(eset)$trt                # +Dox OR -DOX treatment, for the induction of the expression of a given gene
groups <- pData(eset)$Group           # Biological group
MAT <- exprs(eset)                    # Gene expression matrix

# We check that the order of the groups are the name as the colnames in the expression matrix
stopifnot(all(colnames(MAT)==rownames(pData(eset))))

### 2 Exploratory analysis to decide what to do in the contrasts ######
if(!dir.exists("./Imgs")) dir.create("./Imgs")

# Q: Is the KRAS expression already expressed in HPDE cells?
KRASexpr <- data.frame("KRAS_expr"=as.numeric(MAT["KRAS",]),
                       "induced_gene"=indGene,
                       "treatment"=trt,
                       "Group"=groups)
boundaries <- apply(MAT,1, function(z) quantile(z,probs = c(0.25,0.75),na.rm = TRUE))

png("./Imgs/KRASexpr.png",width = 750*3,height = 600*3,res=280)
ggplot(KRASexpr, aes(x=induced_gene,y=KRAS_expr,fill=treatment,colour=treatment)) + 
  geom_boxplot() +   geom_point(position=position_jitterdodge(),alpha=0.25,stroke=2,size=3) +
  # geom_point(alpha=0.25,size=8,stroke=2) +
  # geom_jitter(alpha=0.25,size=4,stroke=2,aes(x=induced_gene,y=KRAS_expr,fill=treatment,colour=treatment)) +
  scale_y_continuous(name = "KRAS expression (log2 intensitiy)",limits=c(0,17)) +
  scale_color_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent") +
  theme_bw() + ggtitle("GSE58055: KRAS-induced transcription analysis in immortalized pancreatic cells") +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=26),
        axis.text.x=element_text(size = 22),
        axis.text.y=element_text(size = 28),
        legend.position = "bottom")
dev.off()


# Q: Take a look at the HALLMARK of KRAS signaling pathway
GSC <- getGmt("./Ext_data/h.all.v6.1.symbols.gmt")
GSC <- GSC[grepl("KRAS",names(GSC))]

ssgsea.es <- gsva(MAT,GSC,method=c("ssgsea"),kcdf="Gaussian")
gsva.es <- gsva(MAT,GSC,method=c("gsva"),mx.diff=FALSE,kcdf="Gaussian")

# tes <- es[1,] - es[2,]

png("./Imgs/GSE58055_explAnalysis.png",width = 1000*3,height = 1000*3,res=330)
par(mar=c(4,4,6,4),mfrow=c(2,2),cex.lab=1)

plot(ssgsea.es[1,],ssgsea.es[2,],bg=as.integer(indGene),pch=c(25,24)[as.integer(trt)],
     xlab=paste0("ssGSEA score - ",rownames(ssgsea.es)[1]),
     ylab=paste0("ssGSEA score - ",rownames(ssgsea.es)[2]),las=1,cex=1.5,
     main="ssGSEA")

plot(gsva.es[1,],gsva.es[2,],bg=as.integer(indGene),pch=c(25,24)[as.integer(trt)],
     xlab=paste0("GSVA score - ",rownames(gsva.es)[1]),
     ylab=paste0("GSVA score - ",rownames(gsva.es)[2]),las=1,cex=1.5,
     main="GSVA")

pca.obj <- prcomp(t(na.omit(MAT)))
plot(PC2 ~ PC1,as.data.frame(pca.obj$x),
     bg=as.integer(indGene),pch=c(25,24)[as.integer(trt)],
     las=1,main="Principal Component Analysis\n(whole transcriptome)",cex=1.5,
     xlab=paste0("PC1 (",round(summary(pca.obj)$importance[2,1]*100,digits=3),"% var. explained)"),
     ylab=paste0("PC2 (",round(summary(pca.obj)$importance[2,2]*100,digits=3),"% var. explained)"))

plot(1:10,1:10,type="n", axes=F, xlab="", ylab="")
legend("top",legend = levels(indGene),pt.bg=c(1:3),pch=21,cex=1.5)
legend("bottom",legend = levels(trt),pt.bg="white",pch=c(25,24),cex=1.5)

title("GSE58055: Exploratory Analysis", outer=TRUE,line=-2)
dev.off()
