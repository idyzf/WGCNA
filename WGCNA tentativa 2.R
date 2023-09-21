library(WGCNA)
library(dplyr)
library(openxlsx)
rm(list = ls())
options(stringsAsFactors = FALSE)

femData<- read.csv("TCGA_GBM_TPM_counts.csv", head = TRUE, sep=",")
rownames(femData)<- femData$Gene
#femData1<- femData
#rownames(femData1) <- femData$Name#[-1:-14, -1:-3]
datExpr0 = femData[,-1];
#datExpr0<- datExpr0[-1:-3, -1:-14]
#datExpr0<- as.matrix(unlist(datExpr0))
#names(datExpr0) = femData$substanceBXH;
datExpr0<- as.data.frame(t(datExpr0))
rownames(datExpr0)


gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 60000, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 60000, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


traitData = read.xlsx("TCGA_Clinical_traits.xlsx");
dim(traitData)
names(traitData)

rownames(traitData) <- traitData$CGGA_ID
allTraits<- traitData

femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$barcode);
datTraits = allTraits[traitRows, -1:-2];
rownames(datTraits) = allTraits[traitRows, 1];
#datExpr<- t(datExpr)
collectGarbage();
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
#save(datExpr, datTraits, file = "FemaleLiver-01-dataInput.RData")

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#se os genes estiverem em int precisam seer passador para num
i=1
for (i in 1:ncol(datExpr)){ 
  datExpr[,i] <- as.numeric(as.character(datExpr[,i]))
  print(i)
}
#

net = blockwiseModules(datExpr, power = 9,
                       TOMType = "unsigned", minModuleSize = 92,
                       reassignThreshold = 0, mergeCutHeight = 0.07,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM", 
                       verbose = 3)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "FemaleLiver-02-networkConstruction-auto.RData")

#########################################################################
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
###################################ideli fazendo arte #############################
##################################################################################
moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
genes_filter <- subset(moduleTraitPvalue, moduleTraitPvalue$Gender <= 0.05 |  moduleTraitPvalue$OS_groups <= 0.05 | 
                         moduleTraitPvalue$Censor..alive.0..dead.1. <= 0.05 | 
                         moduleTraitPvalue$Radio_status..treated.1.un.treated.0. <= 0.05 | 
                         moduleTraitPvalue$Chemo_status..TMZ.treated.1.un.treated.0. <= 0.05|
                         moduleTraitPvalue$IDH.status <= 0.05 | moduleTraitPvalue$MGMT <= 0.05)
genes_filter <- subset(moduleTraitPvalue, moduleTraitPvalue$IDH.status <= 0.05 )

modulenames<- rownames(genes_filter)
MEs0<- as.data.frame(t(MEs0))
MEs0_filter <- subset(MEs0, rownames(MEs0) %in% modulenames)



MEs0 <- MEs0_filter
MEs0<- as.data.frame(t(MEs0))
#go back to line 138
#######################################################################################

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$surv.group);
names(weight) = "surv.group"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

module ="turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for IDH status",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

geneList<-names(datExpr)[moduleColors=="lightcyan"]
paste(geneList,collapse = ", ")

#######################################GO Analysis############################
library(clusterProfiler)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(GO.db)
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

# Entrez gene ID
head(gene)

ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

head(ggo)

ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)

ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(ego2, 3)
ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
ego3$ID
ego3$Description
goplot(ego)


library(readr)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(openxlsx)

ggplot(ego3, aes(x=enrichmentScore,y= Description,size=setSize,color= -pvalue)) +
  geom_point(alpha=0.5) + scale_size(range = c(3, 15), name="Number of genes") + 
  scale_colour_gradientn(colours=c( "lightblue")) + ggtitle("Biological Process") 
dev.off()

