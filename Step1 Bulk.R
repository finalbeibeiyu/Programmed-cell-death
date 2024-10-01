
#####################1------limma package for FPKM  data#################

#library
library(limma)
library(pheatmap)
library(ggplot2)
library(EnhancedVolcano)
library(RColorBrewer)

mydata<-read.table("Input.txt",header = T,sep="\t",row.names = 1)
mydata<-as.matrix(mydata)

colnames(mydata)

# Assuming 'geneExp' is your matrix of gene expression data
# and 'group' is your factor variable as defined earlier
# Create a factor for the groups
group <- factor(rep(c("dpi1",  "dpi3", "dpi7", "Sham"), each = 3))

# Define the comparison groups
comparison_groups <- c("dpi1", "dpi3", "dpi7")

# Load the limma package
library(limma)

# Loop through each comparison group
for (group_name in comparison_groups) {
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  # Fit the linear model
  fit <- lmFit(mydata, design)
  
  # Create contrasts: current group vs Young
  contrast_formula <- paste(group_name, "-Sham", sep = "")
  contrast_name <- paste(group_name, "vs Sham", sep = "_")
  contrasts <- makeContrasts(contrasts = setNames(contrast_formula, contrast_name), levels = design)
  
  # Apply the contrasts
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2)
  
  # The correct coefficient name is the one used in contrast formula
  coef_name <- paste(group_name, "-Sham", sep = "")
  
  # Get the top table with the correct coefficient
  allgene <- topTable(fit2, coef = coef_name, adjust.method = "BH", number = Inf)
  # Ensure row names (gene names) are preserved
  allgene <- cbind(Gene = rownames(allgene), allgene)
  # Create a directory for the group if it doesn't exist
  dir_name <- paste("Results", group_name, sep = "_")
  if (!dir.exists(dir_name)) {
    dir.create(dir_name)
  }
  
  # Filter for differentially expressed genes
  diffgene <- allgene[abs(allgene$logFC) >= 1.5 & allgene$P.Value < 0.05, ]
  upgene <- allgene[allgene$logFC >= 1.5 & allgene$P.Value < 0.05, ]
  downgene <- allgene[allgene$logFC <= -1.5 & allgene$P.Value < 0.05, ]
  
  # Write the results to files in the respective directory
  write.table(allgene, file = paste0(dir_name, "/allgene_", group_name, "_vs_Sham.txt"), sep = "\t", quote = FALSE, row.names = F)
  write.table(diffgene, file = paste0(dir_name, "/diffgene_", group_name, "_vs_Sham.txt"), sep = "\t", quote = FALSE, row.names = F)
  write.table(upgene, file = paste0(dir_name, "/upgene_", group_name, "_vs_Sham.txt"), sep = "\t", quote = FALSE, row.names = F)
  write.table(downgene, file = paste0(dir_name, "/downgene_", group_name, "_vs_Sham.txt"), sep = "\t", quote = FALSE, row.names = F)
}



################################2------Volcano plot###############
library(ggplot2)
library(patchwork)
library(dplyr)
library(plot1cell)
library(scRNAtoolVis)
library(jjPlot)
library(ggsci)
diff <- read.table("merge.txt",header = T,sep="\t")
allcolour <- c("#A64036","#F0C2A2","#4182A4")
pdf("3_上下火山图.pdf",width = 10, height = 6)
markerVocalno(markers = diff,
              topn = 5,
              log2FC=0.58,
              labelCol = ggsci::pal_npg()(9))

dev.off()

pdf("3_差异基因火山图.pdf",width = 8, height = 5)
jjVolcano(diffData = diff,
          log2FC.cutoff =1.5,
          pvalue.cutoff=0.05,
          back.col="white",
          aesCol=c('#305BB0','#D33530'),
          legend.position=c(0.95,0.95),
          #myMarkers=myMarkers,
          tile.col=allcolour)
dev.off()

################################3------GSEA###############
#devtools::install_github("junjunlab/GseaVis")
library(clusterProfiler)
library(org.Mm.eg.db)
library(GseaVis)



# dpi1

dpi1 <- read.table("1dpiallgene.txt",header = T,sep="\t") %>%
  arrange(desc(logFC))

genelist1 <- dpi1$logFC
names(genelist1) <- dpi1$id

# dpi3

dpi3 <- read.table("3dpiallgene.txt",header = T,sep="\t") %>%
  arrange(desc(logFC))

genelist2 <- dpi3$logFC
names(genelist2) <- dpi3$id

# dpi7

dpi7 <- read.table("7dpiallgene.txt",header = T,sep="\t") %>%
  arrange(desc(logFC))

genelist3 <- dpi7$logFC
names(genelist3) <- dpi7$id

# check
head(genelist3,3)
#    Aplnr    Foxf1     Bmp5
# 13.45176 13.35322 12.02845

all_glist <- list(genelist1,genelist2,genelist3)
# load tpm

expr <- read.table("wgcna.txt",header = T,sep="\t")



####################其他通路

lapply(1:3, function(x){
  ego3 <- gseGO(geneList     = all_glist[[x]],
                OrgDb        = org.Mm.eg.db,
                ont          = "BP",
                keyType = "SYMBOL",
                nPerm = 10000,
                minGSSize    = 5,
                maxGSSize    = 500,
                pvalueCutoff = 1,
                verbose      = FALSE)
}) -> m_gsea_list

save(m_gsea_list,file = 'm_gsea_list.RData')
df1 <- data.frame(m_gsea_list[[1]])

df2 <- data.frame(m_gsea_list[[2]])
df3 <- data.frame(m_gsea_list[[3]])

write.table(df1,"dpi1.txt",sep = "\t",quote = F,row.names = F)
write.table(df2,"dpi3.txt",sep = "\t",quote = F,row.names = F)
write.table(df3,"dpi7.txt",sep = "\t",quote = F,row.names = F)

load('m_gsea_list.RData')
library(ggsci)
pdf("WntGSEA.pdf",6,4)
GSEAmultiGP(gsea_list = m_gsea_list,
            geneSetID = "GO:0016055",
            exp_name = c("1dpi","3dpi","7dpi"),
            curve.col = c("#A64036","#F0C2A2","#4182A4"),
            addPval = TRUE,
            lineSize=1,
            rect.bm.col=c('#A5CAE7','white','#EE817F'),
            pvalX = 0.99,pvalY = 0.99,
            legend.position = "right",
            rect.bm.label = c("Adult","Neonatal"))
# )
dev.off()

################################4------WGCNA###############
library(fastcluster,warn.conflicts = F)
library(WGCNA,warn.conflicts = F)

rm()
inputdata1="wgcna.txt" 
data0=read.table(inputdata1,sep="\t",row.names=1,header=T,check.names=F,quote="!")
datavar=apply(data0,1,var)

#选择方差大于25%的基因
data0=data0[which(datavar>quantile(datavar, probs = seq(0, 1, 0.25))[4]),]
datExpr = t(data0)
dim(datExpr)


#选择power值 设置软阈值
powers = c(c(1:10), seq(from = 12, to=20, by=2))
SoftTh= pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
pdf("power.pdf",8,6)
par(mfrow = c(1,2));
cex1 = 0.85;
plot(SoftTh$fitIndices[,1], -sign(SoftTh$fitIndices[,3])*SoftTh$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(SoftTh$fitIndices[,1], -sign(SoftTh$fitIndices[,3])*SoftTh$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.85,col="red")#设定R^2 0.85是最好的
plot(SoftTh$fitIndices[,1], SoftTh$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(SoftTh$fitIndices[,1], SoftTh$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


#查看最佳power值
SoftTh$powerEstimate
mypower =5

#计算邻接矩阵
adjacency = adjacency(datExpr, power = mypower)

# 将邻接矩阵转为 TOM矩阵
TOM = TOMsimilarity(adjacency)

# 计算基因之间的相异度
dissTOM = 1-TOM
hierTOM = hclust(as.dist(dissTOM),method="average")



#划分基因模块（使用动态剪切树）
minModuleSize = 40#可以修改
geneTree = hclust(as.dist(dissTOM), method = "average")
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = F,
                            minClusterSize = minModuleSize)
dynamicColors = labels2colors(dynamicMods)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors(GEO)")

#将相关性系数大于0.8（即相异性系数小于0.2)的模块合并掉
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres=0.1
abline(h=MEDissThres, col = "red")

merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

plotDendroAndColors(geneTree, mergedColors,
                    "Merged",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors(GEO)")

#读取临床信息
clinical=read.table("clinical.txt",header = T,sep = "\t",row.names = 1)
rownames(clinical)
traitrr = match(rownames(datExpr), rownames(clinical))
datcli = clinical[traitrr, ]

#模块与性状的相关性
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

afdir <- paste0(getwd(),"/5.WGCNA")           #路径必须中文
dir.create(afdir)
moduleTraitCor = cor(MEs, datcli, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)
library(RColorBrewer)
#连续型
display.brewer.all(type = "seq")
#离散型
display.brewer.all(type = "div")
#极端型
display.brewer.all(type = "qual")

pdf("Module-trait relationships.pdf",7,8)
par(mar = c(6, 8, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(clinical),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = colorRampPalette(c("#305BB0","white","#D33530"))(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

##基因与性状关系（GS）& 基因与模块关系（MM）
modNames = substring(names(MEs), 3)
traitNames = names(clinical)

#计算 MM
MMcor= as.data.frame(cor(datExpr, MEs, use = "p"))
MMP = as.data.frame(corPvalueStudent(as.matrix(MMcor), 
                                     nSamples))
names(MMcor) = paste("MM", modNames, sep="")
names(MMP) = paste("p.MM", modNames, sep="")

#计算 GS
GScor = as.data.frame(cor(datExpr,datcli, use = "p"))
GSP = as.data.frame(corPvalueStudent(as.matrix(GScor), 
                                     nSamples))
names(GScor) = paste("GS.", traitNames, sep="");
names(GSP) = paste("p.GS.", traitNames, sep="");
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#names of those trait
traitNames=names(datcli)

geneTraitSignificance = as.data.frame(cor(datExpr, datcli, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")


####plot MM vs GS for each trait vs each module


##########example:royalblue and CK
#module="royalblue"
#column = match(module, modNames)
#moduleGenes = moduleColors==module

#trait="CK"
#traitColumn=match(trait,traitNames)

#sizeGrWindow(7, 7)


for (trait in traitNames){
  traitColumn=match(trait,traitNames)
  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      
      #sizeGrWindow(7, 7)
      pdf(file=paste(afdir,"/9_", trait, "_", module,"_Module membership vs gene significance.pdf",sep=""),width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
    }
  }
}

#####



#################export GS and MM############### 

geneInfo0 = data.frame(probes= traitNames,
                       moduleColor = modNames)

for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]

write.table(geneInfo, file = paste0(afdir,"/10_GS_and_MM.xls"),sep="\t",row.names=F)

module="yellow"###修改
column = match(module, modNames)
moduleGenes = moduleColors==module

trait="SCI"#修改
traitColumn=match(trait,traitNames)

par(mfrow = c(1,1))
verboseScatterplot(abs(MMcor[moduleGenes, column]),
                   abs(GScor[moduleGenes, traitColumn]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for ",trait),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
abline(h=0.5,v=0.5,col="red")

#查看模块具体的基因
a=which(merge$colors=="black")#修改
b=datExpr[,a]
c=colnames(b)
write.table(c,"black.txt",quote = F,row.names = F,col.names = F)

################################5------wilcoxon rank sum test###############

library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)

data0=read.table("1.txt",header = T,sep = "\t")
data1=read.table("2.txt",header = T,sep = "\t")

data2=merge(data0,data1,by="id")
rownames(data2)=data2[,1]
write.table(data2,"ann.txt",sep = "\t",quote = F,row.names = F)
rt <- data2
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

exp=data

#提取样品类型
Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))

#基因差异分析
sigVec=c()
sigGeneVec=c()
for(i in row.names(data)){
  test=wilcox.test(data[i,] ~ Type)
  pvalue=test$p.value
  Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
  if(pvalue<0.05){
    sigVec=c(sigVec, paste0(i, Sig))
    sigGeneVec=c(sigGeneVec, i)}
}
#提取差异基因表达量
data=data[sigGeneVec,]
outTab=rbind(ID=colnames(data), data)
write.table(outTab, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)
row.names(data)=sigVec

#对差异基因可视化绘图
rownames(Type)=colnames(data)
Type=as.data.frame(Type)
pdf("heatmap.pdf", width=4, height=15)
pheatmap(data,
         annotation=Type,
         color = colorRampPalette(c("#2066B1", "white", "#A81F24"))(50),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=6,
         fontsize_row=7,
         fontsize_col=6,
         border_color="white")
dev.off()


#把表达数据转换成ggplot2文件输入
exp=as.data.frame(t(exp))
exp=cbind(exp, Type=Type)
data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

#
p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
            xlab="",
            ylab="Gene expression",
            legend.title="Type",
            palette = c("#A81F24", "#2066B1"),
            width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")+coord_flip()

#
pdf(file="boxplot.pdf", width=3,height=22)
print(p1)
dev.off()