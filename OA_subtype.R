library(SingleCellExperiment)
library(scater)
library(SC3)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(grid)
library(ComplexHeatmap)
library(GenomicAlignments)
library(gage)
library(pathview)
library(gageData)
library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGG.db)

setwd("C:/works/OA_Project_submit")
raw_matrix <- read.table("cartilage.count.201905113.txt",header = TRUE , row.names =1)
length(colnames(raw_matrix))
raw_matrix <- raw_matrix[,6:ncol(raw_matrix)]
dim(raw_matrix)
colSums(raw_matrix)[order(colSums(raw_matrix))]
geneNumber <- colSums(raw_matrix>0)
geneNumber[order(geneNumber)]
sample_matrix <- read.table("cart_ann_20190513.txt",header = TRUE , row.names =1)
sample_matrix <- sample_matrix[colnames(raw_matrix),]
dim(sample_matrix)
sampleInfor <- sample_matrix[order(rownames(sample_matrix)),]
treat <- rep("OA",nrow(sampleInfor))
sampleInfor <- cbind(sampleInfor,treat )
sampleInfor$treat <- as.character(sampleInfor$treat)
sampleInfor[c("P2_050","P2_197","P2_201","P2_205","P2_209"),]$treat<-rep("Normal",5)
# filter samples according totoal reads 
matrix <- raw_matrix[,colSums(raw_matrix) > 500000 ]
matrix <- matrix[,colSums(matrix) < 2000000]
dim(matrix)  #[1] 58684   135
matrix.nor <- t(t(matrix)/colSums(matrix)*1000000)    # calculate CPM
matrix  <- matrix.nor  

sce <- SingleCellExperiment(assays = list(counts = as.matrix(matrix),logcounts = log2(as.matrix(matrix) + 1)), colData = sampleInfor[colnames(matrix),] )
rowData(sce)$feature_symbol <- rownames(sce)

# filter genes var lage in normal
sampleInfor.normal <- sampleInfor[c("P2_197","P2_201","P2_205","P2_209"),]
sce.normal <- sce[,rownames(sampleInfor.normal)]
summary ( apply(exprs(sce.normal), 1, var) )
sum(apply(exprs(sce.normal), 1, var) > 3)  #[1] 4423
drop_genes.01 <- apply(exprs(sce.normal), 1, function(x) {var(x) > 3 }) 
matrix.filter01 <- matrix[!drop_genes.01, ]       # filter genes var lage in normal
col_data <- colData(sce)

# filter samples 
matrix.filter02 <- matrix.filter01[,!(colnames(matrix.filter01) %in% c("P2_197","P2_201","P2_205","P2_209")) ]  # delet normal samples
dim(matrix.filter02)  #[1] 54261   131
new_sampleInfor <- col_data[colnames(matrix.filter02),]
new_sampleInfor$sampleIndex <-paste0("S",new_sampleInfor$sampleID)
new_sampleInfor$libID <- rownames(new_sampleInfor) 
colnames(matrix.filter02) <- new_sampleInfor[colnames(matrix.filter02),]$sampleIndex   # change colname of the matrix as sampleName
rownames(new_sampleInfor) <-  new_sampleInfor$sampleIndex       # change infor colenames
dim(matrix.filter02)   # [1] 54261   131
new_sampleInfor.final <- new_sampleInfor[colnames(matrix.filter02),]
dim(new_sampleInfor.final)

matrix.filter03 <- matrix.filter02[rowSums(matrix.filter02>5)>15,]  # filter genes with low count
dim(matrix.filter03)    # [1] 9028  131

#==========================
varGene <- apply(matrix.filter03, 1, var)[order(apply(matrix.filter03, 1, var), decreasing = T)]
matrix.filter.final <- matrix.filter03[names(varGene[1:4000]),]

sce.final <- SingleCellExperiment(assays = list(counts = matrix.filter.final,logcounts = log2(as.matrix(matrix.filter.final) + 1)), colData = new_sampleInfor.final )
rowData(sce.final)$feature_symbol <- rownames(sce.final)
sce.final <- sc3_estimate_k(sce.final)
str(metadata(sce.final)$sc3)       # $ k_estimation: num 3
sce.final <- sc3(sce.final, ks = 2:8, biology = TRUE,n_cores=3,gene_filter = F) #
sc3_interactive(sce.final)
#save(sce.final,file="OA.sec.final.Rdata")
load(file="OA.sec.final.Rdata")
rowData(sce.final)$feature_symbol  <- rownames(sce.final)
col_data <- colData(sce.final)
row_Data <- rowData(sce.final)
table(col_data[,"sc3_4_clusters"])
sc3_export_results_xls(sce.final,file="cart02_sc3_result.xls")
write.table(rowData(sce.final),file="cart02_sc3_result_gene.xls",quote = F, sep = "\t")

####### marker genes 
mkGene.name <- c()
rowData(sce.final)$feature_symbol <- rownames(sce.final)
head(row_Data[ , grep("sc3_4", colnames(row_Data))])
row_Data <- rowData(sce.final)

cart_mk4_c1 <- row_Data[which(row_Data$sc3_4_markers_clusts == 1 & row_Data$sc3_4_markers_padj < 0.05 ),]
cart_mk4_c1 <- cart_mk4_c1[order(cart_mk4_c1$sc3_4_markers_auroc,decreasing = TRUE),]
dim(cart_mk4_c1)  # [1] 211  30  
cart_mk4_c2 <- row_Data[which(row_Data$sc3_4_markers_clusts == 2 & row_Data$sc3_4_markers_padj < 0.05 ),]
cart_mk4_c2 <- cart_mk4_c2[order(cart_mk4_c2$sc3_4_markers_auroc,decreasing = TRUE),]
dim(cart_mk4_c2)  # [1] 37 30  
cart_mk4_c3 <- row_Data[which(row_Data$sc3_4_markers_clusts == 3 & row_Data$sc3_4_markers_padj < 0.05 ),]
cart_mk4_c3 <- cart_mk4_c3[order(cart_mk4_c3$sc3_4_markers_auroc,decreasing = TRUE),]
dim(cart_mk4_c3)  # [1] 86 30  #  0.901652892561983
cart_mk4_c4 <- row_Data[which(row_Data$sc3_4_markers_clusts == 4 & row_Data$sc3_4_markers_padj < 0.05 ),]
cart_mk4_c4 <- cart_mk4_c4[order(cart_mk4_c4$sc3_4_markers_auroc,decreasing = TRUE),]
dim(cart_mk4_c4)  # [1] 230  30    # 0.836413043478261
###***********


##################   syno   
syno_nor <- read.table("../syno/syno_PC_count.matrix.nor.txt",head=T)
syno_ann <- read.table('../syno/syno_PC_sampleInfor.mat.nor.txt',head=T)

dim(syno_nor)    #[1] 15510    85
dim(syno_ann)    #  [1] 85 12
#syno_nor <- syno_nor[rowSums(syno_nor>10)>15,]
dim(syno_nor) 
head(col_data[,1:3])
clusterFromCart <- c()
k<-1
for(i in 1:nrow(syno_ann)){
  for(j in 1:nrow(col_data)){
    if (syno_ann$sampleID[i] == as.character(col_data$sampleID[j])){
      clusterFromCart[k] <- paste("cart",col_data$sc3_4_clusters[j],sep = "_")
      k<-k+1
    }
    next
  }
}
syno_ann.clustered <- cbind(syno_ann[syno_ann$sampleID %in% col_data$sampleID,],clusterFromCart)
write.table(syno_ann.clustered,file="../syno/syno_ann.fromcart.txt",quote = F,sep="\t")
##################   subc   #########
subc_nor <- read.table("../subc/subc_PC_count.matrix.nor.txt",head=T)
subc_ann <- read.table('../subc/subc_PC_sampleInfor.mat.nor.txt',head=T)
dim(subc_nor)     # [1] 13400    79     
dim(subc_ann)
head(col_data)
#subc_ann[subc_ann$sampleID %in% cart_p_data$sampleID,]$sampleID
clusterFromCart <- c()
k<-1
for(i in 1:nrow(subc_ann)){
  for(j in 1:nrow(col_data)){
    if (subc_ann$sampleID[i] == as.character(col_data$sampleID[j])){
      clusterFromCart[k] <- paste("cart",col_data$sc3_4_clusters[j],sep = "_")
      k<-k+1
    }
    next
  }
}
subc_ann.clustered <- cbind(subc_ann[subc_ann$sampleID %in% col_data$sampleID,],clusterFromCart)
write.table(subc_ann.clustered,file="../subc/subc_ann.fromcart.txt",quote = F,sep="\t")
write.table(col_data,file="tissues/in/cart.col_data.txt",quote = F,sep="\t")






