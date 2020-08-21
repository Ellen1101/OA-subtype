#setwd("~/figures")

############## plot expression heatmap
col_data <- colData(sce.final)
table(col_data[,"sc3_4_clusters"])   # 81 24 10 16

library(grid)
library(ComplexHeatmap)
k<-4
dataset <- get_processed_dataset(sce.final)
hc <- metadata(sce.final)$sc3$consensus[[as.character(k)]]$hc

ht_ann <- col_data
ht_ann$age_bin <- sub("^41_50","<60",ht_ann$age_bin)
ht_ann$age_bin <- sub("^51_60","<=60",ht_ann$age_bin)
ht_ann$age_bin <- sub("^81_90",">80",ht_ann$age_bin)
ht_ann$age_bin <- sub("^N","Null",ht_ann$age_bin)
ht_ann$gender <- sub("^N","Null",ht_ann$gender)
ht_ann$KL  <- sub("^0","N",ht_ann$KL)
ht_ann$KL  <- sub("^N","Null",ht_ann$KL)
ht_ann$JSN  <- sub("^0","N",ht_ann$JSN)
ht_ann$JSN  <- sub("^N","Null",ht_ann$JSN)
ht_ann$Osteophytes  <- sub("^0","N",ht_ann$Osteophytes)
ht_ann$Osteophytes  <- sub("^N","Null",ht_ann$Osteophytes)
# change cluster 
ht_ann$sc3_4_clusters  <- sub("^1","C1",ht_ann$sc3_4_clusters)
ht_ann$sc3_4_clusters  <- sub("^2","C2",ht_ann$sc3_4_clusters)
ht_ann$sc3_4_clusters  <- sub("^3","C3",ht_ann$sc3_4_clusters)
ht_ann$sc3_4_clusters  <- sub("^4","C4",ht_ann$sc3_4_clusters)


col.C1=rgb(31, 149, 203, max = 255)
col.C2=rgb(255, 233, 0, max = 255)
col.C4=rgb(236, 95, 92, max = 255)
col.C3=rgb(236,171,85, max = 255)

annotation_col = data.frame(row.names = rownames(ht_ann), cluster=as.factor(ht_ann$sc3_4_clusters), gender=as.factor(ht_ann$gender), age=as.factor(ht_ann$age_bin), KL=as.factor(ht_ann$KL), JSN=as.factor(ht_ann$JSN), Osteophytes=as.factor(ht_ann$Osteophytes) )
ann_colors <- list(
  cluster = c("C1"=col.C1, "C2"=col.C2,"C3"=col.C3,"C4"=col.C4)
  ,gender = c("F"=rgb(255, 240, 0, max = 255), "M"=rgb(0, 150, 211, max = 255), "Null"=rgb(246, 246, 246, max = 255))
  ,age = c("Null"=rgb(246, 246, 246, max = 255),"<=60"=rgb(0, 150, 211, max = 255),"61_70"=rgb(171,197,174, max = 255),"71_80"=rgb(255, 240, 0, max = 255),">80"=rgb(255, 74, 75, max = 255))
  ,KL = c("Null"=rgb(246, 246, 246, max = 255),"2"=rgb(0, 150, 211, max = 255),"3"=rgb(108,173,193, max = 255),"4"=rgb(255, 240, 0, max = 255))
  ,JSN = c("Null"=rgb(246, 246, 246, max = 255),"1"=rgb(0, 150, 211, max = 255),"2"=rgb(0, 150, 211, max = 255),"3"=rgb(108,173,193, max = 255),"4"=rgb(171,197,174, max = 255),"5"=rgb(224, 217, 152, max = 255),"6"=rgb(255, 240, 0, max = 255))
  ,Osteophytes = c("Null"=rgb(246, 246, 246, max = 255),"3"="#0096D3","4"="#1CA0BB","5"="#38AAA4","6"="#55B48C","7"="#71BE75","8"="#8DC85D","9"="#AAD246","10"="#C6DB2E","11"="#E2E517","12"="#FFF000")
)
ha_column = HeatmapAnnotation(df = annotation_col, col = ann_colors,gap = unit(0.5, "mm"),gp = gpar(col = "grey"),annotation_height = 0.35,annotation_name_gp = gpar(col = "black")
                              ,annotation_legend_param = list(JSN=list(ncol=2,title_position="topcenter"),Osteophytes=list(ncol=2,title_position="topcenter")) )

#color.palette =c(colorRampPalette(c("#8e162e", "#f69273"))(30),colorRampPalette(c("#f69273", "#ffffff", "#91cade"))(15),colorRampPalette(c("#91cade", "#265091"))(30))
color.palette =c(colorRampPalette(c("#265091", "#91cade"))(30),colorRampPalette(c("#91cade", "#ffffff", "#f69273"))(10),colorRampPalette(c("#f69273", "#8e162e"))(30))

Dset <- dataset[order(apply(dataset,1,var),decreasing = T),][1:1000,]
#Dset <- Dset[order(apply(Dset,1,var),decreasing = T),][1:900,]
ht1 = Heatmap(t(scale(t(Dset))), cluster_columns= hc, show_row_names = F ,show_row_dend=F,show_column_names=F,
              #              ,column_dend_gp=,column_dend_reorder=
              name="",top_annotation_height = unit(2.5, "cm"),col=color.palette, top_annotation = ha_column)
pdf("figures/expression.heatmap.pdf",7,7)
draw( ht1, annotation_legend_side = "bottom")
dev.off()

#########****************end plot expression heatmap


##***********
### caculate the clinical score in each cluster ###########

pdf("clinical.score_k4.pdf",7,7)

KL_score <- table(ht_ann$KL,as.character(ht_ann[,"sc3_4_clusters"]),exclude="Null")
barplot(prop.table(KL_score,2)*100,legend=T,xlim=c(0,6),col=c("#c6dbef","#4292c6","#a6bddb","#3690c0"),main="KL score",ylab = "Percentage(%)")

JSN_score <- table(ht_ann$JSN,as.character(ht_ann[,"sc3_4_clusters"]),exclude="Null")
barplot(prop.table(JSN_score,2)*100,legend=T,xlim=c(0,6),col=c("#efedf5","#dadaeb","#bcbddc","#9e9ac8","#807dba","#6a51a3","#54278f"),main="JSN score",ylab = "Percentage(%)")

OST_score <- table(ht_ann$Osteophytes,as.character(ht_ann[,"sc3_4_clusters"]),exclude="Null")
OST_score <- rbind(OST_score[4:9,],OST_score[1:3,])
barplot(prop.table(OST_score,2)*100,legend=T,xlim=c(0,6),col=c("#fff5eb","#fee6ce","#fdd0a2","#fdae6b","#fd8d3c","#f16913","#d94801","#a63603","#7f2704"),main="Osteophyte score",ylab = "Percentage(%)")

age_table<-table(ht_ann$age_bin,as.character(ht_ann[,"sc3_4_clusters"]),exclude="Null")
age_table <- rbind(age_table[1,],age_table[3:4,],age_table[2,])
rownames(age_table) <- c("<=60","61_70","71_80",">80")
barplot(prop.table(age_table,2)*100,legend=T,xlim=c(0,6),col=c(brewer.pal(8, "YlGn"))[3:7],main="Age",ylab = "Percentage(%)")

### cacullate significance
library(DescTools)
CochranArmitageTest(t(OST_score[,c(2,1)]),alternative = "decreasing")   # Z = 1.8084, dim = 9, p-value = 0.03527
CochranArmitageTest(t(OST_score[,c(2,3)]),alternative = "decreasing")   # Z = 2.0088, p-value = 0.02228 
CochranArmitageTest(t(OST_score[,c(2,4)]),alternative = "decreasing")   # Z = 0.88814, p-value = 0.1872
p.adjust(c(CochranArmitageTest(t(OST_score[,c(2,1)]),alternative = "decreasing")$p.value,CochranArmitageTest(t(OST_score[,c(2,3)]),alternative = "decreasing")$p.value,CochranArmitageTest(t(OST_score[,c(2,4)]),alternative = "decreasing")$p.value),method = "fdr",n=3)
#[1] 0.05291117 0.05291117 0.18723327

CochranArmitageTest(t(OST_score[,c(1,3)]),alternative = "decreasing")   # Z = 1.5509, dim = 9, p-value = 0.06046
CochranArmitageTest(t(OST_score[,c(1,4)]),alternative = "increasing")   # Z = -0.50297, dim = 9, p-value = 0.3075
CochranArmitageTest(t(OST_score[,c(3,4)]),alternative = "increasing")   # Z = -1.5882, dim = 9, p-value = 0.05613 
p.adjust(c(CochranArmitageTest(t(OST_score[,c(1,3)]),alternative = "decreasing")$p.value,CochranArmitageTest(t(OST_score[,c(1,4)]),alternative = "increasing")$p.value,CochranArmitageTest(t(OST_score[,c(3,4)]),alternative = "increasing")$p.value),method = "fdr",n=3)
#[1] 0.09068595 0.32130460 0.09068595
##
CochranArmitageTest(t(JSN_score[,c(4,1)]),alternative = "decreasing")   # Z = 1.7548, p-value = 0.03964 
CochranArmitageTest(t(JSN_score[,c(4,2)]),alternative = "decreasing")   # Z = 0.28015, p-value = 0.3897 
CochranArmitageTest(t(JSN_score[,c(4,3)]),alternative = "decreasing")   # Z = 1.2372, p-value = 0.108
p.adjust(c(CochranArmitageTest(t(JSN_score[,c(4,1)]),alternative = "decreasing")$p.value,CochranArmitageTest(t(JSN_score[,c(4,2)]),alternative = "decreasing")$p.value,CochranArmitageTest(t(JSN_score[,c(4,3)]),alternative = "decreasing")$p.value),method = "fdr",n=3)
##
CochranArmitageTest(t(age_table[,c(3,1)]),alternative = "increasing")   # Z = 0.3058, dim = 4, p-value = 0.6201
CochranArmitageTest(t(age_table[,c(3,2)]),alternative = "increasing")   # Z = -2, p-value = 0.02275 
CochranArmitageTest(t(age_table[,c(3,4)]),alternative = "increasing")   # Z = -0.78979, p-value = 0.2148
p.adjust(c(CochranArmitageTest(t(age_table[,c(3,1)]),alternative = "increasing")$p.value,CochranArmitageTest(t(age_table[,c(3,2)]),alternative = "increasing")$p.value,CochranArmitageTest(t(age_table[,c(3,4)]),alternative = "increasing")$p.value),method = "fdr",n=3)
##

### end caculate the clinical score in each cluster ###########
##***********



######## plot maker genes heatmap

row_Data <- rowData(sce.final)
rowData(sce.final)$feature_symbol <- rownames(sce.final)
dataset <- get_processed_dataset(sce.final)
hc <- metadata(sce.final)$sc3$consensus[[as.character(4)]]$hc

cart_mk4_c1 <- row_Data[which(row_Data$sc3_4_markers_clusts == 1 & row_Data$sc3_4_markers_padj < 0.05 ),]
cart_mk4_c1 <- cart_mk4_c1[order(cart_mk4_c1$sc3_4_markers_auroc,decreasing = TRUE),]
dim(cart_mk4_c1)  # [1] 211  30  # auroc  211   0.726666666666667
cart_mk4_c2 <- row_Data[which(row_Data$sc3_4_markers_clusts == 2 & row_Data$sc3_4_markers_padj < 0.05 ),]
cart_mk4_c2 <- cart_mk4_c2[order(cart_mk4_c2$sc3_4_markers_auroc,decreasing = TRUE),]
dim(cart_mk4_c2)  # [1] 37 30   #   0.785046728971963
cart_mk4_c3 <- row_Data[which(row_Data$sc3_4_markers_clusts == 3 & row_Data$sc3_4_markers_padj < 0.05 ),]
cart_mk4_c3 <- cart_mk4_c3[order(cart_mk4_c3$sc3_4_markers_auroc,decreasing = TRUE),]
dim(cart_mk4_c3)  # [1] 86 30  #  0.901652892561983
cart_mk4_c4 <- row_Data[which(row_Data$sc3_4_markers_clusts == 4 & row_Data$sc3_4_markers_padj < 0.05 ),]
cart_mk4_c4 <- cart_mk4_c4[order(cart_mk4_c4$sc3_4_markers_auroc,decreasing = TRUE),]
dim(cart_mk4_c4)  # [1] 230  30    # 0.836413043478261

mkGene.name <- c()
mkGene.name <- c(mkGene.name,cart_mk4_c1$feature_symbol[which(cart_mk4_c1$feature_symbol %in% PC_name)][1:20])
mkGene.name <- c(mkGene.name,cart_mk4_c2$feature_symbol[which(cart_mk4_c2$feature_symbol %in% PC_name)][1:20])
mkGene.name <- c(mkGene.name,cart_mk4_c3$feature_symbol[which(cart_mk4_c3$feature_symbol %in% PC_name)][1:20])
mkGene.name <- c(mkGene.name,cart_mk4_c4$feature_symbol[which(cart_mk4_c4$feature_symbol %in% PC_name)][1:20])
dataset <- get_processed_dataset(sce.final)
hc <- metadata(sce.final)$sc3$consensus[[as.character(4)]]$hc
color.palette =c(colorRampPalette(c("#265091", "#91cade"))(30),colorRampPalette(c("#91cade", "#ffffff", "#f69273"))(10),colorRampPalette(c("#f69273", "#8e162e"))(30))
#color.palette =colorRampPalette(c("#8e162e","white", "#265091"))(20)
pheatmap::pheatmap(dataset[mkGene.name,], cluster_cols = hc, cutree_cols = 4, cluster_rows = F,scale = "row",
                   show_rownames = T,show_colnames = F,color = color.palette,border_color = "grey60",cex=0.8 )


#######***********  end maker genes


col.C1=rgb(31, 149, 203, max = 255)
col.C2=rgb(255, 233, 0, max = 255)
col.C4=rgb(236, 95, 92, max = 255)

SF_ann <- read.table("./hydrathrosis/SF_sample_information_150917.txt",sep ="\t",header = T)
SF_ann <- SF_ann[rownames(SF_ann)[!duplicated(rownames(SF_ann))],1:2]
SF_ann$sampleID <- paste0("S",SF_ann$sampleID)
SF_data <- read.table("./hydrathrosis/hydrathrosis_pgml.txt")
SF_data.scaled <- scale(SF_data)
SF_data <- SF_data[apply(SF_data.scaled,1,sum)/ncol(SF_data) < 1,]   # remove outlier

SF_data.filter <- SF_data[,c("MIP.1alpha","SDF.1alpha","RANTES","IL.8","IFN.gamma", "VEGF.A","IL.6","IL.10")]
dim(SF_data.filter)
SF_data_boxplot <- c()
for (i in c(1,2,4)) {
  SF_data_c01<- SF_data.filter[rownames(SF_data.filter) %in% SF_ann[SF_ann$sampleID  %in% rownames(col_data[col_data[,"sc3_4_clusters"]==i,]),]$SFID,]
  condition <- rep(paste0("C",i),nrow(SF_data_c01))
  SF_data_c01 <- cbind(SF_data_c01,condition )
  dim(SF_data_c01)
  SF_data_boxplot <- rbind(SF_data_boxplot,SF_data_c01)
}
pdf("cartilage.hydrathrosis_4.pdf",7,5)
par(mfrow = c(2, 4),mar = rep(2, 4))
for(i in 1:8){
  SF_data_boxplot <- SF_data_boxplot[which(SF_data_boxplot$condition != "C3"),]
  boxplot(SF_data_boxplot[,i] ~ SF_data_boxplot$condition,main=colnames(SF_data_boxplot)[i],cex=0.5,cex.axis=1.3,cex.main=1.5,col=c(col.C1,col.C2,col.C4))
}
dev.off()

p.matrix <- c()
p_adjust.matrix <- c()
for(i in 1:8){
  p_C1C2 <- c()
  p_C1C4 <- c()
  p_C2C4 <- c()
  
  p_C1C2 <-wilcox.test(SF_data_boxplot[which(SF_data_boxplot$condition == "C1"),i],SF_data_boxplot[which(SF_data_boxplot$condition == "C2"),i],exact=FALSE)$p.value
  p_C1C4 <-t.test(SF_data_boxplot[which(SF_data_boxplot$condition == "C1"),i],SF_data_boxplot[which(SF_data_boxplot$condition == "C4"),i])$p.value
  p_C2C4 <-t.test(SF_data_boxplot[which(SF_data_boxplot$condition == "C2"),i],SF_data_boxplot[which(SF_data_boxplot$condition == "C4"),i])$p.value
  
  p.matrix <-  cbind(p.matrix,c(p_C1C2,p_C1C4,p_C2C4))    
  p_adjust <- p.adjust(c(p_C1C2,p_C1C4,p_C2C4),method = "fdr",n=3)
  p_adjust.matrix <- cbind(p_adjust.matrix,p_adjust)
}
colnames(p.matrix) <- colnames(SF_data_boxplot)[1:8]
colnames(p_adjust.matrix) <- colnames(SF_data_boxplot)[1:8]
write.table(p.matrix,"hydrathrosis.p.matrix",sep = "\t")
write.table(p_adjust.matrix,"hydrathrosis.p_adjust.matrix",sep = "\t")

##################  end hydrathrosis ############


##################

temp <- as.data.frame(annotation_col$cluster)
rownames(temp)<-rownames(annotation_col)
colnames(temp)<-c("subgroup")

goGene.names <- go.sets.hs.sym$`GO:0045202 synapse`[go.sets.hs.sym$`GO:0045202 synapse` %in% rownames(dataset)]
outfile = paste0("synapse.k",k,".pdf")
pheatmap(t(scale(t(dataset),scale = F))[goGene.names,],cluster_cols = hc, cutree_cols = k, cluster_rows = T,show_rownames = T,show_colnames = F,
         border_color = "grey60", main="genes in GO synapse term",fontsize_main = 15,fontsize_row = 6, treeheight_col=F,annotation_col = temp,file=outfile)

names(go.sets.hs.sym)[grepl("angiogenic",names(go.sets.hs.sym),perl = T)]
names(go.sets.hs.sym)[grepl("angiogen",names(go.sets.hs.sym),perl = T)]
goGene.names <- go.sets.hs.sym$`GO:0045202 synapse`[go.sets.hs.sym$`GO:0045202 angiogenic` %in% rownames(dataset)]

goGene.names <- go.sets.hs.sym$`GO:0001525 angiogenesis`[go.sets.hs.sym$`GO:0001525 angiogenesis` %in% rownames(dataset)]
length(goGene.names)
outfile = paste0("angiogenesis.k",k,".pdf")
pheatmap(t(scale(t(dataset),scale = F))[goGene.names,],cluster_cols = hc, cutree_cols = k, cluster_rows = T,show_rownames = T,show_colnames = F,
         border_color = "grey60", main="genes in GO angiogenesis term",fontsize_main = 15,fontsize_row = 6, treeheight_col=F,annotation_col = temp,file=outfile)


