
############## marker genes for GO_ALL   
cart_mk4_c4.ego <- enrichGO(gene = cart_mk4_c4$feature_symbol, OrgDb = org.Hs.eg.db, ont = "All",
                            keyType = "SYMBOL",pAdjustMethod = "BH", 
                            pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
                            minGSSize = 30, maxGSSize = 500 )
head(cart_mk4_c4.ego)
write.table(cart_mk4_c4.ego,file="GO/cart_mk4_c4.ego.txt",quote = F,sep = "\t") 

cart_mk4_c3.ego <- enrichGO(gene = cart_mk4_c3$feature_symbol, OrgDb = org.Hs.eg.db, ont = "All",
                            keyType = "SYMBOL",pAdjustMethod = "BH", 
                            pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
                            minGSSize = 30, maxGSSize = 500 )
head(cart_mk4_c3.ego)
write.table(cart_mk4_c3.ego,file="GO/cart_mk4_c3.ego.txt",quote = F,sep = "\t")

cart_mk4_c2.ego <- enrichGO(gene = cart_mk4_c2$feature_symbol, OrgDb = org.Hs.eg.db, ont = "All",
                            keyType = "SYMBOL",pAdjustMethod = "BH", 
                            pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
                            minGSSize = 30, maxGSSize = 500 )
head(cart_mk4_c2.ego)
write.table(cart_mk4_c2.ego,file="GO/cart_mk4_c2.ego.txt",quote = F,sep = "\t") 

cart_mk4_c1.ego <- enrichGO(gene = cart_mk4_c1$feature_symbol, OrgDb = org.Hs.eg.db, ont = "All",
                            keyType = "SYMBOL",pAdjustMethod = "BH", 
                            pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
                            minGSSize = 30, maxGSSize = 500 )
head(cart_mk4_c1.ego)
write.table(cart_mk4_c1.ego,file="GO/cart_mk4_c1.ego.txt",quote = F,sep = "\t") 


cart_mk4_c3.ego<-dropGO(cart_mk4_c3.ego,term=c("GO:0007157","GO:0098742"))
#cart_mk4_c3.ego.matrix<-cart_mk4_c3.ego.matrix[order(cart_mk4_c3.ego.matrix$qvalue),]

#p <- ggplot(cart_mk4_c3.ego.matrix,aes(x=Description,y=Count,colour=qvalue))
#p + geom_bar() + scale_colour_gradient(high='red3',low='dodgerblue4') 

#barplot(cart_mk4_c1.ego, drop=TRUE, showCategory=20,title = "GO enrichment terms in C1",color="qvalue",x="Count")
#barplot(cart_mk4_c2.ego, drop=TRUE, showCategory=20,title = "GO enrichment terms in C2",color="qvalue",x="Count")
barplot(cart_mk4_c3.ego, drop=TRUE, showCategory=15,title = "GO enrichment terms in C3",color="qvalue",x="GeneRatio")
#barplot(cart_mk4_c4.ego, drop=TRUE, showCategory=20,title = "GO enrichment terms in C4",color="qvalue",x="Count")

###***********
############## maker genes for GO_BP

cart_mk4_c4.ego.bp <- enrichGO(gene = cart_mk4_c4$feature_symbol, OrgDb = org.Hs.eg.db, ont = "BP",
                               keyType = "SYMBOL",pAdjustMethod = "BH", 
                               pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
                               minGSSize = 30, maxGSSize = 1000 )
head(cart_mk4_c4.ego.bp,n=20)
write.table(cart_mk4_c4.ego.bp,file="GO/cart_mk4_c4.ego.bp.txt",quote = F,sep = "\t") 
cart_mk4_c4.ego.bp.sp<- clusterProfiler::gofilter(cart_mk4_c4.ego.bp, level=5 )
head(clusterProfiler::gofilter(cart_mk4_c4.ego.bp, level=4),n=20)


cart_mk4_c3.ego.bp <- enrichGO(gene = cart_mk4_c3$feature_symbol, OrgDb = org.Hs.eg.db, ont = "BP",
                               keyType = "SYMBOL",pAdjustMethod = "BH", 
                               pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
                               minGSSize = 30, maxGSSize = 1000 )
head(cart_mk4_c3.ego.bp,n=15)
head(clusterProfiler::gofilter(cart_mk4_c3.ego.bp, level=6 ),n=20)
cart_mk4_c3.ego.bp.sp<-clusterProfiler::gofilter(cart_mk4_c3.ego.bp, level=6 )

cart_mk4_c2.ego.bp <- enrichGO(gene = cart_mk4_c2$feature_symbol, OrgDb = org.Hs.eg.db, ont = "BP",
                               keyType = "SYMBOL",pAdjustMethod = "BH", 
                               pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
                               minGSSize = 30, maxGSSize = 1000 )
head(cart_mk4_c2.ego.bp,n=20)
head(clusterProfiler::gofilter(cart_mk4_c2.ego.bp, level=5 ),n=20)
cart_mk4_c2.ego.bp.sp<-clusterProfiler::gofilter(cart_mk4_c2.ego.bp, level=5 )

cart_mk4_c1.ego.bp <- enrichGO(gene = cart_mk4_c1$feature_symbol, OrgDb = org.Hs.eg.db, ont = "BP",
                               keyType = "SYMBOL",pAdjustMethod = "BH", 
                               pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
                               minGSSize = 30, maxGSSize = 1000 )
head(cart_mk4_c1.ego.bp,n=15)
head(clusterProfiler::gofilter(cart_mk4_c1.ego.bp, level=6 ),n=20)
cart_mk4_c1.ego.bp.sp<-clusterProfiler::gofilter(cart_mk4_c1.ego.bp, level=6 )

cart_mk4_c2.ego.bp.sp<-dropGO(cart_mk4_c2.ego.bp.sp,term=c("GO:0071560","GO:0044259","GO:1900026"))
barplot(cart_mk4_c1.ego.bp.sp, drop=TRUE, showCategory=15,title = "GO enrichment terms in C1",color="qvalue",x="GeneRatio")
barplot(cart_mk4_c2.ego.bp.sp, drop=TRUE, showCategory=15,title = "GO enrichment terms in C2",color="qvalue",x="GeneRatio")
#barplot(cart_mk4_c3.ego.bp.sp, drop=TRUE, showCategory=15,title = "GO enrichment terms in C3",color="qvalue",x="GeneRatio")
barplot(cart_mk4_c4.ego.bp.sp, drop=TRUE, showCategory=15,title = "GO enrichment terms in C4",color="qvalue",x="GeneRatio")

##*******************
###################  cart   kegg    GO #################
data(go.sets.hs)
data(kegg.gs)
data(egSymb)
kegg.gs.sym<-lapply(kegg.gs, eg2sym)
go.sets.hs.sym<-lapply(go.sets.hs, eg2sym)

col_data <- colData(sce.final)
cart.c1.idx <- which(col_data$sc3_4_clusters==1)
cart.c1.other.idx <- which(col_data$sc3_4_clusters!=1)
cart.c1.kegg.t.p <- gage(exprs(sce.final), gsets = kegg.gs.sym, ref = cart.c1.other.idx, samp = cart.c1.idx , compare = "as.group", rank.test = T)
head(cart.c1.kegg.t.p$greater[, 1:5], 10)
head(cart.c1.kegg.t.p$less[, 1:5], 10)
cart.c1.go.t.p <- gage(exprs(sce.final), gsets = go.sets.hs.sym, ref = cart.c1.other.idx, samp = cart.c1.idx , compare = "as.group", rank.test = T)
head(cart.c1.go.t.p$greater[, 1:5], 10)
head(cart.c1.go.t.p$less[, 1:5], 10)

cart.c2.idx <- which(col_data$sc3_4_clusters==2)
cart.c2.other.idx <- which(col_data$sc3_4_clusters!=2)
cart.c2.kegg.t.p <- gage(exprs(sce.final), gsets = kegg.gs.sym, ref = cart.c2.other.idx, samp = cart.c2.idx , compare = "as.group", rank.test = T)
head(cart.c2.kegg.t.p$greater[, 1:5], 10)
head(cart.c2.kegg.t.p$less[, 1:5], 10)
cart.c2.go.t.p <- gage(exprs(sce.final), gsets = go.sets.hs.sym, ref = cart.c2.other.idx, samp = cart.c2.idx , compare = "as.group", rank.test = T)
head(cart.c2.go.t.p$greater[, 1:5], 10)
head(cart.c2.go.t.p$less[, 1:5], 10)

cart.c3.idx <- which(col_data$sc3_4_clusters==3)
cart.c3.other.idx <- which(col_data$sc3_4_clusters!=3)
cart.c3.kegg.t.p <- gage(exprs(sce.final), gsets = kegg.gs.sym, ref = cart.c3.other.idx, samp = cart.c3.idx , compare = "as.group", rank.test = T)
head(cart.c3.kegg.t.p$greater[, 1:5], 10)
head(cart.c3.kegg.t.p$less[, 1:5], 10)
cart.c3.go.t.p <- gage(exprs(sce.final), gsets = go.sets.hs.sym, ref = cart.c3.other.idx, samp = cart.c3.idx , compare = "as.group", rank.test = T)
head(cart.c3.go.t.p$greater[, 1:5], 10)
head(cart.c3.go.t.p$less[, 1:5], 10)

cart.c4.idx <- which(col_data$sc3_4_clusters==4)
cart.c4.other.idx <- which(col_data$sc3_4_clusters!=4)
cart.c4.kegg.t.p <- gage(exprs(sce.final), gsets = kegg.gs.sym, ref = cart.c4.other.idx, samp = cart.c4.idx , compare = "as.group", rank.test = T)
head(cart.c4.kegg.t.p$greater[, 1:5], 10)
head(cart.c4.kegg.t.p$less[, 1:5], 10)
cart.c4.go.t.p <- gage(exprs(sce.final), gsets = go.sets.hs.sym, ref = cart.c4.other.idx, samp = cart.c4.idx , compare = "as.group", rank.test = T)
head(cart.c4.go.t.p$greater[, 1:5], 10)
head(cart.c4.go.t.p$less[, 1:5], 10)

write.table(rbind(as.data.frame(cart.c1.kegg.t.p$greater)[which(as.data.frame(cart.c1.kegg.t.p$greater)$stat.mean > 0),1:5],
                  as.data.frame(cart.c1.kegg.t.p$less)[which(as.data.frame(cart.c1.kegg.t.p$less)$stat.mean < 0),1:5]),
            file ="GO/cart.c1.gage.kegg.xls",sep="\t")
write.table(rbind(as.data.frame(cart.c2.kegg.t.p$greater)[which(as.data.frame(cart.c2.kegg.t.p$greater)$stat.mean > 0),1:5],
                  as.data.frame(cart.c2.kegg.t.p$less)[which(as.data.frame(cart.c2.kegg.t.p$less)$stat.mean < 0),1:5]),
            file ="GO/cart.c2.gage.kegg.xls",sep="\t")
write.table(rbind(as.data.frame(cart.c3.kegg.t.p$greater)[which(as.data.frame(cart.c3.kegg.t.p$greater)$stat.mean > 0),1:5],
                  as.data.frame(cart.c3.kegg.t.p$less)[which(as.data.frame(cart.c3.kegg.t.p$less)$stat.mean < 0),1:5]),
            file ="GO/cart.c3.gage.kegg.xls",sep="\t")
write.table(rbind(as.data.frame(cart.c4.kegg.t.p$greater)[which(as.data.frame(cart.c4.kegg.t.p$greater)$stat.mean > 0),1:5],
                  as.data.frame(cart.c4.kegg.t.p$less)[which(as.data.frame(cart.c4.kegg.t.p$less)$stat.mean < 0),1:5]),
            file ="GO/cart.c4.gage.kegg.xls",sep="\t")

write.table(rbind(as.data.frame(cart.c1.go.t.p$greater)[which(as.data.frame(cart.c1.go.t.p$greater)$stat.mean > 0),1:5],
                  as.data.frame(cart.c1.go.t.p$less)[which(as.data.frame(cart.c1.go.t.p$less)$stat.mean < 0),1:5]),
            file ="GO/cart.c1.gage.go.txt",sep="\t")
write.table(rbind(as.data.frame(cart.c2.go.t.p$greater)[which(as.data.frame(cart.c2.go.t.p$greater)$stat.mean > 0),1:5],
                  as.data.frame(cart.c2.go.t.p$less)[which(as.data.frame(cart.c2.go.t.p$less)$stat.mean < 0),1:5]),
            file ="GO/cart.c2.gage.go.txt",sep="\t")
write.table(rbind(as.data.frame(cart.c3.go.t.p$greater)[which(as.data.frame(cart.c3.go.t.p$greater)$stat.mean > 0),1:5],
                  as.data.frame(cart.c3.go.t.p$less)[which(as.data.frame(cart.c3.go.t.p$less)$stat.mean < 0),1:5]),
            file ="GO/cart.c3.gage.go.txt",sep="\t")
write.table(rbind(as.data.frame(cart.c4.go.t.p$greater)[which(as.data.frame(cart.c4.go.t.p$greater)$stat.mean > 0),1:5],
                  as.data.frame(cart.c4.go.t.p$less)[which(as.data.frame(cart.c4.go.t.p$less)$stat.mean < 0),1:5]),
            file ="GO/cart.c4.gage.go.txt",sep="\t")

###################plot go in clusters 
#rownames(cart.c1.go.t.p$greater)[grepl("carti",rownames(cart.c1.go.t.p$greater))]

go.names <- c()
go.names <- c(go.names,rownames(as.data.frame(cart.c1.go.t.p$greater)[which(as.data.frame(cart.c1.go.t.p$greater)$q.val < 0.05),])[1:8])
go.names <- c(go.names,rownames(as.data.frame(cart.c2.go.t.p$greater)[which(as.data.frame(cart.c2.go.t.p$greater)$q.val < 0.05),])[1:8])
go.names <- c(go.names,rownames(as.data.frame(cart.c3.go.t.p$greater)[which(as.data.frame(cart.c3.go.t.p$greater)$q.val < 0.05),])[1:8])
go.names <- c(go.names,rownames(as.data.frame(cart.c4.go.t.p$greater)[which(as.data.frame(cart.c4.go.t.p$greater)$q.val < 0.05),])[1:8])
go.names <- go.names[!duplicated(go.names)]
go.names <- na.omit(go.names)

go.table <- as.data.frame(cart.c1.go.t.p$greater)[which(as.data.frame(cart.c1.go.t.p$greater)$q.val < 0.05),][ rownames(as.data.frame(cart.c1.go.t.p$greater)[which(as.data.frame(cart.c1.go.t.p$greater)$q.val < 0.05),]) %in% go.names,]
cluster <- rep("C1",dim(as.data.frame(cart.c1.go.t.p$greater)[which(as.data.frame(cart.c1.go.t.p$greater)$q.val < 0.05),][ rownames(as.data.frame(cart.c1.go.t.p$greater)[which(as.data.frame(cart.c1.go.t.p$greater)$q.val < 0.05),]) %in% go.names,])[1])
go.term <- rownames(go.table)
go.martix <- cbind(go.term,cluster,go.table )
go.table <- as.data.frame(cart.c2.go.t.p$greater)[which(as.data.frame(cart.c2.go.t.p$greater)$q.val < 0.05),][ rownames(as.data.frame(cart.c2.go.t.p$greater)[which(as.data.frame(cart.c2.go.t.p$greater)$q.val < 0.05),]) %in% go.names,]
cluster <- rep("C2",dim(as.data.frame(cart.c2.go.t.p$greater)[which(as.data.frame(cart.c2.go.t.p$greater)$q.val < 0.05),][ rownames(as.data.frame(cart.c2.go.t.p$greater)[which(as.data.frame(cart.c2.go.t.p$greater)$q.val < 0.05),]) %in% go.names,])[1])
go.term <- rownames(go.table)
go.martix <- rbind(go.martix,cbind(go.term,cluster,go.table))
go.table <- as.data.frame(cart.c3.go.t.p$greater)[which(as.data.frame(cart.c3.go.t.p$greater)$q.val < 0.05),][ rownames(as.data.frame(cart.c3.go.t.p$greater)[which(as.data.frame(cart.c3.go.t.p$greater)$q.val < 0.05),]) %in% go.names,]
cluster <- rep("C3",dim(as.data.frame(cart.c3.go.t.p$greater)[which(as.data.frame(cart.c3.go.t.p$greater)$q.val < 0.05),][ rownames(as.data.frame(cart.c3.go.t.p$greater)[which(as.data.frame(cart.c3.go.t.p$greater)$q.val < 0.05),]) %in% go.names,])[1])
go.term <- rownames(go.table)
go.martix <- rbind(go.martix,cbind(go.term,cluster,go.table))
go.table <- as.data.frame(cart.c4.go.t.p$greater)[which(as.data.frame(cart.c4.go.t.p$greater)$q.val < 0.05),][ rownames(as.data.frame(cart.c4.go.t.p$greater)[which(as.data.frame(cart.c4.go.t.p$greater)$q.val < 0.05),]) %in% go.names,]
cluster <- rep("C4",dim(as.data.frame(cart.c4.go.t.p$greater)[which(as.data.frame(cart.c4.go.t.p$greater)$q.val < 0.05),][ rownames(as.data.frame(cart.c4.go.t.p$greater)[which(as.data.frame(cart.c4.go.t.p$greater)$q.val < 0.05),]) %in% go.names,])[1])
go.term <- rownames(go.table)
go.martix <- rbind(go.martix,cbind(go.term,cluster,go.table))
go.martix$go.term <- sub("(GO:\\w+ )","",go.martix$go.term)
colnames(go.martix) <- c("GeneOntology","subgroup","p.geomean","stat.mean","p.val","FDR","set.size","exp1")

p <- ggplot(go.martix,aes(x=subgroup,y=GeneOntology,size=stat.mean,colour=-log10(FDR)))
p + geom_point() + scale_colour_gradient(high='red3',low='dodgerblue4') + theme()   #panel.background = element_rect(fill = NA)

## less
go.names <- c()
go.names <- c(go.names,rownames(as.data.frame(cart.c1.go.t.p$less)[which(as.data.frame(cart.c1.go.t.p$less)$q.val < 0.01),])[1:8])
go.names <- c(go.names,rownames(as.data.frame(cart.c2.go.t.p$less)[which(as.data.frame(cart.c2.go.t.p$less)$q.val < 0.01),])[1:8])
go.names <- c(go.names,rownames(as.data.frame(cart.c3.go.t.p$less)[which(as.data.frame(cart.c3.go.t.p$less)$q.val < 0.01),])[1:8])
go.names <- c(go.names,rownames(as.data.frame(cart.c4.go.t.p$less)[which(as.data.frame(cart.c4.go.t.p$less)$q.val < 0.01),])[1:8])
go.names <- go.names[!duplicated(go.names)]
go.names <- na.omit(go.names)

go.table <- as.data.frame(cart.c1.go.t.p$less)[which(as.data.frame(cart.c1.go.t.p$less)$q.val < 0.01),][ rownames(as.data.frame(cart.c1.go.t.p$less)[which(as.data.frame(cart.c1.go.t.p$less)$q.val < 0.01),]) %in% go.names,]
cluster <- rep("C1",dim(as.data.frame(cart.c1.go.t.p$less)[which(as.data.frame(cart.c1.go.t.p$less)$q.val < 0.01),][ rownames(as.data.frame(cart.c1.go.t.p$less)[which(as.data.frame(cart.c1.go.t.p$less)$q.val < 0.01),]) %in% go.names,])[1])
go.term <- rownames(go.table)
go.martix <- cbind(go.term,cluster,go.table )
go.table <- as.data.frame(cart.c2.go.t.p$less)[which(as.data.frame(cart.c2.go.t.p$less)$q.val < 0.01),][ rownames(as.data.frame(cart.c2.go.t.p$less)[which(as.data.frame(cart.c2.go.t.p$less)$q.val < 0.01),]) %in% go.names,]
cluster <- rep("C2",dim(as.data.frame(cart.c2.go.t.p$less)[which(as.data.frame(cart.c2.go.t.p$less)$q.val < 0.01),][ rownames(as.data.frame(cart.c2.go.t.p$less)[which(as.data.frame(cart.c2.go.t.p$less)$q.val < 0.01),]) %in% go.names,])[1])
go.term <- rownames(go.table)
go.martix <- rbind(go.martix,cbind(go.term,cluster,go.table))
go.table <- as.data.frame(cart.c3.go.t.p$less)[which(as.data.frame(cart.c3.go.t.p$less)$q.val < 0.01),][ rownames(as.data.frame(cart.c3.go.t.p$less)[which(as.data.frame(cart.c3.go.t.p$less)$q.val < 0.01),]) %in% go.names,]
cluster <- rep("C3",dim(as.data.frame(cart.c3.go.t.p$less)[which(as.data.frame(cart.c3.go.t.p$less)$q.val < 0.01),][ rownames(as.data.frame(cart.c3.go.t.p$less)[which(as.data.frame(cart.c3.go.t.p$less)$q.val < 0.01),]) %in% go.names,])[1])
go.term <- rownames(go.table)
go.martix <- rbind(go.martix,cbind(go.term,cluster,go.table))
go.table <- as.data.frame(cart.c4.go.t.p$less)[which(as.data.frame(cart.c4.go.t.p$less)$q.val < 0.01),][ rownames(as.data.frame(cart.c4.go.t.p$less)[which(as.data.frame(cart.c4.go.t.p$less)$q.val < 0.01),]) %in% go.names,]
cluster <- rep("C4",dim(as.data.frame(cart.c4.go.t.p$less)[which(as.data.frame(cart.c4.go.t.p$less)$q.val < 0.01),][ rownames(as.data.frame(cart.c4.go.t.p$less)[which(as.data.frame(cart.c4.go.t.p$less)$q.val < 0.01),]) %in% go.names,])[1])
go.term <- rownames(go.table)
go.martix <- rbind(go.martix,cbind(go.term,cluster,go.table))

go.martix$go.term <- sub("(GO:\\w+ )","",go.martix$go.term)
colnames(go.martix) <- c("GeneOntology","subgroup","p.geomean","stat.mean","p.val","FDR","set.size","exp1")

p <- ggplot(go.martix,aes(x=subgroup,y=GeneOntology,size=stat.mean,colour=-log10(FDR)))
p + geom_point() + scale_colour_gradient(high='red3',low='dodgerblue4') + theme()   #panel.background = element_rect(fill = NA)

###********************
