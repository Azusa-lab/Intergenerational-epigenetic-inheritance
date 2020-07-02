# Data normalizaiton, Data visualizaiton in the paper

# hierarchical clustering
library("ape")
library("phytools")
rpkm_bin <- read.table("all_developmental_samples_rpkm.filter.matrix", header=T, sep="\t")
rpkm_bin_matrix <- rpkm_bin[,2:length(rpkm_bin[1,])]
rpkm_bin_matrix <- as.matrix(rpkm_bin_matrix)
hc <- hclust(as.dist(1-cor(rpkm_bin_matrix, method="pearson")),method="complete")
plot(as.phylo(hc),cex=0.5,label.offset=0.4)
hc_out <- as.phylo(hc)
write.tree(hc_out,"samples_hc_phylo_tree.out")

# association between gene density and H2AK119ub1 enrichment level
library(ggplot2)
uH2A_maternal <- read.table("uH2A_1_cell_genome1.out",sep="\t",header=F) # table with columns: uH2A enrichment level; gene counts per 1M
uH2A_paternal <- read.table("uH2A_1_cell_genome2.out",sep="\t",header=F)
maternal_med <- tapply(uH2A_maternal$V1,as.factor(uH2A_maternal$V2),median)
paternal_med <- tapply(uH2A_paternal$V1,as.factor(uH2A_paternal$V2),median)
maternal_sd <- tapply(uH2A_maternal$V1,as.factor(uH2A_maternal$V2),sd)
paternal_sd <- tapply(uH2A_paternal$V1,as.factor(uH2A_paternal$V2),sd)
u_med_df <- data.frame(type=names(maternal_med),med=maternal_med)
u_sd_df <- data.frame(type=names(maternal_sd),med=maternal_sd)
k_med_df <- data.frame(type=names(paternal_med),med=paternal_med)
k_sd_df <- data.frame(type=names(paternal_sd),med=paternal_sd)
ggplot(u_med_df,aes(x=0:15)) + geom_smooth(aes(y=u_med_df$med),se=F,span=0.4,color="#c82f63") + geom_ribbon(aes(ymin=u_med_df$med - u_sd_df$med,ymax=u_med_df$med + u_sd_df$med),fill = "#f9e0ec", alpha=0.4)  + geom_smooth(aes(y=k_med_df$med),se=F,span=0.4,color="#728CA3") + geom_ribbon(aes(ymin=k_med_df$med - k_sd_df$med,ymax=k_med_df$med + k_sd_df$med),fill = "#94c0c2",alpha=0.1) + theme(panel.background=element_rect(color='black',fill="white",size=1.5),axis.text=element_text(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# kmeans clustering
temp <- read.table("uH2A_allelic_average_reads_value.out",sep="\t",header=T)
ipp <- temp[,c(2,3,4,5)]
ipp <- as.matrix(ipp)
rownames(ipp) <- temp$geneID
km <- kmeans(ipp,2)
m.kmeans <- cbind(ipp,km$cluster)
o <- order(m.kmeans[,5])
m.kmeans <- m.kmeans[o,]
write.table(m.kmeans,file="uH2A_allelic_average_reads_value_kmeans.out",sep="\t",quote=FALSE)

# genomic distribution of H2AK119ub1 peaks
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peak <- read.table("uH2A_merge_peaks.merge.filter.bed",sep="\t",header=FALSE)
gr <- GRanges(seqnames = peak$V1,IRanges(start=peak$V2,end=peak$V3))
peakAnno <- annotatePeak(gr,tssRegion=c(-2500,2500),TxDb=txdb)
dfAnno <- as.data.frame(peakAnno)
write.table(dfAnno,"chip_seeker_uH2A_annotation.out",sep="\t",quote=FALSE)


# DEGs identification by edgeR
FGO_count <- read.table("RNA_pcgf_ctr_ko_reads_raw_count_for_edgeR",header=T,sep="\t")
group <- factor(c("ctr","ctr","ko","ko"))
FGO_matrix <- FGO_count[,2:length(FGO_count[1,])]
FGO_matrix <- as.matrix(FGO_matrix)
rownames(FGO_matrix) <- FGO_count$gene_id
FGO_list <- DGEList(counts=FGO_matrix,group=group)
keep <- filterByExpr(FGO_list)
FGO_list <- FGO_list[keep,keep.lib.sizes=FALSE]
FGO_list <- calcNormFactors(FGO_list)
design <- model.matrix(~group)
FGO_list <- estimateDisp(FGO_list,design)
fit <- glmQLFit(FGO_list,design)
lrt <- glmLRT(fit,coef=2)
write.table(lrt$table,file="RNA_Pcgf_ctr_vs_ko_edgeR.out",sep="\t",quote=FALSE)
