#CHECK MSAT GENE EXPRESSION AND TMM NORMALIZE
setwd("/home/nikolay/WABI/O_Hansson/MSAT/")
expr<-read.delim("MSAT_expr_0.000025.txt",header=TRUE,row.names=1,check.names=FALSE,sep="\t")
#expr<-read.delim("MSAT_expr.txt",header=TRUE,row.names=1,check.names=FALSE,sep="\t")
expr[1:5,1:5]
barplot(colSums(expr),las=2,ylab="Library Size")
expr<-expr[rowMeans(expr)>1,]

library("edgeR")
cds<-DGEList(expr)
plotMDS(cds,main="MSAT MDS Plot")

library("tweeDEseq")
cds$counts<-normalizeCounts(cds$counts)
barplot(colSums(cds$counts),las=2,ylab="Library Size", main = "After Normalization")
plotMDS(cds,main="MSAT MDS Plot After TMM Normalization")
write.table(cds$counts,file="MSAT_expr_0.000025_TMM_Normalized.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")


#DECONVOLUTION ANALYSIS
setwd("/home/nikolay/WABI/O_Hansson/scRNAseq/")

cluster0_markers<-c("ABCA5","TRIM63")
cluster1_markers<-c("EMC10","CTD.2545M3.8","TNNT3","ATP2A1","ACTA1","MYH2")
cluster2_markers<-c("MYBPC1","TNNT1","MYH7B","MYH7","ATP2A2")
cluster3_markers<-c("MYH7B","MYH7","ATP2A2")

final_list<-unique(c(cluster0_markers,cluster1_markers,cluster2_markers,cluster3_markers))

#READ SNRNASEQ DATA
library("data.table")
expr <- suppressWarnings(as.data.frame(fread("snRNAseq_Human_Chimp_SKM_SeuratFiltered.txt",sep="\t")))
rownames(expr)<-expr$V1; expr$V1<-NULL;
expr <- expr[rowSums(expr) != 0,]
expr[1:5,1:5]

#READ CELL ASSIGNMENTS
cell_assignment<-read.delim("snRNAseq_cell_assignment.txt",header=TRUE,sep="\t")
cell_assignment$FUNCTION<-"CHIMP_FAST_TWITCH"
cell_assignment$FUNCTION[cell_assignment$CLUSTER==1]<-"HUMAN_FAST_TWITCH"
cell_assignment$FUNCTION[cell_assignment$CLUSTER==2]<-"HUMAN_SLOW_TWITCH"
cell_assignment$FUNCTION[cell_assignment$CLUSTER==3]<-"CHIMP_SLOW_TWITCH"
cell_assignment$FUNCTION[cell_assignment$CLUSTER==4]<-"UNCLEAR_CELLS"
head(cell_assignment)
type<-as.character(cell_assignment$FUNCTION)

#BUILD AVERAGE GENE EXPRESSION
expr_chimp_fast_twitch<-subset(expr,select=colnames(expr)[type=="CHIMP_FAST_TWITCH"])
expr_chimp_fast_twitch[final_list,][1:5,1:5]
average_chimp_fast_twitch<-rowMeans(expr_chimp_fast_twitch[final_list,])
head(average_chimp_fast_twitch)
length(average_chimp_fast_twitch)

expr_human_fast_twitch<-subset(expr,select=colnames(expr)[type=="HUMAN_FAST_TWITCH"])
expr_human_fast_twitch[final_list,][1:5,1:5]
average_human_fast_twitch<-rowMeans(expr_human_fast_twitch[final_list,])
head(average_human_fast_twitch)
length(average_human_fast_twitch)

expr_human_slow_twitch<-subset(expr,select=colnames(expr)[type=="HUMAN_SLOW_TWITCH"])
expr_human_slow_twitch[final_list,][1:5,1:5]
average_human_slow_twitch<-rowMeans(expr_human_slow_twitch[final_list,])
head(average_human_slow_twitch)
length(average_human_slow_twitch)

expr_chimp_slow_twitch<-subset(expr,select=colnames(expr)[type=="CHIMP_SLOW_TWITCH"])
expr_chimp_slow_twitch[final_list,][1:5,1:5]
average_chimp_slow_twitch<-rowMeans(expr_chimp_slow_twitch[final_list,])
head(average_chimp_slow_twitch)
length(average_chimp_slow_twitch)

expr_unclear_cells<-subset(expr,select=colnames(expr)[type=="UNCLEAR_CELLS"])
expr_unclear_cells[final_list,][1:5,1:5]
average_unclear_cells<-rowMeans(expr_unclear_cells[final_list,])
head(average_unclear_cells)
length(average_unclear_cells)

#READ BULK RNASEQ DATA
path<-"/home/nikolay/WABI/O_Hansson/MSAT/"
bulk_expr<-read.delim(paste0(path,"MSAT_expr_0.000025_TMM_Normalized.txt"),header=TRUE,row.names=1,check.names=FALSE,sep="\t")
rownames(bulk_expr)<-substr(rownames(bulk_expr),1,15)
bulk_expr[1:5,1:5]

library("stringr")
annot<-scan("/home/nikolay/WABI/O_Hansson/Annotation/FULL_ANNOT.txt",what="character")
annot<-matrix(unlist(strsplit(annot[str_count(annot,"_")==2],"__")),ncol=2,byrow=TRUE)
head(annot)

intersect_ens<-intersect(rownames(bulk_expr),annot[,1])
length(intersect_ens)
bulk_expr<-bulk_expr[match(intersect_ens,rownames(bulk_expr)),]
annot<-annot[match(intersect_ens,annot[,1]),]
bulk_expr[1:5,1:5]
head(annot)
annot[,2]<-make.names(annot[,2],unique=TRUE)
rownames(bulk_expr)<-annot[,2]
bulk_expr[1:5,1:5]
bulk_expr_human<-bulk_expr
dim(bulk_expr_human)

#BUILD MATRIX OF SIGNATURES
signatures_human<-data.frame(human_fast_twitch=average_human_fast_twitch,human_slow_twitch=average_human_slow_twitch)
head(signatures_human)
dim(signatures_human)

#DECONVOLUTION HUMAN
library("DeconRNASeq")
signatures_human<-na.omit(signatures_human)
signatures_human
result_human<-DeconRNASeq(bulk_expr_human, signatures_human, use.scale = TRUE)
final_result_human<-result_human$out.all
rownames(final_result_human)<-colnames(bulk_expr_human)
final_result_human

#CORRELATION WITH FIBERTYPE DATA
MSAT_meta<-read.delim("/home/nikolay/WABI/O_Hansson/MSAT_samples_RNA_seq.csv",header=TRUE,sep="\t")
head(MSAT_meta)
MSAT_meta_muscle<-MSAT_meta[grepl("Muscle",as.character(MSAT_meta$RNA_Type)),]
head(MSAT_meta_muscle)
final_result_human<-final_result_human[match(as.character(MSAT_meta_muscle$ID_Nr_RNA_Seq),rownames(final_result_human)),]
head(final_result_human)

plot(as.numeric(final_result_human[,"human_slow_twitch"])~as.numeric(MSAT_meta_muscle$Type_1_Fibers_perc),
     xlab="TRUE TYPE1 FRACTION",ylab="PREDICTED TYPE1 FRACTION")
summary(lm(as.numeric(final_result_human[,"human_slow_twitch"])~as.numeric(MSAT_meta_muscle$Type_1_Fibers_perc)))
abline(lm(as.numeric(final_result_human[,"human_slow_twitch"])~as.numeric(MSAT_meta_muscle$Type_1_Fibers_perc)))
cor.test(as.numeric(final_result_human[,"human_slow_twitch"]),as.numeric(MSAT_meta_muscle$Type_1_Fibers_perc),
         method="spearman")
mtext("Spearman correlation: rho = 0.72, p-value = 2e-7")

setwd("/home/nikolay/WABI/O_Hansson/MSAT/")
df<-read.delim("MSAT_deconv_result.txt",header=TRUE,sep="\t")
spearmanCI::spearmanCI(as.numeric(df[,"human_slow_twitch"]),as.numeric(df$Type_1_Fibers_perc), level = 0.9, plot=TRUE)

MSAT_meta_muscle_final<-cbind(MSAT_meta_muscle,final_result_human)
MSAT_meta_muscle_final
write.table(MSAT_meta_muscle_final,file="MSAT_deconv_result.txt",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")





###################################### MEAN SQUAR ERROR PLOT ############################################
my_files<-c("MSAT_expr_TMM_Normalized.txt","MSAT_expr_0.2_TMM_Normalized.txt", "MSAT_expr_0.01_TMM_Normalized.txt", 
            "MSAT_expr_0.001_TMM_Normalized.txt", "MSAT_expr_0.0001_TMM_Normalized.txt","MSAT_expr_0.00005_TMM_Normalized.txt",
            "MSAT_expr_0.000025_TMM_Normalized.txt","MSAT_expr_0.00001_TMM_Normalized.txt")
downsample<-c(1,0.2,0.01,0.001,0.0001,0.00005,0.000025,0.00001)

spearman_rho<-vector()
mse<-vector()
seq_depth<-vector()
for(k in 1:length(downsample))
{
print(paste0("Working with downsample = ",downsample[k]))
path<-"/home/nikolay/WABI/O_Hansson/MSAT/"
#bulk_expr<-read.delim(paste0(path,"MSAT_expr_0.00001_TMM_Normalized.txt"),header=TRUE,row.names=1,check.names=FALSE,sep="\t")
bulk_expr<-read.delim(paste0(path,my_files[k]),header=TRUE,row.names=1,check.names=FALSE,sep="\t")
rownames(bulk_expr)<-substr(rownames(bulk_expr),1,15)
seq_depth<-append(seq_depth,mean(colSums(bulk_expr)))
bulk_expr[1:5,1:5]

library("stringr")
annot<-scan("/home/nikolay/WABI/O_Hansson/Annotation/FULL_ANNOT.txt",what="character")
annot<-matrix(unlist(strsplit(annot[str_count(annot,"_")==2],"__")),ncol=2,byrow=TRUE)
head(annot)

intersect_ens<-intersect(rownames(bulk_expr),annot[,1])
length(intersect_ens)
bulk_expr<-bulk_expr[match(intersect_ens,rownames(bulk_expr)),]
annot<-annot[match(intersect_ens,annot[,1]),]
bulk_expr[1:5,1:5]
head(annot)
annot[,2]<-make.names(annot[,2],unique=TRUE)
rownames(bulk_expr)<-annot[,2]
bulk_expr[1:5,1:5]
bulk_expr_human<-bulk_expr
dim(bulk_expr_human)

#BUILD MATRIX OF SIGNATURES
signatures_human<-data.frame(human_fast_twitch=average_human_fast_twitch,human_slow_twitch=average_human_slow_twitch)
head(signatures_human)
dim(signatures_human)

#DECONVOLUTION HUMAN
library("DeconRNASeq")
signatures_human<-na.omit(signatures_human)
signatures_human
result_human<-DeconRNASeq(bulk_expr_human, signatures_human, use.scale = TRUE)
final_result_human<-result_human$out.all
rownames(final_result_human)<-colnames(bulk_expr_human)
final_result_human

#CORRELATION WITH FIBERTYPE DATA
MSAT_meta<-read.delim("/home/nikolay/WABI/O_Hansson/MSAT_samples_RNA_seq.csv",header=TRUE,sep="\t")
head(MSAT_meta)
MSAT_meta_muscle<-MSAT_meta[grepl("Muscle",as.character(MSAT_meta$RNA_Type)),]
head(MSAT_meta_muscle)
final_result_human<-final_result_human[match(as.character(MSAT_meta_muscle$ID_Nr_RNA_Seq),rownames(final_result_human)),]
head(final_result_human)

spear_rho<-cor.test(as.numeric(final_result_human[,"human_slow_twitch"]),as.numeric(MSAT_meta_muscle$Type_1_Fibers_perc),
                    method="spearman")
spearman_rho<-append(spearman_rho,as.numeric(spear_rho$estimate))
MSAT_meta_muscle_final<-cbind(MSAT_meta_muscle,final_result_human)
my_mse<-sum((MSAT_meta_muscle_final$human_slow_twitch-MSAT_meta_muscle_final$Type_1_Fibers_perc)^2)/39
mse<-append(mse,my_mse)
}
seq_depth
spearman_rho
mse
plot(spearman_rho~log10(seq_depth),type="b",xlab="LOG10(SEQUENCING DEPTH)",ylab="SPEARMAN RHO",
     main="CORRELATION BETWEEN TRUE AND PREDICTED TYPE 1 FRACTION VS. SEQUENCING DEPTH", col="blue")
plot(mse~log10(seq_depth),type="b",xlab="LOG10(SEQUENCING DEPTH)",ylab="MEAN SQUARE ERROR",
     main="MEAN SQUARE DEVIATION OF PREDICTE TYPE 1 FRACTION FROM TRUE VS. SEQUENCING DEPTH", col="red")
  