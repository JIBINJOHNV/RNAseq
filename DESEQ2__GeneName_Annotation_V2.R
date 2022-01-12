set.seed("1234")

suppressMessages(library(argparse))
suppressMessages(library(DESeq2))
suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library(affy))
suppressMessages(library(RColorBrewer))
suppressMessages(library(genefilter ))
suppressMessages(library(dplyr))
suppressMessages(library("stringr"))    # Load stringr
library(dplyr)

parser <- ArgumentParser()

parser$add_argument("--ReadCount",
    help="Provide read count file, it should only contain feature name and counts")

parser$add_argument("--SampleInfo", 
    help="Provide file with sample info; headers [id,condition]; Should be tab separated")

parser$add_argument("--contrast", 
     help="Provide file with contrast info; headers [Reference,Contrast]; Should be tab separated")

parser$add_argument("--OutpuName", 
     help="Provide OutputFile/Folder Prefix")

parser$add_argument("--GeneDetails", 
     help="Gene details to be incorporated in the result/figures; like GeneName,chr, start, end, etc")



# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

# print some progress messages to stderr if "quietly" wasn't requested
if ( length(args$ReadCount)==0 ) { 
    write("Please provide required three files..\n", stderr()) 
    quit()
}

#Loading all files
countdata <-read.delim(args$ReadCount, header = T, row.names = NULL, check.names = F)
metadata <- read.delim(args$SampleInfo, header = T,sep = ",", row.names = NULL)
contrast <- read_delim(args$contrast,delim = ",",show_col_types = FALSE)
metadata2 <- read_delim(args$SampleInfo,delim = ",",show_col_types = FALSE)
GeneDetails<-read_delim(args$GeneDetails,delim = ",",show_col_types = FALSE)
OutputFilename=args$OutpuName


ifelse(!dir.exists("DES"), dir.create("DES"), "Folder exists already")
ifelse(!dir.exists("DES/Datat"), dir.create("DES/Data"), "Folder exists already")
ResultDir='DES'

#Select the row contain only genes
GeneDetails <- GeneDetails[str_detect(GeneDetails$type, "gene"), ]  # Extract matching rows with str_detect

#-------------------------------------Preprocess count data and Sample info file--------------------------------------------------------------------
rownames(countdata)<-countdata$Geneid ;countdata <-  subset(countdata, select=-c(Geneid))

#Adding One to all the samples read couns
#countdata <- countdata + 1
#Remove the genes with total read count <10
countdata <- countdata[rowSums(countdata) >=  10,]
countdata <- data.matrix(countdata)


#Prepare meta data
rownames(metadata)<-metadata$id 


# Checking the sample names in the ReadCount file and SampleInfo files are same or not; if not it will stop 
if ( all(colnames(countdata) != rownames(metadata)) ) { 
    write("\n\nSample names in ReadCount file and SampleInfo file is not same or not in the same order; please correct..\n\n", stderr()) 
    quit()
}
#checking column named condition in the metadata file
if("condition" %in% colnames(metadata2)){ "Column named condition is present in the provided SampleInfo file "
}else{
     print("COlumn named with condition is not in in SampleInfo file provided, may column name is wrong", stderr())
     quit()
     }

##Need to introdue metdata and and cntrast checkin

if(all(c(as.character(contrast$Reference), as.character(contrast$Contrast)) %in% metadata2$condition)){ 
    "Reference and contrast present in the condition column of SampleInfo file  provided"
}else{
     print("contrast file names not present in the condition column of SampleInfo file  provided", stderr())
     quit()
     }




#CREATE DESEQ OBJECTS
dds <- DESeqDataSetFromMatrix(countData=countdata,
                              colData=metadata,
                              design=~condition)




###Preparing separate normalised read for savin to fle and 

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)



### Transform counts for data visualization: https://hbctraining.github.io/DGE_workshop_salmon/lessons/03_DGE_QC_analysis.html; 
        #Transform normalized counts using the rlog transformation

#vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=TRUE)  ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2


normalized_countsN <- cbind(rownames(normalized_counts), data.frame(normalized_counts, row.names=NULL))
names(normalized_countsN)[names(normalized_countsN)=="rownames(normalized_counts)"] <- "Gene_ID"
write.table(normalized_countsN,file.path(ResultDir,paste0(OutputFilename,"_","ReadCountNormalised.xls")),sep = "\t",row.names = F) 
normalized_countsNWithGene<-merge(normalized_countsN, GeneDetails,by.x ="Gene_ID",by.y="gene_id",all.x = TRUE)
write.table(normalized_countsNWithGene,file.path(ResultDir,paste0(OutputFilename,"_","ReadCountNormalised_WithGeneDetails.xls")),sep = "\t",row.names = F) 



RlogNCount<- cbind(rownames(assay(rld)), data.frame(assay(rld), row.names=NULL))
colnames(RlogNCount)[which(names(RlogNCount) == "rownames(assay(rld))")] <- "Gene_ID"
write.table(RlogNCount,file.path(ResultDir,paste0(OutputFilename,"_","ReadCountRlogNorm.xls")),sep = "\t",row.names = F) 
RlogNCountWithGene<-merge(RlogNCount, GeneDetails,by.x ="Gene_ID",by.y="gene_id",all.x = TRUE)
write.table(RlogNCountWithGene,file.path(ResultDir,paste0(OutputFilename,"_","ReadCountRlogNorm_WithGeneDetails.xls")),sep = "\t",row.names = F) 


###----------------------------------------------------sample Level QC------Visualisation------------------------------------------------------------------------------

###Plot PCA
tiff(file=file.path(ResultDir,paste0(OutputFilename,"_PCA.tiff")), units="in", width=12, height=12, res=800)

pcaData <- plotPCA(rld, intgroup=c("condition","id"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=id, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

dev.off()


#Hierarchical Clustering (https://hbctraining.github.io/DGE_workshop_salmon/lessons/03_DGE_QC_analysis.html)

tiff(file=file.path(ResultDir,paste0(OutputFilename,"_H_Clustering.tiff")), units="in", width=8, height=8, res=800)

#Extract the rlog matrix from the object
rld_mat <- assay(rld)    ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2

### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function

### Plot heatmap
pheatmap(rld_cor)

dev.off()

## Plot dispersion estimates
dds <- estimateDispersions(dds)
tiff(file=file.path(ResultDir,paste0(OutputFilename,"_DispersionPlot.tiff")), units="in", width=8, height=8, res=300)
plotDispEsts(dds)
dev.off()


#Density Plots
tiff(file=file.path(ResultDir,paste0(OutputFilename,"_Density.tiff")), units="in", width=12, height=12, res=800)

plotDensity(log2(normalized_counts), 
            xlab="log2(normalized counts)", cex.lab=0.7, panel.first=grid()) 
legend("topright", legend=c(metadata2$id), lwd=2)

dev.off()



## Boxplots
tiff(file=file.path(ResultDir,paste0(OutputFilename,"_Gene_Count_DistributionsBoxplots.tiff")), units="in", width=12, height=12, res=800)
boxplot(log2(normalized_counts), pch=".", 
        horizontal=TRUE, cex.axis=0.5,
        las=1, ylab="Samples", xlab="log2(normalized Counts)")
dev.off()


#Histograms of counts per gene
tiff(file=file.path(ResultDir,paste0(OutputFilename,"_Count_PerGeneHist.tiff")), units="in", width=12, height=12, res=800)

hist(as.matrix(log2(normalized_counts)), breaks=100, col="blue", border="white",
     main="Log2-transformed normalized counts per gene", xlab="log2(counts+1)", ylab="Number of genes", 
     las=1, cex.axis=0.7)

dev.off()

#############--------------------------------------------------####################################################################
GdetailsHeatmap<-GeneDetails[c("gene_name","gene_id")]
rownames(GdetailsHeatmap)<-GdetailsHeatmap$gene_id


#Heatmap 25 variable genes
tiff(file=file.path(ResultDir,paste0(OutputFilename,"_Heatmap_Top25_VariableGene.tiff")), units="in", width=12, height=12, res=800)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 25 )
annotation_col = select(metadata,condition)

Merged<-merge(assay(rld)[ topVarGenes, ],GdetailsHeatmap,by=0,all.x = TRUE, no.dups = TRUE)
Merged$gene_name <- ifelse(is.na(Merged$gene_name), Merged$Row.names, Merged$gene_name)
rownames(Merged)<-Merged$gene_name
Merged = subset(Merged, select = -c(gene_id,Row.names,gene_name) )
pheatmap(Merged, scale="row",annotation_col = annotation_col)

#pheatmap( assay(rld)[ topVarGenes, ], scale="row",annotation_col = annotation_col)

dev.off()



#Heatmap 50 variable genes
tiff(file=file.path(ResultDir,paste0(OutputFilename,"_Heatmap_Top50_VariableGene.tiff")), units="in", width=12, height=12, res=800)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 50 )
annotation_col = select(metadata,condition)

Merged<-merge(assay(rld)[ topVarGenes, ],GdetailsHeatmap,by=0,all.x = TRUE, no.dups = TRUE)
Merged$gene_name <- ifelse(is.na(Merged$gene_name), Merged$Row.names, Merged$gene_name)
rownames(Merged)<-Merged$gene_name
Merged = subset(Merged, select = -c(gene_id,Row.names,gene_name) )
pheatmap(Merged, scale="row",annotation_col = annotation_col)

#pheatmap( assay(rld)[ topVarGenes, ], scale="row",annotation_col = annotation_col)

dev.off()



#Heatmap 100 variable genes
tiff(file=file.path(ResultDir,paste0(OutputFilename,"_Heatmap_Top100_VariableGene.tiff")), units="in", width=12, height=12, res=800)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 100 )
annotation_col = select(metadata,condition)

Merged<-merge(assay(rld)[ topVarGenes, ],GdetailsHeatmap,by=0,all.x = TRUE, no.dups = TRUE)
Merged$gene_name <- ifelse(is.na(Merged$gene_name), Merged$Row.names, Merged$gene_name)
rownames(Merged)<-Merged$gene_name
Merged = subset(Merged, select = -c(gene_id,Row.names,gene_name) )
pheatmap(Merged, scale="row",annotation_col = annotation_col,show_rownames=FALSE)

#pheatmap( assay(rld)[ topVarGenes, ], scale="row",annotation_col = annotation_col,show_rownames=FALSE,)

dev.off()



#Heatmap 500 variable genes
tiff(file=file.path(ResultDir,paste0(OutputFilename,"_Heatmap_Top500_VariableGene.tiff")), units="in", width=12, height=12, res=800)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 500 )
annotation_col = select(metadata,condition)

Merged<-merge(assay(rld)[ topVarGenes, ],GdetailsHeatmap,by=0,all.x = TRUE, no.dups = TRUE)
Merged$gene_name <- ifelse(is.na(Merged$gene_name), Merged$Row.names, Merged$gene_name)
rownames(Merged)<-Merged$gene_name
Merged = subset(Merged, select = -c(gene_id,Row.names,gene_name) )
pheatmap(Merged, scale="row",annotation_col = annotation_col,show_rownames=FALSE)

#pheatmap( assay(rld)[ topVarGenes, ], scale="row",annotation_col = annotation_col,show_rownames=FALSE,)

dev.off()




##------------------------------------------------------Differential Expression analysis-----------------------------------------------------------------------------

dds <- DESeq(dds)


if (all(contrast[1]==contrast['Reference'])) {
}else{
    print("First colun of Contrast file not contain Reference, or column name not correct",stderr())
    quit()
}



#Let’s get output for normal tissue vs primary tumor expression results and view a summary of results. contrast=c("tissueType", "primary colorectal cancer", "normal-looking surrounding colonic epithelium"))
#https://genviz.org/module-04-expression/0004/02/01/DifferentialExpression/


for (i in 1:nrow(contrast)) {
CONDITION="condition"
REFERENCE=toString(contrast[i,1])
TREATED=toString(contrast[i,2])

Name=paste(paste(CONDITION,REFERENCE,sep = "_"),TREATED,sep = "_vs_")

#res <- results(dds, name=Name)
res <- results(dds, contrast=c(CONDITION,TREATED,REFERENCE))


#EXTRACT THE NORMALISED READ COUNTS
Columns=row.names(subset(metadata,  subset=(condition==TREATED|condition==REFERENCE)))
RequiredNCounts=subset(normalized_counts,  select=Columns)
colnames(RequiredNCounts) <- paste("NormalisedRcount", colnames(RequiredNCounts), sep = "_")

res$diffexpressed <- "No_change"
res$diffexpressed[res$log2FoldChange >= 1 & res$padj < 0.05] <- "Up"
res$diffexpressed[res$log2FoldChange <= -1 & res$padj < 0.05] <- "Down"
res$diffexpressed[res$log2FoldChange >= 1 & res$padj >= 0.05] <- "Non_significant"
res$diffexpressed[res$log2FoldChange <= -1 & res$padj >= 0.05] <- "Non_significant"
res_data_NormReadCount<-merge(res,RequiredNCounts, by=0,all = TRUE)
colnames(res_data_NormReadCount)[colnames(res_data_NormReadCount) == 'Row.names'] <- 'Geneid'



ifelse(!dir.exists(paste0("DES/Data/",Name)), dir.create(paste0("DES/Data/",Name)), "Folder exists already")
ResultDir='DES/Data'

write.csv(res_data_NormReadCount,file.path(paste0("DES/Data/",Name),paste0(OutputFilename,"_",Name,"_","DESeq2_with_NormReadCount.csv"))) 
write.csv(res,file.path(paste0("DES/Data/",Name),paste0(OutputFilename,"_",Name,"_","DESeq2.csv"))) 




res_data_NormReadCountWithGene<-merge(res_data_NormReadCount, GeneDetails,by.x ="Geneid",by.y="gene_id",all.x = TRUE)
write.csv(res_data_NormReadCountWithGene,file.path(paste0("DES/Data/",Name),paste0(OutputFilename,"_",Name,"_","DESeq2with_NormReadCount_WithGeneDetails.csv"))) 



row.names(GeneDetails)<-GeneDetails$gene_id
resWithGene<-merge(data.frame(res), GeneDetails,by=0,all.x = TRUE)
resWithGene = subset(resWithGene, select = -c(gene_id) )
colnames(resWithGene)[which(names(resWithGene) == "Row.names")] <- "gene_id"

write.csv(resWithGene,file.path(paste0("DES/Data/",Name),paste0(OutputFilename,"_",Name,"_","DESeq2_WithGeneDetails.csv"))) 




#MA-plot

resLFC <- lfcShrink(dds, contrast=c(CONDITION,TREATED,REFERENCE), type="ashr")
mycolors <- c("blue", "red", "green", "black")
names(mycolors) <- c("Up", "Down", "Non_significant", "No_change")

res=as.data.frame(resLFC)
res$diffexpressed <- "No_change"
res$diffexpressed[res$log2FoldChange >= 1 & res$padj < 0.05] <- "Up"
res$diffexpressed[res$log2FoldChange <= -1 & res$padj < 0.05] <- "Down"
res$diffexpressed[res$log2FoldChange >= 1 & res$padj >= 0.05] <- "Non_significant"
res$diffexpressed[res$log2FoldChange <= -1 & res$padj >= 0.05] <- "Non_significant"


tiff(file=file.path(paste0("DES/Data/",Name),paste0(OutputFilename,"_",Name,"_","LFC_MAPlot_ggplot.tiff")), units="in", width=8, height=8, res=300)
    p<- ggplot(data=as.data.frame(res), aes(x=log10(baseMean), y=log2FoldChange, 
    col=diffexpressed)) + geom_point(size=0.5) + theme_minimal() + geom_hline(yintercept=c(-1,1), col="red") + scale_colour_manual(values = mycolors)
    print(p)

dev.off()



###Volcano Plot
tiff(file=file.path(paste0("DES/Data/",Name),paste0(OutputFilename,"_",Name,"_","LFC_Volcano_ggplot.tiff")), units="in", width=8, height=8, res=300)
    p<- ggplot(data=as.data.frame(res), aes(x=log2FoldChange, y=-log10(pvalue), 
    col=diffexpressed)) + geom_point(size=0.5) + theme_minimal() +  geom_vline(xintercept=c(-1, 1), 
    col="red") + geom_hline(yintercept=-log10(0.05), col="red") + scale_colour_manual(values = mycolors)
    print(p)

dev.off()



### Heatmap Top 25 DE genes
res2<-res[order(res$padj),]
Significant=head(res[order(res$padj),],25)

#Extract the desired column annotation columns
Samples<-metadata[str_detect(metadata$condition, REFERENCE), ]$id
Sam<-metadata[str_detect(metadata$condition, TREATED), ]$id
Samples<-c(Samples,Sam)

annotation_col=select(rbind( metadata[str_detect(metadata$condition, REFERENCE), ],metadata[str_detect(metadata$condition, TREATED), ]),condition)

tiff(file=file.path(paste0("DES/Data/",Name),paste0(OutputFilename,"_",Name,"_Heatmap_Top25_DEGene.tiff")), units="in", width=12, height=12, res=800)

Merged<-merge(assay(rld)[row.names(Significant),Samples],GdetailsHeatmap,by=0,all.x = TRUE, no.dups = TRUE)
Merged$gene_name <- ifelse(is.na(Merged$gene_name), Merged$Row.names, Merged$gene_name)
rownames(Merged)<-Merged$gene_name
Merged = subset(Merged, select = -c(gene_id,Row.names,gene_name) )
pheatmap(Merged, scale="row",annotation_col = annotation_col,show_rownames=TRUE)

#pheatmap(assay(rld)[row.names(Significant),Samples], scale="row",annotation_col = annotation_col)

dev.off()




### Heatmap Top 50 DE genes
res2<-res[order(res$padj),]
Significant=head(res[order(res$padj),],50)

#Extract the desired column annotation columns
Samples<-metadata[str_detect(metadata$condition, REFERENCE), ]$id
Sam<-metadata[str_detect(metadata$condition, TREATED), ]$id
Samples<-c(Samples,Sam)

annotation_col=select(rbind( metadata[str_detect(metadata$condition, REFERENCE), ],metadata[str_detect(metadata$condition, TREATED), ]),condition)

tiff(file=file.path(paste0("DES/Data/",Name),paste0(OutputFilename,"_",Name,"_Heatmap_Top50_DEGene.tiff")), units="in", width=12, height=12, res=800)

Merged<-merge(assay(rld)[row.names(Significant),Samples],GdetailsHeatmap,by=0,all.x = TRUE, no.dups = TRUE)
Merged$gene_name <- ifelse(is.na(Merged$gene_name), Merged$Row.names, Merged$gene_name)
rownames(Merged)<-Merged$gene_name
Merged = subset(Merged, select = -c(gene_id,Row.names,gene_name) )
pheatmap(Merged, scale="row",annotation_col = annotation_col,show_rownames=TRUE)

#pheatmap(assay(rld)[row.names(Significant),Samples], scale="row",annotation_col = annotation_col)

dev.off()


### Heatmap Top 100 DE genes
res2<-res[order(res$padj),]
Significant=head(res[order(res$padj),],100)

#Extract the desired column annotation columns
Samples<-metadata[str_detect(metadata$condition, REFERENCE), ]$id
Sam<-metadata[str_detect(metadata$condition, TREATED), ]$id
Samples<-c(Samples,Sam)

annotation_col=select(rbind( metadata[str_detect(metadata$condition, REFERENCE), ],metadata[str_detect(metadata$condition, TREATED), ]),condition)

tiff(file=file.path(paste0("DES/Data/",Name),paste0(OutputFilename,"_",Name,"_Heatmap_Top100_DEGene.tiff")), units="in", width=12, height=12, res=800)



Merged<-merge(assay(rld)[row.names(Significant),Samples],GdetailsHeatmap,by=0,all.x = TRUE, no.dups = TRUE)
Merged$gene_name <- ifelse(is.na(Merged$gene_name), Merged$Row.names, Merged$gene_name)
rownames(Merged)<-Merged$gene_name
Merged = subset(Merged, select = -c(gene_id,Row.names,gene_name) )
pheatmap(Merged, scale="row",annotation_col = annotation_col,show_rownames=FALSE)

#pheatmap(assay(rld)[row.names(Significant),Samples], scale="row",annotation_col = annotation_col,show_rownames=FALSE)

dev.off()



### Heatmap Top 500 DE genes
res2<-res[order(res$padj),]
Significant=head(res[order(res$padj),],500)

#Extract the desired column annotation columns
Samples<-metadata[str_detect(metadata$condition, REFERENCE), ]$id
Sam<-metadata[str_detect(metadata$condition, TREATED), ]$id
Samples<-c(Samples,Sam)

annotation_col=select(rbind( metadata[str_detect(metadata$condition, REFERENCE), ],metadata[str_detect(metadata$condition, TREATED), ]),condition)

tiff(file=file.path(paste0("DES/Data/",Name),paste0(OutputFilename,"_",Name,"_Heatmap_Top500_DEGene.tiff")), units="in", width=12, height=12, res=800)


Merged<-merge(assay(rld)[row.names(Significant),Samples],GdetailsHeatmap,by=0,all.x = TRUE, no.dups = TRUE)
Merged$gene_name <- ifelse(is.na(Merged$gene_name), Merged$Row.names, Merged$gene_name)
rownames(Merged)<-Merged$gene_name
Merged = subset(Merged, select = -c(gene_id,Row.names,gene_name) )
pheatmap(Merged, scale="row",annotation_col = annotation_col,show_rownames=FALSE)

#pheatmap(assay(rld)[row.names(Significant),Samples], scale="row",annotation_col = annotation_col,show_rownames=FALSE)

dev.off()


}


####Session information#######
sink("DES/DESeq2.session_info.txt")
print("seed is 1234")
sessionInfo()
sink()




