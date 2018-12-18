######################################################
#Vignettes used
#http://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html#10_getting_or_building_ensdb_databasespackages
#http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
#https://www.biostars.org/p/257180/
#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#input-data
######################################################

###Parameters###
# Import the libraries
library(readr)
library(ensembldb)
library(AnnotationHub)
library(tximport)
library(DESeq2)
library(biomaRt)
library(BiocParallel)
library(ashr)
library(pheatmap)
library(ggplot2)
library(vsn)
library(RColorBrewer)

# point to directory
dir <- "GEUVADIS_salmon"

# Load sample meta sheet we will read off of the meta sheet for groupings
samples <- read.table(file.path("sample_GEUVADIS_tximport_unique.txt"), header = TRUE)

# load the files from the directory based on the sample sheet names
files <- file.path(dir, samples$FILE_NAME)

#get the names of the files, set the names to each file
names(files) <- paste0(samples$RUN)

#check if they are present
all(file.exists(files))

######################################################
# Transcripts need to be associated with gene IDs for gene-level summarization. 
# If that information is present in the files, we can skip this step. 
# For Salmon, Sailfish, and kallisto the files only provide the transcript ID. 
# We first make a data.frame called tx2gene with two columns: 
# 1) transcript ID and 2) gene ID. 
# The column names do not matter but this column order must be used. 
# The transcript ID must be the same one used in the abundance files.
######################################################

###Creating the reference Tx2Gene library###

##option from sleuth##

# Use biomart to get the hsapien database
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = 'ensembl.org')

# create tx2gene db
tx2gene <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)

# rename the columns
tx2gene <- dplyr::rename(tx2gene, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

######################################################
#If inferential replicates (Gibbs or bootstrap samples) are present in expected 
#locations relative to the quant.sf file, tximport will import these as well, if 
#txOut=TRUE. tximport will not summarize inferential replicate information to the 
#gene-level. Here we demonstrate using Salmon, run with only 5 Gibbs replicates 
#(usually more Gibbs samples would be useful for estimating variability).
#We first tell salmon that we want to start with the import of transcript 
#level estimates. Downstream we can go to the gene level using the method described later
######################################################

###Tximport Salmon files with inferential replicates###

# create a db of salmon counts
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)

# show names of the imported data 
names(txi.salmon)

# Create the DESEQ2 object from the tximport df and the metadata
dds <- DESeqDataSetFromTximport(txi.salmon, 
								colData = samples,
								design = ~ POP)

#write this dds object as a file or save it in some way
saveRDS(dds, file ="GEUVADIS_POP_dds_FACTOR.rds")


######################################################
#Some advantages of using the above methods for transcript abundance estimation are: 
#(i) this approach corrects for potential changes in gene length across samples 
#(e.g. from differential isoform usage) (Trapnell et al. 2013), 
#(ii) some of these methods (Salmon, Sailfish, kallisto) are substantially faster 
#and require less memory and disk usage compared to alignment-based methods that 
#require creation and storage of BAM files, and (iii) it is possible to avoid 
#discarding those fragments that can align to multiple genes with homologous 
#sequence, thus increasing sensitivity (Robert and Watson 2015).	
######################################################
							
###DEseq2###

#While it is not necessary to pre-filter low count genes before running the DESeq2 
#functions, there are two reasons which make pre-filtering useful: by removing rows 
#in which there are very few reads, we reduce the memory size of the dds data object, 
#and we increase the speed of the transformation and testing functions within DESeq2. 
#Here we perform a minimal pre-filtering to keep only rows that have at least 10 reads 
#total. Note that more strict filtering to increase power is automatically applied via 
#independent filtering on the mean of normalized counts within the results function.
keep <- rowSums(counts(dds)) >= 46
dds <- dds[keep,]
                   
# Produce some results
dds <- DESeq(dds, parallel = TRUE)
saveRDS(dds, file ="GEUVADIS_FACTORS_dds.rds")
write.table(assay(dds), "GEUVADIS_FACTORS_dds.tsv", sep="\t")

res <- results(dds)
saveRDS(res, file = "GEUVADIS_res.rds")
write.csv(as.data.frame(res), 
          file="GEUVADIS_res.csv")
          
res05 <- results(dds, alpha=0.05)
saveRDS(res05, file ="GEUVADIS_POP_res05.rds")
write.csv(as.data.frame(res05), 
          file="GEUVADIS_DEseq2_POP_res05.csv")
          
resOrdered <- res[order(res$pvalue),]
saveRDS(resOrdered, file ="GEUVADIS_POP_resOrdered.rds")
write.csv(as.data.frame(resOrdered), 
          file="GEUVADIS_DEseq2_POP_resOrdered.csv")
          
resNorm <- lfcShrink(dds, coef=5, type="normal")
saveRDS(resNorm, file ="GEUVADIS_POP_resNorm.rds")
write.csv(as.data.frame(resNorm), 
          file="GEUVADIS_DEseq2_POP_resNorm.csv")
          
resAsh <- lfcShrink(dds, coef=5, type="ashr")
saveRDS(resASH, file ="GEUVADIS_POP_ASH.rds")
write.csv(as.data.frame(resAsh), 
          file="GEUVADIS_DEseq2_POP_resAsh.csv")
          
vsd <- vst(dds, blind=FALSE)
saveRDS(vsd, file ="GEVUADIS_vsd.rds")
write.csv(as.data.frame(vsd), 
          file="GEUVADIS_DEseq2_POP_vsd.csv")
          
ntd <- normTransform(dds)
saveRDS(ntd, file ="GEVUADIS_ntd.rds")
write.csv(as.data.frame(ntd), 
          file="GEUVADIS_DEseq2_POP_ntd.csv")
          
rld <- rlog(dds, blind=FALSE)
saveRDS(rld, file ="GEVUADIS_rld.rds")
write.csv(as.data.frame(rld), 
          file="GEUVADIS_DEseq2_POP_rld.csv")
          
resLFC <- lfcShrink(dds, coef=5, type="apeglm")
saveRDS(resLFC, file = "GEUVADIS_POP_resLFC.rds")
write.csv(as.data.frame(resLFC), 
          file="GEUVADIS_DEseq2_POP_resLFC.csv")
          
### Create the image that is going to show MAplot and normalization comparisons
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(res, ylim=ylim, xlim=xlim, main="unnormal")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")

### Create a plot of counts for the minimum differentially expressed gene between groups. 
plotCounts(dds, gene=which.min(res$padj), intgroup="CENTER")

###create a DF of counts based on the set up for use with ggplot
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="CENTER", 
                returnData=TRUE)

ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

# mean SD plots
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

# Create a heatmap to look at clustering based on the top counts
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:100]

df <- as.data.frame(colData(dds)[,c("CENTER","POP")])

pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=FALSE, show_colnames=FALSE,
         cluster_cols=TRUE, annotation_col=df)

# Heatmap of sample to sample distances     
sampleDists <- dist(t(assay(vsd)))
saveRDS(sampleDists, file ="GEVUADIS_sampleDists.rds")

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$CENTER, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup=c("CENTER", "POP"))

plot_group_density(vsd, use_filtered = TRUE, units = "est_counts",
  trans = "log", grouping = setdiff(colnames(so$sample_to_covariates),
  "sample"), offset = 1)