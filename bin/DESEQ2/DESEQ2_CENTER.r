######################################################
#Vignettes used
#http://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html#10_getting_or_building_ensdb_databasespackages
#http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
#https://www.biostars.org/p/257180/
#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#input-data
######################################################

###Parameters###

# install into personal repo if on cluster
install.packages("tidyverse")
install.packages("readr")

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

# point to directory
dir <- "salmon_output"

# Load sample meta sheet we will read off of the meta sheet for groupings
samples <- read.table(file.path("sample_GEUVADIS_tximport_unique.txt"), header = TRUE)

# load the files from the directory based on the sample sheet names
files <- file.path(dir, samples$RUN, "quant.sf")

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

##option from salmon/ensemble##

# open up the annotation hub module
#ah <- AnnotationHub()

# query available databases
#query(ah, "EnsDb")

# Query the species and the version
#edb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 89))
#edb

# add the database version to object
#edb <- edb[[1]]

# get the genes and transcripts from the database.... example of the use
#gns <- genes(edb)
#Tx <- transcripts(edb)

# get an object that is TX by Gns from the database
#TxByGns <- transcriptsBy(edb, by = "gene")

# Create the keys from the TX names in the edb database
#k <- keys(edb, keytype = "TXNAME")

# create the tx2gene database and select the Gene and TX columns
#tx2gene <- select(edb, k, "GENEID", "TXNAME")

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

#show the repnames
names(txi.salmon$infReps)

# Create the DESEQ2 object from the tximport df and the metadata
dds <- DESeqDataSetFromTximport(txi.salmon, 
								colData = samples,
								design = ~ CENTER)

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
res <- results(dds)
res 

plotMA(res, ylim=c(-2,2))
summary(res)


sum(res$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05)

#this is the summarize to gene command tha can be used to look at gene counts

write.csv(as.data.frame(resOrdered), 
          file="GEUVADIS_DEseq2_center.csv")

saveRDS(dds, file ="GEVUADIS_dds.rds")
saveRDS(res, file = "GEUVADIS_res.rds")

### RESULTS ANALYSIS ###
### GEVUADIS_dds is used because there I imported the dataframe from file
res05 <- results(GEVUADIS_dds, alpha=0.05)
summary(res05)

#Find the number of coefficients and then use them for the results analysis
resultsNames(GEVUADIS_dds)

### LFC in two methods on the data
resNorm <- lfcShrink(GEVUADIS_dds, coef=6, type="normal")
resAsh <- lfcShrink(GEVUADIS_dds, coef=6, type="ashr")

### Create the image that is going to show MAplot and normalization comparisons
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(res, ylim=ylim, xlim=xlim, main="unnormal")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")

### Create a plot of counts for the minimum differentially expressed gene between groups. 
plotCounts(GEVUADIS_dds, gene=which.min(res$padj), intgroup="CENTER")

###create a DF of counts based on the set up for use with ggplot
d <- plotCounts(GEVUADIS_dds, gene=which.min(res$padj), intgroup="CENTER", 
                returnData=TRUE)

ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

#normalizing counts and saving the dataframes
vsd <- vst(GEVUADIS_dds, blind=FALSE)
saveRDS(vsd, file ="GEVUADIS_vsd.rds")

#rld <- rlog(GEVUADIS_dds, blind=FALSE)
#saveRDS(rld, file ="GEVUADIS_rld.rds")

#just checking out the table
head(assay(vsd), 3)

# Create the SD plots to look at variance and also look at the different methods for normalization
# this gives log2(n + 1)
ntd <- normTransform(dds)

meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
#meanSdPlot(assay(rld))

# Create a heatmap to look at clustering based on the top counts
library("pheatmap")
select <- order(rowMeans(counts(GEVUADIS_dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]

df <- as.data.frame(colData(GEVUADIS_dds)[,c("CENTER","POP")])

pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=FALSE, show_colnames=FALSE,
         cluster_cols=TRUE, annotation_col=df)

# Heatmap of sample to sample distances     
sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$CENTER, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

