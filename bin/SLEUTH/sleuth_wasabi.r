#load the libraries
library(wasabi)
library(sleuth)
library(BiocParallel)

###Preferences and Directories###

#set the working directories
#setwd("/Volumes/NOSTRADOMUS")

#point to directory
dir <- "salmon_output"

#import metadata file
samples <- read.table(file.path("sample_GEUVADIS_tximport_unique.txt"), header = TRUE)

#import where files are
sfdir <- file.path(dir, samples$RUN)

#use wasabi to convert salmon to kallisto format
#prepare_fish_for_sleuth(sfdir)

#Create a dataframe of the data files that were converted and the corresponding experimental information
s2c <- data.frame(sample = sfdir, 
                  genotype = samples$RUN,
                  population = samples$POP,
                  center= samples$CENTER, 
                  lab= samples$LAB,
                  experiment= samples$EXPERIMENT,
                  path = sfdir)

s2c$path <- as.character(s2c$path)	 	 
s2c	 

###TX to Gns information ###

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

### Use Sleuth ###

# model the experimental question	 	 
model <- ~center	 

# load data and fit the model
so <- sleuth_prep(s2c, target_mapping = tx2gene, aggregation_column = 'ens_gene', extra_bootstrap_summary = FALSE)	 

###PCA###
# 1. Open jpeg file
png("rplot.png", width = 350, height = "350")
# 2. Create the plot
#plot_pca(so, color_by = 'population')
# 3. Close the file
#dev.off()

so <- sleuth_fit(so, ~center, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

#FInally model the question 
models(so)

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 20)

###PCA###
# 1. Open jpeg file
#png("rplot_gene.png", width = 350, height = "350")
#plot_group_density(so, use_filtered = TRUE, units = "est_counts",
#  trans = "log", grouping = setdiff(colnames(so$sample_to_covariates),
 # "sample"), offset = 1)
# 3. Close the file
#dev.off()

write.csv(sleuth_table, file = "Geuvadis_sleuth.csv")
write.csv(sleuth_significant, file = "Geuvadis_sleuth_significant.csv")
saveRDS(so, file ="GEVUADIS_sleuth_so.rds")
