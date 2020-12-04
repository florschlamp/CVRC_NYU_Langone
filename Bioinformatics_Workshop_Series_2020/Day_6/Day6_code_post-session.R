## 
## Friday December 4th
## Day 6: “Data Analysis, Part 3: Downstream Analyses”
##

## RNA-seq data analysis using public 'airway' test data
# ===
# RNA-Seq experiment on four human airway smooth muscle cell lines 
# treated with dexamethasone
# ===

## before session:
# have 'biomaRt', 'clusterProfiler', and 'org.Hs.eg.db' installed

### Code:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)

## 1) set working directory to source file (this script) location
getwd()
# Session > Set Working Directory > To Source File Location

## 2) read in data files
load("DESeq_results.Rdata")

head(results_dex_treated_vs_control)

### translate IDs to gene names ###
input_IDs <- as.character(rownames(results_dex_treated_vs_control))

head(input_IDs)

library(biomaRt)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl") <-- if you have mouse

## alternatives specifying host:
#marts <- useMart("ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org") # host="useast" or "uswest"
#human <- useDataset("hsapiens_gene_ensembl", mart=marts) <-- if you have human
#mouse <- useDataset("mmusculus_gene_ensembl", mart=marts) <-- if you have mouse

# if you want to see all possible attributes and filters you can use:
humanatts <- listAttributes(human)
humanfilters <- listFilters(human)

gene_convert <- getBM(attributes = c("ensembl_gene_id","external_gene_name"),
                      filters="ensembl_gene_id",
                      values=input_IDs, mart=human)
# gene_convert is now your dictionary
head(gene_convert)

# SAVE THIS TRANSLATION DICTIONARY
# names can change

# 1) save as file for your records (.csv)
write.csv(gene_convert, file="geneID_translator_DATE.csv")
# 2) save as Rdata file, because excel LOVES to corrupt gene names (i.e.: Sep1)
save(gene_convert, file="geneID_translator_DATE.Rdata")

load("geneID_translator_DATE.Rdata")

# to create new column with gene names, you can map values using this code:
head(results_dex_treated_vs_control)

results_dex_treated_vs_control$gene_name <- as.character(plyr::mapvalues(x=rownames(results_dex_treated_vs_control),
                                                                         from=gene_convert$ensembl_gene_id,
                                                                         to=gene_convert$external_gene_name))

head(results_dex_treated_vs_control)
# you can get other information (check listAttributes/listFilters options above)

# NOTE: if there is no translation available, this 'plyr::mapvalues' function
# will keep the original name. example:
tail(results_dex_treated_vs_control) # look at gene "ENSG00000283111"



# you can get other information (check listAttributes/listFilters options above)
# (2) For example getting chromosome numbers as well:
gene_convert <- getBM(attributes = c("ensembl_gene_id","chromosome_name"),
                      filters="ensembl_gene_id",
                      values=input_IDs, mart=human)
head(gene_convert)
# we don't want to use the 'plyr::mapvalues' in this case, instead use 'merge'

head(gene_convert)
colnames(gene_convert) <- c("gene_ID","chr")
head(gene_convert)


head(results_dex_treated_vs_control)
results_dex_treated_vs_control$gene_ID <- rownames(results_dex_treated_vs_control)
head(results_dex_treated_vs_control)

results_table <- merge(results_dex_treated_vs_control,gene_convert,all=T)
head(results_table)
tail(results_table)

# IMPORTANT NOTE: the reason we don't want NAs when translating from geneID to gene_name
# is because this translation is very useful to create new data frames 
# (for volcano plots for example!) with gene_name as rownames. And rownames cannot be NAs.

# FINAL NOTE: this translation is not final/static. As the biomart database is updated,
# the translations can change too. So always save the 'gene_convert' dictionary as a file
# ideally with the date when you queried biomart on the name.

# One more example:
### get mouse homolog genes for human genes (or viceversa) ###
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

gene_convert = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                      values = input_IDs, mart = mouse, 
                      attributesL = c("hgnc_symbol"), martL = human)

#####


#####
## Filter genes of interest

## filter NAs
table(complete.cases(results_table$padj))
# 9917 genes have NAs, should be removed

res_filtered <- results_table[complete.cases(results_table$padj),]
dim(results_table) # originally 25,258 genes
dim(res_filtered)  # 15,341 genes

# lets start with significant genes (padj < 0.05), with a logFC cutoff of < -2 or > 2
res_of_interest <- as.data.frame(res_filtered[res_filtered$padj < 0.05 & 
                                                abs(res_filtered$log2FoldChange) > 2,])
head(res_of_interest)
dim(res_of_interest)
# we end up with 169 genes
genes_of_interest <- res_of_interest$gene_name

write.csv(genes_of_interest,file="top_expressed_genes.csv")

### Functional Enrichment analyses:

### Panther
# 1) copy paste IDs to website here: http://www.pantherdb.org/
# 2) choose 'homo sapiens'
# 3) choose 'Statistical annotation set' with the annotation set 'GO biological process complete'

# in next page, it asks for select reference list. 
# you can upload the whole list of genes before the subset,
# or you can just choose on the left:
# "Use a Reference List that includes all genes from a whole genome", and pick human genome

# WARNING: using gene names vs gene IDs can give different results
# show example:
genes_of_interest <- res_of_interest[,c("gene_name","gene_ID")]
write.csv(genes_of_interest,file="top_expressed_genes.csv")

#####

## There are 2 other common enrichment tools, clusterProfiler and GOrilla.
## First we will run clusterProfiler, which is an R package that can do lots of enrichment
## type analyses

library(clusterProfiler)
## For human data, we'll need to install the following database
## BiocManager::install("org.Hs.eg.db")
## For any other organism, you'll need to install the appropriate database.
## For example, here is the database for mouse
## BiocManager::install("org.Mm.eg.db")

## For today, we'll run both enrichGO and enrichKEGG to get GO and KEGG term enrichments
## for our genes of interest

head(genes_of_interest)

## these are in Ensembl gene IDs.  
## For ClusterProfiler, we need them in Entrez IDs (I know, eye roll)
## For the sake of expediency, let's also get external gene names, because we'll need those
## for the last enrichment tool, GOrilla

GOI.convert <- getBM(attributes = c("ensembl_gene_id","external_gene_name", "entrezgene_id"),
                     filters = "ensembl_gene_id",
                     values = genes_of_interest$gene_ID,
                     mart = human)
head(GOI.convert)

## now we can run ClusterProfiler straight from this

?enrichGO

GOI.GO <- enrichGO(GOI.convert$entrezgene_id, OrgDb = "org.Hs.eg.db", ont = "BP")
head(GOI.GO)

## if you want to write out the GO results, you can easily do that

write.csv(as.data.frame(GOI.GO), file = "GO_results_top_expressed_genes.csv")

## there are also some default plotting behaviors in ClusterProfiler

dotplot(GOI.GO)
barplot(GOI.GO)

## we can also run enrichKEGG

?enrichKEGG
GOI.KEGG <- enrichKEGG(GOI.convert$entrezgene_id, organism = "hsa")
dotplot(GOI.KEGG)
barplot(GOI.KEGG)


## now we can also write out results for GOrilla, which is web-based
## it takes the gene name ("external_gene_name") as input

## GOrilla works in 2 modes
## first, it can take a single list of genes ranked by some metric of interest
## let's rank our genes by log2FoldChange, high to low, and run that through GOrilla

res_filtered.rank <- res_filtered[order(-res_filtered$log2FoldChange),]
head(res_filtered.rank)

## we already have the gene name in there, so we can just write that out as a simple file
## we want just one gene per line, and no quotes.

write.table(res_filtered.rank$gene_name, quote = F, sep = "\t", row.names = F, col.names = F,
            file = "genes_ranked_LFC.txt")

## we can take a quick look at this that it looks right by clicking on the file name

## looks good, let's try running that in GOrilla
## http://cbl-gorilla.cs.technion.ac.il/

## As you can see from the home page, you can also run it using 2 unranked lists of genes
## This would be something like genes that meet a cutoff, vs genes that don't meet a cutoff
## In our case, let's try running our genes of interest, vs. all other genes in the filtered
## dataset

write.table(GOI.convert$external_gene_name, quote = F, sep = "\t", row.names = F, col.names = F,
            file = "top_expressed_genes.txt")

## to get the background genes, let's subset the res_filtered gene names to only include those 
## that are NOT in the genes of interest
head(res_filtered)
mybkgrd <- subset(res_filtered, !gene_name %in% GOI.convert$external_gene_name) 
write.table(mybkgrd$gene_name, quote = F, sep = "\t", row.names = F, col.names = F,
            file = "background_genes.txt")