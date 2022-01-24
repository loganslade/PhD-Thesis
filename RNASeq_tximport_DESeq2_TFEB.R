library(tximport)
library(DESeq2)
library(tximportData)
library(rhdf5)
library(dplyr)

######The following code is adapted from the DESeq2 vignette with credit to Michael Love, Simon Anders and Wolfgang Huber######
##Show R where the kallisto files are
dir <- "C:/Program Files/kallisto/MDARNASeq"
list.files(dir)
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
samples
files <- file.path(dir, "output", samples$samples, "abundance.h5")
names(files) <- paste0("MDA_", 559:574)
all(file.exists(files))

##Get the latest ensembl gene annotation files, use the same one the the transcript index was built from

library(ensembldb)
library(AnnotationHub)

ah <- AnnotationHub()
query(ah, "EnsDb.HSapiens")
edb <- ah[["AH69187"]]

##Make the transcript to gene table
txdf <- transcripts(edb, return.type="DataFrame")
tx2gene <- as.data.frame(txdf[,c("tx_id","gene_id")])
head(tx2gene)


##Summarize transcripts by gene
library(readr)
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)
names(txi)
head(txi$counts)
head(txi$abundance)
setwd("C:/Users/Logan Slade/Desktop/Research/Research Project/Cancer projects/48-MDA231RNASeq/Analysis")
write.csv(txi, "kallisto_tximport_genes.csv")

txicounts <- txi$counts
head(txicounts)
write.csv(txicounts, "kallisto_txcounts_genes.csv")

## Get official gene symbols
library(biomaRt)
listMarts()
ensembl <- useEnsembl(biomart = "ensembl")
listDatasets(mart = ensembl)

mart <- biomaRt::useMart(biomart ="ENSEMBL_MART_ENSEMBL", 
                         dataset = "hsapiens_gene_ensembl",
                         host = "www.ensembl.org",
                         ensemblRedirect = FALSE)

t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)

genenames <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                            mart       = mart)

#####Change matrix to data.frame####
library(tibble)
txicounts2 <- as.data.frame(txicounts)
head(txicounts2)

##Change row names to variables 
tibble::has_rownames(txicounts2)
txicounts3 <- tibble::rownames_to_column(txicounts2, "gene_id")

##Match gene symbols to gene_ids and combine
txigene <- dplyr::inner_join(x=txicounts3, y=genenames, by = c("gene_id" = "ensembl_gene_id"))
write.csv(txigene, "kallisto_txcounts_genes.csv")

##Same thing for gene abundances 
txiabun <- txi$abundance
txiabun <- as.data.frame(txiabun)
head(txiabun)
txiabun2 <- tibble::rownames_to_column(txiabun, "gene_id") 
txiabungene <- dplyr::inner_join(x=txiabun2, y=genenames, by = c("gene_id" = "ensembl_gene_id"))
write.csv(txiabungene, "kallisto_txabun_genes.csv")

###BEGIN DESEq2 code####
samples.info <- read.table(file.path(dir, "samples_info.txt"), header = TRUE)
samples.info$samples <- as.factor(samples.info$samples)
samples.info$knockdown <- as.factor(samples.info$knockdown)
samples.info$sirna <- as.factor(samples.info$sirna)
ddstxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples.info,
                                   design = ~ sirna) #change for 1-2 combo

ddstxi$knockdown <- relevel(ddstxi$knockdown, ref = "CTRL")
ddstxi$sirna <- relevel(ddstxi$sirna, ref = "siCTRL")
keep <- rowSums(counts(ddstxi)) >= 5
dds <- ddstxi[keep]


dds <- DESeq(dds)
dds.count <- as.data.frame(counts(dds, normalized=TRUE))
dds.count <- tibble::rownames_to_column(dds.count, "gene_id")
genes.norcount <- 
  dplyr::inner_join(x=genenames, y= dds.count, by = c("ensembl_gene_id" = "gene_id"))

write.csv(genes.norcount, "DESeq2_normalized_Counts_v3.csv")

restf1 <- results(dds, contrast=c("sirna", "siTFEB_1", "siCTRL"))
restf2 <- results(dds, contrast=c("sirna", "siTFEB_2", "siCTRL"))


restf1porder <- restf1[order(restf1$padj),] 
restf2porder <- restf2[order(restf2$padj),]

library(tibble)
restf1porder <- as.data.frame(restf1porder)
restf2porder <- as.data.frame(restf2porder)

restf1porder <- tibble::rownames_to_column(restf1porder, "gene_id_tf1")
restf2porder <- tibble::rownames_to_column(restf2porder, "gene_id_tf2")

res.tfeb.all <- dplyr::inner_join(x=restf1porder, y=restf2porder, by = c("gene_id_tf1" = "gene_id_tf2"))
genes.tfeb.all <- dplyr::inner_join(x=res.tfeb.all, y=genenames, by = c("gene_id_tf1" = "ensembl_gene_id"))
write.csv(genes.tfeb.all, "DESeq2_TFEB_all2.csv")

unique_tfeb <- distinct(genes.tfeb.all)
unique_tfeb <- unique_tfeb[order(unique_tfeb$external_gene_name, unique_tfeb$padj.x),]
unique_tfeb <- test[!duplicated(test$external_gene_name),]

write.csv(unique_tfeb, "DESeq2_TFEB_dedup.csv")
dds.tfeb.unique <- dplyr::inner_join(x=unique_tfeb, y=genes.norcount, by = c("gene_id_tf1" = "ensembl_gene_id"))
write.csv(dds.tfeb.unique, "DESeq2_resTFEB_norcount.csv")

###Using TFEB knockdown as the level not siRNA####
ddstxi2 <- DESeqDataSetFromTximport(txi,
                                   colData = samples.info,
                                   design = ~ knockdown) #change for 1-2 combo

ddstxi2$knockdown <- relevel(ddstxi$knockdown, ref = "CTRL")
keep <- rowSums(counts(ddstxi2)) >= 5
dds2 <- ddstxi2[keep]

dds2 <- DESeq(dds2)

restfeb <- results(dds2, contrast=c("knockdown", "TFEB", "CTRL"))
restfeb <- as.data.frame(restfeb)

res.tfeb <- tibble::rownames_to_column(restfeb, "gene_id_tfeb")

combo.tfeb.all <- dplyr::inner_join(x=genes.tfeb.all, y=res.tfeb, by = c("gene_id_tf1" = "gene_id_tfeb"))
write.csv(combo.tfeb.all, "DESeq2_TFEB_all_combo.csv")

################Get results that are significant in both siRNAs##########

res.sig.tf1 <- subset(restf1porder, padj < 0.01)
res.sig.tf2 <-subset(restf2porder, padj < 0.01)

sig.tf1.genes <- dplyr::inner_join(x=res.sig.tf1, y=genenames, by = c("gene_id_tf1" = "ensembl_gene_id"))
sig.tf2.genes <- dplyr::inner_join(x=res.sig.tf2, y=genenames, by = c("gene_id_tf2" = "ensembl_gene_id"))

write.csv(sig.tf1.genes, "Sig_DESeq2_TFEB_1.csv")
write.csv(sig.tf2.genes, "Sig_DESeq2_TFEB_2.csv")

sig.results.tfeb <- dplyr::inner_join(x=res.sig.tf1, y=res.sig.tf2, by = c("gene_id_tf1" = "gene_id_tf2"))
sig.results.tfeb <- dplyr::inner_join(x=sig.results.tfeb, y=genenames, by = c("gene_id_tf1" = "ensembl_gene_id"))
write.csv(sig.results.tfeb, "Sig_DESeq2_TFEB_siRNA.csv")

resultstfeb <- dplyr::inner_join(x=restf1porder, y=restf2porder, by = c("gene_id_tf1" = "gene_id_tf2"))
resultstfeb <- dplyr::inner_join(x=resultstfeb, y=genenames, by = c("gene_id_tf1" = "ensembl_gene_id"))
setwd("C:/Users/Logan Slade/Desktop/Research/Research Project/Cancer projects/48-MDA231RNASeq/Analysis")
write.csv(resultstfeb, "DESeq2_TFEB_siRNA_all.csv")

###########Graphics###############
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
head(assay(rld), 3)

plotPCA(vsd, intgroup=c("sirna", "knockdown"))
plotPCA(vsd, intgroup=c("knockdown"))
plotPCA(rld, intgroup=c("knockdown"))

rld.tfebsub <- rld[ , rld$sirna %in% c("siCTRL", "siTFEB_1", "siTFEB_2") ]
plotPCA(rld.tfebsub, intgroup=c("sirna"))
rld.kdsub <- rld[ , rld$knockdown %in% c("CTRL", "TFEB") ]
plotPCA(rld.kdsub, intgroup=c("knockdown"))

library(ggplot2)

pca.tfeb <- plotPCA(rld.tfebsub, intgroup=c("sirna"), returnData=TRUE)
percentVar2 <- round(100 * attr(pca.tfeb, "percentVar"))
ggplot(pca.tfeb, aes(PC1, PC2, color=sirna)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar2[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar2[2],"% variance")) + 
  labs(colour = "siRNA")+
  coord_fixed() +
  theme(legend.title = element_text(size=17))+
  theme(legend.text = element_text(size=17))+
  theme(axis.text = element_text(size=17))+
  theme(axis.title = element_text(size=17))

vsd.tfebsub <- vsd[ , vsd$sirna %in% c("siCTRL", "siTFEB_1", "siTFEB_2") ]

sampleDists <- dist(t(assay(vsd.tfebsub)))
library("RColorBrewer")
library("pheatmap")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) 
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
######Enrichment Plots########
install.packages("ggplot2")
library(ggplot2)
library(dplyr)
library(readr)


tfeb.en <- read_csv(file.choose())
a <- tfeb.en

order.a <- a[order(a$logp), ]
order.a$Term <- factor(order.a$'Term', levels = order.a$'Term')

g <- ggplot(order.a, aes(logp, Term)) 
g +  geom_point(aes(colour = `Odds Ratio`), size = 6) +
  scale_colour_gradient2(low = "yellow", mid = "orange", high = "red", midpoint = 2.5) +
  theme_classic() +
  scale_y_discrete() + 
  scale_size(range = c(3,8))+
  labs(title = "Enrichr: ChEA 2016 gene sets",
       subtitle = "Transcription factors associated with downregulated genes",
       y = "",
       x = "-Log10(p-value)",
       colour = "Odds \nRatio"
  ) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, colour = "black", face = "bold", size = 12), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size= 10),
        legend.position = c(.85, .25),
        legend.key.size = unit(0.8, "cm"),
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 10, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

tfeb.gsea <- read.csv(file.choose())

tfeb.gsea <- tfeb.gsea[order(tfeb.gsea$nes), ]
tfeb.gsea$Term <- factor(tfeb.gsea$Term, levels = tfeb.gsea$Term)

head(tfeb.gsea)

h <- ggplot(tfeb.gsea, aes(nes, Term))
h + geom_point(aes(colour = nes, size = logp)) +
  scale_colour_gradient2(low = "yellow", mid = "orange", high = "red", midpoint = 1.77) +
  scale_size(range = c(3,10))+
  scale_y_discrete() +
  theme_classic() +
  labs(title = "GSEA: GO Biological Processes",
       subtitle = "Terms downregulated by TFEB Knockdown",
       y = "",
       x = "Normalized Enrichment Score",
       colour = "NES",
       size = "-Log10 \nNominal \nP-Value") +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, colour = "black", face = "bold"), 
        axis.text.y = element_text(colour = "black", face = "bold", size =14),
        axis.title.x = element_text(face = "bold", size=12, hjust = 0.5),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
        plot.subtitle = element_text(hjust = 0.5, size= 12),
        legend.key.size = unit(0.8, "cm"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank()) 