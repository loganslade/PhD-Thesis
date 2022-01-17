library(ggplot2)
library(readr)
library(tibble)
library(plyr)
library(dplyr)
library(data.table)
library(ggfortify)
library(reshape2)
library(Rtsne)
library(stats)
library(BiocManager)
#BiocManager::install("LMGene")
library(limma)
library(LMGene)

###############Import Data########
dir <- "/Users/lslad/Desktop/Research/Research Project/Cancer projects/43-siRNAviability-Sept418/Cayman Kinase Inhibitor Library"
setwd(dir)
mda <- read.csv(file.choose())
mda <- as.data.frame(mda)
wells <- read.csv(file.choose())
mda_ids <- dplyr::inner_join(wells,mda,by = c("Well_ID" = "Well_ID"))
mda_c <- filter(mda_ids, Con_Flag == F)
mda_c <- select(mda_c, -c(Plate.x,Well.x,Plate.y,Well.y))
mda_c <- select(mda_c, -(Con_Flag))
mda.ratio <- mda_c[,c(1:7)]
mda.log <- mda_c[,c(1:3,11:14)]

######Linear Models#######
Inh <- mda_c[,2]
wells <- mda_c[,1]
mda.lim <- select(mda.log, -c(Well_ID, Product.Name,Target))
rownames(mda.lim) <- wells

names <- c("ctrl", "ctrl", "tf2", "tf2")
groups <- factor(names)
design <- model.matrix(~0 + groups)
colnames(design) <- c( "ctrl","tf2")
con <- makeContrasts(tf2 - ctrl, levels=design)

fit <- lmFit(mda.lim, design)
fit2 <- contrasts.fit(fit, con)
fit2 <- eBayes(fit2)

l <- topTable(fit2, number=Inf)
l <- rownames_to_column(l, "well_id")
l.1 <- topTable(fit2, coef = 1, number=Inf)
l.1 <- rownames_to_column(l.1, "metab")
l.2 <- topTable(fit2, coef = 2, number=Inf)
l.2 <- rownames_to_column(l.2, "metab")
res <- l


res <- dplyr::inner_join(x = res, y= mda_c, by = c("well_id" = "Well_ID"))
write.csv(res, "sitfeb_kinaselibrary_limmares.csv")

