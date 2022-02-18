######Code used for the statistical analysis of MDA-MB-231 metabolimics data presented in 3.12######

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
setwd(dir)
mda <- read.csv(file.choose())
mda <- as.data.frame(mda)


######Normalization#########
samples <- mda$sample_name
treat <- mda$Treatment
mda.sum <- select(mda, -c("sample_name", "Treatment", "Phosphorylcholine"))

mda.sum$sums <- rowSums(mda.sum)
max <- max(mda.sum$sums)
mda.norm <- mda.sum %>% mutate(mda.sum, coef =  max/sums) %>% 
                                      mutate(across(where(is.numeric), ~.*coef)) %>%
                                            select(-c("sums", "coef"))

mda.norm.g <- mda.norm 
mda.norm.t <- as.data.frame(t(mda.norm.g))
colnames(mda.norm.t) <- samples
mda.norm.t <- rownames_to_column(mda.norm.t, "metab")

res.class <- read.csv("metaboclass.csv")
mda.norm.class <-  dplyr::inner_join(x = mda.norm.t, y= res.class, by = c("metab" = "metab"), suffix = c("_x", "_y"))
write.csv(mda.norm.class, "mda_norm_samples.csv")


mda.glog <- mda.norm
mda.glog <- glog(mda.glog, lambda = 0)

#####Calculate Log2 Fold Change####
mda.norm.g$samples <- samples
mda.norm.g$treat <- treat
mda.melt <- melt(mda.norm.g, na.rm = TRUE)
mda.cast <- dcast(mda.melt, treat ~ variable, fun.aggregate = mean)
mda.cast <- column_to_rownames(mda.cast, var = "treat")
mda.cast.t <- as.data.frame(t(mda.cast))

log2 <- function(x,y)(log(x/y,2))
mda.fc <- mda.cast.t %>% rownames_to_column("metabo") %>%
                          mutate(across(where(is.numeric), ~log2(.,siCTRL))) %>%
                          select(-c("siCTRL"))

######Linear Models#######
metabolites <- colnames(mda.glog)
mda.lim <- mda.glog
mda.lim.t <- t(mda.lim)

mda.lim.t.df <- as.data.frame(mda.lim.t)
mda.lim.t.df <- rownames_to_column(mda.lim.t.df, "metabo")

names <- treat
groups <- factor(names)
design <- model.matrix(~0 + groups)
colnames(design) <- c( "siCTRL","siTFEB_1", "siTFEB_2")
con <- makeContrasts(siTFEB_1 - siCTRL, siTFEB_2 - siCTRL, levels=design)

fit <- lmFit(mda.lim.t, design)
fit2 <- contrasts.fit(fit, con)
fit2 <- eBayes(fit2)

l.1 <- topTable(fit2, coef = 1, number=Inf)
l.1 <- rownames_to_column(l.1, "metab")
l.2 <- topTable(fit2, coef = 2, number=Inf)
l.2 <- rownames_to_column(l.2, "metab")

res <- dplyr::inner_join(x = l.1, y= l.2, by = c("metab" = "metab"), suffix = c("_tf_1", "_tf_2"))
res <- dplyr::inner_join(x = res, y= mda.fc, by = c("metab" = "metabo"), suffix = c("_lm", "_tfc"))
write.csv(res, "mda_sitf_limma_normglog.csv")

