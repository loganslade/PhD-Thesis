library(ggplot2)
library(readr)
library(tibble)
library (plyr)
library(dplyr)
library(data.table)
library(ggfortify)
library(reshape2)
library(Rtsne)
library(tidyr)
library(survival)
library(survminer)

dir <- "C:/Users/lslad/Desktop/Research/Research Project/Cancer projects/Bioinformatics/METABRIC-expressiondata"
setwd(dir)

#####Data Import#####
geneex <- read_csv("AllSamples-GeneExpression.csv")
sample <- read_csv("IDs_Clinical.csv")
clinical <- dplyr::rename(sample, 
                          sample_ID = `Sample ID`,
                          age = `Age at Diagnosis`,
                          ihc_er = `ER status measured by IHC`,
                          er = `ER Status`,
                          cluster = `Integrative Cluster`,
                          lymph_nodes = `Lymph nodes examined positive`,
                          index = `Nottingham prognostic index`,
                          OS = `Overall Survival (Months)`,
                          OSstatus= `Overall Survival Status`,
                          three.gene = `3-Gene classifier subtype`,
                          size = `Tumor Size`,
                          stage = `Tumor Stage`,
                          vital = `Patient's Vital Status`,
                          grade = `Neoplasm Histologic Grade`
                          )

g.n <- colnames(geneex)
g.n <- gsub("-","_", g.n)
colnames(geneex) <- g.n

s.n <- clinical$sample_ID
s.n <- gsub("-","_", s.n)
clinical$sample_ID <- s.n

clin.can <- filter(clinical, vital != "Died of Other Causes", 
                    !is.na(vital))


clin.can <- clin.can %>% mutate(surv = ifelse(clin.can$OSstatus == "DECEASED", 1, 0))
clin.can$stage <- factor(clin.can$stage, levels = c(1,2,3,4,0))
clin.can$grade <- factor(clin.can$grade)
#clin.narm <- na.omit(clin.can)
#clin.can <- clin.narm

clin.tnbc <- filter(clin.can, Subtype == "Basal" | Subtype == "claudin-low")

tnbctype <- read.csv("TNBCtype_AO2018.csv")
t.n <- tnbctype$sample
t.n <- gsub("-","_", t.n)
tnbctype$sample <- t.n
clin.tnbc <- dplyr::inner_join(x=clin.can, y=tnbctype, by = c("sample_ID" = "sample"))
ids.type <- clin.tnbc %>% select(sample_ID, TNBCtype)

tnbc.id <- clin.tnbc$sample_ID

ntnbc <- filter(clin.can, !sample_ID %in% tnbc.id)  
ntnbc <- ntnbc %>% select(sample_ID, Subtype)
tnbc <- dplyr::rename(ids.type, Subtype = TNBCtype)

all.types <- rbind(ntnbc, tnbc)
all.types <- all.types %>% rename(pam_tnbc = Subtype)
clin.types <- dplyr::inner_join(x=clin.can, y=all.types, by = c("sample_ID" = "sample_ID"))
clin.type6 <- filter(clin.types, pam_tnbc != "claudin-low", pam_tnbc != "Basal", pam_tnbc != "NC") 
clin.type6 <- filter(clin.type6, !is.na(pam_tnbc))
######TFEB and genesets######
geneset <- c("TFEB", "TFE3", "FLCN", "PPP3CA", "PPP3CB", "PPP3CC", "PPP3R1", "PPP3R2", 
             "RRAGC", "RRAGD", "FNIP1", "FNIP2", "MAP4K3")

bcgeneex <- geneex
bcgeneex <- filter(geneex,
                     Hugo_Symbol %in% geneset) 

n <- bcgeneex$Hugo_Symbol
bcgeneex.w <- bcgeneex[,-1]

bcgene.s <- as.data.frame(t(bcgeneex.w)) 
colnames(bcgene.s) <- n
bc.b <- tibble::rownames_to_column(bcgene.s, "sample")

#######All Samples#####
bc.clin <- dplyr::inner_join(x=bc.b, y=clin.type6, by = c("sample" = "sample_ID"))
bc.surv <- bc.clin
bc.surv$Subtype <- factor(bc.surv$Subtype, levels = c("LumA", "LumB", "Normal", "Her2", "Basal", "claudin-low"))
bc.surv$Subtype <- recode_factor(bc.surv$Subtype, LumA = "Luminal A",
                                                      LumB = "Luminal B", 
                                                      Normal = "Normal-like", 
                                                      Her2 = "HER2-enriched", 
                                                      Basal = "Basal-Like")

bc.surv$pam_tnbc <- factor(bc.surv$pam_tnbc, levels = c("LumA", "LumB", "Normal_like", "Her2", 
                                                            "BL1", "BL2", "IM", "LAR", "M", "MSL", "UNS"))

bc.surv <- bc.surv %>% filter(Subtype != "NC")
bc.surv$Cellularity <- factor(bc.surv$Cellularity, levels = c("Low", "Moderate", "High"))
bc.surv$grade <- as.factor(bc.surv$grade)
bc.surv$cluster <- as.factor(bc.surv$cluster)
bc.surv$cluster <- factor(bc.surv$cluster, levels = c("1", "2", "3", "4ER+", "4ER-", "5", "6", "7", "8", "9", "10"))
bc.surv$three.gene <- factor(bc.surv$three.gene, levels = c("ER+/HER2- Low Prolif", "ER+/HER2- High Prolif", "HER2+", "ER-/HER2-")) 
bc.surv$ihc_er <- as.factor(bc.surv$ihc_er)
bc.surv$ihc_er <- factor(bc.surv$ihc_er, levels = c("Positve", "Negative"))
bc.surv$ihc_er <- recode_factor(bc.surv$ihc_er, Positve = "ER+",
                                  Negative = "ER-")                                                                                                                                 
bc.r <- bc.surv[,-1]

#######COXph models#######
res.cox <- coxph(Surv(OS, surv) ~ TFEB, data = bc.surv)
x <- summary(res.cox)
x
res <- as.data.frame(x$conf.int)
#######KM plots########
q <- quantile(bc.surv$TFEB)
l <- as.numeric(q[2])
m <- as.numeric(q[3])
h <- as.numeric(q[4])

bc.km <- bc.surv %>% mutate(strata = ifelse(TFEB  < l, "LOW", 
                                                ifelse(between(bc.surv$TFEB, l, h), "MEDIAN",
                                                        "HIGH")))


bc.km$strata <- factor(bc.km$strata, levels = c("LOW", "MEDIAN",  "HIGH"))
km <- survfit(Surv(OS, surv) ~ strata, data= bc.km)
k <- autoplot(km, conf.int = F, surv.size = 1, censor.size = 3)
k + labs(title = "METABRIC", x= "Time (Months)", y = "Overall Survival (%)",
         fill = "TFEB expression", color = "TFEB expression")+
  scale_x_continuous(limits = c(0, 240))+
  scale_y_continuous(limits = c(0.3,1), labels = scales::percent)+
  scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D"))+
  theme_classic()+
  theme(text = element_text(size=18),
        axis.text = element_text(colour = "black", size = 20),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8,.7),
        plot.margin = margin(0,1,0,0, "cm"))
ggsave("KM_tfeb_metabric_3s.pdf", device = "pdf", width = 14.94, height = 13.71, units = "cm")
###################
#########Univariate plots########
a <- bc.r %>% drop_na(Subtype)%>% filter(Subtype != "claudin-low") %>% ggplot(aes(Subtype, TFEB)) 
a + geom_jitter(alpha = 0.2, colour = "green4", size = 1, width = 0.2)+
  geom_boxplot(fill = "grey65", colour = "black",  
               outlier.shape = NA, alpha=0.75, 
               fatten = 1, lwd = 0.5, width = 0.5, notch = TRUE)+
  ylim(5.15,6.3)+
  theme_classic()+
  labs(title = "", y=expression('TFEB Log'[2]*' Intensity'), x="")+
  theme(text = element_text(size=18),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "white"))

ggsave("box_tfeb_metabric_pam50_v4.pdf", device = "pdf", width = 5.16, height = 6.61, units = "in")

b <- bc.r %>% drop_na(three.gene)%>% ggplot(aes(three.gene, TFEB)) 
b + geom_jitter(alpha = 0.2, colour = "green4", size = 1, width = 0.2)+
  geom_boxplot(fill = "grey65", colour = "black",  
               outlier.shape = NA, alpha=0.75, 
               fatten = 1, lwd = 0.5, width = 0.5, notch = TRUE)+
  ylim(5.15,6.3)+
  theme_classic()+
  labs(title = "", y=expression('TFEB Log'[2]*' Intensity'), x="")+
  theme(text = element_text(size=18),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "white"),
        plot.margin = margin(0,0,0,1, "cm"))

ggsave("box_tfeb_metabric_3gene_v4.pdf", device = "pdf", width = 4.38, height = 5.61, units = "in")

c <- bc.r %>% drop_na(ihc_er)%>% ggplot(aes(ihc_er, TFEB)) 
c + geom_jitter(alpha = 0.2, colour = "green4", size = 1, width = 0.2)+
  geom_boxplot(fill = "grey65", colour = "black",  
               outlier.shape = NA, alpha=0.75, 
               fatten = 1, lwd = 0.5, width = 0.5, notch = TRUE)+
  ylim(5.15,6.3)+
  theme_classic()+
  labs(title = "", y=expression('TFEB Log'[2]*' Intensity'), x="")+
  theme(text = element_text(size=18),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "white"),
        plot.margin = margin(0,0,0,0, "cm"))

ggsave("box_tfeb_metabric_ihc_er_v4.pdf", device = "pdf", width = 3.77, height = 4.83, units = "in")


fit <- aov(TFEB ~ ihc_er, data = bc.clin)
summary(fit)
resaov <- as.data.frame(TukeyHSD(fit)[1]) %>% write.csv("res_TFEB_ihc_er_v3.csv")

######TCGA Data########
dir <- "C:/Users/lslad/Desktop/Research/Research Project/Cancer projects/Bioinformatics/METABRIC-expressiondata/TCGA/TCGA Broad RSEM Norm"
setwd(dir)

tcgageneex <- read_csv("TCGA_brca_rnaseq_rsem.csv")
tcgasample <- read_csv("TCGA_ids_clin.csv")
ids.pam50 <- read_csv("nature_tcga_id_isubtypes.csv")

#####Generate three gene classifcation based on IHC######
tcgasample$interaction <- interaction(tcgasample$ihc_er, tcgasample$ihc_her2)
tcgasample.2 <- tcgasample %>% mutate(ihc_st = ifelse(interaction == "Positive.Negative", "ER+/HER2-", 
                                              ifelse(interaction == "Positive.Equivocal", "ER+/HER2", 
                                                     ifelse(interaction == "Positive.Positive", "ER+/HER2+",
                                                            ifelse(interaction == "Negative.Equivocal", "HER2",
                                                                   ifelse(interaction == "Negative.Positive", "HER2+",
                                                                          ifelse(interaction == "Negative.Negative", "ER-/HER2-", "NC"))))))) %>%
  drop_na(ihc_st) %>% filter(ihc_st != "NC")

#######Filter and connect genes with IDs######
tcgabcgeneex <- filter(tcgageneex,
                     gene_id %in% geneset) 

n <- tcgabcgeneex$gene_id
tcgabcgeneex.w <- tcgabcgeneex[,-1]

tcgabcgene.s <- as.data.frame(t(tcgabcgeneex.w)) 
colnames(tcgabcgene.s) <- n
tcgabc.b <- tibble::rownames_to_column(tcgabcgene.s, "sample")
tcgabc.u <- tcgabc.b[order(tcgabc.b$sample, decreasing = T),]

sample.id <- tcgabc.u$sample
sample.id.2 <- substr(sample.id, start=1, stop=12)
tcgabc.c <- tcgabc.u
tcgabc.c$sample <- sample.id.2
tcgabc.c <- tcgabc.c[!duplicated(tcgabc.c$sample),]

tcgabc.clin <- dplyr::inner_join(x=tcgabc.c, y=tcgasample.2, by = c("sample" = "patient_id"))
tcgabc.clin <- tcgabc.clin[!duplicated(tcgabc.clin$sample),]
tcgabc.clin$ihc_st <- factor(tcgabc.clin$ihc_st, levels = c("ER+/HER2-", "ER+/HER2","ER+/HER2+","HER2", "HER2+", "ER-/HER2-"))
tcgabc.clin$ihc_er <- factor(tcgabc.clin$ihc_er, levels = c("Positive","Negative"))
tcgabc.clin$ihc_er <- recode_factor(tcgabc.clin$ihc_er, Positive = "ER+", Negative = "ER-")

#####Connect genes with PAM50 ids#######
tcgabc.pam <- dplyr::inner_join(x=tcgabc.c, y=ids.pam50, by = c("sample" = "patient_id"))
tcgabc.pam <- drop_na(tcgabc.pam, subtype)
tcgabc.pam$subtype <- factor(tcgabc.pam$subtype, levels = c("Luminal A","Luminal B","Normal-like","HER2-enriched","Basal-like"))

#####BOX-plots#######
a <- tcgabc.pam %>% drop_na(subtype) %>% ggplot(aes(subtype, TFEB)) 
a + geom_jitter(alpha = 0.2, colour = "green4", size = 1, width = 0.2)+
  geom_boxplot(fill = "grey65", colour = "black",  
               outlier.shape = NA, alpha=0.75, 
               fatten = 1, lwd = 0.5, width = 0.5, notch = TRUE)+
  ylim(0,1600)+
  theme_classic()+
  labs(title = "", y="TFEB Normalized Expression", x="")+
  theme(text = element_text(size=18),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "white"))

ggsave("box_tfeb_tcga_pam50_v5.pdf", device = "pdf", width = 5.16, height = 6.61, units = "in")

b <- tcgabc.clin %>% drop_na(ihc_st) %>% ggplot(aes(ihc_st, TFEB)) 
b + geom_jitter(alpha = 0.2, colour = "green4", size = 1, width = 0.2)+
  geom_boxplot(fill = "grey65", colour = "black",  
               outlier.shape = NA, alpha=0.75, 
               fatten = 1, lwd = 0.5, width = 0.5, notch = TRUE)+
  ylim(0,1600)+
  theme_classic()+
  labs(title = "", y="TFEB Normalized Expression", x="")+
  theme(text = element_text(size=18),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "white"))

ggsave("box_tfeb_tcga_ihc_v5.pdf", device = "pdf", width = 4.15, height = 5.32, units = "in")

c <- tcgabc.clin %>% drop_na(ihc_er) %>% ggplot(aes(ihc_er, MAP4K3)) 
c + geom_jitter(alpha = 0.2, colour = "green4", size = 1, width = 0.2)+
  geom_boxplot(fill = "grey65", colour = "black",  
               outlier.shape = NA, alpha=0.75, 
               fatten = 1, lwd = 0.5, width = 0.5, notch = TRUE)+
  ylim(0,2500)+
  theme_classic()+
  labs(title = "", y="MAP4K3 Normalized Expression", x="")+
  theme(text = element_text(size=18),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "white"))

ggsave("box_MAP4K3_tcga_ihcer_v5.pdf", device = "pdf", width = 4.5, height = 5.77, units = "in")

####Stats####
fit <- aov(TFEB ~ ihc_er, data = tcgabc.clin)
summary(fit)

resaov <- as.data.frame(TukeyHSD(fit)[1]) %>% write.csv("res_TFEB_ihcer.csv")
