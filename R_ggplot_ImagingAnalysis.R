#####Used for creating the Figures in chapter 4.2, 4.4, 4.5, 4.8, and 4.9###### 

setwd(dir)

library(ggplot2)
library(dplyr)
library(tibble)
library(forcats)
library(ggridges)
library(readr)

#####Import and filter data from cell profiler results######
hci.edu <- read_csv(file.choose())
head(hci.edu, n = 5)


hci.edu.e <- filter(hci.edu, !treat %in% c("nedu","nab"))
nedu <- filter(hci.edu, treat == "nedu")
nab <- filter(hci.edu, treat == "nab")

pos <- max(nedu$mean_edu)
max.e <- max(hci.edu.e$mean_edu)
min <- min(hci.edu.e$mean_edu)

hci.edu.e$treat <- as.factor(hci.edu.e$treat)
levels(hci.edu.e$treat) <- c("siCTRL", "siTFEB#1", "siTFEB#2")
hci.edu.e$exptreat <- interaction(hci.edu.e$treat, hci.edu.e$exp)

#######Normalization of intensity value across experiments######
corany <- function(a,b){
  corfac <- a  %>% group_by(exptreat) %>% summarise(mean = mean(.data[[b]])) %>%
    mutate(corfac = mean[1]/mean) %>% select(-mean)  
  
  cor <- left_join(x = a, y = corfac) %>% mutate("cor_{b}" := (.data[[b]]*corfac)) %>% select(-corfac)                   
  return(cor)
}

hci.edu.a <- hci.edu.e %>% corany("mean_edu") %>% corany("int_dna") 

######Downsample######
hci.edu.ds <- hci.edu.a %>% group_by(treat) %>% sample_n(2500)
hci.edu.c <- hci.edu.a

#######Univariate density plots######
qc <- ggplot(hci.edu.c, aes(factor(exptreat),cor_int_dna))
qc + geom_boxplot() +
  scale_x_discrete()


dna.dense <- ggplot(hci.edu.a, aes(cor_int_dna, fill = treat)) 
dna.dense + geom_density(adjust=1.1, alpha=0.3) + 
  theme_classic()+
  labs(title="BT549", x="DNA intensity", y="Density", fill = "siRNA") +  
  scale_x_continuous(limits = c(0, 1500))+
  scale_fill_manual( values = c("red", "turquoise3", "cyan1"))+
  theme(text = element_text(size=16),
        plot.title = element_text(hjust = 0.5))+
  facet_grid(rows = vars(time))

edu.dense <- ggplot(hci.edu.c, aes(cor_mean_edu, fill = treat)) 
edu.dense + geom_density(adjust=1.2, alpha=0.3) + 
  theme_classic()+
  labs(title="BT549", x="MFI EdU", y="Density", fill = "siRNA") +  
  scale_x_continuous(trans = "log10", limits = c(min, max.e))+
  scale_fill_manual( values = c("red", "turquoise3", "turquoise3"))+
  theme(text = element_text(size=16),
        plot.title = element_text(hjust = 0.5))+
  facet_grid(rows = vars(time))

h2a.dense <- ggplot(hci.edu.c, aes(mean_h2ax, fill = treat)) 
h2a.dense + geom_density(adjust=1.2, alpha=0.3) + 
  theme_classic()+
  labs(title="BT549", x="MFI EdU", y="Density", fill = "siRNA") +  
  scale_x_continuous(limits = c(0, 0.2))+
  scale_fill_manual( values = c("red", "turquoise3", "turquoise3"))+
  theme(text = element_text(size=16),
        plot.title = element_text(hjust = 0.5))+
  geom_vline(xintercept = 0.01)

######2d and 3d plots#####
gg <- ggplot(hci.edu.c, aes(x=cor_int_dna, y=mean_edu))
gg + geom_jitter(aes(colour = treat), alpha=0.6, size = 0.5)+
  scale_y_continuous(trans = "log10", limits = c(min, max.e))+
  scale_colour_manual( values = c("red", "turquoise3", "turquoise3"))+
  xlim(0, 1400)+
  labs(title = "hci-MB-231", x="DNA Integrated Intensity", y="MFI EdU", colour = "siRNA")+
  theme(text = element_text(size=16),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "grey88"))+
  facet_grid(cols = vars(treat))

h2ax <- hci.edu.c %>% ggplot(aes(x=cor_int_dna, y=mean_edu))
h2ax + geom_jitter(aes(colour = mean_h2ax), alpha=0.8, size = 0.9)+
  scale_y_continuous(trans = "log10", limits = c(min, max.e))+
  scale_colour_gradient2(low = "white", mid = "deepskyblue1", high = "Red", na.value = "Red", 
                           midpoint = 0.012, limits = c(0,0.32))+
  xlim(0, 1400)+
  labs(title = "BT549", x="DNA Integrated Intensity", y="MFI EdU", colour = expression(paste("MFI ", gamma, "H2AX")))+
  theme_classic()+
  facet_grid(cols = vars(treat))+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 1))

########CC phase intensity values##############
x1 <- 200
x2 <- 700
x3 <- x2+(x2-x1)
y1 <- min*0.9999
y2<- pos
y4 <- max.e

#######Gate validation###########

gg.boxes <- ggplot(hci.edu.c, aes(x=cor_int_dna, y=cor_mean_edu))
gg.boxes + geom_jitter(aes(colour = "red4"), alpha=0.5, size = 0.5)+
  scale_y_continuous(trans = "log10", limits = c(min*0.9999, max.e))+
  xlim(x1, 1800)+
  labs(title = "hci-MB-231", x="DNA Integrated Intensity", y="MFI EdU")+
  theme(text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "grey88"))+
   annotate("rect", xmin= x1, xmax= x2, ymin= y1, ymax= y2, 
           alpha=0, color = "black", linetype = "dotted")+
  annotate("rect", xmin= x1, xmax= x2, ymin= y2, ymax= y4, 
           alpha=0, color = "black", linetype = "dotted")+
  annotate("rect", xmin=x2, xmax=x3, ymin=y2, ymax=y4, 
           alpha=0, color = "black", linetype = "dotted")+
  annotate("rect", xmin=x2, xmax=x3, ymin=y1, ymax=y2, 
           alpha=0, color = "black", linetype = "dotted")


#####Gates to cell cycle phase#####
hci.edu.tib <- hci.edu.c %>% as_tibble 
gated.hci.edu <- hci.edu.tib %>% mutate(phase = ifelse(between(hci.edu.tib$cor_mean_edu,y1,y2) & between(hci.edu.tib$cor_int_dna,x1,x2), "G1",
                                                       ifelse(between(hci.edu.tib$cor_mean_edu,y3,y4) & between(hci.edu.tib$cor_int_dna,x1,x2), "S1",
                                                              ifelse(between(hci.edu.tib$mean_edu,y3,y4) & between(hci.edu.tib$cor_int_dna,x2,x3), "S2",
                                                                     ifelse(between(hci.edu.tib$cor_mean_edu,y1,y2) & between(hci.edu.tib$cor_int_dna,x2,x3), "G2", "OT"))))) 


tgated.hci.edu <- gated.hci.edu %>% mutate(phase2 = ifelse(between(hci.edu.tib$cor_mean_edu,y1,y2) & between(hci.edu.tib$cor_int_dna,x1,x2), "G1",
                                                           ifelse(between(hci.edu.tib$cor_mean_edu,y2,y4) & between(hci.edu.tib$cor_int_dna,x1,x3), "S",
                                                                  ifelse(between(hci.edu.tib$cor_mean_edu,y1,y2) & between(hci.edu.tib$cor_int_dna,x2,x3), "G2", "OT")))) 



tgated.hci.edu <- tgated.hci.edu%>%filter(phase != "OT") %>%filter(phase2 != "OT")


write.csv(tgated.hci.edu, "3pGated_BT5_sitf_EDU_H2AX_v2.csv")

gated.edu.f <- tgated.hci.edu
gated.edu.f$phase <- factor(gated.edu.f$phase2, 
                            levels = c("G1", "S", "G2"))

gated.edu.f$treatphase <- interaction(gated.edu.f$treat, gated.edu.f$phase)

#######Gate validation and intensity by phase######
ww <- ggplot(gated.edu.f, aes(treatphase, mean_h2ax))
ww + 
  geom_jitter(width = 0.2, colour = "deepskyblue4", alpha = 0.2)+
  geom_boxplot(fill = "grey88", colour = "black",  
               outlier.shape = NA, alpha=0.1, 
               fatten = 2, lwd = 0.75, width = 0.75, notch = TRUE)+
  theme_classic()+
  labs(title = "No Primary Antibody", y="MFI H2AX", 
       x="Cell Cycle Gate")+
  ylim(0,0.1)+
  theme(text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black"))

