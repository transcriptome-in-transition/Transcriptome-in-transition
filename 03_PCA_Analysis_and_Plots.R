#######################################################
### PRINCIPAL COMPONENT ANALYSES & PLOTS            ###
#######################################################
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(dplyr)
library(gridExtra)
library(tidyverse)
library(ggthemes)
library(grid)
library(gridExtra)
library(ggpubr)

#Updated from previous version to reflect additional sample exclusions (ie. this is PCA for "vanilla" samples)
### Setting folders:
folderIN <- "."
folderIN_data <- paste0(folderIN,"/Workspaces")
pathIN_data <- paste0(folderIN_data,"/Normalized_data.Rdata")

### Loading data prepared data:
load(pathIN_data)
# FEMALES:
datF  <- exprs(Exgirl)
phenF <- pData(Exgirl)
## Remove outlier samples:
INDF  <- which(!(rownames(phenF) %in% c("7032-07-118.CEL","7032-07-282.CEL")))
datF  <- datF[,INDF]
phenF <- phenF[INDF,]
# MALES:
datM  <- exprs(Exboy)
phenM <- pData(Exboy)
INDM  <- which(!(rownames(phenM) %in% c("7032-07-206.CEL","7032-07-203.CEL","7032-07-114.CEL","7032-07-194.CEL")))
datM  <- datM[,INDM]
phenM <- phenM[INDM,]
## Quick check:
all(rownames(phenF)==colnames(datF))
all(rownames(phenM)==colnames(datM))

pcaModelboy <- prcomp(t(datM))
expVarboy <- round(pcaModelboy$sdev^2/sum(pcaModelboy$sdev^2)*100)

pcaModelgirl <- prcomp(t(datF))
expVargirl <- round(pcaModelgirl$sdev^2/sum(pcaModelgirl$sdev^2)*100)

##Plotting PCA outputs using GGPlot##
dir.create("./PCA_plots")
plots<-c("./PCA_plots")


pcaboy<-as.data.frame(pcaModelboy$x)
pdM<-phenM[,1:5]
MPCA<-cbind(pdM,pcaboy)

pcagirl<-as.data.frame(pcaModelgirl$x)
pdF<-phenF[,1:5]
FPCA<-cbind(pdF,pcagirl)


MPCA <- MPCA %>% 
  mutate(Larval.food.level = recode(Larval.food.level, 
                                    Poor = "0.25x", Control = "1x", Rich = "2.5x"),
         Adult.food.level = recode(Adult.food.level, 
                                   P = "0.25x", C = "1x", R = "2.5x"))
FPCA <- FPCA %>% 
  mutate(Larval.food.level = recode(Larval.food.level, 
                                    Poor = "0.25x", Control = "1x", Rich = "2.5x"),
         Adult.food.level = recode(Adult.food.level, 
                                   P = "0.25x", C = "1x", R = "2.5x"))

bl<-c("#008BBA")
gr<-c("#D0C91F")
re<-c("#DC403B")

opts<-theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.key = element_blank(),
            legend.title = element_text(size=10),
            legend.text = element_text(size = 10)) 

fills<-c(bl, gr, re)
shapes<-c(15,16,17)

### Recreate plots -------------------------------------------------------------
MPCA_plot <- ggplot(data = MPCA,
                    aes(x = PC1,y = PC2,color = Larval.food.level,shape = Adult.food.level)) +
  geom_point(size=3,alpha=0.8) + 
  opts + 
  scale_color_manual(values=fills) + 
  scale_shape_manual(values=shapes) +
  xlab(paste("PC1\n(Variation explained: ",expVarboy[1],"%)")) + 
  ylab(paste("PC2\n(Variation explained: ", expVarboy[2],"%)")) +
  labs(col="Larval diet", shape = "Adult diet") +
  theme(aspect.ratio=1)

ggsave("PCA_plots/Males_PC1_versus_PC2_color_LFL.pdf",MPCA_plot, width=6,height=4,useDingbats=FALSE)

FPCA_plot <- ggplot(data = FPCA,
                    aes(x = PC1,y = PC2,color = Larval.food.level,shape = Adult.food.level)) +
  geom_point(size=3,alpha=0.8) + 
  opts + 
  scale_color_manual(values=fills) + 
  scale_shape_manual(values=shapes) +
  xlab(paste("PC1\n(Variation explained: ",expVargirl[1],"%)")) + 
  ylab(paste("PC2\n(Variation explained: ",expVargirl[2],"%)")) +
  labs(col="Larval diet", shape = "Adult diet") +
  theme(aspect.ratio=1)

ggsave("PCA_plots/Females_PC1_versus_PC2_color_LFL.pdf", FPCA_plot, width=6,height=4,useDingbats=FALSE)

### Create boxplots ------------------------------------------------------------

# Create fn for boxplots
make_pc_boxplots <- function(data, mapping, xlabel, ylabel, var_explained, label_sex){
  grob <- grobTree(textGrob(label_sex, x=0.825,  y=0.95, hjust=0,
                            gp=gpar(col="red", fontsize=13, fontface="italic")))
  ggplot(data, mapping) +
    geom_boxplot() + 
    geom_point() +
    opts +
    xlab(xlabel) +
    ylab(paste(ylabel, "\n(Variation explained: ",var_explained, "%)")) + 
    annotation_custom(grob)
}

# Make all boxplots for PC1 and PC2 by sex and larval and adult diet
FPC1_LFL <- make_pc_boxplots(FPCA, aes(x = Larval.food.level, y = PC1), 
                             "Larval diet", "PC1", expVargirl[1], "Female")
FPC1_AFL <- make_pc_boxplots(FPCA, aes(x = Adult.food.level, y = PC1), 
                             "Adult diet", "PC1", expVargirl[1], "Female")
FPC2_LFL <- make_pc_boxplots(FPCA, aes(x = Larval.food.level, y = PC2), 
                             "Larval diet", "PC2", expVargirl[2], "Female")
FPC2_AFL <- make_pc_boxplots(FPCA, aes(x = Adult.food.level, y = PC2), 
                             "Adult diet", "PC2", expVargirl[2], "Female")
MPC1_LFL <- make_pc_boxplots(MPCA, aes(x = Larval.food.level, y = PC1), 
                             "Larval diet", "PC1", expVarboy[1], "Male")
MPC1_AFL <- make_pc_boxplots(MPCA, aes(x = Adult.food.level, y = PC1), 
                             "Adult diet", "PC1", expVarboy[1], "Male")
MPC2_LFL <- make_pc_boxplots(MPCA, aes(x = Larval.food.level, y = PC2), 
                             "Larval diet", "PC2", expVarboy[2], "Male")
MPC2_AFL <- make_pc_boxplots(MPCA, aes(x = Adult.food.level, y = PC2), 
                             "Adult diet", "PC2", expVarboy[2], "Male")

# Save plots in grid
grid_plot <- grid.arrange(FPC1_LFL,FPC1_AFL,FPC2_LFL,FPC2_AFL, 
                          MPC1_LFL, MPC1_AFL, MPC2_LFL, MPC2_AFL, ncol = 2)
ggsave("PCA_plots/PC_boxplots.pdf",grid_plot, useDingbats=FALSE)

## Combine PCA plots with boxplots ---------------------------------------------
make_pc_boxplots_no_label <- function(data, mapping, xlabel, ylabel, var_explained){
  ggplot(data, mapping) +
    geom_boxplot() + 
    geom_point() +
    opts +
    xlab(xlabel) +
    ylab(ylabel)
}

FPC2_AFL_nl <- make_pc_boxplots_no_label(FPCA, aes(x = Adult.food.level, y = PC2), 
                                         "Adult diet", "PC2", expVar[2,3])

MPC2_AFL_nl <- make_pc_boxplots_no_label(MPCA, aes(x = Adult.food.level, y = PC2), 
                                         "Adult diet", "PC2", expVar[2,2])


FPC1_LFL_nl <- make_pc_boxplots_no_label(FPCA, aes(x = Larval.food.level, y = PC1), 
                                         "Larval diet", "PC1", expVar[1,3])

MPC1_LFL_nl <- make_pc_boxplots_no_label(MPCA, aes(x = Larval.food.level, y = PC1), 
                                         "Larval diet", "PC1", expVar[1,2])


### Try to fix this to include legend
four_panel <- ggarrange(FPCA_plot, MPCA_plot, FPC2_AFL_nl, MPC2_AFL_nl, 
                        heights = c(0.7, 0.3),
                        ncol = 2, nrow = 2, 
                        labels = "auto",
                        align = "hv",
                        common.legend = TRUE, legend = "right")

six_panel <- ggarrange(FPCA_plot, MPCA_plot, FPC1_LFL_nl, MPC1_LFL_nl, FPC2_AFL_nl, MPC2_AFL_nl,
                       heights = c(1.5, 0.7, 0.7),
                       ncol = 2, nrow = 3, 
                       labels = "auto", 
                       align = "hv",
                       common.legend = TRUE, legend = "right")

ggexport(four_panel, filename = "PCA_plots/Figure_1.PCA_potential_update_four_panel.pdf")
ggexport(six_panel, filename = "PCA_plots/Figure_1.PCA_potential_update_six_panel.pdf")