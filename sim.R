#Gene-level model for 53 GTEX transcriptomes
#Prediction performance on the median absolute expression levels across tissues.
#The predicted log2-TPM values for 2,705 genes coded on chromosome 1 were compared with the actual median log2-TPMs across 53 tissues.

setwd("C:/Users/seonyeong/Desktop/SystemBiology/team project/simulation")
#library load
options(stringsAsFactors = FALSE)
#install.packages("see")
library(see)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("ggpubr")
library(ggpubr)
#install.packages("circlize")
library(circlize)
#install.packages("ComplexHeatmap")
library(ComplexHeatmap)
#install.packages("dplyr")
library(dplyr)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("hexbin")
library(hexbin)
source("C:/Users/seonyeong/Desktop/SystemBiology/team project/simulation/test_Chr1/fig_util.R")
# p - Median scatter ####
ds <- readRDS("C:/Users/seonyeong/Desktop/SystemBiology/team project/simulation/test_Chr1/pred.rds")

ds
# gene-level ###
ds %>%
  filter(tissue%in%"MedianExp")%>%
  filter(model=="Promoter & RNA features") -> data

ggplot(data,aes(pred,actual))+
  geom_hex()+
  scale_fill_flat_c(name="") +
  theme_modern()+
  stat_cor(method="spearman")+
  ggtitle("Gene")+
  xlab("Prediction, log-TPM")+
  ylab("Actual, log-TPM") -> p
plot(p)

#Heatmap
# acutal expression

ds%>%
  filter(tissue!="MedianExp")%>%
  dplyr::select(gene,tissue,actual)%>%
  filter(!duplicated(paste0(gene,tissue)))%>%
  spread(tissue,actual) -> mat
mat_actual=mat[,-1]%>%as.matrix()
rownames(mat_actual)=mat$gene

# predicted expression
ds%>%
  filter(tissue!="MedianExp")%>%
  filter(model=="Promoter & RNA features")%>%
  dplyr::select(gene,tissue,pred)%>%
  spread(tissue,pred) -> mat
mat_pred=mat[,-1]%>%as.matrix()
rownames(mat_pred)=mat$gene


jet.colors <- colorRampPalette(c("#00007F", "#007FFF", "cyan",
                                 "#7FFF7F", "yellow", "#FF7F00", "red"))
col_pos = jet.colors(7)

# calc correlation between acutal and prediction
cor_between_pred_actual=cor(mat_actual,mat_pred)
colnames(cor_between_pred_actual)=colnames(cor_between_pred_actual)%>%rename_tissue_bar()
rownames(cor_between_pred_actual)=rownames(cor_between_pred_actual)%>%rename_tissue_bar()

# get predefined tissue color 
alphabet_order <- readRDS("C:/Users/seonyeong/Desktop/SystemBiology/team project/simulation/test_Chr1/cluster_color.rds")
tissue_color=alphabet_order$color
names(tissue_color)=alphabet_order$key

# make column and row annotation
tc = columnAnnotation(Tissue = names(tissue_color),
                      col = list(Tissue = tissue_color),
                      show_annotation_name = FALSE)
ha = rowAnnotation(Tissue = names(tissue_color),
                   col = list(Tissue = tissue_color),
                   show_annotation_name = FALSE)

# cluster correlation matrix
indx=match(names(tissue_color),colnames(cor_between_pred_actual))
cor_between_pred_actual=cor_between_pred_actual[indx,indx]
hc_rows=hclust(d = dist(cor_between_pred_actual),method = "ward.D2")

# generate heatmap
ht = Heatmap(cor_between_pred_actual,
             name="Correlation",
             cluster_rows = rev(as.dendrogram(hc_rows)),
             cluster_columns = hc_rows,
             left_annotation = ha,
             bottom_annotation = tc,
             col = g2_colorRamp2(values =c(-0.5,-0.3,-0.1,
                                           0,
                                           0.1,0.3,0.5 )),
             column_title="Predicted expression",
             column_title_side="bottom",
             row_title="Actual expression",
             show_row_names = TRUE, 
             row_names_side = "left",
             show_row_dend = FALSE,
             show_column_dend = FALSE,
             column_names_max_height = unit(20,"cm"),
             row_names_max_width = unit(20,"cm"))
draw(ht)
# boxplot
# load DeepLIFT score
DeepLIFT_mean <- readRDS("C:/Users/seonyeong/Desktop/SystemBiology/team project/simulation/test_Chr1/DeepLIFT_mean.rds")

DeepLIFT_mean%>%
  filter(tissue_name!="MedianExp")%>%
  filter(Type!="other")%>%
  group_by(feat_name,Type)%>%
  summarise(shap_mean=mean(abs(shap_mean)))%>%
  ungroup()%>%
  arrange(desc(shap_mean))%>%
  mutate(feat_name=paste0(Type,"_",feat_name))%>%
  mutate(Var="Gene")-> ds

ds%>%
  group_by(Var)%>%
  mutate(shap_mean=replace(shap_mean,shap_mean==0,min(shap_mean[shap_mean>0]/2)))%>%
  mutate(Type=factor(Type,levels = c("TF","RBP","miRNA")))%>%
  ggviolin(., x = "Type", y = "shap_mean", fill = "Type",
           add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = list(c(1,2),c(1,3),c(2,3)),size=2.5)+
  scale_y_log10(1e-10, 0)+
  theme_modern(plot.title.size = 20)+
  facet_wrap(~Var,ncol = 1)+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 11))+
  d3_fill(label = NULL)+
  ylab("DeepLIFT score")+
  xlab(NULL)-> p

plot(p)
#second heatmap 
# load DeepLIFT score
a <- readRDS("C:/Users/seonyeong/Desktop/SystemBiology/team project/simulation/test_Chr1/DeepLIFT_mean.rds")
DeepLIFT_mean <- readRDS("C:/Users/seonyeong/Desktop/SystemBiology/team project/simulation/test_Chr1/DeepLIFT_mean.rds")%>%
  filter(tissue_name!="MedianExp",
         Type!="other")

# keep top5 influential regulators
DeepLIFT_mean%>%
  group_by(tissue_name)%>%
  arrange(desc(abs(shap_mean)))%>%
  filter((1:n())%in%c(1:5))%>%.$feat_name%>%unique()-> keep


summary(k)
DeepLIFT_mean%>%
  filter(feat_name%in%keep)%>%
  spread(tissue_name,shap_mean)->DeepLIFT_mean

# make DeepLIFT score matrix for clustering
Importance_mean_rowName=DeepLIFT_mean$feat_name # TF, RBP 이름들 
DeepLIFT_mean%>%
  dplyr::select(-feat_name,-Type)%>%
  as.matrix() -> Importance_mean_matrix #Tissue 들에 대해서만 값. 
rownames(Importance_mean_matrix)=Importance_mean_rowName
colnames(Importance_mean_matrix) = colnames(Importance_mean_matrix)%>%rename_tissue_bar()


# color code
alphabet_order <- readRDS("C:/Users/seonyeong/Desktop/SystemBiology/team project/simulation/test_Chr1/cluster_color.rds")
tissue_color=alphabet_order$color
names(tissue_color)=alphabet_order$key
tissue_color=tissue_color[match(colnames(Importance_mean_matrix),names(tissue_color))]

# make column and row Annotation
tc = columnAnnotation(Tissue = names(tissue_color),
                      col = list(Tissue = tissue_color),
                      show_annotation_name = FALSE)
ha = rowAnnotation(Type = DeepLIFT_mean$Type,
                   col = list(Type = c("TF" = "#24A277", "RBP" = "#EA9423")))

# generate heatmap
Importance_mean_matrix%>%
  Heatmap(cluster_rows = TRUE,cluster_columns = TRUE,
          right_annotation = ha,
          bottom_annotation = tc,
          clustering_method_rows="ward.D2",
          clustering_method_columns="ward.D2",
          col = g2_colorRamp2(values = seq(from=-0.06,to=0.06,0.12/6)),
          column_names_max_height = unit(20,"cm"),
          name = "DeepLIFT score")->p
show(p)
