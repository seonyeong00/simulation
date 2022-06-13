getwd()
#data
actual = read.csv("C:/Users/seonyeong/Desktop/SystemBiology/team project/simulation/new_data/actual.csv")
actual= actual[,-c(1,2)]
pred = read.csv("C:/Users/seonyeong/Desktop/SystemBiology/team project/simulation/new_data/pred.csv")
pred= pred[,-c(1)]
#library
library(tidyr)
actual %>% gather('Adipose_Subcutaneous', 'Adipose_Visceral_Omentum',
                  'Adrenal_Gland', 'Artery_Aorta', 'Artery_Coronary', 'Artery_Tibial',
                  'Bladder', 'Brain_Amygdala', 'Brain_Anterior_cingulate_cortex_BA24',
                  'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere',
                  'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex_BA9',
                  'Brain_Hippocampus', 'Brain_Hypothalamus',
                  'Brain_Nucleus_accumbens_basal_ganglia', 'Brain_Putamen_basal_ganglia',
                  'Brain_Spinal_cord_cervical_c_1', 'Brain_Substantia_nigra',
                  'Breast_Mammary_Tissue', 'Cells_EBV_transformed_lymphocytes',
                  'Cells_Transformed_fibroblasts', 'Cervix_Ectocervix',
                  'Cervix_Endocervix', 'Colon_Sigmoid', 'Colon_Transverse',
                  'Esophagus_Gastroesophageal_Junction', 'Esophagus_Mucosa',
                  'Esophagus_Muscularis', 'Fallopian_Tube', 'Heart_Atrial_Appendage',
                  'Heart_Left_Ventricle', 'Kidney_Cortex', 'Liver', 'Lung',
                  'Minor_Salivary_Gland', 'Muscle_Skeletal', 'Nerve_Tibial', 'Ovary',
                  'Pancreas', 'Pituitary', 'Prostate', 'Skin_Not_Sun_Exposed_Suprapubic',
                  'Skin_Sun_Exposed_Lower_leg', 'Small_Intestine_Terminal_Ileum',
                  'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina',
                  'Whole_Blood', 'MedianExp', key='tissue', value = 'actual') -> actual


pred %>% gather('Adipose_Subcutaneous', 'Adipose_Visceral_Omentum',
                  'Adrenal_Gland', 'Artery_Aorta', 'Artery_Coronary', 'Artery_Tibial',
                  'Bladder', 'Brain_Amygdala', 'Brain_Anterior_cingulate_cortex_BA24',
                  'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere',
                  'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex_BA9',
                  'Brain_Hippocampus', 'Brain_Hypothalamus',
                  'Brain_Nucleus_accumbens_basal_ganglia', 'Brain_Putamen_basal_ganglia',
                  'Brain_Spinal_cord_cervical_c_1', 'Brain_Substantia_nigra',
                  'Breast_Mammary_Tissue', 'Cells_EBV_transformed_lymphocytes',
                  'Cells_Transformed_fibroblasts', 'Cervix_Ectocervix',
                  'Cervix_Endocervix', 'Colon_Sigmoid', 'Colon_Transverse',
                  'Esophagus_Gastroesophageal_Junction', 'Esophagus_Mucosa',
                  'Esophagus_Muscularis', 'Fallopian_Tube', 'Heart_Atrial_Appendage',
                  'Heart_Left_Ventricle', 'Kidney_Cortex', 'Liver', 'Lung',
                  'Minor_Salivary_Gland', 'Muscle_Skeletal', 'Nerve_Tibial', 'Ovary',
                  'Pancreas', 'Pituitary', 'Prostate', 'Skin_Not_Sun_Exposed_Suprapubic',
                  'Skin_Sun_Exposed_Lower_leg', 'Small_Intestine_Terminal_Ileum',
                  'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina',
                  'Whole_Blood', 'MedianExp', key='tissue', value = 'pred') -> pred
merge = cbind(actual,pred[,3])
colnames(merge)
merge = rename(merge,pred = "pred[, 3]")
merge = rename(merge,gene = "geneID")
saveRDS(merge,"merge.rds")

#library load
options(stringsAsFactors = FALSE)
install.packages("see")
library(see)
install.packages("tidyverse")
library(tidyverse)
install.packages("ggpubr")
library(ggpubr)
install.packages("circlize")
library(circlize)
install.packages("ComplexHeatmap")
library(ComplexHeatmap)
install.packages("dplyr")
library(dplyr)
install.packages("ggplot2")
library(ggplot2)
install.packages("hexbin")
library(hexbin)
source("C:/Users/seonyeong/Desktop/SystemBiology/team project/simulation/test_Chr1/fig_util.R")
ds <- readRDS("C:/Users/seonyeong/Desktop/SystemBiology/team project/simulation/new_data/merge.rds")

# gene-level ###
ds %>%
  filter(tissue%in%"MedianExp")  -> data

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
