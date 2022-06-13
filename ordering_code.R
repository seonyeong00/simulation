# data load
DeepLIFT <- readRDS("C:/Users/seonyeong/Desktop/SystemBiology/team project/simulation/test_Chr1/DeepLIFT_mean.rds")
# library load
library("dplyr")
library(ggplot2)

# Type 종류 확인
levels(as.factor(DeepLIFT$Type)) # "miRNA" "other" "RBP"   "TF"  
#Liver 의 TF 개수 파악, 725개
DeepLIFT %>%
  filter(Type == 'TF'& tissue_name == 'Liver') -> TF
#Liver 의 RBP 개수 파악, 170개
DeepLIFT %>%
  filter(Type == 'RBP'& tissue_name == 'Liver') -> RBP
#Liver 의 miRNA 개수 파악, 201개
DeepLIFT %>%
  filter(Type == 'miRNA'& tissue_name == 'Liver') -> miRNA
#Liver 의 other 개수 파악, 1개
DeepLIFT %>%
  filter(Type == 'other'& tissue_name == 'Liver') -> other

# tissue별 df 생성
#Liver
DeepLIFT %>%
  filter(tissue_name == 'Liver') -> liver_df
#Adipose_Subcutaneous
DeepLIFT %>%
  filter(tissue_name == 'Adipose_Subcutaneous') -> Adipose_Subcutaneous_df

# DeepLIFT mean 기준으로 ordering
liver_df <- liver_df[c(order(-liver_df$shap_mean)),]
# top10 추출
top15_liver <- head(liver_df,15)

# ggplot 이용 bar graph 그리기
ggplot(top15_liver,aes(x=reorder(feat_name, -shap_mean),y=shap_mean,fill=Type))+
  geom_bar(stat="identity")+
  labs(title = "DeepLIFT score Top 15 in Liver tissue",x = 'gene', y='DeepLIFT score')

# HNF4A gene에 대해서 .. 
DeepLIFT %>%
  filter(feat_name == 'HNF4A') ->HNF4A_df
# DeepLIFT mean 기준으로 ordering
HNF4A_df <- HNF4A_df[c(order(-HNF4A_df$shap_mean)),]
# 상위 10개 tissue 추출
top10_HNF4A_df <- head(HNF4A_df,10)

# ggplot 
ggplot(top10_HNF4A_df,aes(x=reorder(tissue_name, -shap_mean),y=shap_mean,fill=tissue_name))+
  geom_bar(stat="identity")+
  labs(title = "DeepLIFT score Top 10 in HNF4A gene",x = 'Tissue', y='DeepLIFT score')

# RXRA gene에 대해서 .. 
DeepLIFT %>%
  filter(feat_name == 'RXRA') ->RXRA_df
# DeepLIFT mean 기준으로 ordering
RXRA_df <- RXRA_df[c(order(-RXRA_df$shap_mean)),]
# 상위 10개 tissue 추출
top10_RXRA_df <- head(RXRA_df,10)

# ggplot 
ggplot(top10_RXRA_df,aes(x=reorder(tissue_name, -shap_mean),y=shap_mean,fill=tissue_name))+
  geom_bar(stat="identity")+
  labs(title = "DeepLIFT score Top 10 in RXRA gene",x = 'Tissue', y='DeepLIFT score')

# ATF3 gene에 대해서 .. 
DeepLIFT %>%
  filter(feat_name == 'ATF3') ->ATF3_df
# DeepLIFT mean 기준으로 ordering
ATF3_df <- ATF3_df[c(order(-ATF3_df$shap_mean)),]
# 상위 10개 tissue 추출
top10_ATF3_df <- head(ATF3_df,10)

# ggplot 
ggplot(top10_ATF3_df,aes(x=reorder(tissue_name, -shap_mean),y=shap_mean,fill=tissue_name))+
  geom_bar(stat="identity")+
  labs(title = "DeepLIFT score Top 10 in ATF3 gene",x = 'Tissue', y='DeepLIFT score')
