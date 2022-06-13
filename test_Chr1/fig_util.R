

rename_tissue_dot=function(tissue){
  library(tidyverse)
  tissue_color = read_delim("functions/GTEX_tissue_color.txt",delim = "\t",col_names = F)%>%
    mutate(Name_dot=
             X2%>%make.names())
  
  tissue_color = read_delim("functions/GTEX_tissue_color.txt",delim = "\t",col_names = F)%>%
    mutate(Name_dot=
             X2%>%make.names()%>%gsub("\\.$","",.))%>%
    bind_rows(.,tissue_color)
  
  indx=match(tissue,tissue_color$Name_dot)
  tissue[!is.na(indx)]=tissue_color$X2[indx[!is.na(indx)]]
  tissue=replace(tissue,tissue=="MedianExp","Median expression")
    return(tissue)
}

rename_tissue_bar=function(tissue){
  library(tidyverse)
  
  tissue_color = read_delim("functions/GTEX_tissue_color.txt",delim = "\t",col_names = F)%>%
    mutate(Name_bar=
             X2%>%
             make.names()%>%
             gsub("\\.$","",.)%>%
             gsub("\\.\\.\\.","_",.)%>%
             gsub("\\.\\.","_",.)%>%
             gsub("\\.","_",.))
  indx=match(tissue,tissue_color$Name_bar)
  tissue[!is.na(indx)]=tissue_color$X2[indx[!is.na(indx)]]
  tissue=replace(tissue,tissue=="MedianExp","Median expression")
  return(tissue)
}

get_tissue_color=function(tissue){
  library(tidyverse)
  cor_to_brain <- readRDS("functions/cluster_color.rds")
  #tissue_color = read_delim("~/OneDrive/Project/DEcode/Data/GTEX_phenotype/GTEX_tissue_color.txt",delim = "\t",col_names = F)
  indx=match(tissue,cor_to_brain$key)
  return(cor_to_brain$color[indx])
}


library(ggplot2)
d5_color<-scale_color_manual(values = c("#447cff","#ff9544"))


d5_fill<-function(label=NULL){
  scale_fill_manual(name=label,values = c("#447cff","#ff9544"))%>%
    return()
}


d2_fill<-function(label=NULL){
  scale_fill_manual(name=label,values = c("#EA9423","#24A277","#767171"))%>%
    return()
}


d4_fill<-function(label=NULL){
  scale_fill_manual(name=label,values = c("#E74D3D","#2980BA","#767171"))%>%
    return()
}

d6_fill<-function(label=NULL){
  scale_fill_manual(name=label,values = c("#193065","#EA9423","#24A277","#767171"))%>%
    return()
}

d3_fill<-function(label=NULL){
  scale_fill_manual(name=label,values = c("#24A277","#EA9423","#A15BAB"))%>%
    return()
}

d3_color<-function(label=NULL){
  scale_color_manual(name=label,values = c("#24A277","#EA9423","#A15BAB"))%>%
    return()
}


g2_fill=function(label=NULL,limit=2.8){
  scale_fill_gradientn(name=label,
                       values = seq(0,1,by = 1/7),
                       colors = c("#193065","#00539C","#53B0D8",
                                  "#FFFFFF",
                                  "#FEAB92","#E00016","#730009"),
                       limits=c(-limit,limit))%>%
    return()
}


g2_colorRamp2=function(values=c(-0.05,-0.033,-0.016,0,0.016,0.033,0.05)){
  colorRamp2(values,
             c("#193065","#00539C","#53B0D8",
               "#FFFFFF",
               "#FEAB92","#E00016","#730009"))%>%return()
}



