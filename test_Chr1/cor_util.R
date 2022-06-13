corAndPvalue_to_list<-function(res){
  res$cor%>%t()%>%as.data.frame()%>%
    dplyr::mutate(Var2=rownames(.))%>%
    gather(Var1,cor,-Var2) -> res_gather
  
  for(stat_name in c("p","Z","t","nObs")){
    res[[stat_name]]%>%t()%>%as.data.frame()%>%
      mutate(Var2=rownames(.))%>%
      gather(Var1,value,-Var2) -> res_gather_tmp
    check_match=all.equal(res_gather$Var2,res_gather_tmp$Var2) & all.equal(res_gather$Var1,res_gather_tmp$Var1)
    if (!check_match){stop()}
    res_gather[[stat_name]]=res_gather_tmp$value
  }
  return(res_gather)
}