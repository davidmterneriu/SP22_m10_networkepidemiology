#SIS Functions 

#rm(list=ls())

library(tidyverse)
library(ggplot2)
library(igraph)

infect_share=function(G){
  return(mean(V(G)$infect))
}


compartment_mod=function(G,beta,mu,init_infec=0.1,tmax=100){
  #browser()
  n_nodes=length(V(G))
  G=G%>%set_vertex_attr(name="infect",value=0)
  infect_nodes=round(init_infec*n_nodes)
  node_ids_init_infect=sample(V(G),size = infect_nodes)%>%as.numeric()
  
  G=G%>%set_vertex_attr(name="infect",node_ids_init_infect,value=1)
  G=G%>%set_vertex_attr(name="infect2",value=NA)
  
  
  k_mean=mean(degree(G))
  k2_mean=mean(degree(G)^2)
  ratio1=(k2_mean-k_mean)/k_mean
  r0_1=beta/mu
  r0_real=r0_1*ratio1

  
  infect_share_df=data.frame(t=seq(0,tmax),infect_share=NA)
  infect_share_df$infect_share[1]<-infect_share(G)
  
  for (t in 2:(tmax+1)){
    current_infect=which(V(G)$infect==1)
    current_susp=which(V(G)$infect==0)
    #Simulating Recovery 
    recover_vals=as.numeric((rbinom(n=length(current_infect),size = 1,prob = mu)==0))
    G=G%>%set_vertex_attr(name="infect2",current_infect,value=recover_vals)
    num_infect_neigh=sapply(X=current_susp,FUN = function(x)sum(V(G)$infect[neighbors(G,v=x)]))
    
    no_infect_neigh=num_infect_neigh==0
    new_infect=as.numeric(runif(n=length(current_susp))<(1 - (1 - beta)^num_infect_neigh))
    new_infect[no_infect_neigh]<-0
    
    G=G%>%set_vertex_attr(name="infect2",current_susp,value=new_infect)
    
    #Initializing next iteration
    V(G)$infect=V(G)$infect2
    V(G)$infect2=NA
    
    infect_share_df$infect_share[t]=infect_share(G)
    
    if(infect_share(G)==1){
      break
    }
    
    #print(paste0("Time period: ", t, "|| Share: ",round(infect_share(G),3)))
  }

  return(infect_share_df)
 
}

infect_path=function(G,beta,mu,tmax,init_infec,kmean=NA,t_opt=NA){
  if(is.na(kmean)){
    k_mean=mean(degree(G))
  }else{
    k_mean=kmean
  }
 
  beta_k=k_mean*beta
  C=init_infec/(1-init_infec-mu/beta_k)
  q1=(1-mu/beta_k)
  q2_fun=function(t){
    res=(C*exp((beta_k-mu)*t))/(1+C*exp((beta_k-mu)*t))
    return(res)
  }
  
  if(is.na(t_opt)){
    res_df=data.frame(t=0:tmax,infect_share=NA)
    res_df$infect_share=q1*q2_fun(t=res_df$t)
    return(res_df)
  }else{
    return(q1*q2_fun(t=t_opt))
  }
  

}

equil_outcome=function(beta,mu,avg_k){
  return(1-mu/(beta*avg_k))
}
