#SIS Demo

rm(list=ls())

library(tidyverse)
library(ggplot2)
library(igraph)

source("sis_functions.R")

g=igraph::erdos.renyi.game(n=100,p.or.m = 0.1)

components <- igraph::clusters(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)

# ids
vert_ids <- V(g)[components$membership == biggest_cluster_id]

# subgraph
g1=igraph::induced_subgraph(g, vert_ids)

trial_data=data.frame()
niter=100

for(i in 1:niter ){
  print(i)
  trial_data=rbind(trial_data,compartment_mod(G=g1,beta=0.5,mu=0.2)%>%
    mutate(trial=i))
  
}


trial_data%>%
  group_by(t)%>%
  summarise(infect_share=mean(infect_share))%>%
  ggplot(aes(x=t,y=infect_share))+
  geom_line()

1-(0.2/(0.5*mean(degree(g1))))
