#SIS Demo

rm(list=ls())

library(tidyverse)
library(ggplot2)
library(latex2exp)
library(metR)
library(igraph)

source("sis_functions.R")

# Generate Random Network
g=igraph::erdos.renyi.game(n=100,p.or.m = 0.1)
# Extracting Largest Component 
components <- igraph::clusters(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(g)[components$membership == biggest_cluster_id]
g1=igraph::induced_subgraph(g, vert_ids)

#Simple test of simulation model and infection path
res1=compartment_mod(G=g1,beta=0.05,mu=0.1,init_infec = 0.9)
res1a=infect_path(G=g1,beta=0.05,mu=0.1,init_infec = 0.9,tmax=100)


#Varying Initial Conditions 

trial_data=data.frame()
initial_infect_vect=seq(0.05,0.95,by=0.05)
niter=50

for(k in 1:length(initial_infect_vect)){
  i0=initial_infect_vect[k]
  
  for(i in 1:niter ){
    print(paste0("Initial infect: ",i0, "|| Trial: ",i))
    g=igraph::erdos.renyi.game(n=100,p.or.m = 0.1)
    
    components <- igraph::clusters(g, mode="weak")
    biggest_cluster_id <- which.max(components$csize)
    
    # ids
    vert_ids <- V(g)[components$membership == biggest_cluster_id]
    
    # subgraph
    g1=igraph::induced_subgraph(g, vert_ids)
    trial_data=rbind(trial_data,compartment_mod(G=g1,beta=0.5,mu=0.2,init_infec =  i0)%>%
                       mutate(trial=i,init_infec =  i0))
    
  }
}

trial_data%>%
  group_by(t,init_infec)%>%
  filter(t<11)%>%
  summarise(infect_share=mean(infect_share))%>%
  ggplot(aes(x=t,y=infect_share,color=as.factor(init_infec)))+
  geom_line()+
  scale_color_viridis_d()+
  theme_bw()+
  scale_x_continuous(n.breaks = 10)+
  scale_y_continuous(breaks = seq(0.05,0.95,by=0.1))+
  labs(x=TeX("Time$(t)$"),y=TeX("$i(t)$"),color=TeX("$i_0$"),
       caption = TeX("Note: Each parameterization of $i_0$ is simulated 50 times."),
       title = "SIS Model Dynamics: Infection Path to Equilibrium",
       subtitle = TeX("$\\beta=0.5, \\mu=0.2, n=100,$ and $\\langle k \\rangle \\approx 10$"))+
  theme(legend.position = "bottom")+
  guides(color=guide_legend(ncol=10))


param_grid=expand.grid(beta=seq(0.01,0.99,by=0.01),mu=seq(0.01,0.99,by=0.01))%>%
  as.data.frame()%>%
  mutate(R0=beta*10/mu,
         eq_val=ifelse(R0>=1,1-1/R0,0))

param_grid%>%
  ggplot(aes(x=beta,y=mu,fill=eq_val,z=eq_val))+
  geom_raster(interpolate = T)+
  scale_fill_viridis_c()+
  geom_contour(colour = "white",)+
  geom_text_contour(label.placer = label_placer_fraction(0.8),stroke = 0.1,
                    rotate = FALSE,size=3)+
  theme_bw()+
  scale_x_continuous(n.breaks = 10)+
  scale_y_continuous(n.breaks = 10)+
  labs(x=TeX("$\\beta$"),
       y=TeX("$\\mu$"),
       fill=TeX("$i(\\infty)$"),
       title = TeX("SIS Equilbrium $i(\\infty)$ as a function of $\\beta, \\mu$, and $\\langle k \\rangle =10$"),
       caption =TeX("Note that: $i(\\infty)=1-1/R_0$ when $R_0=\\frac{\\beta \\langle k \\rangle}{\\mu}\\geq 1, 0$ otherwise "))



if(FALSE){

sample_net_size=10^seq(2,4,length.out=20)%>%round()


trial_data=data.frame()
#initial_infect_vect=seq(0.05,0.95,by=0.05)
niter=50

for(k in 1:length(sample_net_size)){
  net_size=sample_net_size[k]
  
  for(i in 1:niter ){
    print(paste0("Net size: ",net_size, "|| Trial: ",i))
    g=igraph::erdos.renyi.game(n=net_size,p.or.m = 0.1)
    
    components <- igraph::clusters(g, mode="weak")
    biggest_cluster_id <- which.max(components$csize)
    
    # ids
    vert_ids <- V(g)[components$membership == biggest_cluster_id]
    
    # subgraph
    g1=igraph::induced_subgraph(g, vert_ids)
    
    avg_k=mean(degree(g1))
    trial_data=rbind(trial_data,compartment_mod(G=g1,beta=0.5,mu=0.2,init_infec =  0.1)%>%
                       mutate(trial=i,avg_k=avg_k,net_size=net_size))
    
  }
}


trial_data%>%mutate(pred_eq=1-0.2/(0.5*avg_k))%>%
  group_by(net_size,trial)%>%
  summarise(final_share=mean(infect_share[75:100]),
            pred_eq=mean(pred_eq))%>%
  ungroup()%>%
  mutate(dif=pred_eq-final_share)%>%
  ggplot(aes(x=as.factor(net_size),y=dif))+
  geom_boxplot()}


