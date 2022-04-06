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
res1=compartment_mod(G=g1,beta=0.1,mu=0.2,init_infec = 0.9)
equil_time(res1$infect_share)

equil_time(res1$infect_share)%>%
  ggplot(aes(x=t))+
  geom_line(aes(y=x,color="Infect"))+
  geom_line(aes(y=cum_mean,color="Average"))+
  scale_color_manual(values=c("Infect"="blue","Average"="red"))

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

#write_csv(trial_data,"trial_data_net_size.csv")
trial_data%>%mutate(pred_eq=1-0.2/(0.5*avg_k))%>%
  group_by(net_size,trial)%>%
  summarise(final_share=mean(infect_share[75:100]),
            pred_eq=mean(pred_eq))%>%
  ungroup()%>%
  mutate(dif=pred_eq-final_share)%>%
  ggplot(aes(x=as.factor(net_size),y=dif))+
  geom_boxplot()


trial_data%>%
  group_by(net_size,t)%>%
  summarise(share=mean(infect_share))%>%
  filter(t<15)%>%
  ggplot(aes(x=t,y=share,color=as.factor(net_size)))+
  geom_line()}


################################################################################
# Vary beta, mu, initial infect
################################################################################
trial_data=data.frame()

vary_parm=expand.grid(beta=c(0.01,0.2,0.85),mu=c(0.05,0.45,0.7),init_frac=c(0.2,0.5,0.9))%>%
  as.data.frame()


niter=50

for(k in 1:nrow(vary_parm)){
  
  beta_k=vary_parm$beta[k]
  mu_k=vary_parm$mu[k]
  init_k=vary_parm$init_frac[k]

  for(i in 1:niter ){
    print(paste0(k, "|| Trial: ",i))
    g=igraph::erdos.renyi.game(n=100,p.or.m = 0.1)
    
    components <- igraph::clusters(g, mode="weak")
    biggest_cluster_id <- which.max(components$csize)
    
    # ids
    vert_ids <- V(g)[components$membership == biggest_cluster_id]
    
    # subgraph
    g1=igraph::induced_subgraph(g, vert_ids)
    adj_1=g1%>%as_adjacency_matrix()%>%as.matrix()
    lambda_1=eigen(adj_1)$values%>%max()
    
    avg_k=mean(degree(g1))
    avg_k2=mean(degree(g1)^2)
    
    
    
    trial_data=rbind(trial_data,compartment_mod(G=g1,beta=  beta_k,mu= mu_k,
                                                init_infec =   init_k)%>%
                       mutate(trial=i,beta=  beta_k,mu= mu_k,
                              init_infec =   init_k, lambda_1= lambda_1))
    
  }
}

write_csv( trial_data,"trial_data_param_vary.csv")

trial_data%>%
  filter(t<25)%>%
  group_by(t,beta,mu,init_infec)%>%
  summarise(avg=mean(infect_share,na.rm=T))%>%
  mutate(inital_frac=paste0("Initial Infection: ",init_infec),
         mu_text=paste0("mu=",mu),
         beta_text=paste0("beta=",beta),
         prama=paste0(beta,",",mu))%>%
  ggplot(aes(x=t,y=avg,color=as.factor(init_infec)))+
  geom_line()+
  facet_wrap(~ beta_text+mu_text)+
  labs(x=TeX("Time $$(t)$$"),y=TeX("Infected Share $$i(t)$$",),
       color=TeX("Initial infection ($i_0$)"),
       title = TeX("Network SIS Dynamics varying $\\beta, \\mu$ and $i_0$"),
       subtitle = TeX("ER network with $N=100$ and $\\langle k \\rangle =10$"))+
  theme_bw()+
  scale_color_viridis_d()

################################################################################
# Time to equilibrium 
################################################################################

delta=10^(-3)

trial_data%>%
  group_by(t,beta,mu,init_infec)%>%
  summarise(avg=mean(infect_share,na.rm=T))%>%
  ungroup()%>%
  arrange(t)%>%
  group_by(beta,mu,init_infec)%>%
  mutate(cum_avg=cummean(avg),
         avg_dif=cum_avg-lag(cum_avg,1))%>%
  filter(is.na( avg_dif)==F)%>%
  summarise(t_final=t[min(which(abs(avg_dif)<delta),na.rm = T)],
            infect=cum_avg[min(which(abs(avg_dif)<delta),na.rm = T)])%>%
  mutate(inital_frac=paste0("Initial Infection: ",init_infec),
         mu_text=paste0("mu=",mu),
         beta_text=paste0("beta=",beta))%>%
  ggplot(aes(x=t_final,y=infect,color=as.factor(init_infec)))+
  geom_point(size=2)+
  facet_wrap(~beta_text+mu_text)+
  labs(x=TeX("Equilibrium Time ($t^*$)"),y=TeX("Equilibrium  Share $m(t^*)$",),
       color=TeX("Initial infection ($i_0$)"),
       title = TeX("Network SIS Equilibrium varying $\\beta, \\mu$ and $i_0$"),
       subtitle = TeX("ER network with $N=100$ and $\\langle k \\rangle =10$"))+
  theme_bw()+
  scale_color_viridis_d()

################################################################################
# Scale Free Graph
################################################################################

scale_free_graph=sample_fitness_pl(
  no.of.nodes=100,
  no.of.edges=1000,
  exponent.out=2.5)

g=igraph::erdos.renyi.game(n=100,p.or.m = 0.2)
# Extracting Largest Component 
components <- igraph::clusters(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(g)[components$membership == biggest_cluster_id]
g1=igraph::induced_subgraph(g, vert_ids)
mean(degree(g1))

lambda_seq=data.frame(lambda=seq(0.01,0.40,length.out=40))%>%
  mutate(beta=0.01,
         mu=beta/lambda)



scale_data=data.frame()
er_data=data.frame()


niter=50
a1=tictoc::tic()
for(k in 1:nrow(lambda_seq)){
  
  beta_k=lambda_seq$beta[k]
  mu_k=lambda_seq$mu[k]
  init_k=0.1
  
  for(i in 1:niter ){
    print(paste0(k, "|| Trial: ",i))
    g=igraph::erdos.renyi.game(n=100,p.or.m = 0.02)
    
    components <- igraph::clusters(g, mode="weak")
    biggest_cluster_id <- which.max(components$csize)
    
    # ids
    vert_ids <- V(g)[components$membership == biggest_cluster_id]
    
    # subgraph
    g1=igraph::induced_subgraph(g, vert_ids)
    adj_1=g1%>%as_adjacency_matrix()%>%as.matrix()
    lambda_1=eigen(adj_1)$values%>%max()
    
    avg_k=mean(degree(g1))
    avg_k2=mean(degree(g1)^2)
    
    #avg_k/avg_k2
    
    
    
    er_data=rbind(er_data,compartment_mod(G=g1,beta=  beta_k,mu= mu_k,
                                                init_infec =   init_k)%>%
                       mutate(trial=i,beta=  beta_k,mu= mu_k,
                              init_infec =   init_k, lambda_1= lambda_1))
    
    
    scale_free_graph=barabasi.game(100, power=2.5, directed=F)
    
    adj_1=scale_free_graph%>%as_adjacency_matrix()%>%as.matrix()
    lambda_1=eigen(adj_1)$values%>%max()
    
    avg_k=mean(degree(scale_free_graph))
    avg_k2=mean(degree(scale_free_graph)^2)
    
    scale_data=rbind(scale_data,compartment_mod(G=scale_free_graph,beta=  beta_k,mu= mu_k,
                    init_infec =   init_k)%>%
      mutate(trial=i,beta=  beta_k,mu= mu_k,
             init_infec =   init_k, lambda_1= lambda_1))
    
  }
}
a2=tictoc::toc()

write_csv(er_data,"lambda_er_data.csv")
write_csv(scale_data,"lambda_scale_data.csv")


compare_df=rbind(er_data%>%mutate(model="ER"),
                 scale_data%>%mutate(model="Scale-Free"))%>%
  group_by(model,beta,mu)%>%
  filter(t>=80)%>%
  summarise(infect=mean(infect_share,na.rm=T))%>%
  ungroup()%>%
  mutate(lambda=beta/mu)

compare_df%>%
  ggplot(aes(x=lambda,y=infect,color=model))+
  geom_line()+
  labs(x=TeX("$\\lambda$"),y=TeX("$i(\\infty)$"),color="Model")+
  ggthemes::scale_color_calc()+
  theme_classic()+
  theme(legend.position = "top")

