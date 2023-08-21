rm(list=ls())
library(ggplot2)
library(patchwork)
###load the procedures
getwd()
source("Exhaustive_Procedures.R")

###Gaussian testing problem for exact p-values
m=2000    #Number of Trials
n=1000    #Number of Hypotheses per Trial
mu_A=4    #Strength of the alternative
mu_N=0    #Conservativeness of null p-values (<0 for conservative null p-values)
#pi_A is defined in the loop below

###Initialise Hyperparameters
seq_1_n=seq(1,n,1)
gamma=(1/((seq_1_n+1)*log(seq_1_n+1)^2))/2.06227      #2.06227 is the approximated value such that the series equals 1
alpha=0.2
tau=0.8
lambda=0.16

w=abs(matrix(1:n-1 , nrow = n, ncol = n, byrow = TRUE) - (1:n-1))
w[w==0]=1
w=matrix(gamma[w],n,n)
w[upper.tri(w)==0]=0


###Predefine vectors for FWER and power of the different procedures
FWER_Alpha=rep(0,9)
power_Alpha=rep(0,9)
FWER_ADDIS_Graph=rep(0,9)
power_ADDIS_Graph=rep(0,9)
FWER_I_ADDIS_Graph=rep(0,9)
power_I_ADDIS_Graph=rep(0,9)
###Set seed to make the results reproducible
set.seed(12345)

###Generate p-values and compute FWER and power for the desired procedures
for(l in 1:9){
  pi_A=l/10
  p=matrix(,nrow=n,ncol=m)
  hypo=matrix(,nrow=n,ncol=m)
  for(j in 1:m){
    hypo[,j]=rbinom(n,1,pi_A)
    X=rnorm(n)
    Z=mu_N*(hypo[,j]-1)*(-1)+mu_A*hypo[,j]+X
    p[,j]=pnorm(-Z)
  }
  
  ##Alpha-Spending
  V=rep(0,m)      #Indicates, whether there was atleast one type 1 error in a trial
  power=rep(0,m)  #Power within each trial
  for(j in 1:m){
    alpha_ind=alpha_spending(alpha, gamma, p[,j],n)
    hypo_est=alpha_ind>=p[,j]
    V[j]=max((hypo[,j]==0 & hypo_est==1))
    D=(hypo[,j]==1 & hypo[,j]==hypo_est)
    power[j]=sum(D)/sum(hypo[,j])
  }
  FWER_Alpha[l]=mean(V)
  power_Alpha[l]=mean(power)
  
  ##ADDIS-Graph
  V=rep(0,m)      #Indicates, whether there was atleast one type 1 error in a trial
  power=rep(0,m)  #Power within each trial
  
  for(j in 1:m){
    alpha_ind=ADDIS_Graph(alpha, gamma,w, tau, lambda,0, p[,j], n)
    hypo_est=alpha_ind>=p[,j]
    V[j]=max((hypo[,j]==0 & hypo_est==1))
    D=(hypo[,j]==1 & hypo[,j]==hypo_est)
    power[j]=sum(D)/sum(hypo[,j])
  }
  FWER_ADDIS_Graph[l]=mean(V)
  power_ADDIS_Graph[l]=mean(power)
  
  
  ##Evenly-Improved-ADDIS-Graph
  V=rep(0,m)      #Indicates, whether there was atleast one type 1 error in a trial
  power=rep(0,m)  #Power within each trial
  
  for(j in 1:m){
    alpha_ind=EI_ADDIS_Graph(alpha, gamma,w, tau,lambda, p[,j], n)
    hypo_est=alpha_ind>=p[,j]
    V[j]=max((hypo[,j]==0 & hypo_est==1))
    D=(hypo[,j]==1 & hypo[,j]==hypo_est)
    power[j]=sum(D)/sum(hypo[,j])
  }
  FWER_I_ADDIS_Graph[l]=mean(V)
  power_I_ADDIS_Graph[l]=mean(power)
  
  
}



###Create Plot for conservative p-values
cols <- c("limegreen",  "#f84f4f","cornflowerblue")
Sys.setlocale('LC_CTYPE', 'greek')

results_df=data.frame(seq(0.1,0.9,0.1), power_Alpha, FWER_Alpha, power_ADDIS_Graph,FWER_ADDIS_Graph,power_I_ADDIS_Graph,FWER_I_ADDIS_Graph)

p1=ggplot(results_df, aes(seq.0.1..0.9..0.1.)) + 
  geom_line(aes(y = power_Alpha, colour = "Alpha-Spending")) + 
  geom_point(aes(y = power_Alpha, colour = "Alpha-Spending", shape = "Alpha-Spending")) +
  geom_line(aes(y = FWER_Alpha, colour = "Alpha-Spending"))+
  geom_point(aes(y = FWER_Alpha, colour = "Alpha-Spending", shape = "Alpha-Spending"))+
  geom_line(aes(y = power_ADDIS_Graph, colour = "ADDIS-Graph")) + 
  geom_point(aes(y = power_ADDIS_Graph, colour = "ADDIS-Graph", shape = "ADDIS-Graph")) +
  geom_line(aes(y = FWER_ADDIS_Graph, colour = "ADDIS-Graph"))+
  geom_point(aes(y = FWER_ADDIS_Graph, colour = "ADDIS-Graph", shape = "ADDIS-Graph"))+
  geom_line(aes(y = power_I_ADDIS_Graph, colour = "EI-ADDIS-Graph")) + 
  geom_point(aes(y = power_I_ADDIS_Graph, colour = "EI-ADDIS-Graph", shape = "EI-ADDIS-Graph")) +
  geom_line(aes(y = FWER_I_ADDIS_Graph, colour = "EI-ADDIS-Graph"))+
  geom_point(aes(y = FWER_I_ADDIS_Graph, colour = "EI-ADDIS-Graph", shape = "EI-ADDIS-Graph"))+
  geom_hline(yintercept=alpha)+
  scale_shape_manual(name  ="Procedure",values=c("Alpha-Spending"=8,"ADDIS-Graph"=16,"EI-ADDIS-Graph"=17))+
  scale_color_manual(name  ="Procedure", values = c("Alpha-Spending"=cols[1],"ADDIS-Graph"=cols[2],"EI-ADDIS-Graph"=cols[3]))+
  xlab(expression(pi[A]))+
  ylab("FWER / Power")+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2), limits=c(0.1,0.9),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0,1.2,0.2), limits=c(0,1),expand = c(0, 0))+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5,linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

###Gaussian testing problem for conservative p-values

mu_N=-2    #Conservativeness of null p-values (<0 for conservative null p-values)

###Predefine vectors for FWER and power of the different procedures
FWER_Alpha=rep(0,9)
power_Alpha=rep(0,9)
FWER_ADDIS_Graph=rep(0,9)
power_ADDIS_Graph=rep(0,9)
FWER_I_ADDIS_Graph=rep(0,9)
power_I_ADDIS_Graph=rep(0,9)
###Set seed to make the results reproducible
set.seed(12345)

###Generate p-values and compute FWER and power for the desired procedures
for(l in 1:9){
  pi_A=l/10
  p=matrix(,nrow=n,ncol=m)
  hypo=matrix(,nrow=n,ncol=m)
  for(j in 1:m){
    hypo[,j]=rbinom(n,1,pi_A)
    X=rnorm(n)
    Z=mu_N*(hypo[,j]-1)*(-1)+mu_A*hypo[,j]+X
    p[,j]=pnorm(-Z)
  }
   
   ##Alpha-Spending
   V=rep(0,m)      #Indicates, whether there was atleast one type 1 error in a trial
   power=rep(0,m)  #Power within each trial
   for(j in 1:m){
     alpha_ind=alpha_spending(alpha, gamma, p[,j],n)
     hypo_est=alpha_ind>=p[,j]
     V[j]=max((hypo[,j]==0 & hypo_est==1))
     D=(hypo[,j]==1 & hypo[,j]==hypo_est)
     power[j]=sum(D)/sum(hypo[,j])
   }
   FWER_Alpha[l]=mean(V)
   power_Alpha[l]=mean(power)
  
  ##ADDIS-Graph
  V=rep(0,m)      #Indicates, whether there was atleast one type 1 error in a trial
  power=rep(0,m)  #Power within each trial
  
  for(j in 1:m){
    alpha_ind=ADDIS_Graph(alpha, gamma,w, tau, lambda,0, p[,j], n)
    hypo_est=alpha_ind>=p[,j]
    V[j]=max((hypo[,j]==0 & hypo_est==1))
    D=(hypo[,j]==1 & hypo[,j]==hypo_est)
    power[j]=sum(D)/sum(hypo[,j])
  }
  FWER_ADDIS_Graph[l]=mean(V)
  power_ADDIS_Graph[l]=mean(power)
  
  
  ##Evenly-Improved-ADDIS-Graph
  V=rep(0,m)      #Indicates, whether there was atleast one type 1 error in a trial
  power=rep(0,m)  #Power within each trial
  
  for(j in 1:m){
    alpha_ind=EI_ADDIS_Graph(alpha, gamma,w, tau,lambda, p[,j], n)
    hypo_est=alpha_ind>=p[,j]
    V[j]=max((hypo[,j]==0 & hypo_est==1))
    D=(hypo[,j]==1 & hypo[,j]==hypo_est)
    power[j]=sum(D)/sum(hypo[,j])
  }
  FWER_I_ADDIS_Graph[l]=mean(V)
  power_I_ADDIS_Graph[l]=mean(power)
  
  
}



###Create Plot for conservative p-values
cols <- c("limegreen",  "#f84f4f","cornflowerblue")
Sys.setlocale('LC_CTYPE', 'greek')

results_df=data.frame(seq(0.1,0.9,0.1), power_Alpha, FWER_Alpha, power_ADDIS_Graph,FWER_ADDIS_Graph,power_I_ADDIS_Graph,FWER_I_ADDIS_Graph)

p2=ggplot(results_df, aes(seq.0.1..0.9..0.1.)) + 
  geom_line(aes(y = power_Alpha, colour = "Alpha-Spending")) + 
  geom_point(aes(y = power_Alpha, colour = "Alpha-Spending", shape = "Alpha-Spending")) +
  geom_line(aes(y = FWER_Alpha, colour = "Alpha-Spending"))+
  geom_point(aes(y = FWER_Alpha, colour = "Alpha-Spending", shape = "Alpha-Spending"))+
  geom_line(aes(y = power_ADDIS_Graph, colour = "ADDIS-Graph")) + 
  geom_point(aes(y = power_ADDIS_Graph, colour = "ADDIS-Graph", shape = "ADDIS-Graph")) +
  geom_line(aes(y = FWER_ADDIS_Graph, colour = "ADDIS-Graph"))+
  geom_point(aes(y = FWER_ADDIS_Graph, colour = "ADDIS-Graph", shape = "ADDIS-Graph"))+
  geom_line(aes(y = power_I_ADDIS_Graph, colour = "EI-ADDIS-Graph")) + 
  geom_point(aes(y = power_I_ADDIS_Graph, colour = "EI-ADDIS-Graph", shape = "EI-ADDIS-Graph")) +
  geom_line(aes(y = FWER_I_ADDIS_Graph, colour = "EI-ADDIS-Graph"))+
  geom_point(aes(y = FWER_I_ADDIS_Graph, colour = "EI-ADDIS-Graph", shape = "EI-ADDIS-Graph"))+
  geom_hline(yintercept=alpha)+
  scale_shape_manual(name  ="Procedure",values=c("Alpha-Spending"=8,"ADDIS-Graph"=16,"EI-ADDIS-Graph"=17))+
  scale_color_manual(name  ="Procedure", values = c("Alpha-Spending"=cols[1],"ADDIS-Graph"=cols[2],"EI-ADDIS-Graph"=cols[3]))+
  xlab(expression(pi[A]))+
  ylab("FWER / Power")+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2), limits=c(0.1,0.9),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0,1.2,0.2), limits=c(0,1),expand = c(0, 0))+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5,linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

###Create combined plot
combined <- p1 + p2 & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")
ggsave("Plot_EI_ADDIS_Graph_logq.pdf")