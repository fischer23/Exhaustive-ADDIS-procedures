###This R-file contains the procedures that were used for the simulations and
###real data application in the paper "Exhaustive ADDIS procedures for online FWER control"

#Input (at most): alpha, gamma, w, tau, lambda,lags, p, n. 

#alpha: overall significance level. Number between 0 and 1
#gamma: Initial weights for all procedures. Non-negative n-dimensional vector with sum less than 1.
#w:     Weights for graphical procedures. Upper triangle n x n matrix with row sum less than 1
#tau:   Used for ADDIS procedures. n-dimensional vector with values between 0 and 1. Besides, it can be chosen as a fixed number between 0 and 1.
#lambda:Used for the ADDIS procedures. n-dimensional vector with values in [0, tau_i) for ADDIS procedures or [tau_i*alpha^(i), tau_i) for exhaustive adaptive procedures. Besides, it can be chosen as a fixed number fulfilling these conditions.
#lags:  Representing the given local dependence structure. n-dimensional vector of natural numbers (including 0) or a fix natural number (including 0). Choose lags=0 for independence.
#p:     n-dimensional vector of p-values. 
#n:     Number of hypotheses.

#Output: Vector with length n of individual significance levels.

#ADDIS-Spending under local dependence
ADDIS_Spending=function(alpha, gamma, tau, lambda,lags, p, n){
  if(length(lags)== 1){lags=rep(lags,n)}
  if(length(tau)== 1){tau=rep(tau,n)}
  if(length(lambda)== 1){lambda=rep(lambda,n)}
  if(length(gamma)!=n | length(tau)!= n | length(lambda)!=n| length(p)!=n| length(lags)!=n){
    warning("mismatching length")
  }
  S_C=rep(0,n)
  t=rep(1,n)
  alpha_ind=rep(0,n)
  for(i in 1:n){
    if(lags[i]>=i-1){
      t[i]=i
    }else{
      t[i]=1+lags[i]+sum(S_C[1:(i-lags[i]-1)])
    }
    alpha_ind[i]=alpha*gamma[t[i]]*(tau[i]-lambda[i])
    if(p[i]<=tau[i] & p[i]>lambda[i]){
      S_C[i]=1
    }
  }
  return(alpha_ind)
}


 #ADDIS-Graph under local dependence
 ADDIS_Graph=function(alpha, gamma, w, tau, lambda,lags, p, n){
   if(length(lags)== 1){lags=rep(lags,n)}
   if(length(tau)== 1){tau=rep(tau,n)}
   if(length(lambda)== 1){lambda=rep(lambda,n)}
   if(length(gamma)!=n | length(tau)!= n | length(lambda)!=n| length(p)!=n| length(lags)!=n){
     warning("mismatching length")
   }
   C_S=rep(0,n)
   alpha_ind=rep(0,n)
   alpha_ind[1]=alpha*gamma[1]*(tau[1]-lambda[1])
   if(p[1]>tau[1] | p[1]<=lambda[1]){
     C_S[1]=1
   }
   for(i in 2:n){
     if(lags[i]>=i-1){
       alpha_ind[i]=alpha*gamma[i]*(tau[i]-lambda[i])
     }else{
       alpha_ind[i]=(alpha*gamma[i]+sum(C_S[1:(i-lags[i]-1)]*alpha_ind[1:(i-lags[i]-1)]*w[1:(i-lags[i]-1),i]/(tau[1:(i-lags[i]-1)]-lambda[1:(i-lags[i]-1)])))*(tau[i]-lambda[i]) 
     }
     if(p[i]>tau[i] | p[i]<=lambda[i]){
       C_S[i]=1
     }
   }
   return(alpha_ind)
 }

#Exhaustive ADDIS-Spending
E_ADDIS_Spending=function(alpha, gamma, tau, lambda, p, n){
  if(length(tau)== 1){tau=rep(tau,n)}
  if(length(lambda)== 1){lambda=rep(lambda,n)}
  if(length(gamma)!=n | length(tau)!= n | length(lambda)!=n| length(p)!=n){
    warning("mismatching length")
  }
  t=1
  alpha_k=alpha
  alpha_ind=rep(0,n)
  for(i in 1:n){
    alpha_ind[i]=alpha*gamma[t]*(tau[i]-lambda[i])/(1-alpha_k)
    if(p[i]<=tau[i] & p[i]>lambda[i]){
      alpha_k=alpha_k-(alpha_ind[i]*(1-alpha_k)/(tau[i]-lambda[i]))
      t=t+1
    }
    }
    return(alpha_ind)
  }
  
#Exhaustive ADDIS-Graph
E_ADDIS_Graph=function(alpha, gamma, w, tau, lambda, p, n){
  if(length(tau)== 1){tau=rep(tau,n)}
  if(length(lambda)== 1){lambda=rep(lambda,n)}
  if(length(gamma)!=n | length(tau)!= n | length(lambda)!=n| length(p)!=n){
    warning("mismatching length")
  }
  C_S=rep(0,n)
  alpha_ind=rep(0,n)
  alpha_k=rep(0,n)
  alpha_k[1]=alpha
  alpha_ind[1]=alpha*gamma[1]*(tau[1]-lambda[1])/(1-alpha_k[1])
  if(p[1]>tau[1] | p[1]<=lambda[1]){
    C_S[1]=1
    alpha_k[2]=alpha_k[1]
  }else{
    alpha_k[2]=alpha_k[1]-(alpha_ind[1]*(1-alpha_k[1])/(tau[1]-lambda[1]))
  }
  for(i in 2:n){
      alpha_ind[i]=(alpha*gamma[i]+sum(C_S[1:(i-1)]*alpha_ind[1:(i-1)]*w[1:(i-1),i]*(1-alpha_k[1:(i-1)])/(tau[1:(i-1)]-lambda[1:(i-1)])))*(tau[i]-lambda[i])/(1-alpha_k[i]) 
      if(p[i]>tau[i] | p[i]<=lambda[i]){
        C_S[i]=1
        alpha_k[i+1]=alpha_k[i]
      }else{
        alpha_k[i+1]=alpha_k[i]-(alpha_ind[i]*(1-alpha_k[i])/(tau[i]-lambda[i]))
      }
  }
  return(alpha_ind)
}

#Alpha-Spending
alpha_spending=function(alpha, gamma, p, n){
  if(length(gamma)!=n | length(p)!=n){
    warning("mismatching length")
  }
  alpha_ind=rep(0,n)
  for(i in 1:n){
    alpha_ind[i]=alpha*gamma[i]
  }
  return(alpha_ind)
} 


#EI ADDIS-Graph
EI_ADDIS_Graph=function(alpha, gamma, w, tau, lambda, p, n){
  if(length(tau)== 1){tau=rep(tau,n)}
  if(length(lambda)== 1){lambda=rep(lambda,n)}
  if(length(gamma)!=n | length(tau)!= n | length(lambda)!=n| length(p)!=n){
    warning("mismatching length")
  }
  C_S=rep(0,n)
  alpha_ind=rep(0,n)
  alpha_k=rep(0,n)
  alpha_k[1]=alpha
  alpha_ind[1]=alpha*gamma[1]*(tau[1]-lambda[1])
  if(p[1]>tau[1] | p[1]<=lambda[1]){
    C_S[1]=1
    alpha_k[2]=alpha_k[1]
  }else{
    alpha_k[2]=alpha_k[1]-(alpha_ind[1]*(1-alpha_k[1])/(tau[1]-lambda[1]))
  }
  for(i in 2:n){
    alpha_ind[i]=(alpha*gamma[i]+sum(C_S[1:(i-1)]*alpha_ind[1:(i-1)]*w[1:(i-1),i]*1/(tau[1:(i-1)]-lambda[1:(i-1)]))+sum((1-C_S[1:(i-1)])*alpha_ind[1:(i-1)]*w[1:(i-1),i]*alpha_k[1:(i-1)]/(tau[1:(i-1)]-lambda[1:(i-1)])))*(tau[i]-lambda[i])
    if(p[i]>tau[i] | p[i]<=lambda[i]){
      C_S[i]=1
      alpha_k[i+1]=alpha_k[i]
    }else{
      alpha_k[i+1]=alpha_k[i]-(alpha_ind[i]*(1-alpha_k[i])/(tau[i]-lambda[i]))
    }
  }
  return(alpha_ind)
}
