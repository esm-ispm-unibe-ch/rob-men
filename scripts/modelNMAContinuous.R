modelNMAContinuous <- function(){
  
  for(i in 1:ns){
    w[i,1] <-0
    delta[i,t[i,1]]<-0
    u[i] ~ dnorm(0,.0001)
    for (k in 1:na[i])  {
      #normal likelihood
      y[i,t[i,k]]~dnorm(phi[i,t[i,k]],prec[i,t[i,k]])
      phi[i,t[i,k]]<-(u[i]+delta[i,t[i,k]])*pooled.sd[i]
    }
    # model
    
    for (k in 2:na[i]) {
      delta[i,t[i,k]] ~ dnorm(md[i,t[i,k]],taud[i,t[i,k]])             # trial-specific SMD distributions
      md[i,t[i,k]] <-  d[t[i,k]] - d[t[i,1]]  + sw[i,k]                   # mean of SMD distributions
      taud[i,t[i,k]] <- PREC *2*(k-1)/k                                    #precision of SMD distributions
      w[i,k] <- (delta[i,t[i,k]]  - d[t[i,k]] + d[t[i,1]])          #adjustment, multi-arm RCTs
      sw[i,k] <-sum(w[i,1:(k-1)])/(k-1) }                 # cumulative adjustment for multi-arm trials
  }
  
  
  d[ref]<-0
  for (k in 1:(ref-1)){d[k] ~ dnorm(0,.0001) }
  for (k in (ref+1):nt){d[k] ~ dnorm(0,.0001) }
  
  tau~dunif(0,5)                                      #  vague prior for random effects standard deviation
  PREC<-1/pow(tau,2)
  
  # Collection of results###########
  # pairwise SMDs
  # for all comparisons
  for (c in 1:(nt-1)) {  for (k in (c+1):nt)  { SMD[c,k] <- d[c] - d[k]  }  }
  #compared to baseline
  for (c in 1:nt) {SMD.ref[c] <-d[c] - d[ref] }
  #predictions
  for (c in 1:(ref-1)) { X[c]<-d[c] - d[ref]
  predSMD.ref[c] ~dnorm( X[c],PREC)}
  for (c in (ref+1):nt) { X[c]<-d[c] - d[ref]
  predSMD.ref[c] ~dnorm( X[c],PREC)}
  for (c in 1:(nt-1)) {  for (k in (c+1):nt)  { predSMD[c,k] ~ dnorm( SMD[c,k],PREC)  }  }
  
  # #Treatment hierarchy
  order[1:nt]<- nt + 1 - rank(d[1:nt])
    # this is when the outcome is positive - omit 'nt+1-' when the outcome is negative
  for(k in 1:nt) {
    most.effective[k]<-equals(order[k],1)
    for(j in 1:nt) {
      effectiveness[k,j]<- equals(order[k],j)
    }
  }

  for(k in 1:nt) {
    for(j in 1:nt) {
      cumeffectiveness[k,j]<- sum(effectiveness[k,1:j])
    }
  }

  #SUCRAS#

  for(k in 1:nt) {
    SUCRA[k]<- sum(cumeffectiveness[k,1:(nt-1)]) /(nt-1)
  }
  #Fit of the Model#

  for(i in 1:ns) {
    for(k in 1:na[i]) {
      Darm[i,k]<-(y[i,t[i,k]]-phi[i,t[i,k]])*(y[i,t[i,k]]-phi[i,t[i,k]])*prec[i,t[i,k]]
    }
    D[i]<- sum(Darm[i,1:na[i]])
  }
  D.bar<- sum(D[])
  
}