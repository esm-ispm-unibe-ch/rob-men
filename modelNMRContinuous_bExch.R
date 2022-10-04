modelNMRContinuous_bExch <- function () 
{
  for (i in 1:ns) {
    w[i, 1] <- 0
    delta[i, t[i, 1]] <- 0
    u[i] ~ dnorm(0, 1e-04)
    for (k in 1:na[i]) {
      y[i, t[i, k]] ~ dnorm(phi[i, t[i, k]], prec[i, t[i, 
                                                       k]])
      phi[i, t[i, k]] <- (u[i] + delta1[i, t[i, k]]) * 
        pooled.sd[i]
      delta1[i, t[i, k]] <- delta[i, t[i, k]] + beta[t[i, 
                                                       1], t[i, k]] * variab[i]
    }
    for (k in 2:na[i]) {
      delta[i, t[i, k]] ~ dnorm(md[i, t[i, k]], taud[i, 
                                                     t[i, k]])
      md[i, t[i, k]] <- d[t[i, k]] - d[t[i, 1]] + sw[i, 
                                                     k]
      taud[i, t[i, k]] <- PREC * 2 * (k - 1)/k
      w[i, k] <- (delta[i, t[i, k]] - d[t[i, k]] + d[t[i, 
                                                       1]])
      sw[i, k] <- sum(w[i, 1:(k - 1)])/(k - 1)
    }
  }
  d[ref] <- 0
  for (k in 1:(ref - 1)) {
    d[k] ~ dnorm(0, 1e-04)
  }
  for (k in (ref + 1):nt) {
    d[k] ~ dnorm(0, 1e-04)
  }
  tau ~ dunif(0, 5)
  PREC <- 1/pow(tau, 2)
  for (c in 1:(nt - 1)) {
    for (k in (c + 1):nt) {
      SMD[c, k] <- d[c] - d[k]
    }
  }
  for (c in 1:nt) {
    SMD.ref[c] <- d[c] - d[ref]
  }
  for (c in 1:(ref - 1)) {
    X[c] <- d[c] - d[ref]
    predSMD.ref[c] ~ dnorm(X[c], PREC)
  }
  for (c in (ref + 1):nt) {
    X[c] <- d[c] - d[ref]
    predSMD.ref[c] ~ dnorm(X[c], PREC)
  }
  for (c in 1:(nt - 1)) {
    for (k in (c + 1):nt) {
      predSMD[c, k] ~ dnorm(SMD[c, k], PREC)
    }
  }
  order[1:nt] <- rank(d[1:nt])
  for (k in 1:nt) {
    most.effective[k] <- equals(order[k], 1)
    for (j in 1:nt) {
      effectiveness[k, j] <- equals(order[k], j)
    }
  }
  for (k in 1:nt) {
    for (j in 1:nt) {
      cumeffectiveness[k, j] <- sum(effectiveness[k, 1:j])
    }
  }
  for (k in 1:nt) {
    SUCRA[k] <- sum(cumeffectiveness[k, 1:(nt - 1)])/(nt - 
                                                        1)
  }
  for (i in 1:nt) {
    for (j in 1:nt) {
      beta[i, j] <- b[j] - b[i]
    }
  }
  b[ref] <- 0
  for (k in 1:(ref - 1)) {
    b[k] ~ dnorm(B, 1/pow(tau_b,2))
  }
  for (k in (ref + 1):nt) {
    b[k] ~ dnorm(B, 1/pow(tau_b,2))
  }
  B ~ dnorm(0, 1e-04)
  tau_b ~ dunif(0, 5)
}




