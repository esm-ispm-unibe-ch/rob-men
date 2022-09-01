nmafunnel <- function(nma, small.values="good") {
  # This is a function that takes a netmeta object as argument, produce contour-enhanced funnel plots and results from the test for funnel plot asymmetry.
  #             The argument small.values specify whether small treatment effects indicate a beneficial ("good", default) or harmful ("bad") effect. So, if the outcome of the nma is positive, small.value must be set to "bad".
  # It returns: a dataframe containing slope, p-values and their interpretation for each comparison with at least 10 studies;
  #             the method used for the test for funnel plot asymmetry; the number of comparisons with at least 10 studies;
  #             the effect measure used (same as in the netmeta object); the matrix reporting the number of studies per comparison.
  # If there are no studies with at least 10 studies, the function stops.

  t <- nma$A.matrix
  comp.10 <- sum(t[upper.tri(t)]>9)
  if(comp.10==0) {
    print("There are no comparisons with at least 10 studies")
    return(NULL)
  }
  else {
    print(paste("There are ", comp.10, "comparisons with at least 10 studies"))

    wo <- nma$data$treat1 > nma$data$treat2
    if (any(wo)) {
      nma$data$TE[wo] <- -nma$data$TE[wo]
      ttreat1 <- nma$data$treat1
      nma$data$treat1[wo] <- nma$data$treat2[wo]
      nma$data$treat2[wo] <- ttreat1[wo]
      if (!is.null(nma$data$n1) & !is.null(nma$data$n2)) {
        tn1 <- nma$data$n1
        nma$data$n1[wo] <- nma$data$n2[wo]
        nma$data$n2[wo] <- tn1[wo]
      }
      if (!is.null(nma$data$event1) & !is.null(nma$data$event2)) {
        tevent1 <- nma$data$event1
        nma$data$event1[wo] <- nma$data$event2[wo]
        nma$data$event2[wo] <- tevent1[wo]
      }
    }

    tests <- NULL
    for (i in 1:nrow(t)) {
      for (j in 1:ncol(t)) {
        if(t[i,j]>9 & i<j) {
          if (nma$sm=="OR" | nma$sm=="RR") {
            ma <- metabin(event1,n1,event2,n2, subset = nma$treat1==rownames(t)[i] & nma$treat2==colnames(t)[j], data=nma$data)
            funnel(ma, contour = c(0.9, 0.95, 0.99), col.contour = c("darkred", "red", "lightcoral"),
                   cex=1.5, col="darkblue", bg="blue", cex.lab=1.2, xlab = paste(ma$sm, "effect reported as", rownames(t)[i], "over", colnames(t)[j]))
            legend("topright", c("0.1 > p > 0.05", "0.05 > p > 0.01", "< 0.01"), fill = c("darkred", "red", "lightcoral"), bg="white")
            mtext(paste(rownames(t)[i], "vs", colnames(t)[j]), cex=2)
            mb <- metabias(ma, method.bias = "score")
          }
          else {
            ma <- metagen(nma$TE,nma$seTE,subset = nma$treat1==rownames(t)[i] & nma$treat2==colnames(t)[j], sm=nma$sm)
            funnel(ma, contour = c(0.9, 0.95, 0.99), col.contour = c("darkred", "red", "lightcoral"),
                   cex=1.5, col="darkblue", bg="blue", cex.lab=1.2, xlab = paste(ma$sm, "effect reported as", rownames(t)[i], "-", colnames(t)[j]))
            legend("topright", c("0.1 > p > 0.05", "0.05 > p > 0.01", "< 0.01"), fill = c("darkred", "red", "lightcoral"), bg="white")
            mtext(paste(rownames(t)[i], "vs", colnames(t)[j]), cex=2)
            mb <- metabias(ma)
          }
          small.2nd <- "Small studies favour 2nd intervention"
          small.1st <- "Small studies favour 1st intervention"
          tests <- rbind(tests, data.frame(comparison=paste(rownames(t)[i], "vs", colnames(t)[j]),
                                           bias=round(mb$estimate["bias"], digits = 2),
                                           p.value=round(mb$p.value, digits = 2),
                                           interpretation=ifelse(small.values=="good",
                                                                 ifelse(mb$estimate["bias"]>0, small.2nd, small.1st),
                                                                 ifelse(mb$estimate["bias"]>0, small.1st, small.2nd))))
        }
      }
    }
    rownames(tests) <- NULL
    res <- list(tests=tests, test.method=mb$method, num.comp.10studies=paste(comp.10, "comparisons with more than 10 studies"), effect.measure=ma$sm, comp.matrix=t)
  }
}
