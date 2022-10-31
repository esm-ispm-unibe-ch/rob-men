outJagsNMAmedian <- function (JAGSobject, parameter = parameter, forestplot = F, 
          treatnames = NA, rounding = 2) 
{
  resultstable = JAGSobject$BUGSoutput$summary
  allvariablenames = rownames(resultstable)
  rowsmatching = substr(allvariablenames, 1, nchar(parameter))
  rowstokeep = startsWith(rowsmatching, parameter)
  resultstabletokeep = resultstable[rowstokeep, c(5, 3, 7, 
                                                  2)]
  rowstokeep2 = startsWith(dimnames(resultstabletokeep)[[1]], 
                           paste(parameter, "[", sep = ""))
  resultstabletokeep = resultstabletokeep[rowstokeep2, ]
  library(stringr)
  tosplit = unlist(strsplit(rownames(resultstabletokeep), 
                            ","))
  tosplit2 <- as.numeric(str_extract(tosplit, "[0-9]+"))
  nroftreatments = max(tosplit2)
  location = matrix(tosplit2, ncol = 2, byrow = T)
  medianmat = CImat = sdmat = matrix(NA, nrow = nroftreatments, 
                                   ncol = nroftreatments)
  for (i in 1:nrow(location)) {
    medianmat[location[i, 1], location[i, 2]] = resultstabletokeep[i, 
                                                                 1]
    sdmat[location[i, 1], location[i, 2]] = resultstabletokeep[i, 
                                                               4]
    CImat[location[i, 1], location[i, 2]] = resultstabletokeep[i, 
                                                               3]
    CImat[location[i, 2], location[i, 1]] = resultstabletokeep[i, 
                                                               2]
  }
  if (forestplot) {
    library(metafor)
    slab1 = rep(1:(nroftreatments - 1), (nroftreatments - 
                                           1):1)
    a = t(sapply(1:nroftreatments, rep, nroftreatments))
    slab2 = a[lower.tri(a, F)]
    slab = paste(slab1, "vs", slab2, sep = "")
    forest(x = medianmat[upper.tri(medianmat)], ci.lb = CImat[lower.tri(CImat)], 
           ci.ub = CImat[upper.tri(CImat)], slab = slab, xlab = "Network meta-analysis results")
  }
  medianmat = round(medianmat, rounding)
  CImat = round(CImat, rounding)
  Ttreat = dim(medianmat)[1]
  toprintmat = matrix(nrow = Ttreat, ncol = Ttreat)
  for (i in c(1:(Ttreat - 1))) {
    for (j in c((i + 1):Ttreat)) {
      toprintmat[i, j] = paste(medianmat[i, j], "(", CImat[j, 
                                                         i], ",", CImat[i, j], ")", sep = "")
      toprintmat[j, i] = paste(c(-medianmat[i, j]), "(", 
                               c(-CImat[i, j]), ",", c(-CImat[j, i]), ")", 
                               sep = "")
    }
  }
  if (!missing(treatnames)) {
    diag(medianmat) = treatnames
    diag(CImat) = treatnames
    diag(toprintmat) = treatnames
  }
  list(Medians = medianmat, CI = CImat, leaguetable = toprintmat)
}
