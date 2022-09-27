read.data <- function(inFile) {
  # Reads .csv file, ensures numeric fields are cast as such, then checks if
  # the file conforms to _either_ our binary or continuous file spec. If yes,
  # it returns: a binary indicator; the set of direct links; the set of other
  # links; and the entire data.frame, in a list.
  #
  # If the file does not conform to either file spec, the function throws an
  # error.
  tryCatch({
    alldata <- data.table::fread(inFile$datapath, data.table=FALSE)  # automatically checks .csv delimiters
    binaryColumns <- c("id", "study", "t", "n", "r")
    continuousColumns <- c("id", "study", "t", "n", "mean", "sd")
    # check whether the  data contains all of _either_ the binary, or continuous column names
    isBinary <- all(unlist(lapply(binaryColumns, function(x) { x %in% colnames(alldata) })))
    isContinuous <- all(unlist(lapply(continuousColumns, function(x) { x %in% colnames(alldata) })))
    if (isBinary) {
      alldata <- alldata %>%
        mutate(n = as.integer(n)) %>%
        mutate(r = as.integer(r))
      directs <- alldata %>%
        filter(!is.na(r))
      otherOutcomes <- alldata %>%
        filter(is.na(r))
    } else { # continuous
      if (!isContinuous) {
        stop("missing columns or missmatching column names. Please refresh page")
      }
      alldata <- alldata %>%
        mutate(n = as.integer(n)) %>%
        mutate(mean = as.numeric(mean)) %>%
        mutate(sd = as.numeric(sd))
      directs <- alldata %>%
        filter(!is.na(mean))
      otherOutcomes <- alldata %>%
        filter(is.na(mean))
    }
    out <- list(isBinary = isBinary,
                directs = directs,
                otherOutcomes = otherOutcomes,
                alldata = alldata)
    print("Data read successful")
    return(out)
  },
  error = function(e) {
    print(paste(e))
    state$error <- paste(e)
    return(NULL)
  })
}

nma <- function(data, sm, modelFixed, modelRandom) {
  tryCatch({
    if (sm == "OR" | sm == "RR") {
      pw <- pairwise(treat=t, event=r, n=n, data = data, studlab = id)
    } else {
      pw <- pairwise(treat=t, n=n, mean=mean, sd=sd, data = data, studlab = id)
    }
    out <- netmeta(TE, seTE, treat1, treat2, studlab,
                   data = pw, sm = sm, n1 = n1, n2 = n2,
                   fixed = modelFixed, random = modelRandom)
    return(out)
    },
    error = function(e) {
      print(paste(e))
      state$error <- paste(e)
      return(NULL)
  })
}

bnma <- function(sm, bdata, ref, eff){
  if (sm == "OR" | sm == "RR") {
    model <- BUGSnet::nma.model(data=bdata,
      outcome="r",
      N="n",
      reference=ref,
      family="binomial",
      link=ifelse(sm=="OR","logit", "log"),
      effects= eff)
    
    results <- BUGSnet::nma.run(model,
                                n.burnin=2000,
                                n.iter=10000,
                                n.chains = 2)
  } else if (sm == "SMD") {
    # The model file used here is loaded directly from the `NMAJags` library
    results <- jags(data = bdata, inits = NULL,
                    parameters.to.save = c("SMD", "SMD.ref", "tau"), n.chains = 2, n.iter = 10000,
                    n.burnin = 2000, DIC=F, n.thin=1,
                    model.file = modelNMAContinuous)
  } else if (sm == "MD") {
   model <- BUGSnet::nma.model(data=bdata,
                               outcome="mean",
                               sd="sd",
                               N="n",
                               reference=ref,
                               family="normal",
                               link="identity",
                               effects= eff)
   
   results <- BUGSnet::nma.run(model,
                                n.burnin=2000,
                                n.iter=10000,
                                n.chains = 2)
  }
  return(results)
}


bData <- function(x, sm, ref) {
  # Convenience function to prepare data for BUGSnet analysis.
  if(sm == "SMD") {
    res <- make.jagsNMA.data(id, n=n, y=mean, sd=sd, t=t, data=x, reference = ref)
  } else {
    res <- BUGSnet::data.prep(arm.data = x,
                              varname.t = "t",
                              varname.s = "id")     
  }
  return(res)
}

calculate_pooled_variance <- function(nma, study_id) {
  x <- nma$data[nma$data$studlab == study_id, c("studlab","n1","n2","treat1","treat2","var", ".narms")]
  num <- (x$n1 + x$n2-1) * x$var
  den <- (x$.narms[1] - 1) * (sum(x$n1 + x$n2)) - (x$.narms[1] * (x$.narms[1] - 1))/2
  round(sum(num) / den, 3)
}

pool_variances <- function(nma, data) {
  nma$data$var <- round(nma$seTE.adj^2, 3)  # handle the case of only two arms
  to.calculate <- which(nma$narms > 2) # note that narms and n.arms have different dimensions
  pooled <- sapply(nma$studies[to.calculate], function(study_id) {
    calculate_pooled_variance(nma, study_id)
  })  # returns a named vector
  to.join <- nma$data %>%         # get id and var only for two-arm studies to merge with data
    select(studlab, var) %>%
    mutate(studlab = as.numeric(studlab)) %>%
    group_by(studlab) %>%
    filter(n() == 1)

  if (any(nma$n.arms>2)) {
    pooled <- data.frame(id_multi = as.numeric(names(pooled)), pooled_var = pooled)
    data <- data %>% dplyr::left_join(pooled, by = c("id" = "id_multi")) %>%
      left_join(to.join, by = c("id" = "studlab")) %>%
      mutate(pooled_var = coalesce(pooled_var, var)) %>%  # TODO: give a meaningful name to new column
      select(-c(var))
  } else {
    data <- data %>% dplyr::left_join(to.join, by= c("id" = "studlab")) %>%
      mutate(pooled_var = var) %>%  # TODO: give a meaningful name to new column
      select(-c(var))
  }
  return(data)
}

#-------------------------------------------------------------------------------
# Pairwise comparison table
#-------------------------------------------------------------------------------
get_pairwise_table_rows <- function(long.data, is.binary) {
  if (is.binary) {
    pw <- as_tibble(pairwise(treat=t, event=r, n=n, data = long.data, studlab = id))
  } else {
    pw <- as_tibble(pairwise(treat=t, mean=mean, sd=sd, n=n, data = long.data, studlab = id))
  }
  pw <- pw %>%
    mutate(comparison = mapply(function(x,y){comp = sort(c(x,y)); return(paste(comp[1],comp[2],sep=":"))}, treat1, treat2))
  dirspw <- pw %>%
    group_by(comparison) %>%
    summarise(ns=length(study),n=n1+n2) %>%
    group_by(comparison,ns) %>%
    summarise(n = sum(n))
  return(dirspw)
}

build_pairwise_table <- function(treatments, directs, other.outcomes, is.binary) {
    alldata <- rbind(directs, other.outcomes)
    combinations <- t(combn(treatments,2))
    colnames(combinations) <- c("treat1","treat2")
    allcomparisons <- as_tibble(combinations) %>%
      mutate(comparison = mapply(function(x,y){comp = sort(c(x,y)); return(paste(comp[1],comp[2],sep=":"))}, treat1, treat2)) %>%
      select(comparison,treat1,treat2)

    allDirectComparisons <- get_pairwise_table_rows(alldata, is.binary) %>% rename(total.numstudies=ns, total.samplesize=n)
    # Group A: all direct comparisons
    groupAcolumn1 <- get_pairwise_table_rows(directs, is.binary) %>%
      rename(numstudies=ns, samplesize=n)
    groupA <- groupAcolumn1 %>%
      left_join(allDirectComparisons, by="comparison") %>%
      mutate(compgroup="groupA")

    # Group B: observed for other outcomes
    groupB <- allDirectComparisons %>%
      anti_join(groupAcolumn1, by="comparison") %>%
      filter(comparison %in% allcomparisons$comparison) %>%
      mutate(compgroup="groupB")

    # Group C: unobserved
    groupC <- allcomparisons %>%
      select(comparison) %>%
      anti_join(groupA, by="comparison") %>%
      anti_join(groupB, by="comparison") %>%
      mutate(compgroup="groupC")

    pairwise.table <- rbind(groupA,groupB,groupC) %>%
      left_join(allcomparisons, by="comparison") %>%
      mutate(known_unknowns = 0) %>%
      mutate(unknown_unknowns = 0) %>%
      mutate(proposed = 0) %>%
      mutate(overall_bias = 0)
    column_to_rownames(pairwise.table, var = "comparison")
    return(pairwise.table)
}

FormatPairwiseTableDropdown <- function(column,treat1,treat2,comparison,groupLabel,level,proposed) {
  groupCLabel <- "Group C: unobserved"
  choices <- NULL
  if (column == "known_unknowns") {
    if (groupLabel != groupCLabel) {
      choices <- c( "", "No bias detected", paste0("Suspected bias favouring ", treat1), paste0("Suspected bias favouring ",treat2))
    }
  } else {
    choices <- c( "", "No bias detected", paste0("Suspected bias favouring ", treat1), paste0("Suspected bias favouring ", treat2))
  }

  if (length(choices) > 0) {
    chs <- tibble::rowid_to_column(as.data.frame(choices), "n")
    chs$optiontag <- mapply(function(l, choice) {
        if (l == level + 1) {
          selectedstring <- " selected"
        } else {
          selectedstring <- " "
        }
        proposedstring <- " "
        if (column == "overall_bias") {
          if (l == proposed + 1){
            proposedstring <- " style='color:grey;font-style:italic;font-weight:bold'"
          }
        }
        optag <- paste0("<option ", "value=", l-1, selectedstring, proposedstring, ">", choice, "</option>")
        return(optag)
        }, chs$n,chs$choices)

    #check if changed from proposed
    changedstring = " "
    if(column=="overall_bias") {
      if(proposed != level){
        changedstring = " style='background-color: #ffecbf;'"
      }
    }

    # assign onchange callback and class
    res <- paste0("<select onchange=selectPairwiseTableItem(id) class='custom-select'",
                  changedstring, "id='", paste(column, make.names(comparison), sep="-vs-"),"'>",
                  unite(chs,optiontag,sep=""), "</select>")
  } else {
    res <- ""
  }
  return(res)
}

labelGroup <- function(comp.group) {
  group.labels <- c("Group A: observed for this outcome",
                    "Group B: observed for other outcomes",
                    "Group C: unobserved")
  out <- NA
  if (comp.group == "groupA"){
   out <- group.labels[1]
  } else if (comp.group == "groupB") {
    out <- group.labels[2]
  } else if(comp.group == "groupC"){
    out <- group.labels[3]
  }
  return(out)
}

proposeOverallJudgement <- function(known, unknown) {
  res <- 0
  if(known %in% c(0,1,4)){
    res <- unknown
  }else{
    res <- known
  }
  return(res)
}
#-------------------------------------------------------------------------------
# RoB table
#-------------------------------------------------------------------------------
comparisonToTreatments = function(comparison){
  ctrs <- strsplit(comparison,":")[[1]]
  return(ctrs)
}

iscomparison = function(treat1, treat2, comparison) {
  ctrs <- comparisonToTreatments
  res <- (ctrs[1]==treat1 & ctrs[2]==treat2) |
         (ctrs[2]==treat1 & ctrs[1]==treat2)
  return(res)
}

iscompared = function(treat1, treat2, directs) {
  res <- any(unlist(lapply(directs,function(x){iscomparison(treat1,treat2,x)})))
 return(res)
}

reversecomparison = function(comparison){
  ctrs <- comparisonToTreatments(comparison)
  return(paste(ctrs[2],ctrs[1],sep=":"))
}

buildRobTable <- function(table1, contribution.matrix, nmaleague, nmrleague, treat1, treat2) {
  namedMixed <- function(compgroup){
    if(compgroup == "groupA"){
      res <- "mixed/only direct"
    }else{
      res <- "indirect"
    }
    res
  }

  getBiasContribution <- function(compr, treat) {
    cmrow <- contribution.matrix[contribution.matrix$comparison == reversecomparison(compr)
                                | contribution.matrix$comparison == compr,]
    biasedComparisons <- c()
    if (missing(treat)){
      biasedComparisons <- as.vector(unlist(table1 %>%
        filter( compgroup == "groupA" ) %>%
        filter( (overall_bias == 2 | overall_bias == 3)
              ) %>%
        select(comparison)
      ))
    }else{
      biasedComparisons <- as.vector(unlist(table1 %>%
        filter( compgroup == "groupA" ) %>%
        filter( (treat==treat1 & overall_bias == 2)
              | (treat==treat2 & overall_bias == 3)
              ) %>%
        select(comparison)
      ))
    }

    if(identical(biasedComparisons, character(0))){
      return(0)
    }else{
      revbiasedComparisons <- unlist(lapply(biasedComparisons, reversecomparison))
      biasedComparisons <- c(biasedComparisons,revbiasedComparisons)
      bias <- select(cmrow, -comparison) %>%
        select(one_of(biasedComparisons))
      return(sum(bias))
    }
  }

  formatEffects <- function(efstr) {
    return(str_replace(efstr, "\\(","<br>\\("))
  }

  robTable <- table1 %>%
    rename(table1_overall_bias = overall_bias) %>%
    mutate(mixed = namedMixed(compgroup)) %>%
    mutate(contrTreat1 = suppressWarnings(getBiasContribution(comparison, treat1))) %>%
    mutate(contrTreat2 = suppressWarnings(getBiasContribution(comparison, treat2))) %>%
    mutate(contrTreat3 = suppressWarnings(getBiasContribution(comparison))) %>%
    mutate(contrEvaluation = 0) %>%
    mutate(nmaEffect = nmaleague[treat2,treat1]) %>%
    mutate(nmaEffect = formatEffects(nmaEffect)) %>%
    mutate(nmrEffect = nmrleague[treat2,treat1]) %>%
    mutate(nmrEffect = formatEffects(nmrEffect)) %>%
    mutate(effectsEvaluation = 0) %>%
    mutate(proposedOverall = 0) %>%
    mutate(overallRob = 0) %>%
    select(comparison
          ,treat1
          ,treat2
          ,contrTreat1
          ,contrTreat2
          ,contrTreat3
          ,contrEvaluation
          ,table1_overall_bias
          ,nmaEffect
          ,nmrEffect
          ,effectsEvaluation
          ,proposedOverall
          ,overallRob
          ,mixed)
  return(robTable)
}

rebuildRobTable <- function(table1, contribution.matrix, nmaleague, nmrleague, robTable){
  res <- buildRobTable(table1, contribution.matrix, nmaleague, nmrleague)
  res$contrEvaluation <- robTable$contrEvaluation
  res$effectsEvaluation <- robTable$effectsEvaluation
  res$proposedOverall <- ProposeOverallRob(res)
  res$overallRob <- robTable$overallRob
  return(res)
}

ProposeOverallRob <- function(robTable){
  mixed <- robTable$mixed
  cls <- as.numeric(robTable$contrEvaluation)
  ols <- as.numeric(robTable$table1_overall_bias)
  els <- as.numeric(robTable$effectsEvaluation)

  propose <- function(cl, ol, el, mix){
    if((cl+el)==0){
      out <- 0
    }else{
      out <- 2 # default some concerns
    }
    if(((cl == 1) | (cl == 4))){
      if((el == 1) | (mix=="indirect" & ol == 1)){
        out <- 1
      }
    }else{
      if(((cl==2)|(cl==3))){
        if(((cl == el) | ((cl == ol) & (mix == "indirect")))){
          out <- 3
        }
      }
    }
    out
  }
  res <- mapply(propose, cls, ols, els, mixed)
}

# Functions for formatting RoB table columns -----------------------------------
FormatColumnContrEval <- function(comparison, treat1, treat2, level) {
  # Formats the four possible options for each dropdown
  choices <- c("",
    "No substantial contribution from bias",
    paste0("Substantial contribution from bias favouring ", treat1),
    paste0("Substantial contribution from bias favouring ", treat2),
    "Substantial contribution from bias balanced")

  chs <- tibble::rowid_to_column(as.data.frame(choices), "n")

  chs$optiontag <- mapply(function(l, choice) {
    selected.string <- " "
    if (l == level + 1) {
      selected.string <- " selected"
    }
    optag <- paste0("<option ", "value=", l-1, selected.string, ">", choice, "</option>")
    return(optag)
    }, chs$n, chs$choices)

  # Binds an onchange listener and formats the dropdown element's ID
  res <- paste0("<select onchange=changeContrEval(id) class='custom-select'",
                "id='",
                paste("contrEvaluation", make.names(comparison), sep="-vs-"),"'>",
                unite(chs,optiontag,sep=""),
                "</select>")
  return(res)
}

# TODO: refactor a single column-formatting function that accepts the col number
FormatColumnStudyEffects <- function(comparison, treat1, treat2, level) {
  choices <- c("",
    "No evidence of small-study effects",
     paste0("Small-study effects favouring ", treat1),
     paste0("Small-study effects favouring ", treat2))

  chs <- tibble::rowid_to_column(as.data.frame(choices), "n")

  chs$optiontag <- mapply(function(l,choice){
    selected.string <- " "
    if (l==level+1) {
      selected.string = " selected"
    }
    optag <- paste0("<option ", "value=", l-1, selected.string, ">", choice, "</option>")
    return(optag)
  }, chs$n,chs$choices)

  res <- paste0("<select onchange=changeStudyEffects(id) class='custom-select'",
               "id='",
               paste("effectsEvaluation", make.names(comparison), sep="-vs-"),"'>",
               unite(chs,optiontag,sep=""),
               "</select>")
  return(res)
}

CalculateRobTableOverall <- function(comparison, treat1, treat2, proposed, level){
      choices <- c( "", "Low risk", "Some concerns", "High risk")
      chs <- tibble::rowid_to_column(as.data.frame(choices), "n")

      chs$optiontag = mapply(function(l,choice){
          if (l == level + 1) {
            selectedstring <- " selected"
          } else {
            selectedstring <- " "
          }
          if ( l == proposed + 1) {
            proposedstring <- " style='color:grey;font-style:italic;font-weight:bold'"
          } else {
              proposedstring <- " "
          }
          optag <- paste0("<option ", "value=", l - 1, selectedstring,
                          proposedstring, ">", choice, "</option>")
          return(optag)
       }, chs$n, chs$choices)

       if (proposed != level) {
         changedstring <- " style='background-color: #ffecbf;'"
       } else {
         changedstring <- " "
       }

      res <- paste0("<select onchange=changeOverallRob(id) class='custom-select'",
                 changedstring,
                 "id='",
                 paste("robTableOverall", make.names(comparison), sep="-vs-"),"'>",
                 unite(chs,optiontag,sep=""),
                 "</select>")
      return(res)
    }

overallWeb <- function(comparison, treat1, treat2, table1_overall_bias, mixed){
  choices = c("",
    "No bias detected",
    paste0("Suspected bias favouring ", treat1),
    paste0("Suspected bias favouring ", treat2))
  if (mixed == "indirect") {
    res <- choices[table1_overall_bias+1]
  } else {
    res <- paste0("<span style='color:lightgrey'>", choices[table1_overall_bias+1],"</span>")
  }
  return(res)
}
