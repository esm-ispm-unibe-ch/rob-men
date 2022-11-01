library(shiny)
library(dplyr)
library(data.table)
library(DT)
library(netmeta)
library(BUGSnet)
library(bslib)  # used to override default  bootstrap version (4 instead of 3)
require(tidyverse)
source("scripts/util.R")
source("scripts/netcontrib.R")
source("scripts/nmafunnel.R")
source("scripts/makeJagsNMAdata.R")
source("scripts/modelNMAContinuous.R")
source("scripts/modelNMRContinuous.R")
source("scripts/modelNMRContinuous_bExch.R")
source("scripts/outJagsNMAmedian.R")
library(R2jags)

server <- function(input, output, session) {
  #-----------------------------------------------------------------------------
  # State variables
  #-----------------------------------------------------------------------------
  state <- reactiveValues(data = {},  # list of state variables
    parametersSet = FALSE,
    analysisStarted = FALSE,
    # analysis options
    inputSM = "",
    inputMod = "",
    inputRef = "",
    inputBH = "",
    inputBeta = "",
    treatments = "",
    modelFixed = FALSE,
    modelRandom = FALSE,
    burnIn = 1000,
    numIter = 10000,
    # analysis outputs and state
    error = "",
    nma = "",
    nmaDone = FALSE,
    bData = "",
    bnma = "",
    bnmaDone = FALSE,
    bnmr = "",
    bnmrDone = FALSE,
    nmrData = NULL,
    nmrLeague = "",
    # tables calculated from results
    contributionMatrix = "",
    pairwiseTable = tibble(),
    robTable = tibble()
    )

  #-----------------------------------------------------------------------------
  # Reactives
  #-----------------------------------------------------------------------------
  dataset <- reactive({state$data$alldata})  # reactive value to store the data
  directs <- reactive({state$data$directs})  # just the direct links

  btab <- reactive({
    validate(need(state$bData != "", "bData not ready"),
             need(state$parametersSet, "parameters not set"))
    net.tab(data = state$bData,
            outcome = ifelse(state$inputSM =="OR" | state$inputSM == "RR", "r", "mean"),
            N = "n",
            type.outcome = ifelse(state$inputSM=="OR" | state$inputSM=="RR", "binomial", "continuous"))
  })

  #-----------------------------------------------------------------------------
  # Observers
  #-----------------------------------------------------------------------------
  observeEvent(input$file, {
    state$data <- read.data(input$file)
  }, ignoreInit = TRUE)

  observeEvent(state$data$directs, {
     res <- unique(state$data$directs$t) %>% sort()
     state$treatments <- res
  })

  # Analysis hook --------------------------------------------------------------
  observeEvent(input$startAnalysis, {
    print("starting analysis")
    state$analysisStarted <- TRUE
  }, ignoreInit = TRUE, once = TRUE)

  observeEvent(state$analysisStarted, {
    if (state$analysisStarted) {
      print("starting NMA")
      state$nma <- nma(state$data$directs, state$inputSM, state$modelFixed, state$modelRandom)
      state$nmaDone <- TRUE
      print("NMA done")
    } else {
      state$nma <- ""
      state$nmaDone <- FALSE
    }
  }, ignoreNULL = TRUE)

  # Analysis parameters --------------------------------------------------------
  observeEvent(input$inputSM,{
      state$inputSM <- input$inputSM
  })

  observeEvent(state$inputMod,{
    state$modelFixed <- input$inputMod == "fixed"
    state$modelRandom <- input$inputMod == "random"
  }, ignoreNULL = TRUE)


  observeEvent(input$numIter, {
    state$numIter <- input$numIter
    # AH: this enforces a minimum burn in, but the input does not necessarily reflect this
    updateNumericInput(session, "burnIn", min = state$numIter * 0.1)
    state$burnIn <- max(state$burnIn, state$numIter * 0.1)
  })

  observeEvent(input$burnIn, {
     state$burnIn <- max(input$burnIn, state$numIter * 0.1)
  })

  observeEvent(input$inputMod,
     state$inputMod <- input$inputMod
  )

  observeEvent(input$inputRef,
     state$inputRef <- input$inputRef
  )

  observeEvent(input$inputBH,
     state$inputBH <- input$inputBH
  )

  observeEvent(input$inputBeta,
    state$inputBeta <- input$inputBeta
  )

  # Check that all parameters are set
  observe({
     validate(
       need(state$inputSM != "", "SM not selected"),
       need(state$inputMod != "", "Model not selected"),
       need(state$inputRef != "", "Reference not selected"),
       need(state$inputBH != "", "Value direction not selected"),
       need(state$inputBeta != "", "Assumption coefficients not selected"),
       need(state$burnIn >= (0.1 * state$numIter), "Burn-in must be at least 10% of iterations")
     )
    print("parameters set")
    state$parametersSet <- TRUE
  })
  
  observe({
    validate(
      need(state$inputSM != "", "no sm"),
      need(state$inputRef != "", "no ref")
    )
    state$bData <- bData(state$data$directs, state$inputSM, state$inputRef)
    print("bData calculated")
    print(state$bData)    
  })

  # NMA completion status
  observeEvent(state$nmaDone, {
    validate(need(state$bData !="", "bData not ready"),
             need(state$nmaDone, "nma not started"))
    print("starting bnma")
    state$bnma <- bnma(state$inputSM, state$bData, state$inputRef, state$inputMod)
    state$bnmaDone <- TRUE
    print("bnma done")
    # calculate data for Bayesian NMR
    print("Calculating pooled variance")
    bnmrData <- pool_variances(state$nma, directs())
    print(head(bnmrData))
    if(state$inputSM == "SMD"){
      state$nmrData <- makeJagsNMAdata(id, n=n, y=mean, sd=sd, t=t, data=bnmrData, reference = state$inputRef, othervar = pooled_var)
      state$nmrData$orig <- state$nmrData$variab 
      state$nmrData$variab <- state$nmrData$orig - min(state$nmrData$orig)
    } else {
      state$nmrData <- data.prep(arm.data = bnmrData, varname.t = "t", varname.s = "id")
    }
    print("calculated nmrData")
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  # Bayesian NMA completion status
  observeEvent(state$bnmaDone, {
    validate(need(state$bnmaDone, "bnma not ready"),
             need(state$nmrData != "", "nmrData not ready"))
    print("calculating bnma league")
    larger_better <- ifelse(input$inputBH == "good", FALSE, TRUE)
    if (state$inputSM=="SMD") {
      league <- outJagsNMAmedian(state$bnma, parameter = "SMD", treatnames = state$treatments)$leaguetable
      colnames(league) <- state$treatments
      rownames(league) <- state$treatments
      state$bleague <- league
    } else {
      state$bleague <- BUGSnet::nma.league(state$bnma,
                                           central.tdcy = "median",
                                           order = nma.rank(state$bnma, largerbetter = larger_better)$order,
                                           log.scale = FALSE)$table
    }
    print(state$bleague)
    print(state$nmrData)
    #output$league <- renderTable({ state$bleague })
    
    # start Bayesian NMR
    if (state$inputSM == "OR" | state$inputSM == "RR") {
      model <- BUGSnet::nma.model(data=state$nmrData,
                                  outcome="r",
                                  N="n",
                                  reference=input$inputRef,
                                  family="binomial",
                                  link=ifelse(state$inputSM=="OR","logit", "log"),
                                  effects= input$inputMod,
                                  covariate = "pooled_var",
                                  prior.beta = input$inputBeta)
      
      state$bnmr <- BUGSnet::nma.run(model, 
                                     n.burnin = state$burnIn,
                                     n.iter=state$numIter, 
                                     n.chains = 2, 
                                     DIC = F)
      state$bnmrDone <- TRUE      
    } else if (state$inputSM == "SMD") {
      # The model file used here is loaded directly from the `NMAJags` library or sourced from the directory (for modelNMRContinuous_bExch)
      print("inputSM else")
      if (input$inputBeta=="UNRELATED") {
        state$bnmr <- jags(data = state$nmrData, inits = NULL,
                           parameters.to.save = c("SMD", "SMD.ref", "b", "tau"), n.chains = 2, 
                           n.iter = state$numIter, n.burnin = state$burnIn, 
                           DIC=F, n.thin=1, model.file = modelNMRContinuous)
      } else {
        state$bnmr <- jags(data = state$nmrData, inits = NULL,
                           parameters.to.save = c("SMD", "SMD.ref", "B", "tau"), n.chains = 2, 
                           n.iter = state$numIter, n.burnin = state$burnIn, 
                           DIC=F, n.thin=1, model.file = modelNMRContinuous_bExch)
      }
      state$bnmrDone <- TRUE
    } else if (state$inputSM== "MD") {
      model <- BUGSnet::nma.model(data=state$nmrData,
                                  outcome="mean",
                                  sd="sd",
                                  N="n",
                                  reference=input$inputRef,
                                  family="normal",
                                  link="identity",
                                  effects= input$inputMod,
                                  covariate = "pooled_var",
                                  prior.beta = input$inputBeta)
      
      state$bnmr <- BUGSnet::nma.run(model, 
                                   DIC=F,
                                   n.burnin = state$burnIn,
                                   n.iter=state$numIter, 
                                   n.chains = 2)

    state$bnmrDone <- TRUE
    }

    print("bnmr is Done")

    # AH: this output must be inside an observer, as it uses a reactive value (inputBeta)
    output$coefficients <- renderTable({
      #coef <- NULL
      if (state$inputBeta == "UNRELATED") {
        if (state$inputSM== "SMD") {
          print(state$bnmr$BUGSoutput$summary[grep("b", rownames(state$bnmr$BUGSoutput$summary)),"mean"])
          coef <- state$bnmr$BUGSoutput$summary[grep("b", rownames(state$bnmr$BUGSoutput$summary)),"mean"]
        } else {
          c1 <- state$bnmr$samples[1][[1]][,grep("beta",colnames(state$bnmr$samples[1][[1]]))]
          c2 <- state$bnmr$samples[2][[1]][,grep("beta",colnames(state$bnmr$samples[2][[1]]))]
          coef <- colMeans(rbind(c1, c2))
        }
        coef
      }
    }, rownames = T, colnames = F)
    
    output$coefficientMean <- renderText({
      if (state$inputBeta == "EXCHANGEABLE") {
        if (state$inputSM== "SMD") {
          print(state$bnmr$BUGSoutput$summary[grep("B", rownames(state$bnmr$BUGSoutput$summary)),"mean"])
          coef <- round(state$bnmr$BUGSoutput$mean$B, digits=3)
        } else {
          c1 <- state$bnmr$samples[1][[1]][,grep("beta",colnames(state$bnmr$samples[1][[1]]))]
          c2 <- state$bnmr$samples[2][[1]][,grep("beta",colnames(state$bnmr$samples[2][[1]]))]
          coef <- round(mean(rbind(c1, c2), digits=3))
        }
        coef
      }
    })

    # build pairwise comparison table
    state$pairwiseTable <- build_pairwise_table(state$treatments, state$data$directs, state$data$otherOutcomes, state$data$isBinary)
    
    # create funnel plots
    state$hasFunnels <- !is.null(fp())
    
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  
  # Bayesian NMR completion status
  observeEvent(state$bnmrDone, {
    validate(need(state$analysisStarted, "Analysis not started"),
             need(state$bnmrDone, "bnmr not Done"))
    # calculate contribution matrix
    print("calculating contribution")
    cm <- netcontrib(state$nma, state$inputMod)
    #print(attributes(cm))
    state$contributionMatrix <- round(cm, digits = 1) %>%
      mutate(comparison=rownames(cm)) %>%
      relocate(comparison)
    print("calculated contribution matrix")
    output$contributionMatrix <- shiny::renderTable(state$contributionMatrix, digits=1)

    ## calculate NMR league table
    print("Calculating NMR league table")
    #print(head(state$nmrData))
    larger_better <- ifelse(input$inputBH == "good", FALSE, TRUE)
    if (state$inputSM == "SMD") {
      nmrLeague <- outJagsNMAmedian(state$bnmr, parameter = "SMD", treatnames = state$treatments)$leaguetable
      colnames(nmrLeague) <- state$treatments
      rownames(nmrLeague) <- state$treatments
    } else  {
      cov_value <- min(state$nmrData$arm.data$pooled_var)
      rank <- nma.rank(state$bnmr, largerbetter = larger_better, cov.value = cov_value)
      nmrLeague <- BUGSnet::nma.league(state$bnmr,
                                     central.tdcy = "median",
                                     order = rank$order,
                                     log.scale = FALSE,
                                     cov.value = cov_value)$table
      
    }
    state$nmrLeague <- nmrLeague
    output$nmr <- renderTable({
      state$nmrLeague
    }, rownames = TRUE)
    
  }, ignoreInit = TRUE, ignoreNULL = T)

  # Observers for pairwise table -----------------------------------------------
  observeEvent(input$setNoBiasWithin, {
    print("Setting within-study assessment as no bias")
    state$pairwiseTable <- mutate(state$pairwiseTable, known_unknowns = 1)
    # calculate the overall judgement according to known_unknowns and unknown_unknowns
    state$pairwiseTable <- state$pairwiseTable %>%
      mutate(overall_bias = proposed) %>%
      mutate(proposed = mapply(proposeOverallJudgement, known_unknowns, unknown_unknowns))
  })

  observeEvent(input$setNoBiasAcross, {
    print("Setting across-study assessment as no bias")
    state$pairwiseTable <- mutate(state$pairwiseTable, unknown_unknowns = 1)
    # calculate the overall judgement according to known_unknowns and unknown_unknowns
    state$pairwiseTable <- state$pairwiseTable %>%
      mutate(overall_bias = proposed) %>%
      mutate(proposed = mapply(proposeOverallJudgement, known_unknowns, unknown_unknowns))
  })

  observeEvent(input$calculateOverallJudgement, {
    print("Applying proposed to overall")
    state$pairwiseTable <- state$pairwiseTable %>%
      mutate(overall_bias = proposed) %>%
      mutate(proposed = mapply(proposeOverallJudgement, known_unknowns, unknown_unknowns))
    # AH: we only calculate the bias contributions when the calculate button is
    # clicked because it is an expensive operation
    print("Rebuilding RoB table")
    state$robTable <- rebuildRobTable(state$pairwiseTable,
                                      state$contributionMatrix,
                                      state$bleague,
                                      state$nmrLeague,
                                      state$robTable)
  })

  observeEvent(input$pairwiseSelect, {
    selected <- input$pairwiseSelect
    sel <-  unlist(strsplit(selected$id,"-vs-",fixed = TRUE))
    icolumn <- sel[[1]]
    icomparison <- sel[[2]]
    chr <- state$pairwiseTable %>%
      filter(make.names(comparison) == icomparison) %>%
      mutate("{icolumn}" := as.integer(selected$value)) %>%
      mutate(proposed = mapply(proposeOverallJudgement, known_unknowns, unknown_unknowns))
    state$pairwiseTable <- rows_update(state$pairwiseTable, chr)
  })

  # observer to construct RoB table when all other tables exist
  observe({
    validate(need(state$analysisStarted, "Analysis in progress..."),
      need(nrow(state$pairwiseTable) > "0", "pairwise comparison table not ready"),
      need(state$contributionMatrix != "", "contribution matrix not ready"),
      need(state$nmrLeague != "", "NMR league table not ready"))

    isolate({
      if (nrow(state$robTable) == 0) {
        print(c("building not rebuilding",state$state2))
        state$robTable <- buildRobTable(state$pairwiseTable, state$contributionMatrix, state$bleague, state$nmrLeague)
      }
    })
  })

  # Observers for RoB table ----------------------------------------------------
  observeEvent(input$setSSEUndetected, {
    print("Setting SSE to No bias detected")
    state$robTable <- mutate(state$robTable, effectsEvaluation = 1)
    state$robTable$proposedOverall <- ProposeOverallRob(state$robTable)
  })

  observeEvent(input$setNoContribution, {
    print("Setting Evaluation of contribution to No contribution")
    state$robTable <- mutate(state$robTable, contrEvaluation = 1)
    state$robTable$proposedOverall <- ProposeOverallRob(state$robTable)
  })

  observeEvent(input$applyProposedRobTable, {
    print("Applying proposed to Final")
    state$robTable <- mutate(state$robTable, overallRob = proposedOverall)
  })

  observeEvent(input$contrEval, {
      selected <- input$contrEval
      sel <-  unlist(strsplit(selected$id, "-vs-", fixed = TRUE))
      icomparison <- sel[[2]]
      chr <- filter(state$robTable, icomparison == make.names(comparison)) %>%
             mutate(contrEvaluation = as.integer(selected$value))
      state$robTable <- rows_update(state$robTable, chr)
      state$robTable$proposedOverall <- ProposeOverallRob(state$robTable)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  observeEvent(input$studyEffects, {
      selected <- input$studyEffects
      sel <-  unlist(strsplit(selected$id,"-vs-", fixed = TRUE))
      icomparison <- sel[[2]]
      chr <- filter(state$robTable, icomparison == make.names(comparison)) %>%
             mutate(effectsEvaluation = as.integer(selected$value))
      state$robTable <- rows_update(state$robTable, chr)
      state$robTable$proposedOverall <- ProposeOverallRob(state$robTable)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  observeEvent(input$overallRob, {
      selected <- input$overallRob
      sel <-  unlist(strsplit(selected$id ,"-vs-", fixed = TRUE))
      icomparison <- sel[[2]]
      chr <- filter(state$robTable, icomparison == make.names(comparison)) %>%
        mutate(overallRob = as.integer(selected$value))
      print(chr)
      state$robTable <- rows_update(state$robTable, chr)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  #-------------------------------------------------------------------------------
  # Load tab UI components
  #-------------------------------------------------------------------------------
  output$contents <- DT::renderDataTable({ DT::datatable(dataset()) })

  output$tabLoad <- renderUI({
    fluidPage(theme = bs_theme(version = 4),
      sidebarLayout(
        sidebarPanel(
          fileInput("file", "Choose CSV file",
                     accept = c("text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv"))
        ),
        mainPanel(
          DT::dataTableOutput("contents")
        )
      )
    )
  })

  #-----------------------------------------------------------------------------
  # Analysis tab UI components
  #-----------------------------------------------------------------------------
  output$smOptions <- renderUI({
    chs <- state$inputSM
    if (state$data$isBinary) {
      if (!state$analysisStarted) {
        chs <- c("Odds Ratio" = "OR", "Risk Ratio" = "RR")
      }
    } else {
      if (!state$analysisStarted) {
        chs <- c("Standardized mean difference" = "SMD",
                 "Mean difference" = "MD")
      }
    }
    radioButtons(inputId = "inputSM",
                 label = "Summary measure",
                 choices = chs,
                 selected = state$inputSM)
  })

  output$bhOptions <- renderUI({
    if (!state$analysisStarted) {
      chs = c("Desirable" = "good", "Undesirable" = "bad")
    } else {
      if (state$inputBH == "good") {
        chs = c("Desirable" = "good")
      } else {
        chs = c("Undesirable" = "bad")
      }
    }
    radioButtons(inputId = "inputBH",
                label = "Smaller outcome values are",
                selected = state$inputBH,
                choices = chs)
  })

  output$ModelOptions <- renderUI({
    if (!state$analysisStarted){
      chs = c("Random effects" = "random", "Fixed effects" = "fixed")
    } else {
      if (state$inputMod == "fixed") {
        chs = c("Fixed effects" = "fixed")
      } else {
        chs = c("Random effects" = "random")
      }
    }
    radioButtons(inputId = "inputMod",
                label = "Synthesis model",
                selected = state$inputMod,
                choices = chs)
  })

  output$ref <- renderUI({
    if (!state$analysisStarted) {
      chs <- state$treatments
    } else {
      chs <- state$inputRef
    }
    radioButtons(inputId = "inputRef",
                 label = "Reference treatment",
                 selected = state$inputRef,
                 choices = chs)
  })

  output$bugsnetOptions <- renderUI({
    if (!state$analysisStarted) {
      bin <- state$burnIn
      iter <- state$numIter
      tags$div(
        numericInput(
          inputId = "burnIn",
          label = "Burn In",
          value = bin,
          min = (0.1 * state$numIter),
          max = NA,
          step = NA,
          width = NULL
        ),
        numericInput(
          inputId = "numIter",
          label = "Iterations",
          value = iter,
          min = 10000,
          max = NA,
          step = 10,
          width = NULL
        )
      )
    } else {
      tags$div(
        tags$div(class = "control-label", "Burn In"),
        tags$p(state$burnIn),
        tags$div(class = "control-label", "Iterations"),
        tags$p(state$numIter)
      )
    }
  })

  output$priorbeta <- renderUI({
    if (!state$analysisStarted ){
      chs <- c("Unrelated treatment-specific interactions" = "UNRELATED",
              "Exchangeable/related treatment-specific interactions" = "EXCHANGEABLE")
    } else {
      if (state$inputBeta == "UNRELATED") {
        chs <- c("Unrelated treatment-specific interactions" = "UNRELATED")
      }else{
        chs <- c("Exchangeable/related treatment-specific interactions" = "EXCHANGEABLE")
      }
    }
    radioButtons(inputId = "inputBeta",
                 label = "Assumption for treatment-specific interactions",
                 selected = state$inputBeta,
                 choices = chs
    )
  })


  output$args <- renderUI({
    tags$div(
    uiOutput("smOptions"),
    # if(state$data$isBinary==F) {
    #   strong("NOTE: currently ROB-MEN cannot calculate standardised mean differences.", style = "color:blue")
    # },
    uiOutput("bhOptions"),
    uiOutput("ModelOptions"),
    uiOutput("ref"),
    uiOutput("bugsnetOptions"),
    uiOutput("priorbeta")
    )
  })

  output$run <- reactive({
    res <- ""
    if (state$parametersSet & !state$analysisStarted) {
      res <- paste0("<button class='btn btn-primary btn-block' ",
                   "onclick='startAnalysis()'>",
                   "<i class='fa-solid fa-play'></i> Start Analysis</button>")
    } else {
      res <- ""
    }
    return(res)
  })

  output$tabAnalysis <- renderUI({
    if (!is.null(state$data$directs)) {
    fluidPage(theme = bs_theme(version = 4),
      sidebarLayout(
        sidebarPanel(
          uiOutput("args"),
          uiOutput("run"),
          width = 3
        ),
        mainPanel(
          uiOutput("mainAnalysis"),
          width = 9
        )
      )
    )
    } else {
      fluidPage(theme = bs_theme(version = 4), tags$h4("Dataset not present"))
    }
  })

  output$mainAnalysis <- renderUI({
    if (state$analysisStarted) {
      tabsetPanel(
        tabPanel("Data Summary",uiOutput("summary")),
        tabPanel("Bayesian network meta-analysis", uiOutput("bayesianNMA")),
        tabPanel("Bayesian network meta-regression", uiOutput("bayesianNMR")),
        tabPanel("Funnel plots and test for small-study effects",
                 tabPanel("Contour-enhanced funnel plots",
                          uiOutput("funnelplots")
                 )
        ),
        tabPanel("Contribution matrix", uiOutput("tabContributionMatrix"))
      )
    } else {
      tags$h4("Analysis not started")
    }
  })

  output$summary <- renderUI({
    validate(need(state$analysisStarted, "Analysis parameters not set"))
    if (state$inputSM=="SMD") {
      fluidPage(theme = bs_theme(version = 4),
                h4("Network graph", align = "center"),
                plotOutput("netgraph", width = "100%", height = "500px"),
                tags$br(),
                h6(paste("There are", state$nma$k, "studies reporting the outcome of interest. 
                        Below are the total number of participants in each of the included interventions.")),
                tableOutput("netchar")
      )
    }
    else {
      fluidPage(theme = bs_theme(version = 4),
                h4("Network graph", align = "center"),
                plotOutput("netgraph", width = "100%", height = "500px"),
                tags$br(),
                h4("Network characteristics", align = "left"),
                tableOutput("netinfo"),
                h4("Interventions characteristics", align = "left"),
                tableOutput("intinfo"),
                h4("Direct comparisons characteristics", align = "left"),
                tableOutput("compinfo"))
    }
  })
  
  output$netchar <- renderTable({
    tapply(state$data$directs$n, state$data$directs$t, sum, na.rm=T)
  }, colnames = FALSE, rownames = TRUE)

  output$netgraph <- renderPlot({
    netgraph(state$nma, col = "black", plastic=FALSE,
             points = TRUE, col.points = "darkgreen", cex.points =10*sqrt(n.trts/max(n.trts)),
             thickness="number.of.studies", lwd.max = 12, lwd.min = 1, multiarm=F)
  })

  output$netinfo <- renderTable({
    btab()$network
  }, colnames = FALSE)

  output$intinfo <- renderTable({
    out <- btab()$intervention
    if (state$data$isBinary) {
      colnames(out) <- c("Intervention","Total no. of studies", "Total no. of events",
                         "Total no. of patients","Min observed event rate", "Max observed event rate", "Average event rate")
    } else {
      colnames(out) <- c("Intervention","Total no. of studies", "Total no. of patients",
                         "Min outcome value","Max outcome value","Average outcome value")
    }
    out
  })

  output$compinfo <- renderTable({
    if (state$data$isBinary) {
      out <- btab()$comparison[,-5]
      colnames(out) <- c("Comparison","Total no. of studies", "Total no. of patients","Total no. of events")
    } else {
      out <- btab()$comparison
      colnames(out) <- c("Comparison","Total no. of studies", "Total no. of patients")
    }
    out
  })

  #-----------------------------------------------------------------------------
  # Outputs for Bayesian NMA tab
  #-----------------------------------------------------------------------------
  output$plot_forest<- renderPlot({
    if (state$inputSM!="SMD") {
      BUGSnet::nma.forest(state$bnma, comparator = state$inputRef) + 
        ylab(paste(input$inputSM, "relative to", input$inputRef )) + 
        theme(axis.text = element_text(size=15))
    }
  })
  
  output$tau <- renderText({    
    if (state$inputSM=="SMD") {
      het <- round(state$bnma$BUGSoutput$mean$tau, 3)
    }
    else{
      het <- round(mean(c(mean(state$bnma$samples[1][[1]][,"sigma"]), mean(state$bnma$samples[2][[1]][,"sigma"]))),3)
    }
    het
  })

  output$tab_league <- renderTable({
    state$bleague
  }, rownames = TRUE)

  output$bayesianNMA <- renderUI({
    validate(need(state$analysisStarted, "analysis not started"),
             need(state$bnmaDone, "waiting for analysis"))
    if (state$inputSM=="SMD") {
      tags$div(
        br(),
        h6("The heterogeneity (tau) is estimated at ", textOutput("tau", inline = T), align = "center"),
        br(),
        h5("League table", align = "center"),
        div(tableOutput("tab_league"), style = "font-size:80%", align = "center"),
        p(paste(state$inputSM, "and 95% credible intervals of treatment in the row versus treatment in the column"), align = "center")
      )
    }
    else {
      tags$div(
        h4("Posterior medians and 95% Cr.I.", align = "center"),
        div(plotOutput("plot_forest", height = "500px", width = "800px"), align = "center"),
        br(),
        p("The heterogeneity (tau) is estimated at ", textOutput("tau", inline = T), align = "center"),
        br(),
        h5("League table", align = "center"),
        div(tableOutput("tab_league"), style = "font-size:80%", align = "center"),
        p(paste(state$inputSM, "and 95% credible intervals of treatment in the column versus treatment in the row"), align = "center")
      )
    }
  })

  #-----------------------------------------------------------------------------
  # Outputs for Bayesian NMR tab
  #-----------------------------------------------------------------------------

  output$bayesianNMR <- renderUI({
    validate(need(state$bnmrDone, "waiting for analysis"))
    isolate({
      if(state$inputSM=="SMD") {
        tags$div (
          h5("Checks for convergence of network meta-regression model", align = "center"),
          p("Check the trace plots (download) and the Gelman-Rubin diagnostic values (table below) being close to 1 for convergence. If needed, increase number of iterations and burn-in or change the assumption for treatment-specific interactions to 'Exchangeable' and rerun analysis"),
          conditionalPanel(condition = "$('html').hasClass('shiny-busy')",
                           tags$div(class = "loading", tags$img(src = "./loading.gif"))),
          div(tableOutput("rhat"), align= "center"),
          downloadButton("downloadTrace", "Download Trace Plots as PDF", class = "btn-primary"),
          h4("Network meta-regression for variance of the (linear) treatment effect", align = "center"),
          p(ifelse(state$inputBeta=="UNRELATED", "Values of the coefficients (betas) in the regression model between relative treatment effects and study variance", 
                   "The average value of the regression coefficients (betas) of the interaction between relative treatment effects and study variance is "), 
            textOutput("coefficientMean", inline = T), align="center"),
          conditionalPanel(condition = "$('html').hasClass('shiny-busy')",
                           tags$div(class = "loading", tags$img(src = "./loading.gif"))),
          div(tableOutput("coefficients"), align= "center"),
          h5("League table", align = "center"),
          p("League table showing results for the minimum observed variance of", textOutput("minvar", inline = T), align= "center"),
          conditionalPanel(condition = "$('html').hasClass('shiny-busy')",
                           tags$div(class = "loading", tags$img(src = "./loading.gif"))),
          div(tableOutput("nmr"), style = "font-size:80%", align = "center"),
          p(paste(state$inputSM, "and 95% credible intervals of treatment in the row versus treatment in the column"), align = "center")
        )
      }
      else {
        tags$div (
          h5("Checks for convergence of network meta-regression model", align = "center"),
          p("Check the trace plots (download) and the Gelman-Rubin diagnostic values (table below) being close to 1 for convergence. If needed, increase number of iterations and burn-in or change the assumption for treatment-specific interactions to 'Exchangeable' and rerun analysis"),
          conditionalPanel(condition = "$('html').hasClass('shiny-busy')",
                           tags$div(class = "loading", tags$img(src = "./loading.gif"))),
          div(tableOutput("rhat"), align= "center"),
          downloadButton("downloadTrace", "Download Trace Plots as PDF", class = "btn-primary"),
          h4("Network meta-regression for variance of the (linear) treatment effect", align = "center"),
          p(ifelse(state$inputBeta=="UNRELATED", "Values of the coefficients (betas) in the regression model between relative treatment effects and study variance", 
                   "The average value of the regression coefficients (betas) of the interaction between relative treatment effects and study variance is"), 
            textOutput("coefficientMean", inline = T), align="center"),
          conditionalPanel(condition = "$('html').hasClass('shiny-busy')",
                           tags$div(class = "loading", tags$img(src = "./loading.gif"))),
          div(tableOutput("coefficients"), align= "center"),
          h6("Press the button below to download the network meta-regression plot as PDF."), 
          p("Each line shows how the linear effect of each treatment versus reference changes for different study variances. 
            The value at variance 0 are the extrapolated linear effects of each treatment versus reference for an imaginary study with 0 variance."),
          downloadButton("downloadNmr", "Download Regression Plot as PDF", class = "btn-primary"),
          h5("League table", align = "center"),
          p("League table showing results for the minimum observed variance of", textOutput("minvar", inline = T), align= "center"),
          conditionalPanel(condition = "$('html').hasClass('shiny-busy')",
                           tags$div(class = "loading", tags$img(src = "./loading.gif"))),
          div(tableOutput("nmr"), style = "font-size:80%", align = "center"),
          p(paste(state$inputSM, "and 95% credible intervals of treatment in the column versus treatment in the row"), align = "center")
        )
      }
    })
  })
  
  output$rhat <- renderTable({
    if (state$inputSM=="SMD") {
      nmaDiag <- state$bnmr$BUGSoutput$summary[grep("ref|tau", rownames(state$bnmr$BUGSoutput$summary)),"Rhat"]
      print(nmaDiag)
    }
    else {
      nmaDiag <- nma.diag(state$bnmr, plot_prompt = FALSE)
      nmaDiag$gelman.rubin$psrf[, -2]
    }
  }, rownames = T, colnames = F)
  
  output$downloadTrace <- downloadHandler(
    filename = "traceplots.pdf",
    content = function(file)
      if (state$inputSM=="SMD") {
        pdf(file)
        R2jags::traceplot(state$bnmr, varname=c("SMD.ref","tau"), ask=FALSE)
        dev.off()
      } else {
        pdf(file)
        nma.diag(state$bnmr, plot_prompt = FALSE)  # This has to be repeated on plot because the plots
        dev.off()                                  # are created as a side-effect.
      },
    contentType = 'pdf')
  
  output$downloadNmr <- downloadHandler(
    filename <- "nmrPlot.pdf",
    content = function(file) {
      plot <- nma.regplot(state$bnmr) +
        xlab("Study variance of the (linear) treatment effect") +
        ylab(paste("Treatment effect (linear scale) versus", input$inputRef))
      ggsave(file, plot = plot, device = "pdf")
    },
    contentType = 'pdf')
  
  output$minvar <- renderText({
    if (state$inputSM=="SMD") {
      minvar <- min(min(state$nmrData$orig))
    } else {
      minvar <- min(state$nmrData$arm.data$pooled_var)
    }
    minvar
  })
  
  #-----------------------------------------------------------------------------
  # Outputs for  funnel plots tab
  #-----------------------------------------------------------------------------
  output$funnelplots <- renderUI({
    validate(need(state$nmaDone, "netmeta not ready"))
    if (state$hasFunnels) {
          fluidPage(theme = bs_theme(version = 4),
            fluidRow(
              # verbatimTextOutput("fpprint"),
              dataTableOutput("fptable")),
            hr(),
            fluidRow(
              tags$h6("Only some funnel plots are shown here. To view all plots, press the button below to download them as PDF."),
              plotOutput("plot_funnel")
            ),
            downloadButton("downloadFunnel", "Download funnel plots as PDF", class = "btn-primary")
          )
    } else {
      tags$h6("All comparisons have fewer than 10 studies")
    }
  })

  output$plot_funnel <- renderPlot({
    validate(need(state$nmaDone, "netmeta not ready"))
    par(mfrow=c(2,3))
    fp()
  })

  fp <- function() {
    nmafunnel(state$nma, small.values = state$inputBH)
  }

  output$fpprint <- renderPrint({
      fp()
  })

  output$fptable <- DT::renderDataTable(fp()$tests)


  output$downloadFunnel <- downloadHandler(
    filename = "funnelPlots.pdf",
    content = function(file) {
      pdf(file)
      fp()
      dev.off()
    },
    contentType = "pdf")

  
  #-----------------------------------------------------------------------------
  # Outputs for contribution matrix tab
  #-----------------------------------------------------------------------------
  output$tabContributionMatrix <- renderUI({
    validate(need(state$contributionMatrix != "", "contribution matrix not calculated"))
    fluidPage(theme = bs_theme(version = 4),
              tags$br(),
              p("Each cell entry provides the percentage contribution that the direct comparison (column) makes to the calculation of the corresponding NMA relative treatment effect (row)."),
              downloadButton("downloadContributionMatrix", "Download Contribution Matrix", class = "btn-primary"),
              tags$br(),
              tableOutput("contributionMatrix")
    )
  })
  
  output$downloadContributionMatrix <- downloadHandler(
    filename <- "contribution_matrix.csv",
    content <- function(file) {
      write.csv(state$contributonMatrix, file, row.names = FALSE)
    }
  )
  
  
  #-------------------------------------------------------------------------------
  # Output for pairwise comparison table
  #-------------------------------------------------------------------------------
  output$tabPairwise <- renderUI({
    if (state$analysisStarted) {
      fluidPage(theme = bs_theme(version = 4),
        downloadButton("pairwiseTableDownload", "Download Pairwise Comparison Table", class = "btn-primary"),
        tags$div(class = "with-overflow",
        uiOutput("pairwiseTableHeader"),
        DT::dataTableOutput("pairwiseTable"))
        )
    } else {
      if (is.null(state$data$directs)) {
        fluidPage(theme = bs_theme(version = 4), tags$h4("Dataset not present"))
      } else {
        fluidPage(theme = bs_theme(version = 4), tags$h4("Analysis not started"))
      }
    }
  })

  output$pairwiseTable <- DT::renderDataTable({
    validate(need(state$nmaDone, "netmeta not ready"))

    # format the HTML header for the table
    pairwiseTableHeader <- htmltools::withTags(table(
      class = 'display',
      thead(
        tr(
          th(rowspan = 2, colspan = 1, ' '),
          th(rowspan = 2, colspan = 1, 'Pairwise comparison'),
          th(rowspan = 2, colspan = 1, 'group'),
          th(rowspan = 1, colspan = 2, 'Number of studies in each comparison'),
          th(rowspan = 2, colspan = 1, 'Within-study assessment of bias', br(),
             actionButton("setNoBiasWithin", 'set all to "No bias"',
               onclick = "Shiny.setInputValue(\'setNoBiasWithin\', this.id, {priority: \'event\'})",
               class = "btn-secondary")),
          th(rowspan = 2, colspan = 1, 'Across-study assessment of bias', br(),
             actionButton("setNoBiasAcross", 'set all to "No bias"',
              onclick = "Shiny.setInputValue(\'setNoBiasAcross\', this.id, {priority: \'event\'})",
              class = "btn-secondary")),
          th(rowspan = 2, colspan = 1, 'Overall judgement', br(),
             actionButton("calculateOverallJudgement", "calculate overall judgement",
               onclick = "Shiny.setInputValue(\'calculateOverallJudgement\', this.id, {priority: \'event\'})",
               class = "btn-secondary"))
        ),
        tr(
          th(rowspan = 1, colspan = 1, 'Reporting this outcome (sample size)'),
          th(rowspan = 1, colspan = 1, 'Total identified in the SR (total sample size)')
        )
      )
    ))
    # format the columns that make up the table
    pairwiseTableWeb <- tibble()
    pairwiseTableWeb <- state$pairwiseTable %>%
      mutate(groupLabel = mapply(labelGroup, compgroup)) %>%
      mutate(known_unknownsWeb = mapply(FormatPairwiseTableDropdown, "known_unknowns", treat1, treat2, comparison, groupLabel, known_unknowns, proposed)) %>%
      mutate(unknown_unknownsWeb = mapply(FormatPairwiseTableDropdown, "unknown_unknowns", treat1, treat2, comparison, groupLabel, unknown_unknowns, proposed)) %>%
      mutate(overall_biasWeb = mapply(FormatPairwiseTableDropdown, "overall_bias", treat1, treat2, comparison, groupLabel, overall_bias, proposed)) %>%
      mutate(num_studies = paste0(numstudies," (", samplesize, ")")) %>%
      mutate(total_studies = paste0(total.numstudies," (", total.samplesize, ")")) %>%
      select(comparison, groupLabel, num_studies, total_studies, known_unknownsWeb, unknown_unknownsWeb, overall_biasWeb)

    # output the js DataTable
    datatable(pairwiseTableWeb,
              container = pairwiseTableHeader,
              escape = FALSE,
              extensions = c('RowGroup'),
              options = list(rowGroup = list(dataSrc = 2),
                             paging = FALSE,
                             ordering = FALSE,
                             columnDefs = list(list(visible=FALSE, targets=c(2))),
                             dom = 'Bfrtip'),
              selection = 'none')
  })

  output$pairwiseTableDownload <- downloadHandler(
    filename <- "pairwise_comparison_table.csv",
    content <- function(file) {
      write.csv(state$pairwiseTable, file, row.names = FALSE)
    }
  )


  #-----------------------------------------------------------------------------
  # Output for RoB-MEN table
  #-----------------------------------------------------------------------------
  output$robTableHeader<- renderUI({
   validate(need(state$nma != "", "netmeta not ready")
           , need(nrow(state$pairwiseTable)!="0","pairwise comparison table empty"))
  })

  output$tabRob <- renderUI({
    if (state$analysisStarted) {
      fluidPage(theme = bs_theme(version = 4),
        uiOutput("robTableHeader"),
        downloadButton("robTableDownload", "Download ROB-MEN Table", class = "btn-primary"),
        tags$div(class = "with-overflow",
          tabPanel("View data", DT::dataTableOutput('robTable')))
        )
    } else {
      if (is.null(state$data$directs)) {
        fluidPage(theme = bs_theme(version = 4), tags$h4("Dataset not present"))
      } else {
        fluidPage(theme = bs_theme(version = 4), tags$h4("Analysis not started"))
      }
    }
  })

  output$robTable <- DT::renderDataTable({
    validate(need(state$analysisStarted, "analysis not started"),
             need(nrow(state$pairwiseTable) != "0", "pairwise comparison table not ready"),
             need(state$contributionMatrix != "", "contribution matrix not ready"))

    robTableWeb <- state$robTable %>%
      mutate(table1_overall_bias_web = mapply(overallWeb,comparison,treat1,treat2,table1_overall_bias,mixed)) %>%
      mutate(contrTreat1Web = mapply(function(x){round(x, digits=0)},contrTreat1)) %>%
      mutate(contrTreat2Web = mapply(function(x){round(x, digits=0)},contrTreat2)) %>%
      mutate(contrTreat3Web = mapply(function(x){round(x, digits=0)},contrTreat3)) %>%
      mutate(contrEvaluationWeb = mapply(FormatColumnContrEval, comparison,treat1, treat2, contrEvaluation)) %>%
      mutate(effectsEvaluationWeb = mapply(FormatColumnStudyEffects, comparison,treat1, treat2, effectsEvaluation)) %>%
      mutate(overallRobWeb = mapply(CalculateRobTableOverall, comparison, treat1, treat2, proposedOverall, overallRob)) %>%
      select(mixed, comparison, contrTreat1Web, contrTreat2Web,
             contrEvaluationWeb, table1_overall_bias_web, nmaEffect,
             nmrEffect, effectsEvaluationWeb, overallRobWeb)

  table2Header <- htmltools::withTags(table(
    class = 'display',
    thead(
      tr(
        th(rowspan = 2, colspan=2, ''),
        th(rowspan = 2, colspan = 1, 'NMA estimate'),
        th(rowspan = 1, colspan =2, '% contribution of evidence from pairwise comparisons with suspected bias'),
        th(rowspan = 2, colspan = 1, 'Evaluation of contribution from evidence with suspected bias', br(),
          # For details about onclick binding see: https://shiny.rstudio.com/articles/communicating-with-js.html
          actionButton("setSSEUndetected", 'set all to "No substantial contribution"',
                       onclick = "Shiny.setInputValue(\'setNoContribution\', this.id, {priority: \'event\'})",
                       class = "btn-secondary")),
        th(rowspan = 2, colspan = 1, 'Bias assessment for indirect evidence'),
        th(rowspan = 2, colspan = 1, "NMA treatment effect"),
        th(rowspan = 2, colspan = 1, 'NMR treatment effect at the smallest observed variance'),
        th(rowspan = 2, colspan = 1, 'Evaluation of small-study effects', br(),
          actionButton("setSSEUndetected", 'set all to "No evidence"',
                       onclick = "Shiny.setInputValue(\'setSSEUndetected\', this.id, {priority: \'event\'})",
                       class = "btn-secondary")),
        th(rowspan = 2, colspan = 1, 'Overall risk of bias', br(),
          actionButton("applyProposedRobTable", "calculate overall RoB",
                       onclick = "Shiny.setInputValue(\'applyProposedRobTable\', this.id, {priority: \'event\'})",
                       class = "btn-secondary")),
        tr(
          th(colspan = 1, 'Favouring first treatment'),
          th(colspan = 1, 'Favouring second treatment')
        )
      )
    )
  ))

  datatable(robTableWeb,
    container = table2Header,
    escape = F,
    extensions = c('RowGroup'),
    options = list(rowGroup = list(dataSrc = 1)
                                 , paging = F
                                 , columnDefs = list(list(visible=FALSE, targets=c(1)))
                                 , dom = 'Bfrtip', ordering = FALSE),
    selection = 'none')
  }, server=F
  )

  output$robTableDownload <- downloadHandler(
    filename <- "rob_table.csv",
    content <- function(file) {
      write.csv(state$robTable, file, row.names = FALSE)
    }
  )
}
