#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)
library(mcmapper)
library(bayespmtools)
library(knitr)
library(reshape2)


source("include.R")


global <- list()


# Define UI for application that draws a histogram
ui <- fluidPage(
  tags$head(
    tags$style(HTML(
      "
            .footer {
              position: fixed;
              left: 0;
              bottom: 0;
              width: 100%;
              background-color: #f8f9fa;
              color: black;
              text-align: center;
              padding: 10px;
              border-top: 1px solid #e7e7e7;
            }
          "
    ))
  ),
  # Application title
  titlePanel(
    "Bayesain precision, assurance, and VoI calculator (CONFIDENTIAL - please do not share the link)"
  ),
  tabsetPanel(
    id = "input_tabset",
    tabPanel(
      "Introduction",
      HTML(
        "
               <H3><B>Welcome to the precsion, assurance, and VoI calculator for Bayesian and decision-theoretic approach for sample size considerations in risk prediction models</B></H3><BR/>
               <P>This Web app helps you to determine the precision, sample size, or expected value of learning for the external validation of your risk prediction model based on uncertainty around its net benefit (NB).</P>
               <P><B>To use this app, you need the following categories of information: </B><BR/>
               
                  1. Outcome prevalence in the target population. <BR/>
                  
                  2. The performance of the model in terms of discrimination and calibration.<BR/>

                  3. Metrics of interest and precision / assurance targets.<BR/>

                  4. General setup for calculations (range of target sample size, number of Monte Carlo simulations).<BR/></P
                  
                  </P>

                  <HR/>
      "
      ),
      selectInput(
        "purpose",
        "How do you want to use this tool?",
        choices = c(
          "I know the sample size and want to estimate precision" = "prec",
          "I want to calculate sample size" = "samp"
        ),
        width = "300px"
      )
    ),
    
    tabPanel(
      "Evidence on model performance",
      
      HTML(
        "Here we solicit your assessment of the model performance and your uncertainty around this assessment."
      ),
      hr(),
      
      fluidRow(
        
        ## ===== LEFT AREA (8 columns) =====
        column(
          8,
          
          ## ---- Top row inside left column: prevalence + c-stat ----
          fluidRow(
            ## Prevalence (left half of left column)
            column(
              6,
              wellPanel(
                h4("Evidence on outcome prevalence"),
                
                selectInput(
                  "evidence_prev_dist_type",
                  "Distribution type",
                  choices = c("logitnorm", "beta", "probitnorm"),
                  selected = "beta"
                ),
                selectInput(
                  "evidence_prev_input_type",
                  "How do you want to parameterize",
                  choices = unname(choices$evidence_prev_input_type),
                  selected = "Mean and SD"
                ),
                
                fluidRow(
                  column(
                    6,
                    numericInput(
                      "evidence_prev_input_parm1",
                      "Parm 1",
                      value = 0.427967,
                      min = 0,
                      max = 1
                    )
                  ),
                  column(
                    6,
                    numericInput(
                      "evidence_prev_input_parm2",
                      "Parm 2",
                      value = 0.0295,
                      min = 0,
                      max = 1
                    )
                  )
                )
              )
            ),
            
            ## C-statistic (right half of left column)
            column(
              6,
              wellPanel(
                h4("Evidence on c-statistic"),
                
                selectInput(
                  "evidence_cstat_dist_type",
                  "Distribution type",
                  choices = c("logitnorm", "beta", "probitnorm"),
                  selected = "beta"
                ),
                selectInput(
                  "evidence_cstat_input_type",
                  "How do you want to parameterize",
                  choices = unname(choices$evidence_cstat_input_type),
                  selected = unname(choices$evidence_cstat_input_type)[2]
                ),
                
                fluidRow(
                  column(
                    6,
                    numericInput(
                      "evidence_cstat_input_parm1",
                      "Parm 1",
                      value = 0.7607,
                      min = 0,
                      max = 1
                    )
                  ),
                  column(
                    6,
                    numericInput(
                      "evidence_cstat_input_parm2",
                      "Parm 2",
                      value = 0.0062,
                      min = 0,
                      max = 1
                    )
                  )
                )
              )
            )
          ), # end nested top row
          
          ## ---- Test evidence (below prevalence + c-stat) ----
          wellPanel(
            h4("Test evidence"),
            
            actionButton(
              "btn_test_evidence",
              label = "Test evidence",
              class = "btn-primary"
            ),
            br(), br(),
            pre(uiOutput("test_evidence_results"))
          )
        ), # end left column
        
        ## ===== RIGHT AREA (4 columns): calibration (tall) =====
        column(
          4,
          wellPanel(
            h4("Evidence on calibration"),
            
            selectInput(
              "evidence_cal_other_type",
              "Calibration metric",
              choices = unname(choices$evidence_cal_other_type),
              selected = "Mean calibration"
            ),
            
            conditionalPanel(
              "input.evidence_cal_other_type != ''",
              
              selectInput(
                "evidence_cal_other_dist_type",
                "Distribution type",
                choices = c("norm", "lognorm")
              ),
              selectInput(
                "evidence_cal_other_input_type",
                "How do you want to parameterize",
                choices = unname(choices$evidence_cal_other_input_type),
                selected = unname(choices$evidence_cal_other_input_type)[2]
              ),
              
              fluidRow(
                column(
                  6,
                  numericInput(
                    "evidence_cal_other_input_parm1",
                    "Parm 1",
                    value = -0.00934
                  )
                ),
                column(
                  6,
                  numericInput(
                    "evidence_cal_other_input_parm2",
                    "Parm 2",
                    value = 0.1245
                  )
                )
              )
            ),
            
            hr(),
            
            h5("Calibration slope"),
            
            selectInput(
              "evidence_cal_slp_dist_type",
              "Distribution type",
              choices = c("norm")
            ),
            selectInput(
              "evidence_cal_slp_input_type",
              "How do you want to parameterize",
              choices = unname(choices$evidence_cal_slp_input_type),
              selected = unname(choices$evidence_cal_slp_input_type)[2]
            ),
            
            fluidRow(
              column(
                6,
                numericInput(
                  "evidence_cal_slp_input_parm1",
                  "Parm 1",
                  value = 0.9950,
                  min = 0.5,
                  max = 2
                )
              ),
              column(
                6,
                numericInput(
                  "evidence_cal_slp_input_parm2",
                  "Parm 2",
                  value = 0.0237,
                  min = 0
                )
              )
            )
          )
        ) # end right column
      ) # end main row
    ),

    #3#TARGETS
    tabPanel(
      "Targets",
      HTML("Please choose outcome types, metrics, and target values."),
      hr(),
      fluidRow(
        column(
          12,
          sidebarPanel(
            checkboxInput(
              "b_eciw",
              "I want to investigate the expected value of targets for statistical metrics (discrimination and calibration)",
              value = TRUE
            ),
            conditionalPanel(
              condition = "input.b_eciw == 1",
              
              tags$h5("Please specify metrics of interest"),
              
              ## ---- c-statistic ----
              div(
                style = "
      display: grid;
      grid-template-columns: 1fr auto;
      align-items: center;
      column-gap: 16px;
      margin-bottom: 8px;
    ",
                
                checkboxInput("b_eciw_cstat", "c-statistic", value = TRUE),
                
                conditionalPanel(
                  "input.purpose == 'samp' && input.b_eciw_cstat == 1",
                  div(
                    style = "display: flex; align-items: center; gap: 6px;",
                    tags$span("Target value:"),
                    numericInput(
                      "eciw_cstat",
                      label = NULL,
                      value = 0.75,
                      min = 0.51,
                      max = 0.99,
                      width = "90px"
                    )
                  )
                )
              ),
              
              ## ---- calibration slope ----
              div(
                style = "
      display: grid;
      grid-template-columns: 1fr auto;
      align-items: center;
      column-gap: 16px;
      margin-bottom: 8px;
    ",
                
                checkboxInput("b_eciw_cal_slp", "Calibration slope", value = TRUE),
                
                conditionalPanel(
                  "input.purpose == 'samp' && input.b_eciw_cal_slp == 1",
                  div(
                    style = "display: flex; align-items: center; gap: 6px;",
                    tags$span("Target value:"),
                    numericInput(
                      "eciw_cal_slp",
                      label = NULL,
                      value = 1,
                      min = 0.5,
                      max = 1.5,
                      width = "90px"
                    )
                  )
                )
              ),
              
              ## ---- O/E ratio ----
              div(
                style = "
      display: grid;
      grid-template-columns: 1fr auto;
      align-items: center;
      column-gap: 16px;
      margin-bottom: 8px;
    ",
                
                checkboxInput("b_eciw_cal_oe", "O/E ratio", value = TRUE),
                
                conditionalPanel(
                  "input.purpose == 'samp' && input.b_eciw_cal_oe == 1",
                  div(
                    style = "display: flex; align-items: center; gap: 6px;",
                    tags$span("Target value:"),
                    numericInput(
                      "eciw_cal_oe",
                      label = NULL,
                      value = 1,
                      min = 0.5,
                      max = 1.5,
                      width = "90px"
                    )
                  )
                )
              ),
              
              ## ---- calibration intercept ----
              div(
                style = "
      display: grid;
      grid-template-columns: 1fr auto;
      align-items: center;
      column-gap: 16px;
      margin-bottom: 8px;
    ",
                
                checkboxInput("b_eciw_cal_int", "Calibration intercept", value = TRUE),
                
                conditionalPanel(
                  "input.purpose == 'samp' && input.b_eciw_cal_int == 1",
                  div(
                    style = "display: flex; align-items: center; gap: 6px;",
                    tags$span("Target value:"),
                    numericInput(
                      "eciw_cal_int",
                      label = NULL,
                      value = 1,
                      min = 0.5,
                      max = 1.5,
                      width = "90px"
                    )
                  )
                )
              ),
              
              ## ---- calibration mean ----
              div(
                style = "
      display: grid;
      grid-template-columns: 1fr auto;
      align-items: center;
      column-gap: 16px;
    ",
                
                checkboxInput("b_eciw_cal_mean", "Calibration mean", value = TRUE),
                
                conditionalPanel(
                  "input.purpose == 'samp' && input.b_eciw_cal_mean == 1",
                  div(
                    style = "display: flex; align-items: center; gap: 6px;",
                    tags$span("Target value:"),
                    numericInput(
                      "eciw_cal_mean",
                      label = NULL,
                      value = 1,
                      min = 0.5,
                      max = 1.5,
                      width = "90px"
                    )
                  )
                )
              )
            )
          ),
          
          sidebarPanel(
            checkboxInput(
              "b_qciw",
              "I want to investigate the assurance around targets for statistical metrics (discrimination and calibration)",
              value = TRUE
            ),
            conditionalPanel(
              condition = "input.b_qciw == 1",
              
              tags$h5("Please specify metrics of interest"),
              
              ## ---- c-statistic ----
              div(
                style = "
      display: grid;
      grid-template-columns: 1fr auto;
      align-items: center;
      column-gap: 16px;
      margin-bottom: 8px;
    ",
                
                checkboxInput("b_qciw_cstat", "c-statistic", value = TRUE),
                
                conditionalPanel(
                  "input.purpose == 'samp' && input.b_qciw_cstat == 1",
                  div(
                    style = "display: flex; align-items: center; gap: 6px;",
                    tags$span("Target value:"),
                    numericInput(
                      "qciw_cstat",
                      label = NULL,
                      value = 0.75,
                      min = 0.51,
                      max = 0.99,
                      width = "90px"
                    )
                  )
                )
              ),
              
              ## ---- calibration slope ----
              div(
                style = "
      display: grid;
      grid-template-columns: 1fr auto;
      align-items: center;
      column-gap: 16px;
      margin-bottom: 8px;
    ",
                
                checkboxInput("b_qciw_cal_slp", "Calibration slope", value = TRUE),
                
                conditionalPanel(
                  "input.purpose == 'samp' && input.b_qciw_cal_slp == 1",
                  div(
                    style = "display: flex; align-items: center; gap: 6px;",
                    tags$span("Target value:"),
                    numericInput(
                      "qciw_cal_slp",
                      label = NULL,
                      value = 1,
                      min = 0.5,
                      max = 1.5,
                      width = "90px"
                    )
                  )
                )
              ),
              
              ## ---- O/E ratio ----
              div(
                style = "
      display: grid;
      grid-template-columns: 1fr auto;
      align-items: center;
      column-gap: 16px;
      margin-bottom: 8px;
    ",
                
                checkboxInput("b_qciw_cal_oe", "O/E ratio", value = TRUE),
                
                conditionalPanel(
                  "input.purpose == 'samp' && input.b_qciw_cal_oe == 1",
                  div(
                    style = "display: flex; align-items: center; gap: 6px;",
                    tags$span("Target value:"),
                    numericInput(
                      "qciw_cal_oe",
                      label = NULL,
                      value = 1,
                      min = 0.5,
                      max = 1.5,
                      width = "90px"
                    )
                  )
                )
              ),
              
              ## ---- calibration intercept ----
              div(
                style = "
      display: grid;
      grid-template-columns: 1fr auto;
      align-items: center;
      column-gap: 16px;
      margin-bottom: 8px;
    ",
                
                checkboxInput("b_qciw_cal_int", "Calibration intercept", value = TRUE),
                
                conditionalPanel(
                  "input.purpose == 'samp' && input.b_qciw_cal_int == 1",
                  div(
                    style = "display: flex; align-items: center; gap: 6px;",
                    tags$span("Target value:"),
                    numericInput(
                      "qciw_cal_int",
                      label = NULL,
                      value = 1,
                      min = 0.5,
                      max = 1.5,
                      width = "90px"
                    )
                  )
                )
              ),
              
              ## ---- calibration mean ----
              div(
                style = "
      display: grid;
      grid-template-columns: 1fr auto;
      align-items: center;
      column-gap: 16px;
    ",
                
                checkboxInput("b_qciw_cal_mean", "Calibration mean", value = TRUE),
                
                conditionalPanel(
                  "input.purpose == 'samp' && input.b_qciw_cal_mean == 1",
                  div(
                    style = "display: flex; align-items: center; gap: 6px;",
                    tags$span("Target value:"),
                    numericInput(
                      "qciw_cal_mean",
                      label = NULL,
                      value = 1,
                      min = 0.5,
                      max = 1.5,
                      width = "90px"
                    )
                  )
                )
              ),
              sliderInput("qciw", "Quantile (Assurance) value", 0.90, min=0.01, max=0.99),
            )
          ),
          
          sidebarPanel(
            checkboxInput(
              "b_nb",
              "I want to specify decision-theoretic targets for net benefit (NB)",
              value = TRUE
            ),
            conditionalPanel(
              condition = "input.b_nb == 1",
              sliderInput(
                "threshold",
                "Risk threshold",
                min = 0,
                max = 1,
                value = 0.2
              ),
              checkboxInput("b_voi_nb", "VoI (EVPI and EVSI)", value = TRUE),
              checkboxInput("b_assurance_nb", "Assurance", value = 1),
              conditionalPanel(
                condition = "input.b_assurance_nb == 1",
                sliderInput("assurance_nb", "Value", 0.0, 1, value = 0.9)
              )
            )
          )
        ),
        actionButton("btn_test_targets", "Test targets"),
        pre(uiOutput("test_targets_results"))
      )
    ),
    tabPanel(
      "Setup & run",
      
      sidebarPanel(
        
        inlineInput(
          "Sample sizes",
          textAreaInput(
            "N",
            label = NULL,
            "250,500,1000,2000,4000,8000",
            rows = 2
          )
        ),
        
        inlineInput(
          "Monte Carlo sample size",
          numericInput(
            "n_sim",
            label = NULL,
            min = 100,
            max = 10^6,
            value = 1000,
            width = "120px"
          )
        ),
        
        inlineInput(
          "Random seed",
          numericInput(
            "seed",
            label = NULL,
            value = NULL,
            width = "80px"
          )
        ),
        
        inlineInput(
          "Computation method",
          selectInput(
            "method",
            label = NULL,
            choices = unname(choices$method_type),
            selected = unname(choices$method_type)[2]
          )
        ),
        
        div(
          style = "margin-bottom: 8px;",
          checkboxInput(
            "b_impute_cor",
            "Impute correlation",
            value = TRUE
          )
        ),
        
        inlineInput(
          "Risk distribution",
          selectInput(
            "dist_type",
            label = NULL,
            choices = c("logitnorm", "beta", "probitnorm")
          )
        )
      ),
      
      actionButton("btn_run", "Run"),
      actionButton("btn_show_args", "Show args"),
      actionButton("btn_clear_console", "Clear output"),
      
      pre(uiOutput("console")),
      
      p(
        "For long computations, please use the R package or the local version of this app on your computer"
      )
    )
    ,
    tabPanel(
      "Report",
      actionButton("btn_gen_report", "Generate report"),
      actionButton("btn_show_all", "Show results"),
      uiOutput("report1")
    )
  ),
  # Footer
  tags$div(
    class = "footer",
    "App version 2025.12.31. For questions and bug reports contact msafavi@mail.ubc.ca"
  )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
  require(dplyr)

  observeEvent(input$btn_test_evidence, {
    isolate({
      e1 <- collect_evidence(input)
      e2 <- summary(process_evidence(e1))
    })
    output$test_evidence_results <- renderTable(e2)
  })

  observeEvent(input$btn_test_targets, {
    isolate(targets <- collect_targets(input))
    output$test_targets_results <- renderText(paste(
      deparse(targets),
      collapse = ""
    ))
  })

  observeEvent(input$btn_run, {
    args <- gen_args(input)

    global$args <<- args

    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = "Computing", value = 0)
    args$ex_args$f_progress <- function(txt) {
      progress$inc(message = txt)
    }

    if (is.numeric(input$seed)) {
      set.seed(input$seed)
    }

    console <- capture.output(
      {
        res <- do.call(
          ifelse(input$purpose == "prec", bpm_valprec, bpm_valsamp),
          args = args
        )
      },
      type = "message"
    )

    global$results <<- res

    output$console <- renderUI(HTML(paste(
      "Run complete. Check here for any error messages and then proceed to generate the report on the next panel. \n",
      paste(console, collapse = "\n")
    )))

    .GlobalEnv$output <- global #For debugging
  })

  #   output$console <- reactive({
  #
  #
  #
  #   }) %>%  bindEvent(input$btn_run)
  #
  #

  observeEvent(input$btn_clear_console, {
    output$console <- renderText("")
  })

  observeEvent(input$btn_show_args, {
    output$console <- renderText(paste(
      deparse(isolate(gen_args(input))),
      collapse = "\n"
    ))
  })

  observeEvent(input$btn_gen_report, {
    #rmarkdown::render(input="Report.rmd", params=list(data=global$results), clean=T, runtime="shiny")
    #output$report1 <- renderUI(includeHTML("Report.html"))
    fl <- paste0(tempdir(), "/report.html")
    browser()
    rmarkdown::render(
      input = ifelse(
        input$purpose == "prec",
        "Report_prec.Rmd",
        "Report_samp.Rmd"
      ),
      output_file = fl,
      params = list(data = global),
      clean = T
    )
    output$report1 <- renderUI(tags$iframe(
      src = base64enc::dataURI(file = fl, mime = "text/html"),
      style = 'width:80vw;height:80vh;'
    ))
  })

  observeEvent(input$btn_show_all, {
    output$report1 <- renderUI(pre(paste(deparse(global), collapse = "\n")))
  })
}


shinyApp(ui = ui, server = server)
