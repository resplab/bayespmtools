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
library(bayescpm)
library(knitr)
library(reshape2)


source("include.R")



global <- list()



# Define UI for application that draws a histogram
ui <- fluidPage(
        tags$head(
          tags$style(HTML("
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
          "))
        ),
    # Application title
    titlePanel("Bayesain power, assurance, and VoI calculator (CONFIDENTIAL - please do not share the link)"),
    tabsetPanel(id="input_tabset",
      tabPanel("Introduction",HTML("
               <H3><B>Welcome to the power, assurance, and VoI calculator for Bayesian and decision-theoretic approach for sample size considerations in risk prediction models</B></H3><BR/>
               <P>This Web app helps you to determine the power, sample size, or expected value of learning for the external validation of your risk prediction model based on uncertainty around its net benefit (NB).</P>
               <P><B>To use this app, you need the following categories of information: </B><BR/>
               
                  1. Outcome prevalence in the target population. <BR/>
                  
                  2. The performance of the model in terms of discrimination and calibration.<BR/>

                  3. Metrics of interest and precision / assurance targets.<BR/>

                  4. General setup for calculations (range of target sample size, number of Monte Carlo simulations).<BR/></P
                  
                  </P>

                  <HR/>
      "), selectInput("purpose", "How do you want to use this tool?", choices=c("I know the sample size and want to estimate precision"="pow", "I want to calculate sample size"="samp"), width="300px")
      ),
      tabPanel("Evidence on model performance",
        HTML("Here we solicit your assessment of the model performance and your uncertainty around this asessment."),
        hr(),
        
        sidebarPanel("Evidence on outcome prevalence",
             selectInput("evidence_prev_dist_type", "Distribution type", choices=c("logitnorm", "beta", "probitnorm"), selected="beta"),
             selectInput("evidence_prev_input_type", "How do you want to parameterize", choices=unname(choices$evidence_prev_input_type), selected="Mean and SD"),
             
             fluidRow(
                column(4,
                  numericInput("evidence_prev_input_parm1","Parm 1",value=0.427967,min=0,max=1,width="75px")),
                column(4,
                  numericInput("evidence_prev_input_parm2","Parm 2",value=0.0295,min=0,max=1,width="75px")))
        ),
        sidebarPanel("Evidence on c-statistic",
          selectInput("evidence_cstat_dist_type", "Distribution type", choices=c("logitnorm", "beta", "probitnorm"), selected="beta"),
          selectInput("evidence_cstat_input_type", "How do you want to parameterize", choices=unname(choices$evidence_cstat_input_type), selected=unname(choices$evidence_cstat_input_type)[2] ),
          fluidRow(
            column(4,
              numericInput("evidence_cstat_input_parm1","Parm 1",value=0.7607,min=0,max=1,width="75px")),
            column(4,
              numericInput("evidence_cstat_input_parm2","Parm 2",value=0.0062,min=0,max=1,width="75px")))
        ),
        
        sidebarPanel("Evidence on calibration slope",
          selectInput("evidence_cal_slp_dist_type", "Distribution type", choices=c("norm")),
          selectInput("evidence_cal_slp_input_type", "How do you want to parameterize", choices=unname(choices$evidence_cal_slp_input_type), selected=unname(choices$evidence_cal_slp_input_type)[2] 
                         ),
          fluidRow(
            column(4,
                   numericInput("evidence_cal_slp_input_parm1","Parm 1",value=0.9950,min=0.5,max=2,width="75px")),
            column(4,
                   numericInput("evidence_cal_slp_input_parm2","Parm 2",value=0.0237,min=0,max=1,width="75px"))),
        ), 
        sidebarPanel("Evidence on O/E or equivalent",
          selectInput("evidence_cal_other_type", "And one of:", choices=unname(choices$evidence_cal_other_type), selected="Mean calibration"),
          conditionalPanel("input.evidence_cal_other_type!=''",
                           selectInput("evidence_cal_other_dist_type", "Distribution type", choices=c("norm","lognorm")),
                           selectInput("evidence_cal_other_input_type", "How do you want to parameterize", choices=unname(choices$evidence_cal_other_input_type), selected=unname(choices$evidence_cal_other_input_type)[2]),
                           fluidRow(
                             column(4,
                                    numericInput("evidence_cal_other_input_parm1","Parm 1",value=-0.00934,min=0,max=1,width="75px")),
                             column(4,
                                    numericInput("evidence_cal_other_input_parm2","Parm 2",value=0.1245,min=0,max=1,width="75px")))
          )
        ),
        # sidebarPanel("Evidence on outcome prevalence in target",
        #   checkboxInput("b_dprev","Prevalence in the target population is different"),
        #   conditionalPanel("input.b_dprev==1",
        #                selectInput("evidence_dprev_dist_type", "Distribution type", choices=c("logitnorm", "beta", "probitnorm"), selected="beta"),
        #                selectInput("evidence_dprev_input_type", "How do you want to parameterize", choices=unname(choices$evidence_prev_input_type), selected="Mean and SD"),
        #                
        #                fluidRow(
        #                  column(4,
        #                         numericInput("evidence_dprev_input_parm1","Parm 1",value=0.427967,min=0,max=1,width="75px")),
        #                  column(4,
        #                         numericInput("evidence_dprev_input_parm2","Parm 2",value=0.0295,min=0,max=1,width="75px")))
        #   )
        # ),
        sidebarPanel(
          checkboxInput("b_impute_cor","Impute correlation", value=1),
          selectInput("dist_type", "Distribution of calibrated risks", choices=c("logitnorm","beta","probitnorm")),
        ),
        actionButton("btn_test_evidence", label="Test evidence"), pre(uiOutput("test_evidence_results"))
      ),
      
      
      
      
      
#3#TARGETS
      tabPanel("Targets",
        HTML("Please choose outcome types, metrics, and target values."),
        hr(),
        fluidRow(
             column(12,
               sidebarPanel(
                 checkboxInput("b_stats","I want to investigate precision targets for statistical metrics (discrimination and calibration)", value=TRUE),
                 hr(),
                 conditionalPanel(
                    HTML("Please specify the type of outcomes"),
                    condition = "input.b_stats==1",
                    checkboxInput("b_feciw","Expected CI width", value=TRUE),
                    conditionalPanel("input.b_feciw==1",
                        selectInput("feciw_type", "Conventional or Bayesian?", choices=c("Point estimate of CI width (conventional)"="fciw", "Expected CI width (Bayesian)"="eciw")),
                    ),
                    checkboxInput("b_qciw","Quantile of CI width", value = TRUE),
                    conditionalPanel(
                     condition = "input.b_qciw==1",
                     sliderInput("qciw", "Quantile value", 0.0, 1, value=0.9),
                    ),hr(),
                     HTML("Please specify metrics of interest"),
                     conditionalPanel(condition="input.b_fciw+input.b_eciw+input.b_qciw==0", HTML("<font color='red'>At least one item needs to be selected.</font>")),
                     checkboxInput("b_target_cstat","c-statistc", value=TRUE),
                     checkboxInput("b_target_cal_slp","Calibration slope", value = TRUE), 
                     checkboxInput("b_target_cal_oe","Calibration O/E", value=TRUE), 
                     checkboxInput("b_target_cal_int","Calibration intercept"), 
                     checkboxInput("b_target_cal_mean","Calibration in the large (mean calibration)")
                 )
             ),
             conditionalPanel("input.purpose=='samp'",
                conditionalPanel("input.b_target_cstat==1",
                      sidebarPanel(
                                HTML("c-statistic:"),
                                conditionalPanel("input.b_feciw==1",
                                  sliderInput("eciw_cstat", "Expected CI width", min=0.0, max=0.5, value=0.1)
                                ),
                                conditionalPanel("input.b_qciw==1",
                                  sliderInput("qciw_cstat", "Assurance CI width", min=0.0, max=0.5, value=0.1)
                                )
                        )
                ),
                conditionalPanel("input.b_target_cal_slp==1",
                        sidebarPanel(
                                 HTML("Calibration slope:"),
                                 conditionalPanel("input.b_feciw==1",
                                                  sliderInput("eciw_cal_slp", "Expected CI width", min=0.0, max=0.5, value=0.1)
                                 ),
                                 conditionalPanel("input.b_qciw==1",
                                                  sliderInput("qciw_cal_slp", "Assurance CI width", min=0.0, max=0.5, value=0.1)
                                 )
                        )
                ),
                conditionalPanel("input.b_target_cal_oe==1",
                        sidebarPanel(
                                 HTML("Calibration O/E:"),
                                 conditionalPanel("input.b_feciw==1",
                                                  sliderInput("eciw_cal_oe", "Expected CI width", min=0.0, max=0.5, value=0.1)
                                 ),
                                 conditionalPanel("input.b_qciw==1",
                                                  sliderInput("qciw_cal_oe", "Assurance CI width", min=0.0, max=0.5, value=0.1)
                                 )
                        )
                ),
                conditionalPanel("input.b_target_cal_mean==1",
                        sidebarPanel(
                                HTML("Mean calibration:"),
                                conditionalPanel("input.b_feciw==1",
                                                 sliderInput("eciw_cal_mean", "Expected CI width", min=0.0, max=0.5, value=0.1)
                                ),
                                conditionalPanel("input.b_qciw==1",
                                                 sliderInput("qciw_cal_mean", "Assurance CI width", min=0.0, max=0.5, value=0.1)
                                )
                        )
                ),
                conditionalPanel("input.b_target_cal_int==1",
                        sidebarPanel(
                                 HTML("Calibration intercept:"),
                                 conditionalPanel("input.b_feciw==1",
                                                  sliderInput("eciw_cal_int", "Expected CI width", min=0.0, max=0.5, value=0.1)
                                 ),
                                 conditionalPanel("input.b_qciw==1",
                                                  sliderInput("qciw_cal_int", "Assurance CI width", min=0.0, max=0.5, value=0.1)
                                 )
                        )
                )
             ),
             sidebarPanel(
               checkboxInput("b_nb","I want to specify decision-theoretic targets for net benefit (NB)", value=TRUE), 
               conditionalPanel(
                 condition = "input.b_nb == 1",
                 sliderInput("threshold", "Risk threshold", min=0, max=1, value=0.2),
                 checkboxInput("b_voi_nb","VoI (EVPI and EVSI)", value=TRUE),
                 checkboxInput("b_assurance_nb","Assurance", value = 1),
                 conditionalPanel(
                   condition = "input.b_assurance_nb == 1",
                   sliderInput("assurance_nb", "Value", 0.0, 1, value=0.9)
                 )
               )
             )
          ), actionButton("btn_test_targets","Test targets"), pre(uiOutput("test_targets_results"))
        )
      ),
      tabPanel("Setup & run",
        sidebarPanel(
          textAreaInput("N", "Please enter sample sizes of interest (comma separated)","250,500,1000,2000,4000,8000"),
          numericInput("n_sim", "Monte Carlo sample size", min=100, max=10^6, value=1000),
          numericInput("seed", "Random number seed (optional)", value = NULL, width="50px"),
          selectInput("method","Computation method", choices = unname(choices$method_type), selected=unname(choices$method_type)[2])
        ),
        actionButton("btn_run", "Run"), 
        actionButton("btn_show_args", "Show args"),
        actionButton("btn_clear_console", "Clear output"),
        pre(uiOutput("console")), #getwd(), paste(dir(getwd()),collapse=","),
        p("For long computations, please use the R package or the local version of this app on your computer")
      ),
      tabPanel("Report",
               actionButton("btn_gen_report", "Generate report"),
               actionButton("btn_show_all", "Show results"),
               uiOutput("report1")
      )
    ),
    # Footer
    tags$div(class = "footer", "App version 2025.01.28. For questions and bug reports contact msafavi@mail.ubc.ca")
)





















# Define server logic required to draw a histogram
server <- function(input, output)
{
  require(dplyr)
  
  observeEvent(input$btn_test_evidence,
    {
      isolate({e1 <- collect_evidence(input); e2 <- express_evidence(process_evidence(e1));})
      output$test_evidence_results <- renderText(paste0(paste(deparse(e1), collapse=""),"\n",e2))
    })

  observeEvent(input$btn_test_targets,
               {
                 isolate(targets <- collect_targets(input))
                 output$test_targets_results <- renderText(paste(deparse(targets), collapse=""))
               })
  
  observeEvent(input$btn_run,
   {
     args <- gen_args(input)
     
     global$args <<- args
     
     progress <- shiny::Progress$new()
     # Make sure it closes when we exit this reactive, even if there's an error
     on.exit(progress$close())
     progress$set(message = "Computing", value = 0)
     args$ex_args$f_progress <- function(txt){progress$inc(message=txt)}
     
     if(is.numeric(input$seed)) {set.seed(input$seed)}
     
     console <- capture.output({res <- do.call(ifelse(input$purpose=="pow", bpm_valpow, bpm_valsamp), args=args)}, type="message")
     
     global$results <<- res
     
     output$console <- renderUI(HTML(paste("Run complete. Check here for any error messages and then proceed to generate the report on the next panel. \n", paste(console, collapse = "\n"))))
     
     .GlobalEnv$output <- global #For debugging
     
   })
  
  
#   output$console <- reactive({
#     
#     
#     
#   }) %>%  bindEvent(input$btn_run)
#   
#   
  
  observeEvent(input$btn_clear_console,
              {
                 output$console <- renderText("")
              })
  
  observeEvent(input$btn_show_args,
               {
                 output$console <- renderText(paste(deparse(isolate(gen_args(input))), collapse="\n"))
               })
  
  observeEvent(input$btn_gen_report, {
    #rmarkdown::render(input="Report.rmd", params=list(data=global$results), clean=T, runtime="shiny")
    #output$report1 <- renderUI(includeHTML("Report.html"))
    fl <- paste0(tempdir(),"/report.html")
    rmarkdown::render(input=ifelse(input$purpose=="pow", "Report_pow.Rmd", "Report_samp.Rmd"), output_file=fl, params=list(data=global), clean=T)
    output$report1 <- renderUI(tags$iframe(src=base64enc::dataURI(file=fl, mime="text/html"), style='width:80vw;height:80vh;'))
  })
  
  observeEvent(input$btn_show_all, {
    output$report1 <- renderUI(pre(paste(deparse(global),collapse = "\n")))
  })
  
}




shinyApp(ui = ui, server = server)
