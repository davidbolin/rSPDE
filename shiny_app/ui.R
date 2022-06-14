library(shiny)
library(plotly)
library(shinythemes)


nu_val_loglik = c("0.1","0.4","0.6","0.8","1.2", "1.8", "2.2","2.8","3.1")

ui <- navbarPage(
  "Rational SPDE Approach",
  tabPanel("Covariance error",
           sidebarLayout(position = "left",
                         sidebarPanel(width = 3,
                                      checkboxGroupInput(inputId = "whichApprox", label = "Type of rational approximation",choices = c("Covariance-based", "Operator-based"), 
                                                         selected = c("Covariance-based", "Operator-based"),
                                                         inline=TRUE),
                                      checkboxGroupInput(inputId = "orderRat", label = "Order of the approximation",choices = c("1", "2","3","4"), selected = "2",
                                                         inline=TRUE),
                                      radioButtons(inputId = "meshSize", label = "Mesh size",
                                                   choices = c("10", "20", "50", "100"),
                                                   selected = "10",
                                                   inline=TRUE),
                                      radioButtons(inputId = "rangeParameter", label = "Range Parameter",
                                                   choices = c("0.1", "0.5", "1"),
                                                   selected = "0.1",
                                                   inline=TRUE),
                                      radioButtons(inputId = "approxNorm", label = "Which norm?",
                                                   choices = c("L2", "Sup"),
                                                   selected = "L2",
                                                   inline=TRUE),
                                      checkboxInput("includeINLAcoverror", "Should INLA errors be included?", value = TRUE),
                                      checkboxInput("logScaleCoverror", "Use log scale?", value = TRUE),
                                      checkboxInput("useTrick", "Use decomposition?", value = TRUE),
                                      downloadButton('downloadPlotCov', 'Download Plot'),
                                      radioButtons(inputId = "fileExtensionCov", label = "File extension",
                                                   choices = c("png", "pdf"),
                                                   selected = "png",
                                                   inline=TRUE),
                                      downloadButton('downloadDataFrameCov', 'Download Data Frame'),),
                         mainPanel(
                           width = 12 - 3,
                           plotlyOutput("raterrors")
                         )
           ),
           
  ),
  
  tabPanel("Likelihood error", sidebarLayout(position = "left",
                                              sidebarPanel(width = 3,
                                                           radioButtons(inputId = "orderRatLoglik", label = "Order of the approximation",
                                                                        choices = c("1", "2", "3", "4"),
                                                                        selected = "2",
                                                                        inline=TRUE),
                                                           radioButtons(inputId = "rangeParameterLoglik", label = "Range Parameter",
                                                                        choices = c("0.1", "0.5", "1"),
                                                                        selected = "0.1",
                                                                        inline=TRUE),
                                                           radioButtons(inputId = "measError", label = "Measurment noise std. dev.",
                                                                        choices = c("0.01", "0.1", "1"),
                                                                        selected = "0.1",
                                                                        inline=TRUE),
                                                           checkboxGroupInput(inputId = "smoothloglik", label = "Smoothness parameters",
                                                                              choices = nu_val_loglik, selected = nu_val_loglik,
                                                                              inline=TRUE),
                                                           checkboxInput("includeINLA", "Should INLA errors be included?", value = TRUE),
                                                           radioButtons(inputId = "whichComparisonsLoglik", label = "Which approximations should be included?",
                                                                        choices = c("Operator-based", "Covariance-based", "Both"),
                                                                        selected = "Both",
                                                                        inline=TRUE),
                                                           checkboxInput("useAbsRelError", "Use absolute relative errors?", value = TRUE),
                                                           downloadButton('downloadPlotLogLik', 'Download Plot'),
                                                           radioButtons(inputId = "fileExtensionLogLik", label = "File extension",
                                                                        choices = c("png", "pdf"),
                                                                        selected = "png",
                                                                        inline=TRUE),
                                                           downloadButton('downloadDataFrameViolin', 'Download Data Frame')),
                                              
                                              mainPanel(
                                                width = 12 - 3,
                                                plotOutput("loglikerrorsfig")
                                                
                                              )
  ),),
  
  tabPanel("Likelihood error - Mean", sidebarLayout(position = "left",
                                             sidebarPanel(width = 3,
                                                          checkboxGroupInput(inputId = "orderRatLoglikMean", label = "Order of the approximation",choices = c(1, 2,3,4), selected = 2,
                                                                             inline=TRUE),
                                                          radioButtons(inputId = "rangeParameterLoglikMean", label = "Range Parameter",
                                                                       choices = c("0.1", "0.5", "1"),
                                                                       selected = "0.1",
                                                                       inline=TRUE),
                                                          radioButtons(inputId = "measErrorMean", label = "Measurment noise std. dev.",
                                                                       choices = c("0.01", "0.1", "1"),
                                                                       selected = "0.1",
                                                                       inline=TRUE),
                                                          checkboxGroupInput(inputId = "smoothloglikMean", label = "Smoothness parameters",
                                                                             choices = nu_val_loglik, selected = nu_val_loglik,
                                                                             inline=TRUE),
                                                          checkboxInput("includeINLAmean", "Should INLA errors be included?", value = TRUE),
                                                          checkboxInput("logScaleMean", "Use log scale?", value = TRUE),
                                                          radioButtons(inputId = "whichComparisonsLoglikMean", label = "Which approximations should be included?",
                                                                       choices = c("Operator-based", "Covariance-based", "Both"),
                                                                       selected = "Both",
                                                                       inline=TRUE),
                                                          checkboxInput("useAbsRelErrorMean", "Use absolute relative errors?", value = TRUE),
                                                          downloadButton('downloadPlotLogLikMean', 'Download Plot'),
                                                          radioButtons(inputId = "fileExtensionLogLikMean", label = "File extension",
                                                                       choices = c("png", "pdf"),
                                                                       selected = "png",
                                                                       inline=TRUE),
                                                          downloadButton('downloadDataFrameMean', 'Download Data Frame')),
                                             
                                             mainPanel(
                                               width = 12 - 3,
                                               plotlyOutput("loglikerrorsfigMean")
                                               
                                             )
  ),),
  
  tags$style(type = "text/css", "body {padding-top: 70px;}"),
  theme = shinytheme("cosmo"),
  position = "fixed-top"
  
)
