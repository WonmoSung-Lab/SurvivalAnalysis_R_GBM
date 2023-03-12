library(shiny)
library(data.table)
library(randomForest)
library(DT) # Must include for deploying step. 
library(plotly) 
library(randomForestSRC)
library(shinyalert) # To add alert message
library(dplyr) # To rename data.frame 
library(reshape) # To make data frame into melted_df
library(reshape2)

## 1) Read-in the RFC model
feature_description <- read.csv(file="feature_description.csv")

Rs = c(10, 42, 50, 100, 250, 567, 750, 1000, 1234, 2022,
       5678, 8765, 123456, 234567, 345678, 456789, 567890, 678901, 789012, 890123)
Ks = c(1,2,3,4,5)

OS_rsfs = list()
for (R in Rs){
  for (K in Ks){
    fPath = paste0('./finalmodels/OS/RSF_numFeat11_R',as.character(R), '_K', as.character(K), '.rds')
    OS_rsfs = append(OS_rsfs, list(readRDS(fPath)))
  }
}
PFS_rsfs = list()
for (R in Rs){
  for (K in Ks){
    fPath = paste0('./finalmodels/PFS/RSF_numFeat11_R',as.character(R), '_K', as.character(K), '.rds')
    PFS_rsfs = append(PFS_rsfs, list(readRDS(fPath)))
  }
}

#### Global variables. 
real_time = seq(0, 30*12*5, length=12*5+1)
interp_times_probs = function(rsfs, acc_pt_df){
  interp_y = list()
  for (rsf in rsfs){
    pred_result <- predict.rfsrc(rsf, acc_pt_df) #, outcome = "test" 
    times_raw = pred_result$time.interest
    pred_result$survival[1,] <- pred_result$survival[1,]*100
    
    interp_1 = approx(x = times_raw, y= pred_result$survival[1,], xout = real_time)
    interp_1$y = ifelse( !is.na(interp_1$y), interp_1$y, 100.0)
    interp_y = append(interp_y, list(interp_1$y))
  }
  plotly_df = data.frame(matrix(unlist(interp_y), nrow=25, byrow=TRUE))
  return(plotly_df)
}

accumulate_df <- function(orig_pt_df, new_pt_df){
  return(rbind(orig_pt_df,new_pt_df))
}

ui <- fluidPage(
  useShinyalert(),
  
  #0. Web page name
  titlePanel(title=div(img(height=90, src="Catholic_logo.png"), img(height=110, src="Yonsei_cc_logo.png"),  '  Glioblastoma prognosis prediction calculator'), windowTitle = "RSF is all you need."), 
  # 0.Logo
  tags$p(img(height=50, src="Yesiri.png"),
         "Made by ",
         tags$a(href="https://sites.google.com/view/sunglab", tags$strong("MO Lab")),
         ".",
         align='right'),
  
  sidebarPanel(style = "height: 90vh; overflow-y: auto;background-color: #ffffff;",
               fluidRow(
                 column(12,
                        wellPanel(
                          h3("1. Biometric information", style='color:#1B4F72'),
                          
                          sliderInput(inputId = "Age",
                                      label = h4('1) Age'),
                                      value =58, min = 16, max=79),
                          radioButtons(inputId = "Gender",
                                       h4('2) Gender'),
                                       choices = list("Female"= 0,"Male" = 1),
                                       selected = 1),
                          radioButtons(inputId = "KPSG",
                                       h4('3) Preoperative KPS '),
                                       choices = list("<70"= 1,"80" = 2, ">90" = 3),
                                       selected = 2)),
                        wellPanel(
                          h3("2. Tumor-related factors", style='color:#1B4F72'),
                          radioButtons(inputId = "MGMT",
                                       h4('1) MGMT promoter methylation'), 
                                       choices = list("Unmethylation"= 0,"Methylation" = 1),
                                       selected = 0),
                          radioButtons(inputId = "ExtResc",
                                       h4('2) Extent of resection'),
                                       choices = list("Total"= 2,"Sub total" = 1, "Biopsy" = 0),
                                       selected = 2),
                          radioButtons(inputId = "IDH",
                                       h4('3) IDH1 mutations'),
                                       choices = list("negative"= -1,"mutation" = 1, "unknown"=0), #??????????
                                       selected = -1),
                          radioButtons(inputId = "Subvent",
                                       h4('4) Subventricular zone involvement'),
                                       choices = list("Uninvolved"=0, "Involved"=1),
                                       selected = 0),
                          radioButtons(inputId = "Edma",
                                       h4('5) Edma included plan'),
                                       choices = list("Yes"= 1,"No" = 2),
                                       selected = 2)),
                        
                        wellPanel(
                          h3("3. Radiotherapy details", style='color:#1B4F72'),
                          
                          textOutput("cdPTV_value_to_be_unknown"),
                          textOutput("cdPTV_value_to_be_unknown2"),
                          
                          sliderInput(inputId = "CDPTV",
                                      label = h4('1) Cone down Planning target volume (PTV) in cm3'),
                                      value = 121, min = 6, max=760), 
                          actionButton(inputId = "cdPTV_unknown_bt",label = "Unknown", icon = icon("-")), 
                          actionButton(inputId = "cdPTV_unknown_undo_bt",label = "Undo unknown"),
                          
                          sliderInput(inputId = "TotdoseCDPTV",
                                      label = h4('2) Total dose (Gy) to cone down PTV'),
                                      value = 60, min = 40, max=75),
                          actionButton(inputId = "TotDose_unknown_bt",label = "Unknown", icon = icon("-")), 
                          actionButton(inputId = "TotDose_unknown_undo_bt",label = "Undo unknown"),
                          
                          sliderInput(inputId = "Fraction",
                                      label = h4('3) Fractions (#) to cone down PTV'),
                                      value = 30, min = 10, max=40),
                          actionButton(inputId = "Fx_unknown_bt",label = "Unknown", icon = icon("-")), 
                          actionButton(inputId = "Fx_unknown_undo_bt",label = "Undo unknown"),
                        ),
                        
                        wellPanel(
                          h4("< Single patient prediction >",style = "font-weight: bold;color:#1B4F72;"),
                          actionButton(inputId = "whole_cal_button",
                                       label = "Calculate")),
                        wellPanel(
                          h4("< Multiple patients to compare >",style = "font-weight: bold;color:#1B4F72;"),
                          actionButton(inputId = "cum_button",
                                       label = "Cumulate"),
                          actionButton(inputId = "erase_button",
                                       label = "Erase"),
                          actionButton(inputId = "cum_cal_button",
                                       label = "Calculate",
                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                          textOutput("cum_pt_list"), tags$style(type="text/css", "#cum_pt_list {white-space: pre-wrap;}"))
                 )
               )
  ),
  mainPanel(tabsetPanel(
    tabPanel("Descriptions",
             tableOutput("feature_description"), 
             textOutput("about_unknown_button"),
             downloadButton("download_UNKNOWN_detail_pdf", "Details about the UNKNOWN option")
    ),
    
    tabPanel("Overall Survival",
             h5("* Description: This is a predicted survival curve of the patient."),
             h3("\n "),
             plotlyOutput(outputId = "rsf_plot", height = 600, width = 800), 
             textOutput("about_unknown_button_OS"), tags$style(type="text/css", "#about_unknown_button_OS {white-space: pre-wrap;}")
    ),
    tabPanel("Progression-free Survival",
             h5("* Description: This is a predicted progression-free survival curve of the patient."),
             h3("\n "),
             plotlyOutput(outputId = "rsf_plot_PFS", height = 600, width = 800), 
             textOutput("about_unknown_button_PFS"), tags$style(type="text/css", "#about_unknown_button_PFS {white-space: pre-wrap;}")
    ),
    
    tabPanel("Compare multiple patients",
             wellPanel(
               h3("< Patients list >",style = "font-weight: bold;color:#1B4F72;font-family: sans-serif;"),
               tableOutput("table_cum_pt_df"), style="width:auto; overflow-y: scroll;overflow-x: scroll;", 
               textOutput("unknown_selected_info"), tags$style(type="text/css", "#unknown_selected_info {white-space: pre-wrap;}")),
             wellPanel(
               h3("< Overall survival curve >",style = "font-weight: bold;color:#1B4F72;font-family: sans-serif;"),
               plotlyOutput(outputId = "cum_rsf_plot", height = 500, width = 800), style="width:auto; overflow-y: scroll;overflow-x: scroll;"),
             wellPanel(
               h3("< Progression-free survival curve >",style = "font-weight: bold;color:#1B4F72;font-family: sans-serif;"),
               plotlyOutput(outputId = "cum_rsf_plot_PFS", height = 500, width = 800), style="width:auto; overflow-y: scroll;overflow-x: scroll;")
    )
  )
  )
)


server <- function(input, output, session) {
  ################### 1. Alert messege "Academic use only!" ###################
  observe({ shinyalert(title = "Warning!", text = "Academic use only.", type = "warning") })
  
  output$feature_description <- renderTable({return(feature_description)})
  
  
  ################### ETC. Unknown options for radiotherapy details. ###################
  output$cdPTV_value_to_be_unknown = renderText({"* You can choose UNKNOWN options for radiotherapy details."})
  output$about_unknown_button <- renderText({"* The UNKNOWN buttons allow the model to preidct survival curves by taking a range of values - covering Min to Max - from the population. \n"})
  output$download_UNKNOWN_detail_pdf <- downloadHandler(
    filename = "unknown_supple.pdf",
    content = function(file) {
      file.copy("www/unknown_supple.pdf", file)
    }
  )
  
  ## Show the "unknown" button has been clicked. 
  
  unknown_val_check <- reactiveValues()
  unknown_val_check$a = 0 # For OS prediction; 0: default OS draw/ 1: cdPTV unknown OS draw / 2: TotDosePTV unknown OS draw / 3: FxPTV unknown OS draw
  unknown_val_check$b = 0 # For OS prediction; 0: default OS draw/ 1: FxPTV unknown OS draw / 2: cdPTV unknown OS draw / 3: TotDosePTV unknown OS draw
  unknown_val_check$c = 0 # Unknwon type control value.
  
  ###### 1. cdPTV ######
  is_unknown_cdPTV <- reactiveValues()
  is_unknown_cdPTV$a = FALSE
  do_update_cdPTV = function(){
    if (is_unknown_cdPTV$a == TRUE){
      updateActionButton(session, "cdPTV_unknown_bt", icon = icon("check"))
    }else{
      updateActionButton(session, "cdPTV_unknown_bt", icon = icon("-"))}}
  
  observeEvent(input$cdPTV_unknown_bt,{is_unknown_cdPTV$a = TRUE
  # For OS prediction.
  unknown_val_check$a = 1
  
  # For PFS prediction.
  if (is_unknown_Fx$a == TRUE){
    unknown_val_check$b = 1
  }else{
    unknown_val_check$b = 2
  }
  
  # For UNKNOWN type.
  unknown_val_check$c = is_unknown_cdPTV$a*4 + is_unknown_TotDose$a*2 + is_unknown_Fx$a
  
  do_update_cdPTV()  })
  
  observeEvent(input$cdPTV_unknown_undo_bt,{is_unknown_cdPTV$a = FALSE
  
  #### change "unknown_val_check" value.
  # For OS prediction.
  if (is_unknown_TotDose$a == TRUE){
    unknown_val_check$a = 2
  }else if (is_unknown_Fx$a == TRUE){
    unknown_val_check$a = 3
  }else{
    unknown_val_check$a = 0
  }
  
  # For PFS prediction.
  if (is_unknown_Fx$a == TRUE){
    unknown_val_check$b = 1
  }else if (is_unknown_TotDose$a == TRUE){
    unknown_val_check$b = 3
  }else{
    unknown_val_check$b = 0
  }
  
  # For UNKNOWN type.
  unknown_val_check$c = is_unknown_cdPTV$a*4 + is_unknown_TotDose$a*2 + is_unknown_Fx$a
  
  do_update_cdPTV()})
  
  ###### 2. Total dose ######
  is_unknown_TotDose <- reactiveValues()
  is_unknown_TotDose$a = FALSE
  do_update_TotDose = function(){
    if (is_unknown_TotDose$a == TRUE){
      updateActionButton(session, "TotDose_unknown_bt",
                         icon = icon("check"))
    }else{
      updateActionButton(session, "TotDose_unknown_bt",
                         icon = icon("-"))
    }
  }
  
  observeEvent(input$TotDose_unknown_bt,{is_unknown_TotDose$a = TRUE
  # For OS prediction.
  #### change "unknown_val_check" value. 
  if (is_unknown_cdPTV$a == TRUE){
    unknown_val_check$a = 1
  }else{
    unknown_val_check$a = 2}
  
  # For PFS prediction.
  if (is_unknown_Fx$a == TRUE){
    unknown_val_check$b = 1
  }else if (is_unknown_cdPTV$a == TRUE){
    unknown_val_check$b = 2
  }else{
    unknown_val_check$b = 3
  }
  
  # For UNKNOWN type.
  unknown_val_check$c = is_unknown_cdPTV$a*4 + is_unknown_TotDose$a*2 + is_unknown_Fx$a
  
  do_update_TotDose()
  })
  
  
  observeEvent(input$TotDose_unknown_undo_bt,{is_unknown_TotDose$a = FALSE
  #### change "unknown_val_check" value. 
  if (is_unknown_cdPTV$a == TRUE){
    unknown_val_check$a = 1
  }else if (is_unknown_Fx$a == TRUE){
    unknown_val_check$a = 3}
  else{unknown_val_check$a = 0}
  
  # For PFS prediction.
  if (is_unknown_Fx$a == TRUE){
    unknown_val_check$b = 1
  }else if (is_unknown_cdPTV$a == TRUE){
    unknown_val_check$b = 2
  }else{
    unknown_val_check$b = 0
  }
  
  # For UNKNOWN type.
  unknown_val_check$c = is_unknown_cdPTV$a*4 + is_unknown_TotDose$a*2 + is_unknown_Fx$a
  
  do_update_TotDose()})
  
  ###### 3. Fractions ######
  is_unknown_Fx <- reactiveValues()
  is_unknown_Fx$a = FALSE
  do_update_Fx = function(){
    if (is_unknown_Fx$a == TRUE){
      updateActionButton(session, "Fx_unknown_bt",
                         icon = icon("check"))
    }else{
      updateActionButton(session, "Fx_unknown_bt",
                         icon = icon("-"))
    }
  }
  
  observeEvent(input$Fx_unknown_bt,{is_unknown_Fx$a = TRUE
  
  #### change "unknown_val_check" value. 
  # For OS prediction.
  if (is_unknown_cdPTV$a == TRUE){
    unknown_val_check$a = 1
  }else if (is_unknown_TotDose$a == TRUE){
    unknown_val_check$a = 2}
  else{unknown_val_check$a = 3}
  
  # For PFS prediction.
  unknown_val_check$b = 1
  
  # For UNKNOWN type.
  unknown_val_check$c = is_unknown_cdPTV$a*4 + is_unknown_TotDose$a*2 + is_unknown_Fx$a
  
  do_update_Fx()})
  
  observeEvent(input$Fx_unknown_undo_bt,{is_unknown_Fx$a = FALSE
  
  #### change "unknown_val_check" value. 
  # For OS prediction.
  if (is_unknown_cdPTV$a == TRUE){
    unknown_val_check$a = 1
  }else if (is_unknown_TotDose$a == TRUE){
    unknown_val_check$a = 2}
  else{unknown_val_check$a = 0}
  
  # For PFS prediction.
  if (is_unknown_cdPTV$a == TRUE){
    unknown_val_check$b = 2}
  else if (is_unknown_TotDose$a == TRUE){
    unknown_val_check$b = 3}
  else{unknown_val_check$b = 0}
  
  # For UNKNOWN type.
  unknown_val_check$c = is_unknown_cdPTV$a*4 + is_unknown_TotDose$a*2 + is_unknown_Fx$a
  
  do_update_Fx()})
  
  
  
  ################### 2. Draw predicted OS survival curve of a single patient ###################
  draw_plotly_func = function(real_time, interp_y_mean, interp_y_min,interp_y_max){
    
    x_label <- list(title = "Time (month)")
    y_label <- list( title = "Survival probability (%)")
    
    plotly_df <- data.frame("time"=real_time, "survival_probability"=unlist(interp_y_mean),
                            "interp_y_min"=interp_y_min, "interp_y_max"=interp_y_max)
    
    
    results_plot = ggplot(data = plotly_df, aes(x=time, y=survival_probability))+
      geom_smooth(aes(ymin = interp_y_min, 
                      ymax = interp_y_max),
                  stat = "identity") +
      theme(panel.background = element_rect(fill = "white"),
            axis.line = element_line(colour = 'black', size=2, linetype='solid'),
            panel.grid.major = element_line(color='gray91', size=0.1, linetype='solid'))
    
    results_plot = ggplotly(results_plot)
    
    results_plot <- results_plot %>% layout(xaxis = x_label, yaxis = y_label,  title = 'Predicted Overall Survival Curve')
    results_plot  <- results_plot %>%
      layout(
        xaxis = list(
          ticktext = seq(from = 0, to = 114, by = 6),
          tickvals = seq(from = 0, to = 3420, by = 180), 
          tickmode = "array"
        ),
        yaxis = list(
          ticktext = seq(from = 10, to = 100, by = 10), 
          tickvals = seq(from = 10, to = 100, by = 10),
          tickmode = "array"
        ))
    
    return(results_plot)
  }
  draw_plotly_with_unknown = function(real_time, interp_y_1, interp_y_2, interp_y_3, interp_y_4, interp_y_5){
    x_label <- list(title = "Time (month)")
    y_label <- list( title = "Survival probability (%)")
    
    plotly_df <- data.frame("interp_y_1"=unlist(interp_y_1), "interp_y_2"=unlist(interp_y_2),
                            "interp_y_3"=unlist(interp_y_3), "interp_y_4"=unlist(interp_y_4), "interp_y_5"=unlist(interp_y_5))
    interp_y_min = apply(plotly_df, 1, min)
    interp_y_max = apply(plotly_df, 1, max)
    interp_y_mean = apply(plotly_df, 1, mean)
    
    plotly_minmax_df <- data.frame("time"=real_time, "interp_y_min"=interp_y_min, "interp_y_max"=interp_y_max, "interp_y_mean"=interp_y_mean)
    
    results_plot = ggplot()+ 
      geom_line(data = plotly_minmax_df, aes(x=time, y=interp_y_mean), colour = 'blue')+ 
      geom_ribbon(data = plotly_minmax_df, aes(x=time, ymin = interp_y_min, ymax = interp_y_max), alpha=0.2, colour='blue')+
      theme(panel.background = element_rect(fill = "white"),
            axis.line = element_line(colour = 'black', size=2, linetype='solid'),
            panel.grid.major = element_line(color='gray91', size=0.1, linetype='solid'), 
            legend.position = c(0.85, 0.85))
    
    results_plot = ggplotly(results_plot)
    results_plot <- results_plot %>% layout(xaxis = x_label, yaxis = y_label,  title = 'Predicted Overall Survival Curve')
    results_plot  <- results_plot %>%layout(xaxis = list(ticktext = seq(from = 0, to = 114, by = 6),
                                                         tickvals = seq(from = 0, to = 3420, by = 180),
                                                         tickmode = "array"),
                                            yaxis = list(ticktext = seq(from = 10, to = 100, by = 10),
                                                         tickvals = seq(from = 10, to = 100, by = 10),
                                                         tickmode = "array"))
    return(results_plot)
  }
  get_avgPreds_withConstFact = function(input_df, OS_rsfs, constFactor, constVal){
    input_df[[constFactor]] = as.numeric(constVal)
    interp_ys = interp_times_probs(OS_rsfs, input_df)
    interp_ys_avg = apply(interp_ys, 2, mean)
    return(interp_ys_avg)
  }
  
  draw_pt_rsf <- reactive({
    input_df <- data.frame("MGMT"= as.numeric(input$MGMT), "Age"= as.numeric(input$Age), "Extent_of_resection"= as.numeric(input$ExtResc),
                           "IDH"= as.numeric(input$IDH), "Subvent_involve"= as.numeric(input$Subvent), "Preop_KPSGP"= as.numeric(input$KPSG),
                           "TotaldoseforconedownPTV"= as.numeric(input$TotdoseCDPTV), "conedownPTV"= as.numeric(input$CDPTV), "Edema_include_plan"= as.numeric(input$Edma),
                           "FractionforconedownPTV"= as.numeric(input$Fraction), "Gender_0_female"= as.numeric(input$Gender))
    
    interp_y = interp_times_probs(OS_rsfs, input_df)
    
    interp_y_mean = apply(interp_y, 2, mean)
    interp_y_sd = apply(interp_y, 2, sd)
    interp_y_min = interp_y_mean - 2.57* interp_y_sd/5
    interp_y_max = interp_y_mean + 2.57* interp_y_sd/5
    
    results_plot = draw_plotly_func(real_time, interp_y_mean, interp_y_min,interp_y_max)
    
    return(results_plot)
    
  })
  draw_pt_rsf_cdPTVunknown <- reactive({
    input_temp_df <- data.frame("MGMT"= as.numeric(input$MGMT), "Age"= as.numeric(input$Age), "Extent_of_resection"= as.numeric(input$ExtResc),
                                "IDH"= as.numeric(input$IDH), "Subvent_involve"= as.numeric(input$Subvent), "Preop_KPSGP"= as.numeric(input$KPSG),
                                "TotaldoseforconedownPTV"= as.numeric(input$TotdoseCDPTV), "conedownPTV"= as.numeric(input$CDPTV), "Edema_include_plan"= as.numeric(input$Edma),
                                "FractionforconedownPTV"= as.numeric(input$Fraction), "Gender_0_female"= as.numeric(input$Gender))
    
    ## conedownPTV: min, 25% quantile, 50% quantile, 75% quantile, max
    interp_y_min_avg = get_avgPreds_withConstFact(input_temp_df, OS_rsfs, constFactor = "conedownPTV", constVal = 6)
    interp_y_25quant_avg = get_avgPreds_withConstFact(input_temp_df, OS_rsfs, constFactor = "conedownPTV", constVal = 74)
    interp_y_50quant_avg = get_avgPreds_withConstFact(input_temp_df, OS_rsfs, constFactor = "conedownPTV", constVal = 121)
    interp_y_75quant_avg = get_avgPreds_withConstFact(input_temp_df, OS_rsfs, constFactor = "conedownPTV", constVal = 181)
    interp_y_max_avg = get_avgPreds_withConstFact(input_temp_df, OS_rsfs, constFactor = "conedownPTV", constVal = 760)
    
    results_plot = draw_plotly_with_unknown(real_time, interp_y_min_avg,interp_y_25quant_avg, interp_y_50quant_avg, interp_y_75quant_avg, interp_y_max_avg)
    return(results_plot)
  })
  draw_pt_rsf_TotDoseunknown <- reactive({
    input_temp_df <- data.frame("MGMT"= as.numeric(input$MGMT), "Age"= as.numeric(input$Age), "Extent_of_resection"= as.numeric(input$ExtResc),
                                "IDH"= as.numeric(input$IDH), "Subvent_involve"= as.numeric(input$Subvent), "Preop_KPSGP"= as.numeric(input$KPSG),
                                "TotaldoseforconedownPTV"= as.numeric(input$TotdoseCDPTV), "conedownPTV"= as.numeric(input$CDPTV), "Edema_include_plan"= as.numeric(input$Edma),
                                "FractionforconedownPTV"= as.numeric(input$Fraction), "Gender_0_female"= as.numeric(input$Gender))
    
    ## TotaldoseforconedownPTV: min, 25% quantile, 50% quantile, 75% quantile, max
    interp_y_min_avg = get_avgPreds_withConstFact(input_temp_df, OS_rsfs, constFactor = "TotaldoseforconedownPTV", constVal = 54)
    interp_y_25quant_avg = get_avgPreds_withConstFact(input_temp_df, OS_rsfs, constFactor = "TotaldoseforconedownPTV", constVal = 60)
    interp_y_50quant_avg = get_avgPreds_withConstFact(input_temp_df, OS_rsfs, constFactor = "TotaldoseforconedownPTV", constVal = 60)
    interp_y_75quant_avg = get_avgPreds_withConstFact(input_temp_df, OS_rsfs, constFactor = "TotaldoseforconedownPTV", constVal = 60)
    interp_y_max_avg = get_avgPreds_withConstFact(input_temp_df, OS_rsfs, constFactor = "TotaldoseforconedownPTV", constVal = 70)
    
    results_plot = draw_plotly_with_unknown(real_time, interp_y_min_avg,interp_y_25quant_avg, interp_y_50quant_avg, interp_y_75quant_avg, interp_y_max_avg)
    return(results_plot)
  })
  
  output$rsf_plot <- renderPlotly({
    if (input$whole_cal_button > 0 & unknown_val_check$a == 0){
      isolate(draw_pt_rsf())
    }else if(input$whole_cal_button > 0 & unknown_val_check$a == 1){
      isolate(draw_pt_rsf_cdPTVunknown())
    }else if(input$whole_cal_button > 0 & unknown_val_check$a == 2){
      isolate(draw_pt_rsf_TotDoseunknown())
    }else if(input$whole_cal_button > 0 & unknown_val_check$a == 3){
      isolate(draw_pt_rsf())
    }
  })
  
  
  ## Print unknown options have been selected.
  output$about_unknown_button_OS <- renderText({
    if (unknown_val_check$c == 0){ return(NULL)
    }else{
      lines = "* Selected unknown options are replaced with the median / a range of values (Min ~ Max) from the patient population."
      if (unknown_val_check$c  ==1){ #unknown_val_check$a = 3
        new_line = paste0("- Fractions to cone down PTV: It is not used for OS prediction. \n")
        lines = paste(lines, new_line, sep = "\n")
      }else if (unknown_val_check$c  ==2){ #unknown_val_check$a = 2
        new_line = paste0("- Total dose to cone down PTV: Based on the patient population, a range of values has been assigned (min: 54 ~ max: 70 Gy). \n")
        lines = paste(lines, new_line, sep = "\n")
      }else if (unknown_val_check$c ==3){ #unknown_val_check$a = 2
        new_line = paste0("- Total dose to cone down PTV: Based on the patient population, a range of values has been assigned (min: 54 ~ max: 70 Gy). \n",
                          "- Fractions to cone down PTV: It is not used for OS prediction.  \n")
        lines = paste(lines, new_line, sep = "\n")
      }else if (unknown_val_check$c ==4){ #unknown_val_check$a = 1
        new_line = paste0("- Cone down PTV: Based on the patient population, a range of values has been assigned (min: 6 ~ max: 760 cm2). \n")
        lines = paste(lines, new_line, sep = "\n")
      }else if (unknown_val_check$c ==5){ #unknown_val_check$a = 1
        new_line = paste0("- Cone down PTV: Based on the patient population, a range of values has been assigned (min: 6 ~ max: 760 cm2). \n",
                          "- Fractions to cone down PTV: It is not used for OS prediction.  \n")
        lines = paste(lines, new_line, sep = "\n")
      }else if (unknown_val_check$c ==6){ #unknown_val_check$a = 1
        new_line = paste0("- Cone down PTV: Based on the patient population, a range of values has been assigned (min: 6 ~ max: 760 cm2). \n",
                          "- Total dose to cone down PTV: Based on the patient population, the median value (60 Gy) has been assigned. \n")
        lines = paste(lines, new_line, sep = "\n")
      }else if (unknown_val_check$c ==7){
        new_line = paste0("- Cone down PTV, Total dose and Fractions for cone down PTV are unknown. \n",
                          "- Total dose to cone down PTV: Based on the patient population, the median value (60 Gy) has been assigned. \n",                          
                          "- Fractions to cone down PTV: It is not used for OS prediction.  \n")
        lines = paste(lines, new_line, sep = "\n")
      }
      return(lines)}
  })
  
  
  
  
  ################### 3. Draw predicted PFS survival curve of a single patient ###################
  draw_pt_rsf_PFS <- reactive({
    input_df <- data.frame("MGMT"= as.numeric(input$MGMT), "Age"= as.numeric(input$Age), "Extent_of_resection"= as.numeric(input$ExtResc),
                           "IDH"= as.numeric(input$IDH), "Subvent_involve"= as.numeric(input$Subvent), "Preop_KPSGP"= as.numeric(input$KPSG),
                           "TotaldoseforconedownPTV"= as.numeric(input$TotdoseCDPTV), "conedownPTV"= as.numeric(input$CDPTV), "Edema_include_plan"= as.numeric(input$Edma),
                           "FractionforconedownPTV"= as.numeric(input$Fraction), "Gender_0_female"= as.numeric(input$Gender))
    
    x_label <- list(title = "Time (month)")
    y_label <- list( title = "Progression-free survival probability (%)")
    
    interp_y = interp_times_probs(PFS_rsfs, input_df)
    
    interp_y_mean = apply(interp_y, 2, mean)
    interp_y_sd = apply(interp_y, 2, sd)
    interp_y_min = interp_y_mean - 2.57* interp_y_sd/5
    interp_y_max = interp_y_mean + 2.57* interp_y_sd/5
    
    results_plot = draw_plotly_func(real_time, interp_y_mean, interp_y_min, interp_y_max)
    
    return(results_plot)
    
  })
  
  draw_pt_rsf_cdPTVunknown_PFS <- reactive({
    input_temp_df <- data.frame("MGMT"= as.numeric(input$MGMT), "Age"= as.numeric(input$Age), "Extent_of_resection"= as.numeric(input$ExtResc),
                                "IDH"= as.numeric(input$IDH), "Subvent_involve"= as.numeric(input$Subvent), "Preop_KPSGP"= as.numeric(input$KPSG),
                                "TotaldoseforconedownPTV"= as.numeric(input$TotdoseCDPTV), "conedownPTV"= as.numeric(input$CDPTV), "Edema_include_plan"= as.numeric(input$Edma),
                                "FractionforconedownPTV"= as.numeric(input$Fraction), "Gender_0_female"= as.numeric(input$Gender))
    
    ## conedownPTV: min, 25% quantile, 50% quantile, 75% quantile, max
    interp_y_min_avg = get_avgPreds_withConstFact(input_temp_df, PFS_rsfs, constFactor = "conedownPTV", constVal = 6)
    interp_y_25quant_avg = get_avgPreds_withConstFact(input_temp_df, PFS_rsfs, constFactor = "conedownPTV", constVal = 74)
    interp_y_50quant_avg = get_avgPreds_withConstFact(input_temp_df, PFS_rsfs, constFactor = "conedownPTV", constVal = 121)
    interp_y_75quant_avg = get_avgPreds_withConstFact(input_temp_df, PFS_rsfs, constFactor = "conedownPTV", constVal = 181)
    interp_y_max_avg = get_avgPreds_withConstFact(input_temp_df, PFS_rsfs, constFactor = "conedownPTV", constVal = 760)
    
    results_plot = draw_plotly_with_unknown(real_time, interp_y_min_avg,interp_y_25quant_avg, interp_y_50quant_avg, interp_y_75quant_avg, interp_y_max_avg)
    return(results_plot)
    
  })
  draw_pt_rsf_TotDoseunknown_PFS <- reactive({
    input_temp_df <- data.frame("MGMT"= as.numeric(input$MGMT), "Age"= as.numeric(input$Age), "Extent_of_resection"= as.numeric(input$ExtResc),
                                "IDH"= as.numeric(input$IDH), "Subvent_involve"= as.numeric(input$Subvent), "Preop_KPSGP"= as.numeric(input$KPSG),
                                "TotaldoseforconedownPTV"= as.numeric(input$TotdoseCDPTV), "conedownPTV"= as.numeric(input$CDPTV), "Edema_include_plan"= as.numeric(input$Edma),
                                "FractionforconedownPTV"= as.numeric(input$Fraction), "Gender_0_female"= as.numeric(input$Gender))
    
    ## TotaldoseforconedownPTV: min, 25% quantile, 50% quantile, 75% quantile, max
    interp_y_min_avg = get_avgPreds_withConstFact(input_temp_df, PFS_rsfs, constFactor = "TotaldoseforconedownPTV", constVal = 54)
    interp_y_25quant_avg = get_avgPreds_withConstFact(input_temp_df, PFS_rsfs, constFactor = "TotaldoseforconedownPTV", constVal = 60)
    interp_y_50quant_avg = get_avgPreds_withConstFact(input_temp_df, PFS_rsfs, constFactor = "TotaldoseforconedownPTV", constVal = 60)
    interp_y_75quant_avg = get_avgPreds_withConstFact(input_temp_df, PFS_rsfs, constFactor = "TotaldoseforconedownPTV", constVal = 60)
    interp_y_max_avg = get_avgPreds_withConstFact(input_temp_df, PFS_rsfs, constFactor = "TotaldoseforconedownPTV", constVal = 70)
    
    results_plot = draw_plotly_with_unknown(real_time, interp_y_min_avg,interp_y_25quant_avg, interp_y_50quant_avg, interp_y_75quant_avg, interp_y_max_avg)
    return(results_plot)
  })
  draw_pt_rsf_Fxunknown_PFS <- reactive({
    input_temp_df <- data.frame("MGMT"= as.numeric(input$MGMT), "Age"= as.numeric(input$Age), "Extent_of_resection"= as.numeric(input$ExtResc),
                                "IDH"= as.numeric(input$IDH), "Subvent_involve"= as.numeric(input$Subvent), "Preop_KPSGP"= as.numeric(input$KPSG),
                                "TotaldoseforconedownPTV"= as.numeric(input$TotdoseCDPTV), "conedownPTV"= as.numeric(input$CDPTV), "Edema_include_plan"= as.numeric(input$Edma),
                                "FractionforconedownPTV"= as.numeric(input$Fraction), "Gender_0_female"= as.numeric(input$Gender))
    
    ## TotaldoseforconedownPTV: min ~ max
    interp_y_min_avg = get_avgPreds_withConstFact(input_temp_df, PFS_rsfs, constFactor = "FractionforconedownPTV", constVal = 30)
    interp_y_25quant_avg = get_avgPreds_withConstFact(input_temp_df, PFS_rsfs, constFactor = "FractionforconedownPTV", constVal = 31)
    interp_y_50quant_avg = get_avgPreds_withConstFact(input_temp_df, PFS_rsfs, constFactor = "FractionforconedownPTV", constVal = 32)
    interp_y_75quant_avg = get_avgPreds_withConstFact(input_temp_df, PFS_rsfs, constFactor = "FractionforconedownPTV", constVal = 33)
    interp_y_max_avg = get_avgPreds_withConstFact(input_temp_df, PFS_rsfs, constFactor = "FractionforconedownPTV", constVal = 35)
    
    results_plot = draw_plotly_with_unknown(real_time, interp_y_min_avg,interp_y_25quant_avg, interp_y_50quant_avg, interp_y_75quant_avg, interp_y_max_avg)
    return(results_plot)
  })
  
  output$rsf_plot_PFS <- renderPlotly({
    if (input$whole_cal_button > 0 & unknown_val_check$b == 0){
      isolate(draw_pt_rsf_PFS())
    }else if(input$whole_cal_button > 0 & unknown_val_check$b == 1){
      isolate(draw_pt_rsf_Fxunknown_PFS())
    }else if(input$whole_cal_button > 0 & unknown_val_check$b == 2){
      isolate(draw_pt_rsf_cdPTVunknown_PFS())
    }else if(input$whole_cal_button > 0 & unknown_val_check$b == 3){
      isolate(draw_pt_rsf_TotDoseunknown_PFS())
    }
  })
  
  ## Print unknown options have been selected.
  output$about_unknown_button_PFS <- renderText({
    if (unknown_val_check$c == 0){ return(NULL)
    }else{
      lines = "* Selected unknown options are replaced with the median / a range of values (Min ~ Max) from the patient population."
      if (unknown_val_check$c  ==1){ #unknown_val_check$b = 1
        new_line = paste0("- Fractions to cone down PTV: Based on the patient population, a range of values has been assigned (min: 30 ~ max: 35 #). \n")
        lines = paste(lines, new_line, sep = "\n")
      }else if (unknown_val_check$c  ==2){ #unknown_val_check$b = 3
        new_line = paste0("- Total dose to cone down PTV: Based on the patient population, a range of values has been assigned (min: 54 ~ max: 70 Gy). \n")
        lines = paste(lines, new_line, sep = "\n")
      }else if (unknown_val_check$c ==3){ #unknown_val_check$b = 1
        new_line = paste0("- Total dose to cone down PTV: Based on the patient population, the median value (60 Gy) has been assigned. \n",
                          "- Fractions to cone down PTV: Based on the patient population, a range of values has been assigned (min: 30 ~ max: 35 #). \n")
        lines = paste(lines, new_line, sep = "\n")
      }else if (unknown_val_check$c ==4){ #unknown_val_check$b = 2
        new_line = paste0("- Cone down PTV: Based on the patient population, a range of values has been assigned (min: 6 ~ max: 760 cm2). \n")
        lines = paste(lines, new_line, sep = "\n")
      }else if (unknown_val_check$c ==5){ #unknown_val_check$b = 1
        new_line = paste0("- Cone down PTV: Based on the patient population, the median value (121 cm2) has been assigned.  \n",
                          "- Fractions to cone down PTV: Based on the patient population, a range of values has been assigned (min: 30 ~ max: 35 #). \n")
        lines = paste(lines, new_line, sep = "\n")
      }else if (unknown_val_check$c ==6){ #unknown_val_check$b = 2
        new_line = paste0("- Cone down PTV: Based on the patient population, a range of values has been assigned (min: 6 ~ max: 760 cm2). \n",
                          "- Total dose to cone down PTV: Based on the patient population, the median value (60 Gy) has been assigned. \n")
        lines = paste(lines, new_line, sep = "\n")
      }else if (unknown_val_check$c ==7){ #unknown_val_check$b = 1
        new_line = paste0("- Cone down PTV: Based on the patient population, the median value (121 cm2) has been assigned.  \n",
                          "- Total dose to cone down PTV: Based on the patient population, the median value (60 Gy) has been assigned. \n",                          
                          "- Fractions to cone down PTV: Based on the patient population, a range of values has been assigned (min: 30 ~ max: 35 #). \n")
        lines = paste(lines, new_line, sep = "\n")
      }
      return(lines)}
  })
  ################### 4. Accumulate multiple patients data & compare prediction results ###################
  
  unknown_selected_pt_info = reactiveValues()
  
  # 1) Accumulate multiple patients' data
  cumulate_pt_df <- reactiveValues()
  cumulate_pt_df$df <-data.frame("MGMT"= NULL, "Age"= NULL, "Extent_of_resection"= NULL,
                                 "IDH"= NULL, "Subvent_involve"= NULL, "Preop_KPSGP"= NULL,
                                 "TotaldoseforconedownPTV"= NULL, "conedownPTV"= NULL, "Edema_include_plan"= NULL,
                                 "FractionforconedownPTV"= NULL, "Gender_0_female"= NULL)
  ## Accumulate each patient.
  observeEvent(input$cum_button,{
    cumu_input_df <- data.frame("MGMT"= as.numeric(input$MGMT), "Age"= as.numeric(input$Age), "Extent_of_resection"= as.numeric(input$ExtResc),
                                "IDH"= as.numeric(input$IDH), "Subvent_involve"= as.numeric(input$Subvent), "Preop_KPSGP"= as.numeric(input$KPSG),
                                "TotaldoseforconedownPTV"= as.numeric(input$TotdoseCDPTV), "conedownPTV"= as.numeric(input$CDPTV), "Edema_include_plan"= as.numeric(input$Edma),
                                "FractionforconedownPTV"= as.numeric(input$Fraction), "Gender_0_female"= as.numeric(input$Gender))
    
    ## For the UNKNOWN option selected, 
    if (is_unknown_cdPTV$a == TRUE){
      cumu_input_df$conedownPTV = 121}
    if (is_unknown_TotDose$a == TRUE){
      cumu_input_df$TotaldoseforconedownPTV = 60}
    if (is_unknown_Fx$a == TRUE){
      cumu_input_df$FractionforconedownPTV = 30}
    
    cumulate_pt_df$df <- rbind(cumulate_pt_df$df, cumu_input_df)
    
    num_pt = nrow(cumulate_pt_df$df)
    unknown_selected_pt_info$a[[num_pt]] = is_unknown_cdPTV$a*4 + is_unknown_TotDose$a*2 + is_unknown_Fx$a
    
  })
  
  ## Erase last inserted patient.
  observeEvent(input$erase_button,{
    cumulate_pt_df$df <- head(cumulate_pt_df$df, -1)
  })
  ## Show how many patients were added.
  output$cum_pt_list <- renderText({
    if(nrow(cumulate_pt_df$df)== 0){
      return(NULL)
    }else {
      print_text <- NULL
      for(i in 1:nrow(cumulate_pt_df$df)){
        print_text <- paste0(print_text, "Patient no.", as.character(i), " added!", "\n")
      }
      return(print_text)
    }
  })
  
  # 2) Check multiple patients' data
  cum_pt_datainput <- reactive({  
    pt_num_df <- data.frame("pt_number"=c(1:nrow(cumulate_pt_df$df)))
    print_df <- cbind(pt_num_df, cumulate_pt_df$df)
  })
  output$table_cum_pt_df <- renderTable({
    if (input$cum_cal_button > 0){
      isolate(cum_pt_datainput())
    }
  }, width=400)
  
  output$unknown_selected_info = renderText({
    if (sum(unlist(unlist(unlist(unknown_selected_pt_info$a)))) == 0){
      return("ALL no unknown.")
    }else{
      lines = " Selected unknown options are replaced with the median from the patient population.\n (Cone down PTV: 121 cm2 / Total dose to cone down PTV: 60 Gy / Fractions to cone down PTV: 30 #)"
      for (pt_i in 1:length(unknown_selected_pt_info$a)){
        if (unknown_selected_pt_info$a[[pt_i]] ==0){
          new_line = paste0("* Pt No.", as.character(pt_i), ":  No UNKNOWN option.")
          lines = paste(lines, new_line, sep = "\n")
        }else if (unknown_selected_pt_info$a[[pt_i]] ==1){
          new_line = paste0("* Pt No.", as.character(pt_i), ":  Fractions to cone down PTV is unknown.")
          lines = paste(lines, new_line, sep = "\n")
        }else if (unknown_selected_pt_info$a[[pt_i]] ==2){
          new_line = paste0("* Pt No.", as.character(pt_i), ":  Total dose to cone down PTV is unknown.")
          lines = paste(lines, new_line, sep = "\n")
        }else if (unknown_selected_pt_info$a[[pt_i]] ==3){
          new_line = paste0("* Pt No.", as.character(pt_i), ":  Total dose and Fractions to cone down PTV are unknown.")
          lines = paste(lines, new_line, sep = "\n")
        }else if (unknown_selected_pt_info$a[[pt_i]] ==4){
          new_line = paste0("* Pt No.", as.character(pt_i), ":  Cone down PTV is unknown.")
          lines = paste(lines, new_line, sep = "\n")
        }else if (unknown_selected_pt_info$a[[pt_i]] ==5){
          new_line = paste0("* Pt No.", as.character(pt_i), ":  Cone down PTV and Fractions to cone down PTV are unknown.")
          lines = paste(lines, new_line, sep = "\n")
        }else if (unknown_selected_pt_info$a[[pt_i]] ==6){
          new_line = paste0("* Pt No.", as.character(pt_i), ":  Cone down PTV and Total dose to cone down PTV are unknown.")
          lines = paste(lines, new_line, sep = "\n")
        }else if (unknown_selected_pt_info$a[[pt_i]] ==7){
          new_line = paste0("* Pt No.", as.character(pt_i), ":  Cone down PTV, Total dose and Fractions to cone down PTV are unknown.")
          lines = paste(lines, new_line, sep = "\n")
        }
      }
      return(lines)}
  })
  
  
  # 3) Draw predicted survival curve of multiple patients
  draw_cum_pt_rsf <- reactive({
    input_df <- cumulate_pt_df$df
    
    
    plotly_df <- data.frame("time"=real_time)
    for (i in 1:nrow(input_df)){
      input_row = input_df[i,]
      
      interp_y = interp_times_probs(OS_rsfs, input_row)
      interp_y_mean = apply(interp_y, 2, mean)
      
      plotly_df[[i+1]] = interp_y_mean
    }
    
    x_label <- list(title = "Time (month)")
    y_label <- list(title = "Survival probability (%)")
    
    results_plot <- plot_ly()
    
    results_plot <- results_plot %>% layout(xaxis = x_label, yaxis = y_label, title = 'Predicted Overall Survival Curve')
    for (i in 2:ncol(plotly_df)){
      results_plot <- results_plot %>% add_trace(x=plotly_df[["time"]], y=plotly_df[[i]],type='scatter',name=paste0("pt_",i-1), mode='line') 
    }
    results_plot  <- results_plot %>%
      layout(
        xaxis = list(
          ticktext = seq(from = 6, to = 114, by = 6),
          tickvals = seq(from = 180, to = 3420, by = 180), 
          tickmode = "array"
        ),
        yaxis = list(
          ticktext = seq(from = 10, to = 100, by = 10),
          tickvals = seq(from = 10, to = 100, by = 10), 
          tickmode = "array"
        ))
    return(results_plot)
  })
  
  output$cum_rsf_plot <- renderPlotly({
    if (input$cum_cal_button > 0){
      isolate(draw_cum_pt_rsf())}})
  
  #### PFS curve. 
  draw_cum_pt_rsf_PFS <- reactive({
    input_df <- cumulate_pt_df$df
    
    plotly_df <- data.frame("time"=real_time)
    for (i in 1:nrow(input_df)){
      input_row = input_df[i,]
      
      interp_y = interp_times_probs(PFS_rsfs, input_row)
      interp_y_mean = apply(interp_y, 2, mean)
      
      plotly_df[[i+1]] = interp_y_mean
    }
    x_label <- list(title = "Time (month)")
    y_label <- list(title = "Progression-free survival probability (%)")
    
    results_plot <- plot_ly()
    
    results_plot <- results_plot %>% layout(xaxis = x_label, yaxis = y_label, title = 'Predicted Overall Survival Curve')
    for (i in 2:ncol(plotly_df)){
      results_plot <- results_plot %>% add_trace(x=plotly_df[["time"]], y=plotly_df[[i]],type='scatter',name=paste0("pt_",i-1), mode='line') 
    }
    results_plot  <- results_plot %>%
      layout(
        xaxis = list(
          ticktext = seq(from = 6, to = 114, by = 6),
          tickvals = seq(from = 180, to = 3420, by = 180), 
          tickmode = "array"
        ),
        yaxis = list(
          ticktext = seq(from = 10, to = 100, by = 10),
          tickvals = seq(from = 10, to = 100, by = 10), 
          tickmode = "array"
        ))
    
    return(results_plot)
  })
  
  output$cum_rsf_plot_PFS <- renderPlotly({
    if (input$cum_cal_button > 0){
      isolate(draw_cum_pt_rsf_PFS())
    }})}

shinyApp(ui=ui, server=server)
