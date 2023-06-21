library(caret)
library(ranger)
library(e1071)
library(pROC)
################################################ Interface ######################################

library(shiny)
library(magrittr)

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    tags$head(tags$style("#loading{color: red;}")),
    
    tabsetPanel(tabPanel("Model",
    # Application title
    titlePanel("Utility for Randomforest predicting resistance to ceftriaxone"),
    # tags$head(
    #     tags$style(type="text/css", "label{ display: table-cell; text-align: center; vertical-align: middle; } .form-group { display: table-row;}")
    # ),
    # Sidebar with a slider input for number of bins 
    fluidRow(
        ########################################## Inputs ############################################
        column(4,
            wellPanel(
                fluidRow(
                    column(6, h4("Parameters")),
                    column(6, align="right",
                        actionButton("button", "Refresh", value = 0, style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                    )
                ),
                splitLayout(cellWidths = c("45%", "55%"), cellArgs = list(style = "padding: 6px"),
                            verticalLayout(h5("Drug cost (USD)"),
                                           fluidRow(
                                               column(6, align = "right", offset = 0, style='padding:2px;',
                                                      p("Imipenem:")),
                                               column(5, align = "left", offset = 0, style='padding:0px;',
                                                      numericInput("imicost", "X:", label = NULL, value = 7.5, width = 80))
                                           ),
                                           fluidRow(
                                               column(6, align = "right", offset = 0, style='padding:2px;',
                                                      p("Ceftriaxone:")),
                                               column(5, align = "left", offset = 0, style='padding:0px;',
                                                      numericInput("ceftbacost", label = NULL, value = 0.45, width = 80))
                                           ),
                                           fluidRow(
                                               column(6, align = "right", offset = 0, style='padding:2px;',
                                                      HTML("Treatment <br />days:")),
                                               column(5, align = "left", offset = 0, style='padding:0px;',
                                                      numericInput("treatdays", label = NULL, value = 3, width = 80))
                                           )
                            ),
                            verticalLayout(h5("Risk of death"),
                                           fluidRow(
                                               column(7, align = "right", offset = 0, style='padding:2px;',
                                                      HTML("Correct <br />treatment:")),
                                               column(4, align = "left", offset = 0, style='padding:0px;',
                                                      numericInput("riskcorrect", "X:", label = NULL, value = 0.0700, width = 90))
                                           ),
                                           fluidRow(
                                               column(7, align = "right", offset = 0, style='padding:2px;',
                                                      HTML("Incorrect <br />treatment:")),
                                               column(4, align = "left", offset = 0, style='padding:0px;',
                                                      numericInput("riskwrong", label = NULL, value = 0.1050, width = 90))
                                           ),
                                           fluidRow(
                                               column(7, align = "right", offset = 0, style='padding:2px;',
                                                      HTML("Years of life lost <br />for each death:")),
                                               column(4, align = "left", offset = 0, style='padding:0px;',
                                                      numericInput("lyl", label = NULL, value = 53.5, width = 90))
                                           )
                            )
                ),
                splitLayout(cellWidths = c("45%", "55%"), cellArgs = list(style = "padding: 6px"),
                            verticalLayout(HTML("<h5>Willingness to pay <br />(USD)</h>"),
                                           fluidRow(
                                               column(5, align = "right", offset = 0, style='padding:2px;',
                                                      HTML("GDP per <br />capita:")),
                                               column(6, align = "left", offset = 0, style='padding:0px;',
                                                      numericInput("wtpgdp", "X:", label = NULL, value = 1269.9, width = 90))
                                           ),
                                           fluidRow(
                                               column(7, align = "right", offset = 0, style='padding:2px;',
                                                      HTML("Avoid <br />carbapenem:")),
                                               column(4, align = "left", offset = 0, style='padding:0px;',
                                                      numericInput("wtpcarb", label = NULL, value = 200, width = 90))
                                           )
                            ),
                            verticalLayout(HTML("<h5>Impact of delayed <br />appropriate treatment</h>"),
                                           fluidRow(
                                               column(7, align = "right", offset = 0, style='padding:2px;',
                                                      HTML("Added hosp. <br />stay days:")),
                                               column(4, align = "left", offset = 0, style='padding:0px;',
                                                      numericInput("delaydays", "X:", label = NULL, value = 2, width = 75))
                                           ),
                                           fluidRow(
                                               column(7, align = "right", offset = 0, style='padding:2px;',
                                                      HTML("Cost per <br />bed day (USD):")),
                                               column(4, align = "left", offset = 0, style='padding:0px;',
                                                      numericInput("delaycost", label = NULL, value = 7.65, width = 75))
                                           )
                            )
                ),
                ########################################## Stats ############################################
                verticalLayout(h5("Cost weight (counts)"),
                               #verbatimTextOutput("debug"),
                               verbatimTextOutput("stat"),
                               verbatimTextOutput("stat_general"),
                               verbatimTextOutput("stat_condi")
                )
            )
        ),
        
        # Show a plot of the generated distribution
        # https://stackoverflow.com/questions/46735135/how-to-hide-a-conditional-panel-in-shiny
        column(8,
            conditionalPanel("output.hide_panel",
                textOutput("loading")
            ),
            fluidRow(
                splitLayout(cellWidths = c("50%", "50%"), plotOutput("disttest"), plotOutput("disttrain")),
                plotOutput("cost")
            )
        )
    )
    ),
    tabPanel("Explain", 
            uiOutput("pdfview")
        )
    )
)

######################################## Internal ################################################

#Import dataset
load("shiny_dat.Rda")

# https://github.com/joyofdata/joyofdata-articles/blob/master/roc-auc/plot_pred_type_distribution.R
plot_pred_type_distribution <- function(df, name, threshold, positive, negative) {
    v <- rep(NA, nrow(df))
    v <- ifelse(df$pred >= threshold & df$obs == positive, "TP", v)
    v <- ifelse(df$pred >= threshold & df$obs == negative, "FP", v)
    v <- ifelse(df$pred < threshold & df$obs == positive, "FN", v)
    v <- ifelse(df$pred < threshold & df$obs == negative, "TN", v)
    
    df$pred_type <- v
    
    p <- ggplot(data=df, aes(x=obs, y=pred)) + 
        geom_violin(fill=rgb(1,1,1,alpha=0.6), color=NA) + 
        geom_jitter(aes(color=pred_type), alpha=0.6) +
        geom_hline(yintercept=threshold, color="red", alpha=0.6) +
        scale_color_discrete(name = "type") +
        labs(title=paste(name)) +
        theme(aspect.ratio=1) +
        xlab("observed") + ylab("prediction (probability not susceptible)") + scale_x_discrete(labels= c('not susceptible', 'susceptible'))
    
    return(p)
}

# https://github.com/joyofdata/joyofdata-articles/blob/master/roc-auc/calculate_roc.R
calculate_roc <- function(df, cost_fp, cost_fn, cost_tp, cost_tn, positive, negative, n = 1000) {
    tpr_func <- function(df, threshold) {
        sum(df$pred >= threshold & df$obs == positive) / sum(df$obs == positive)
    }
    
    fpr_func <- function(df, threshold) {
        sum(df$pred >= threshold & df$obs == negative) / sum(df$obs == negative)
    }
    
    cost <- function(df, threshold, cost_fp, cost_fn, cost_tp, cost_tn) {
        sum(df$pred >= threshold & df$obs == negative) * cost_fp + 
            sum(df$pred < threshold & df$obs == positive) * cost_fn +
            sum(df$pred >= threshold & df$obs == positive) * cost_tp + 
            sum(df$pred < threshold & df$obs == negative) * cost_tn
    }
    
    roc <- data.frame(threshold = seq(0, 1 ,length.out = n), tpr=NA, fpr=NA)
    roc$tpr <- sapply(roc$threshold, function(th) tpr_func(df, th))
    roc$fpr <- sapply(roc$threshold, function(th) fpr_func(df, th))
    roc$cost <- sapply(roc$threshold, function(th) cost(df, th, cost_fp, cost_fn, cost_tp, cost_tn))
    
    return(roc)
}

require(grid)
require(gridExtra)
require(RColorBrewer)

norm_vec <- function(v){
    
    if(diff(range(v)) != 0){
        v <- (v - min(v))/diff(range(v))
    }else{
        v <- rep(0, length(v))
    }
    
    return(v)
}

######################################### Server ##########################################

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    scale <- 100
    
    imicost <- eventReactive(c(input$button), {imicost <- input$imicost})
    ceftbacost <- eventReactive(c(input$button), {ceftbacost <- input$ceftbacost})
    treatdays <- eventReactive(c(input$button), {treatdays <- input$treatdays})
    
    riskcorrect <- eventReactive(c(input$button), {riskcorrect <- input$riskcorrect})
    riskwrong <- eventReactive(c(input$button), {riskwrong <- input$riskwrong})
    lyl <- eventReactive(c(input$button), {lyl <- input$lyl})
    
    wtpgdp <- eventReactive(c(input$button), {wtpgdp <- input$wtpgdp})
    wtpcarb<- eventReactive(c(input$button), {wtpcarb <- input$wtpcarb})
    
    delaydays <- eventReactive(c(input$button), {delaydays <- input$delaydays})
    delaycost <- eventReactive(c(input$button), {delaycost <- input$delaycost})
    
    ### Show warning parameter values have changed
    old_eval <- 0
    
    output$hide_panel <- eventReactive(c(input$imicost, input$ceftbacost, input$treatdays,
                   input$riskcorrect, input$riskwrong, input$lyl,
                   input$wtpgpd, input$wtpcarb,
                   input$delaydays, input$dwlaycost, 
                   input$button), {
                       if(input$button > old_eval){
                           old_eval <<- old_eval+1
                           #print("button")
                           FALSE
                       }else{
                           #print("input")
                           TRUE
                       }
                   }, ignoreInit=TRUE)
    outputOptions(output, "hide_panel", suspendWhenHidden = FALSE)
    
    ### Show warning values exceed
    isAboveZero <- reactive({
        input$imicost>=0 &&
            input$ceftbacost>=0 &&
            input$treatdays>=0 &&
            input$riskcorrect>=0 && 
            input$riskwrong>=0 &&
            input$lyl>=0 &&
            input$wtpgdp>0 &&
            input$wtpcarg>0 &&
            input$delaydays>0 &&
            input$delaycost>0
    })
    
    isProb <- reactive({input$riskcorrect<=1 && input$riskwrong<=1})
    
    # imicost <- imicost_d %>% debounce(1000)
    # ceftbacost <- ceftbacost_d %>% debounce(1000)
    # treatdays <- treatdays_d %>% debounce(1000)
    # 
    # riskcorrect <- riskcorrect_d %>% debounce(1000)
    # riskwrong <- riskwrong_d %>% debounce(1000)
    # lyl <- lyl_d %>% debounce(1000)
    # 
    # wtpgdp <- wtpgdp_d %>% debounce(1000)
    # wtpcarb <- wtpcarb_d %>% debounce(1000)
    # 
    # delaydays <- delaydays_d %>% debounce(1000)
    # delaycost <- delaycost_d %>% debounce(1000)
    
    # Calculate cost weights
    cost_tp <- reactive({cost_tp <- (imicost()*treatdays() + wtpcarb() + 
            wtpgdp()*riskcorrect()*lyl())/scale})
    
    cost_tn <- reactive({(ceftbacost()*treatdays() + 
            wtpgdp()*riskcorrect()*lyl())/scale})
    
    cost_fp <- reactive({(imicost()*treatdays() + wtpcarb() + 
            wtpgdp()*riskcorrect()*lyl())/scale})
    
    cost_fn <- reactive({(ceftbacost()*treatdays() + delaycost()*delaydays() +
            wtpgdp()*riskwrong()*lyl())/scale})
    
    #output$debug <- renderText({wtpcarb()})
    
    # Variabels for stats
    total_pop <- nrow(cro_test_dat)
    condi_positive <- sum(cro_test_dat$outcome == 'N')
    condi_negative <- sum(cro_test_dat$outcome == 'Y')
    
    pred_test <- predict(cro_model, newdata = cro_test_dat, type = "prob")
    test_roc_df <- data.frame(pred_test[, 'N'], cro_test_dat[, "outcome"])
    names(test_roc_df) <- c('pred', 'obs')
    
    threshold <- reactive({
        roc_df_test <- calculate_roc(test_roc_df, cost_fp(), cost_fn(),
                                     cost_tp(), cost_tn(), 'N', 'Y')
        
        threshold <- roc_df_test$threshold[which(roc_df_test$cost  == min(roc_df_test$cost))[1]]
    })   
    
    
    tp <- reactive({tp <- sum(test_roc_df$pred >= threshold() & test_roc_df$obs == 'N')})
    fp <- reactive({fp <- sum(test_roc_df$pred < threshold() & test_roc_df$obs == 'N')})
    tn <- reactive({tn <- sum(test_roc_df$pred < threshold() & test_roc_df$obs == 'Y')})
    fn <- reactive({fn <- sum(test_roc_df$pred >= threshold() & test_roc_df$obs == 'Y')})
    

    output$stat <- renderText({
        
        paste0("TP = ", format(cost_tp(), digits = 4, trim = 1), "(", tp(), ")", "\t\t", "FP = ", format(cost_fp(), digits = 4, trim = 1), "(", fp(), ")",
        "\n",
        "TN = ", format(cost_tn(), digits = 4, trim = 1), "(", tn(), ")", "\t\t", "FN = ", format(cost_fn(), digits = 4, trim = 1), "(", fn(), ")")})
    
    output$stat_general <- renderText({
        
        acc <- (tp()+tn())/total_pop
        prev <- condi_positive/total_pop
        
        paste0("Prevalence:", "\t\t", "Accuracy:",
        "\n",
        format(prev, digits=2, trim=1), "\t\t\t", format(acc, digits = 2, trim = 1))})
    
    output$stat_condi <- renderText({
        
        tpr_out <- tp()/condi_positive
        fpr_out <- fp()/condi_negative
        fnr <- fp()/condi_positive
        tnr <- tn()/condi_negative
        
        paste0("Sensitivity:", "\t", "False Positive Rate:",
        "\n",
        format(tpr_out, digits = 2, trim = 1), "\t\t\t", format(fpr_out, digits = 2, trim = 1),
        "\n\n",
        "False negative rate:", "\t", "Specificity:",
        "\n",
        format(fnr, digits = 2, trim = 1), "\t\t\t", format(tnr, digits = 2, trim = 1))})
    
    # Move this up to above format
        #ppv <- tp/sum(test_roc_df$obs == 'N')
        #fdr <- fp/sum(test_roc_df$obs == 'N')
        #For <- fn/sum(test_roc_df$obs == 'Y')
        #npv <- tn/sum(test_roc_df$obs == 'Y')
    
    output$disttest <- renderPlot({
        
        pred_test <- predict(cro_model, newdata = cro_test_dat, type = "prob")
        cal_dat_test <- data.frame(cro_test_dat[, "outcome"], pred_test[, 'N'])
        names(cal_dat_test) <- c('obs', 'pred')
        pred_dist_plot_test <- plot_pred_type_distribution(cal_dat_test, "Test", threshold(), 'N', 'Y')
        print(pred_dist_plot_test)
        
    })
    
    output$disttrain <- renderPlot({
        
        pred_train <- predict(cro_model, newdata = cro_train_dat, type= "prob")
        cal_dat_train <- data.frame(cro_train_dat[, "outcome"], pred_train[, 'N'])
        names(cal_dat_train) <- c('obs', 'pred')
        pred_dist_plot_train <- plot_pred_type_distribution(cal_dat_train, "Train", threshold(), 'N', 'Y')
        print(pred_dist_plot_train)
    })
    
    output$cost <- renderPlot({
        
        roc_df_test <- calculate_roc(test_roc_df, cost_fp(), cost_fn(),
                                     cost_tp(), cost_tn(), 'N', 'Y')
        
        idx_thres_test = which.min(abs(roc_df_test$threshold-threshold()))
        
        ramp_test <- colorRampPalette(c("green","orange","red","black"))(100)
        
        test_by_cost <- ramp_test[ceiling(norm_vec(roc_df_test$cost)*99)+1]
        p_roc <- ggplot(roc_df_test, aes(fpr,tpr)) + 
            #geom_line(color=rgb(0,0,1,alpha=0.3)) +
            geom_point(data = roc_df_test, color=test_by_cost, size=3, alpha=0.5) +
            #geom_point(data = roc_df_train, color=train_by_cost, size=3, alpha=0.5) +
            coord_fixed() +
            geom_line(aes(threshold, threshold), color=rgb(0,0,1,alpha=0.5)) +
            labs(title = paste("ROC")) + xlab("False positive rate (1- specificity)") + ylab("True positive rate (sensitivity)") +
            # geom_hline(yintercept=roc[idx_threshold,"tpr"], alpha=0.5, linetype="dotdash") +
            # geom_vline(xintercept=roc[idx_threshold,"fpr"], alpha=0.5, linetype="dotdash") +
            geom_hline(yintercept=roc_df_test[idx_thres_test,"tpr"], linetype="longdash", alpha=0.5) +
            geom_vline(xintercept=roc_df_test[idx_thres_test,"fpr"], linetype="longdash", alpha=0.5)
        # geom_hline(yintercept=roc_df_train[idx_thres_train,"tpr"], alpha=0.5, linetype="dotted") +
        # geom_vline(xintercept=roc_df_train[idx_thres_train,"fpr"], alpha=0.5, linetype="dotted")
        
        p_cost <- ggplot(roc_df_test, aes(threshold, cost, color=cost)) +
            geom_line(color=rgb(0,0,1,alpha=0.3)) +
            geom_point() +
            scale_color_gradientn(colours=c("green","orange","red","black"), 
                                  breaks=c(min(roc_df_test$cost),max(roc_df_test$cost)),labels=c("min","max")) +
            #geom_point(data = roc_df_train, color=train_by_cost, size=3, alpha=0.5) +
            labs(title = sprintf("Cost function")) +
            #geom_vline(xintercept=threshold, linetype="longdash") +
            geom_vline(xintercept=roc_df_test$threshold[which(roc_df_test$cost  == min(roc_df_test$cost))[1]], 
                       alpha=0.5, linetype="longdash") +
            ylab("Net monetary value")
        
        grid.arrange(p_roc, p_cost, ncol=2)
        
    })
    
    #https://stackoverflow.com/questions/4785657/r-how-to-draw-an-empty-plot
    output$loading <- renderText("*Press Refresh to update plot")
    
    #https://stackoverflow.com/questions/19469978/displaying-a-pdf-from-a-local-drive-in-shiny
    output$pdfview <- renderUI({
        tags$iframe(style="position:absolute; height:100%; width:100%", src="S4_Appendix.pdf")
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

