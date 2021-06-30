library(shiny);library(shinycustomloader);library(ggplot2);library(survival);library(jsmodule)
source("global.R")
nfactor.limit <- 20

ui <- navbarPage("Lithium",
                 theme = bslib::bs_theme(version = 4),
                 tabPanel("Data", icon = icon("table"),
                          sidebarLayout(
                            sidebarPanel(
                              radioButtons("data_main", "Data", c("All", "Remove eGFR <60 within 1yrs", "Remove eGFR <60 within 3yrs"), "All", inline = T)
                            ),
                            mainPanel(
                              tabsetPanel(type = "pills",
                                          tabPanel("Data", withLoader(DTOutput("data"), type="html", loader="loader6")),
                                          tabPanel("Label", withLoader(DTOutput("data_label", width = "100%"), type="html", loader="loader6"))
                              )
                              
                            )
                          )
                 ),
                 tabPanel("N Profile",
                          sidebarLayout(
                            sidebarPanel(
                              
                            ),
                            mainPanel(
                              tabsetPanel(type = "pills",
                                          tabPanel("Inclusion & Exclusion", withLoader(DTOutput("NProfile"), type="html", loader="loader6")),
                                          tabPanel("year: N", withLoader(DTOutput("YearN"), type="html", loader="loader6")),
                                          tabPanel("eGFR<60", withLoader(DTOutput("eGFRbelow60"), type="html", loader="loader6")))
                            )
                          )
                 ),
                 tabPanel("Table 1", icon = icon("percentage"),
                          sidebarLayout(
                            sidebarPanel(
                              radioButtons("sub_tb1", "Subgroup", c("All", "Lithium", "Valproate"), "All", inline = T),
                              tb1moduleUI("tb1")
                            ),
                            mainPanel(
                              withLoader(DTOutput("table1"), type="html", loader="loader6"),
                              wellPanel(
                                h5("Normal continuous variables  are summarized with Mean (SD) and t-test(2 groups) or ANOVA(> 2 groups)"),
                                h5("Non-normal continuous variables are summarized with median [IQR or min,max] and kruskal-wallis test"),
                                h5("Categorical variables  are summarized with table")
                              )
                            )
                          )
                          
                 ),
                 tabPanel("Logistic regression",
                          sidebarLayout(
                            sidebarPanel(
                              regressModuleUI("logistic")
                            ),
                            mainPanel(
                              withLoader(DTOutput("logistictable"), type="html", loader="loader6")
                            )
                          )
                 ),
                 tabPanel("Figure 1", icon = icon("percentage"),
                          sidebarLayout(
                            sidebarPanel(
                              checkboxInput("onefig1", "with 1 plot", F)
                            ),
                            mainPanel(
                              withLoader(plotOutput("fig1"), type="html", loader="loader6"),
                              h3("Download options"),
                              wellPanel(
                                uiOutput("downloadControls_fig1"),
                                downloadButton("downloadButton_fig1", label = "Download the plot")
                              )
                            )
                          )
                          
                 ),
                 tabPanel("Kaplan-meier plot",
                          sidebarLayout(
                            sidebarPanel(
                              kaplanUI("kaplan")
                            ),
                            mainPanel(
                              optionUI("kaplan"),
                              withLoader(plotOutput("kaplan_plot"), type="html", loader="loader6"),
                              ggplotdownUI("kaplan")
                            )
                          )
                 ),
                 tabPanel("Cox model",
                          sidebarLayout(
                            sidebarPanel(
                              coxUI("cox")
                            ),
                            mainPanel(
                              withLoader(DTOutput("coxtable"), type="html", loader="loader6")
                            )
                          )
                 )
)



server <- function(input, output, session) {
  
  data.main1 <- reactive({
    switch(input$data_main, 
           "All"= data.main,
           "Remove eGFR <60 within 1yrs" = data.main[!(eGFRbelow60 == 1 & year_FU < 1)],
           "Remove eGFR <60 within 3yrs" = data.main[!(eGFRbelow60 == 1 & year_FU < 3)])
  })
  
  output$data <- renderDT({
    datatable(data.main1(), rownames = F, caption = "Data",
              options = c(opt.data("data"), list(scrollX = T)))
  })
  
  output$data_label <- renderDT({
    datatable(label.main, rownames=F, extensions= "Buttons", caption = "Label of data",
              options = c(opt.data("label"), list(scrollX = TRUE))
    )
  })
  
  output$NProfile <- renderDT({
    datatable(N_profile)
  })
  
  output$YearN <- renderDT({
    datatable(Year_N, rownames=FALSE)
  })
  
  output$eGFRbelow60 <- renderDT({
    datatable(eGFRbelow60ratio)
  })
  
  data.tb1 <- reactive({
    switch(input$sub_tb1, 
           "All" = data.main1(),
           "Lithium" = data.main1()[drug == 1],
           "Valproate" = data.main1()[drug == 0])
  })
  
  out_tb1 <- callModule(tb1module2, "tb1", data = data.tb1, data_label = reactive(label.main), data_varStruct = NULL, nfactor.limit = nfactor.limit)
  
  output$table1 <- renderDT({
    tb = out_tb1()$table
    cap = out_tb1()$caption
    out.tb1 = datatable(tb, rownames = T, extensions= "Buttons", caption = cap,
                        options = c(opt.tb1("tb1"),
                                    list(columnDefs = list(list(visible=FALSE, targets= which(colnames(tb) %in% c("test","sig"))))
                                    ),
                                    list(scrollX = TRUE)
                        )
    )
    if ("sig" %in% colnames(tb)){
      out.tb1 = out.tb1 %>% formatStyle("sig", target = 'row' ,backgroundColor = styleEqual("**", 'yellow'))
    }
    return(out.tb1)
  })
  
  out_logistic <- callModule(logisticModule2, "logistic", data = data.main1, data_label = reactive(label.main), data_varStruct = reactive(varlist_kmcox), nfactor.limit = nfactor.limit, default.unires = F)
  
  output$logistictable <- renderDT({
    hide = which(colnames(out_logistic()$table) == "sig")
    datatable(out_logistic()$table, rownames=T, extensions = "Buttons", caption = out_logistic()$caption,
              options = c(jstable::opt.tbreg(out_logistic()$caption),
                          list(columnDefs = list(list(visible=FALSE, targets =hide))
                          ),
                          list(scrollX = TRUE)
              )
    ) %>% formatStyle("sig", target = 'row',backgroundColor = styleEqual("**", 'yellow'))
  })
  
  
  
  obj.fig1 <- reactive({
    zz <- data.f1[, .SD]
    zz$drug <- ifelse(zz$drug == 1, "Lithium", "Valproate")
    if (input$onefig1 == F){
      ggplot(zz, aes(x=cumulativePrescriptionDay/365.25,y=eGFR))+
        geom_point(color="coral2",size=0.1)+
        geom_smooth()+ theme_bw() + xlab("Cumulative years") + facet_wrap(~ drug)
    } else{
      ggplot(zz, aes(x=cumulativePrescriptionDay/365.25,y=eGFR, color = drug))+
        geom_point(aes(fill=drug), size=0.1)+
        geom_smooth()+ theme_bw() + xlab("Cumulative years")
    }
    
  })
  
  
  output$fig1 <- renderPlot({
    obj.fig1()
  })
  
  
  output$downloadControls_fig1 <- renderUI({
    fluidRow(
      column(4,
             selectizeInput("fig1_file_ext", "File extension (dpi = 300)", 
                            choices = c("jpg","pdf", "tiff", "svg", "emf"), multiple = F, 
                            selected = "jpg"
             )
      ),
      column(4,
             sliderInput("fig_width_fig1", "Width (in):",
                         min = 5, max = 20, value = 8
             )
      ),
      column(4,
             sliderInput("fig_height_fig1", "Height (in):",
                         min = 5, max = 20, value = 6
             )
      )
    )
  })
  
  output$downloadButton_fig1 <- downloadHandler(
    filename =  function() {
      paste("fig1.", input$fig1_file_ext ,sep="")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      withProgress(message = 'Download in progress',
                   detail = 'This may take a while...', value = 0, {
                     for (i in 1:15) {
                       incProgress(1/15)
                       Sys.sleep(0.01)
                     }
                     
                     if (input$fig1_file_ext == "emf"){
                       devEMF::emf(file, width = input$fig_width_fig1, height =input$fig_height_fig1, coordDPI = 300, emfPlus = F)
                       plot(obj.fig1())
                       dev.off()
                       
                     } else{
                       ggsave(file, obj.fig1(), dpi = 300, units = "in", width = input$fig_width_fig1, height =input$fig_height_fig1)
                     }
                     
                   })
      
    })
  
  out_kaplan <- callModule(kaplanModule, "kaplan", data = data.main1, data_label = reactive(label.main), data_varStruct = reactive(varlist_kmcox), nfactor.limit = nfactor.limit)
  
  output$kaplan_plot <- renderPlot({
    print(out_kaplan())
  })
  
  out_cox <- callModule(coxModule, "cox", data = data.main1, data_label = reactive(label.main), data_varStruct = reactive(varlist_kmcox), default.unires = F, nfactor.limit = nfactor.limit)
  
  output$coxtable <- renderDT({
    hide <- which(colnames(out_cox()$table) == c("sig"))
    datatable(out_cox()$table, rownames=T, extensions= "Buttons", caption = out_cox()$caption,
              options = c(opt.tbreg(out_cox()$caption),
                          list(columnDefs = list(list(visible=FALSE, targets= hide))
                          )
              )
    )  %>% formatStyle("sig", target = 'row',backgroundColor = styleEqual("**", 'yellow'))
  })
  
  
}


shinyApp(ui, server)