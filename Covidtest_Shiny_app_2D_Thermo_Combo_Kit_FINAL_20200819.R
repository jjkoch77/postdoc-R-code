#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# Authors: Joe Kochmanski & Nathan Kuhn

library(shiny)
library(tidyverse)
library(readxl)

# Define UI for data upload app ----
ui <- fluidPage(
    
    # App title ----
    titlePanel("Sars Covid Results"),
    
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
            
            img(src = "Spartan-helmet-Green-150-pxls.png", height = 140, width = 140),
            
            # Input: Select a file ----
            fileInput("file1", "Plate File",
                      multiple = TRUE,
                      accept = c("xlsx")),
            
            # Horizontal line ----
            tags$hr(),
            
            # Input: Select a file ----
            fileInput("file2", "QPCR Columns",
                      multiple = TRUE,
                      accept = c("xlsx")),
            
            # Horizontal line ----
            tags$hr(),
            
            # Input: Select a file ----
            fileInput("file3", "QPCR Rows",
                      multiple = TRUE,
                      accept = c("xlsx")),
            
            # Horizontal line ----
            tags$hr(),
            
            # Input: Choose dataset ----
            selectInput("dataset", "Choose a dataset:",
                        choices = c("Covid_Results")),
            
            # Button
            downloadButton("downloadData", "Download")
            
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
            
            # Output: Data file ----
            tableOutput("contents")
            
        )
        
    )
    
)

# Define server logic to read selected file ----
server <- function(input, output) {
    
    output$contents <- renderTable({
        
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, head of that data file by default,
        # or all rows if selected, will be shown.
        
        req(input$file1)
        req(input$file2)
        req(input$file3)

            #bring in excel files. once file location is know. R scripte can be place there to make pathway easy to write in.
            Master <- read_excel(input$file1$datapath, sheet = "Master")
            Rows_Index <- read_excel(input$file1$datapath, sheet = "Rows")
            Columns_Index <- read_excel(input$file1$datapath, sheet = "Columns")
            Columns_qPCR <-  read_excel(input$file2$datapath, sheet = "Results", skip = 42)
            Rows_qPCR <-  read_excel(input$file3$datapath, sheet = "Results", skip = 42)
            
            ##Clean up Rows qpcr data, make calls for posistive based on CT   
            Clean_Rows_qPCR <- Rows_qPCR %>%
                mutate(CT = replace_na(as.numeric(CT), 40)) %>%
                select("Well Position","Target Name", "Reporter", "CT", "Ct Mean", "Amp Status") %>%
                pivot_wider(names_from = "Target Name", values_from = c("Reporter", "CT", "Ct Mean", "Amp Status")) %>%
                mutate(N_gene_positive = ifelse(CT_N <= 37, 1, 0)) %>%
                mutate(ORF1ab_positive = ifelse(CT_ORF1ab <= 37, 1, 0)) %>%
                mutate(S_gene_positive = ifelse(CT_S <= 37, 1, 0)) %>%
                mutate(Row_Positive = ifelse(select(., N_gene_positive:S_gene_positive) %>%
                                                 rowSums() >= 2, TRUE, FALSE))
            
            #Merge Rows qPCR data with Rows location data
            Rows_QR <- Rows_Index %>%
                left_join(Clean_Rows_qPCR, by = "Well Position")
            
            #Clean up columns qpcr data, make calls for posistive based on CT
            Clean_Columns_qPCR <- Columns_qPCR %>%
                mutate(CT = replace_na(as.numeric(CT), 40)) %>%
                select("Well Position","Target Name", "Reporter", "CT", "Ct Mean", "Amp Status") %>%
                pivot_wider(names_from = "Target Name", values_from = c("Reporter", "CT", "Ct Mean", "Amp Status")) %>%
                mutate(N_gene_positive = ifelse(CT_N <= 37, 1, 0)) %>%
                mutate(ORF1ab_positive = ifelse(CT_ORF1ab <= 37, 1, 0)) %>%
                mutate(S_gene_positive = ifelse(CT_S <= 37, 1, 0)) %>%
                mutate(Columns_Positive = ifelse(select(., N_gene_positive:S_gene_positive) %>%
                                                     rowSums() >= 2, TRUE, FALSE))
            
            #Merge Columns qPCR data with columns location data
            Columns_QR <- Columns_Index %>%
                left_join(Clean_Columns_qPCR, by = "Well Position")
           
            #merge rows_qr to master
            Master_Rows_QR <- Master %>%
                left_join(Rows_QR, by = "Rows")
            #merge columns_qr to master_rows_qr
            Master_QR <- Master_Rows_QR %>%
                left_join(Columns_QR, by = "Columns")
            
            
            #check for positive 
            Master_QR_positive <<- Master_QR %>%
                mutate(Results = ifelse(Row_Positive == TRUE & Columns_Positive == TRUE, "Positive", "Negative")) %>%
                select("Barcode", "Results") %>%
                filter(str_length(Barcode)>1)%>%
                arrange(desc(Results))
            
            Master_QR_positive
          
        
    })
    # Reactive value for selected dataset ----
    datasetInput <- reactive({
        switch(input$dataset,
               "Covid_Results" =  Master_QR_positive
               )
    })
    
    # Table of selected dataset ----
    output$table <- renderTable({
        datasetInput()
    })
    
    # Downloadable csv of selected dataset ----
    output$downloadData <- downloadHandler(
        filename = function() {
            paste(input$dataset, ".csv", sep = "")
        },
        content = function(file) {
            write_excel_csv(datasetInput(), file)
        }
    )
}

# Create Shiny app ----
shinyApp(ui, server)
