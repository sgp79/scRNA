# Run this from a `R` commandline (in tmux to keep it running):  
# 
# shiny::runApp("~/scRNA/shiny_apps/cell_filter_plots.R", port = 3838, host = "0.0.0.0")
# 

library(shiny)
library(Seurat)
library(ggplot2)
library(plotly)
library(scales)
options(Seurat.object.assay.version = "v5")
options(browser = "false")

setwd(file.path("..", "..", "working", "single_batch"))

## --------------------------
## UI COMPONENTS
## --------------------------

# File selection UI
file_selection_ui <- function() {
    tagList(
            textInput("dir_path", "Directory path:", value = getwd()),
            actionButton("load_dir", "Load Directory"),
            uiOutput("file_selector")
           )
}

# Plot controls UI
plot_controls_ui <- function() {
    tagList(
            selectInput("plot_type", "Plot Type:", choices = c("VlnPlot", "ScatterPlot")),
    
            # Dynamic UI elements will be inserted here
            uiOutput("plot_specific_controls")
            )
}

# Main UI
ui <- fluidPage(
            titlePanel("Seurat Data Visualization"),
            sidebarLayout(
                    sidebarPanel(
                            file_selection_ui(),
                            actionButton("plot_button", "Generate Plot"),
                            plot_controls_ui()
                            ),
                    mainPanel(
                        tabsetPanel(
                            tabPanel("Static Plot", plotOutput("seurat_plot")),
                            tabPanel("Interactive Plot", plotlyOutput("interactive_plot")),
                            tabPanel("Metadata Summary", 
                                     div(style = 'overflow-x: auto;',
                                         tableOutput("metadata_table")
                                        )
                                    ),
                            tabPanel("Object Info", verbatimTextOutput("object_info"))
                        )))
)

## --------------------------
## SERVER MODULES
## --------------------------

# Handle file loading
file_handling <- function(input, output, session) {
    # initialise reactive variables
    seurat_obj <- reactiveVal(NULL)
    available_files <- reactiveVal(NULL)

    # list available files
    observeEvent(input$load_dir, {
                 req(input$dir_path)
                 if (!dir.exists(input$dir_path)) {
                     showNotification("Path does not exist", type = "error")
                     return()
                 }
                 # list files with .rds extension only
                 files <- list.files(input$dir_path, pattern='\\.rds', recursive = FALSE, full.names = FALSE)
                 # files <- list.files(data_dir, pattern='\\.h5seurat', recursive = FALSE, full.names = FALSE)
                 available_files(files)
    })

    # select the file to load from the file_selector options (available_files)
    output$file_selector <- renderUI({
        req(available_files())  # wait for available_files (set above
        selectInput("selected_file", "Select file:", choices = available_files())
    })

    # load the file
    observeEvent(
                input$selected_file,  # wait for a file to be selected 
                {
                    req(input$selected_file, input$dir_path)  # wait for selections to be made
                    file_path <- file.path(input$dir_path, input$selected_file)
                
                    tryCatch({
                                obj <- readRDS(file_path)
                                #obj <- SeuratDisk::LoadH5Seurat(file_path)  # load the file
                                seurat_obj(obj)  # add to server-defined reactive variable 
                                showNotification("File loaded", type = "message")
                              },
                              error = function(e) {
                                    showNotification(paste("Error loading file:", e$message), type = "error")
                                    seurat_obj(NULL)  # set reactive variable to NULL if there was an error
                              }
                            )
                }
    )
  
    return(seurat_obj)
}

                 
# Handle plot controls
plot_controls <- function(input, output, session, s_obj) {
    # Dynamic UI for plot-specific controls
    output$plot_specific_controls <- renderUI({
        req(s_obj())
        metadata_cols <- names(s_obj()@meta.data)
        # only numeric columns make sense as plot axes
        numeric_metadata_cols <- metadata_cols[sapply(s_obj()@meta.data, is.numeric)]
        switch(input$plot_type,
               "VlnPlot" = tagList(
                                selectInput("metadata_col",
                                            "Metadata Column:",
                                            choices = numeric_metadata_cols
                                           ),
                                checkboxInput("vln_log_y", "Log Y Axis"),
                                numericInput("vln_ymin", "Y Min", value = NA),
                                numericInput("vln_ymax", "Y Max", value = NA)
                                ),
               "ScatterPlot" = tagList(
                                     selectInput("metadata_col_x",
                                               "Metadata Column (x axis):",
                                               choices = numeric_metadata_cols
                                              ),
                                     selectInput("metadata_col_y",
                                                 "Metadata Column (y axis):",
                                                 choices = numeric_metadata_cols
                                                ),
                                     checkboxInput("scatter_log_x", "Log X Axis"),
                                     numericInput("scatter_xmin", "X Min", value = NA),
                                     numericInput("scatter_xmax", "X Max", value = NA),
                                     checkboxInput("scatter_log_y", "Log Y Axis"),
                                     numericInput("scatter_ymin", "Y Min", value = NA),
                                     numericInput("scatter_ymax", "Y Max", value = NA)
                                     )
                )
    })
}

# chose plot type to generate
plot_generation <- function(input, output, session, seurat_obj) {
  # Use eventReactive for plot generation
    plot_data <- eventReactive(input$plot_button, {
        req(seurat_obj())
        generate_plot(seurat_obj(), input)
    })
  
    # Render all plots as static ggplot outputs
    output$seurat_plot <- renderPlot({
        req(plot_data())
        plot_data()
    })
  
    # Add interactive version
    output$interactive_plot <- renderPlotly({
        req(plot_data())
        ggplotly(plot_data()) %>% 
          config(displayModeBar = TRUE)
      })
}


# plot generation function
generate_plot <- function(obj, input) {
    req(input$plot_type)
    switch(input$plot_type,  # picks the type of plot from below using input$plot_type
           "VlnPlot" = {
                         req(input$metadata_col)
                         p <- VlnPlot(obj, features = input$metadata_col) 
                         p <- p + theme(legend.position = "none")
                         # Apply log transformation if requested
                         if(input$vln_log_y) {
                                p <- p + scale_y_continuous(trans = "log10", breaks = pretty_breaks(n = 10))
                         }else{
                                p <- p + scale_y_continuous(breaks = pretty_breaks(n = 10))
                         }
                         # Apply axis limits if provided
                         if(!is.na(input$vln_ymin) && !is.na(input$vln_ymax)) {
                                p <- p + coord_cartesian(ylim = c(input$vln_ymin, input$vln_ymax))
                         }
                         p
                       },
           "ScatterPlot" = {
                         req(input$metadata_col_x, input$metadata_col_y)
                         p <- FeatureScatter(obj, input$metadata_col_x, input$metadata_col_y,
                                        shuffle=TRUE, pt.size = 1.5, slot = 'counts'
                                        )
                         p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
                         # Apply log transformations if requested
                         if(input$scatter_log_x) {
                                p <- p + scale_x_continuous(trans = "log10", breaks = pretty_breaks(n = 10))
                         }else{
                                p <- p + scale_x_continuous(breaks = pretty_breaks(n = 10))
                         }
                         if(input$scatter_log_y) {
                                p <- p + scale_y_continuous(trans = "log10", breaks = pretty_breaks(n = 10))
                         }else{
                                p <- p + scale_y_continuous(breaks = pretty_breaks(n = 10))
                         }
                         # Apply axis limits if provided
                         if(!is.na(input$scatter_xmin) && !is.na(input$scatter_xmax)) {
                                p <- p + coord_cartesian(xlim = c(input$scatter_xmin, input$scatter_xmax))
                         }
                         if(!is.na(input$scatter_ymin) && !is.na(input$scatter_ymax)) {
                                p <- p + coord_cartesian(ylim = c(input$scatter_ymin, input$scatter_ymax))
                         }
                         p
                       }
            )
}

# data display - the two text info tabs
data_display <- function(input, output, session, seurat_obj) {
    # Metadata table - shows the head of the cell metadata df
    output$metadata_table <- renderTable({
                                            req(seurat_obj())
                                            head(seurat_obj()@meta.data, 50)  # Show first 50 rows
                                         }, 
                                         rownames = TRUE, striped = TRUE, hover = TRUE, width = "100%"
                                        )
  
    # Object info
    output$object_info <- renderPrint({
                                        req(seurat_obj())
                                        cat("Seurat Object Information\n")
                                        cat("=========================\n\n")
                                        cat("Number of cells:", ncol(seurat_obj()), "\n")
                                        cat("Number of features:", nrow(seurat_obj()), "\n\n")
                                        cat("Assays:", paste(Seurat::Assays(seurat_obj()), collapse = ", "), "\n")
    })
}

# SERVER FUNCTION
server <- function(input, output, session) {
    seurat_obj <- file_handling(input, output, session)
    plot_controls(input, output, session, seurat_obj)
    plot_generation(input, output, session, seurat_obj)
    data_display(input, output, session, seurat_obj)
}

shinyApp(ui, server)