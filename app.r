# Install required packages if missing
required_packages <- c("shiny", "ggplot2", "DT")
packages_to_install <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(packages_to_install) > 0){
  install.packages(packages_to_install)
  print(paste("Installing missing package(s):", paste(packages_to_install, collapse = ", ")) )
}

# Load the packages
library(shiny)
library(ggplot2)
library(DT)

# Define the user interface (UI)
ui <- fluidPage(
  
  # App title
  titlePanel("TaqMan qPCR Standard Curve Analysis"),
  
  # Sidebar layout with input and output sections
  sidebarLayout(
    
    # Sidebar panel for inputs
    sidebarPanel(
      h4("Step I: Input Your Standard Curve Data"),
      
      # An editable table for user input
      DTOutput("dataTable"),
      
      # Instructions for the user
      br(),
      helpText("Edit the Cq and SQ values directly in the table above."),
      helpText("Ensure SQ values are greater than 0."),
      
      # Move Step II to the sidebar
      hr(),
      h4("Step II: Predict Sample Quantity (SQ)"),
      helpText("Enter a Cq value from your unknown sample to predict its starting quantity (SQ)."),
      numericInput("predictCq", "Cq value:", value = 20),
      textOutput("predictedSq"),
      
      # Move Step III to the sidebar
      hr(),
      h4("Step III: Calculate Virus Titer"),
      helpText("Provide the viral genome length and dilution factor to calculate the titer."),
      div(
        HTML("<b>Formula:</b>"),
        HTML("Titer = (SQ * 10<sup>-9</sup> * 6.022e<sup>23</sup> * dilution facto) / (length * 650 r)")
      ),
      br(),
      fluidRow(
        column(6, numericInput("genomeLength", "Viral Genome Length (bp):", value = 33000, min = 1)),
        column(6, numericInput("dilutionFactor", "Dilution Factor:", value = 400, min = 1))
      ),
      textOutput("virusTiter")
    ),
    
    # Main panel for displaying outputs
    mainPanel(
      h4("Standard Curve Results"),
      
      # Display the calculated equation and metrics
      verbatimTextOutput("results"),
      
      # Display the standard curve plot
      plotOutput("standardCurvePlot")
    )
  )
)

# Define the server logic
server <- function(input, output, session) {
  
  # Reactive values to store the data and update it
  rv <- reactiveValues(
    data = data.frame(
      Cq = c(22.9, 20.51, 17.8, 16.1),
      SQ = c(0.00001, 0.01, 0.1, 1)
    )
  )
  
  # Render the editable data table
  output$dataTable <- renderDT({
    datatable(rv$data, editable = 'cell')
  })
  
  # Observer to update reactive values when the table is edited
  observeEvent(input$dataTable_cell_edit, {
    info <- input$dataTable_cell_edit
    i <- info$row
    j <- info$col
    v <- info$value
    rv$data[i, j] <- as.numeric(v)
  })
  
  # Reactive expression to perform the analysis
  analysis <- reactive({
    data <- rv$data
    
    # Check for valid data (SQ must be > 0)
    if (any(data$SQ <= 0) || any(is.na(data$Cq)) || any(is.na(data$SQ))) {
      return(NULL)
    }
    
    # Calculate log10(SQ)
    data$logSQ <- log10(data$SQ)
    
    # Perform linear regression
    lm_model <- lm(Cq ~ logSQ, data = data)
    
    # Extract key metrics
    slope <- coef(lm_model)[2]
    intercept <- coef(lm_model)[1]
    r_squared <- summary(lm_model)$r.squared
    efficiency <- 10^(-1/slope) - 1
    
    # Return a list of results
    list(
      model = lm_model,
      slope = slope,
      intercept = intercept,
      r_squared = r_squared,
      efficiency = efficiency,
      data = data
    )
  })
  
  # Render the results text
  output$results <- renderPrint({
    results <- analysis()
    if (is.null(results)) {
      cat("Please ensure all Cq and SQ values are numeric and SQ > 0.\n")
      return()
    }
    
    cat("Standard Curve Equation:\n")
    cat("Cq =", round(results$slope, 4), " * log10(SQ) +", round(results$intercept, 4), "\n\n")
    
    cat("Regression Metrics:\n")
    cat("Slope (m):", round(results$slope, 4), "\n")
    cat("Y-intercept (b):", round(results$intercept, 4), "\n")
    cat("R-squared (R²):", round(results$r_squared, 4), "\n\n")
    
    cat("PCR Efficiency:\n")
    cat("Efficiency:", round(results$efficiency * 100, 2), "%\n")
    cat("(Ideal efficiency is 90% - 110%)\n")
  })
  
  # Render the standard curve plot
  output$standardCurvePlot <- renderPlot({
    results <- analysis()
    if (is.null(results)) {
      return()
    }
    
    ggplot(results$data, aes(x = logSQ, y = Cq)) +
      geom_point(color = "blue", size = 3) +
      geom_smooth(method = "lm", se = FALSE, color = "red") +
      labs(
        title = "Standard Curve",
        x = "Log10(Starting Quantity)",
        y = "Quantification Cycle (Cq)"
      ) +
      theme_minimal(base_size = 14)
  })
  
  # Reactive expression to predict SQ from a given Cq
  predicted_sq <- reactive({
    req(analysis())
    req(input$predictCq)
    
    m <- analysis()$slope
    b <- analysis()$intercept
    cq_val <- input$predictCq
    
    log_sq <- (cq_val - b) / m
    sq <- 10^log_sq
    
    return(sq)
  })
  
  # Render the predicted SQ to the UI
  output$predictedSq <- renderText({
    req(predicted_sq())
    
    paste("Predicted SQ:", format(predicted_sq(), scientific = TRUE, digits = 4))
  })
  
  # Reactive expression to calculate virus titer
  virus_titer <- reactive({
    req(predicted_sq())
    req(input$genomeLength)
    req(input$dilutionFactor)
    
    # Get inputs
    sq <- predicted_sq()
    length <- input$genomeLength
    dilution_factor <- input$dilutionFactor
    
    # Calculate titer using the user's specified formula
    # Note: 650 is the average molecular weight of a single base pair in g/mol.
    # 6.022e23 is Avogadro's number (molecules per mole).
    # The calculation assumes specific unit conversions.
    titer_val <- (sq * 1e-9 * 6.022e23 * dilution_factor) / (length * 650)
    
    return(titer_val)
  })
  
  # Render the virus titer to the UI
  output$virusTiter <- renderText({
    req(virus_titer())
    
    paste("Virus Titer:", format(virus_titer(), scientific = TRUE, digits = 4), "genomic copies/mL")
  })
}

# Run the app
shinyApp(ui, server)