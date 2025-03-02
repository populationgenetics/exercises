library(shiny)
library(ggplot2)

# Define UI for application
ui <- fluidPage(
  titlePanel("Allele Frequency Trajectory Simulation"),
  sidebarLayout(
    sidebarPanel(
      numericInput("init_freq", "Initial Allele Frequency:", 0.5, min = 0, max = 1),
      numericInput("selection_coeff", "Selection Coefficient (s):", 0.1, min = 0, max = 1),
      numericInput("pop_size", "Population Size (N):", 1000, min = 100),
      numericInput("generations", "Number of Generations:", 50, min = 10),
      numericInput("num_simulations", "Number of Simulations:", 5, min = 1),
      numericInput("y_min", "Y-axis Minimum:", 0, min = 0, max = 1),
      numericInput("y_max", "Y-axis Maximum:", 1, min = 0, max = 1),
      actionButton("simulate", "Simulate")
    ),
    mainPanel(
      plotOutput("freqPlot")
    )
  )
)

# Define server logic
server <- function(input, output) {
  simulate_data <- eventReactive(input$simulate, {
    init_freq <- input$init_freq
    selection_coeff <- input$selection_coeff
    pop_size <- input$pop_size
    generations <- input$generations
    num_simulations <- input$num_simulations

    simulations <- lapply(1:num_simulations, function(i) {
      allele_freq <- numeric(generations)
      allele_freq[1] <- init_freq

      for (j in 2:generations) {
        allele_count <- rbinom(1, pop_size, allele_freq[j-1])
        allele_freq[j] <- allele_count / pop_size
        allele_freq[j] <- allele_freq[j] + selection_coeff * allele_freq[j] * (1 - allele_freq[j])
      }

      data.frame(Generation = 1:generations, Frequency = allele_freq, Simulation = i)
    })

    do.call(rbind, simulations)
  })

  output$freqPlot <- renderPlot({
    req(simulate_data())
    ggplot(simulate_data(), aes(x = Generation, y = Frequency, group = Simulation, color = as.factor(Simulation))) +
      geom_line(alpha = 0.7) +
      labs(title = "Allele Frequency Trajectory", x = "Generation", y = "Allele Frequency", color = "Simulation") +
      scale_y_continuous(limits = c(input$y_min, input$y_max)) +
      theme_minimal()
  })
}

# Run the application
shinyApp(ui = ui, server = server)
