# app.R
# TRUENORTH Shiny app entry point

# Source global setup
source("global.R")

# Source UI and server
source("ui.R")
source("server.R")

# Run the app
shinyApp(ui = ui, server = server)
