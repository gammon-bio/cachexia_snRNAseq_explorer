# app.R — Entry point

source("global.R")
source("ui.R")
source("server.R")

shinyApp(ui = ui, server = server)
