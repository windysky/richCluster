
homeTabUI <- function(id, tabName) {
  ns <- NS(id)
  tabItem(tabName = tabName,
          h1("richCluster Shiny App"),
          p("Developed by Sarah Hong")
  )
}


homeTabServer <- function(id) {

  moduleServer(id, function(input, output, session) {


  })

}
