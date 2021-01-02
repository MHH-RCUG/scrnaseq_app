tabItem(tabName = "help",
        fluidRow(

          box(title = "User Guide", width = 11, solidHeader = T, status = "primary",
              column(12,
                     includeMarkdown("HELP.md")
              )
          )
        )
)
