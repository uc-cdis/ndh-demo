#if (!require('shiny')) install.packages('shiny', repos='http://cran.us.r-project.org')
#if (!require('googleVis')) install.packages('googleVis', repos='http://cran.us.r-project.org')
require('shiny')
require('googleVis')

# Read data
data <- read.table("/srv/shiny-server/ndh/data.tsv", sep='\t', header=TRUE)
columns <- colnames(data)
#columns <- columns[-7]

# Get values for each field
get_options <- function(filter_group){
    
  options <- sort(unique(data[filter_group])[,filter_group])
  return(options)
}

# Define UI for application
ui <- fluidPage(
  
  # Application title
  #div(img(src="logo.jpg", width=300), align="center"),
  h3('NIAID Data Portal - Exploration', align = "left"),
  
  # Sidebar with a slider input for number of observations
  sidebarPanel(
        width = 2,
        checkboxGroupInput("project", 
                           label = "Project:", 
                           choices = get_options("Project")),
        
        checkboxGroupInput("study", 
                        label = "Study Name:", 
                        choices = get_options("Study")),
        
        checkboxGroupInput("gender", 
                           label = "Gender:", 
                           choices = get_options("Gender")),
                 
        checkboxGroupInput("race", 
                           label = "Race:", 
                           choices = get_options("Race")),

        checkboxGroupInput("ethnicity", 
                           label = "Ethnicity:", 
                           choices = get_options("Ethnicity")),
        
        checkboxGroupInput("vital", 
                           label = "Vital Status:", 
                           choices = get_options("Vital.Status"))         

  ),
  
  # Show output in charts and table
  mainPanel(
     fluidRow(
        column(4, htmlOutput("ChartCase")),
        column(4, htmlOutput("ChartRace")),
        column(4, htmlOutput("ChartVisits"))
     ),
     fluidRow(     
        column(6, htmlOutput("distLab")),
        column(6, htmlOutput("distDrug"))
     ), 
     column(12, dataTableOutput("ATab"))
  )
)

# Define server logic app
server <- function(input, output, session) {

    # Create filtered dataframe
    caseTab <- reactive({
      filterTab <- data
      if(length(input$project) > 0){
        filterTab <- subset(filterTab, Project %in% input$project)
      }        
      if(length(input$study) > 0){
        filterTab <- subset(filterTab, Study %in% input$study)
      }  
      if(length(input$gender) > 0){
        filterTab <- subset(filterTab, Gender %in% input$gender)
      }
      if(length(input$vital) > 0){
        filterTab <- subset(filterTab, Vital.Status %in% input$vital)
      }
      if(length(input$race) > 0){
        filterTab <- subset(filterTab, Race %in% input$race)
      }
      if(length(input$ethnicity) > 0){
        filterTab <- subset(filterTab, Ethnicity %in% input$ethnicity)
      }      
      filterTab
    })
 
    # Render table
    output$ATab <- renderDataTable({caseTab()[, columns]},  
                                   options = list(autoWidth = TRUE,
                                                  striped = TRUE, 
                                                  bordered = TRUE,  
                                                  hover = TRUE,  
                                                  width = '100%',
                                                  digits = 2, 
                                                  na = '--', 
                                                  spacing = 'xs',
                                                  pageLength = 20))  

    # Render pie charts
    output$ChartCase <- renderGvis({
        filterTab <- caseTab()
        if (nrow(filterTab) > 0){
          cases <- as.data.frame(table(unique(filterTab[, c("ID", "Study")])[, "Study"]))
          if(!all(cases$Freq == 0)){
            gvisPieChart(cases, options = list(height = '200px', width = '100%', legend="none", title='Cases', titleTextStyle="{fontSize:16}"))
          }
          else {
            cases <- as.data.frame(table(c('No Data')))
            gvisPieChart(cases, options = list(pieSliceTextStyle="{color: 'black', fontSize:14}", pieSliceText='label', colors="['#D0D0D0']", height = '200px', width = '100%', legend="none", title='Cases', titleTextStyle="{fontSize:16}"))
          }
          
        }
        else{
          cases <- as.data.frame(table(c('No Data')))
          gvisPieChart(cases, options = list(pieSliceTextStyle="{color: 'black', fontSize:14}", pieSliceText='label', colors="['#D0D0D0']", height = '200px', width = '100%', legend="none", title='Cases', titleTextStyle="{fontSize:16}"))
        }  
        
    })
    
    output$ChartRace <- renderGvis({
        filterTab <- caseTab()
        if (nrow(filterTab) > 0){
          race <- as.data.frame(table(filterTab[, "Race"]))
          if(!all(race$Freq == 0)){
            gvisPieChart(race, options = list(height = '200px', width = '100%', legend="none", title='Race', titleTextStyle="{fontSize:16}"))
          }
          else {
            race <- as.data.frame(table(c('No Data')))
            gvisPieChart(race, options = list(pieSliceTextStyle="{color: 'black', fontSize:14}", pieSliceText='label', colors="['#D0D0D0']", height = '200px', width = '100%', legend="none", title='Race', titleTextStyle="{fontSize:16}"))
          }
          
        }
        else{
          race <- as.data.frame(table(c('No Data')))
          gvisPieChart(race, options = list(pieSliceTextStyle="{color: 'black', fontSize:14}", pieSliceText='label', colors="['#D0D0D0']", height = '200px', width = '100%', legend="none", title='Race', titleTextStyle="{fontSize:16}"))
        }                  
        
    })

    output$ChartVisits <- renderGvis({
        filterTab <- caseTab()
        if (nrow(filterTab) > 0){
          visits <- aggregate(Number.Visits ~ Study, filterTab[,c("Number.Visits", "Study")], sum)
          if(!all(visits$Number.Visits == 0)){
            gvisPieChart(visits, options = list(height = '200px', width = '100%', legend="none", title='Number Visits', titleTextStyle="{fontSize:16}"))
          }
          else {
            visits <- as.data.frame(table(c('No Data')))
            gvisPieChart(visits, options = list(pieSliceTextStyle="{color: 'black', fontSize:14}", pieSliceText='label', colors="['#D0D0D0']", height = '200px', width = '100%', legend="none", title='Number Visits', titleTextStyle="{fontSize:16}"))
          }
          
        }
        else{
          visits <- as.data.frame(table(c('No Data')))
          gvisPieChart(visits, options = list(pieSliceTextStyle="{color: 'black', fontSize:14}", pieSliceText='label', colors="['#D0D0D0']", height = '200px', width = '100%', legend="none", title='Number Visits', titleTextStyle="{fontSize:16}"))
        }              
    })
    
    output$distLab <- renderGvis({
      filterTab <- caseTab()
      if (nrow(filterTab) > 0){
          histvisit <- filterTab[,c("Number.Lab.Records")]
          histvisit <- sort(as.numeric(as.character(histvisit[histvisit!="None"])))
          histvisit <- data.frame(histvisit)
          gvisHistogram(histvisit, options=list(
            legend="none",
            colors= "['#75AADB']", 
            border = "white",
            width='500px',
            height='350px',
            autoWidth = TRUE,
            title="# Lab Records per case",
            titleTextStyle="{fontSize:16, textPosition: 'center'}"
          ))
      }
    })    
     
    output$distDrug <- renderGvis({
      filterTab <- caseTab()
      if (nrow(filterTab) > 0){
          histdrugs <- filterTab[,c("Number.Drug.Records")]
          histdrugs <- sort(as.numeric(as.character(histdrugs[histdrugs!="None"])))
          #histdrugs <- histdrugs[histdrugs!=0]
          histdrugs <- data.frame(histdrugs)
          gvisHistogram(histdrugs, options=list(
            legend="none",
            colors= "['#75AADB']", 
            border = "white",
            width='500px',
            height='350px',
            autoWidth = TRUE,
            title="# AIDS Drug Records per case",
            titleTextStyle="{fontSize:16, textPosition: 'center'}"
            ))
       }
    })
}    

## Run app
shinyApp(ui, server)
