#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyWidgets) #za nastavljanje slike ozadja
library(shinythemes)
library(XML)

Sys.setlocale("LC_CTYPE", "Slovenian") # da pravilno izpisuje šumnike

pod_dnevni <-
    "http://www.arso.gov.si/xml/zrak/ones_zrak_dnevni_podatki_zadnji.xml"
pod <- xmlParse(pod_dnevni)
tabela_dnevni <-
    xmlToDataFrame(pod, nodes = getNodeSet(pod, "//postaja"))



shinyUI(
    fluidPage(
        #dodamo sliko iz interneta za ozadje aplikacije
        setBackgroundImage(src = "https://images.unsplash.com/photo-1483825366482-1265f6ea9bc9?ixid=MXwxMjA3fDB8MHxwaG90by1wYWdlfHx8fGVufDB8fHw%3D&ixlib=rb-1.2.1&auto=format&fit=crop&w=1050&q=80"),
        
        
        titlePanel("Onesnaženost zraka v Sloveniji"),
        theme = shinythemes::shinytheme("superhero"),
        
        tabsetPanel(
            tabPanel("Uporaba aplikacije",
                     hr(),
                     textOutput("uporaba"),
                     hr(),
                     style="background-color: #000000"
                     ),
            
            tabPanel("Napoved",
                     sidebarLayout(
                         sidebarPanel(
                             selectInput(
                                 inputId = "napoved_cas",
                                 label = "Izberi napoved",
                                 choices = c("danes", "jutri")
                             )
                         ),
                         
                         mainPanel(tabsetPanel(
                             tabPanel(
                                 title = "Napoved",
                                 hr(),
                                 textOutput("napoved"),
                                 hr(),
                                 uiOutput("napoved_slika")
                             ),
                             
                             tabPanel(title = "Priporočila",
                                      hr(),
                                      textOutput("priporocila"),
                                      hr(),
                                      style="background-color: black")
                            
                                      
                         ))
                     )),
            
            tabPanel("Trenutno stanje",
                     sidebarLayout(
                         sidebarPanel(
                             selectInput(
                                 inputId = "cas",
                                 label = "Izberi obdobje",
                                 choices = c("povprecje vceraj", "povprecje zadnje ure")
                             ),
                             
                             
                             checkboxGroupInput(
                                 inputId = "merilno_mesto",
                                 label = "Izberi merilno mesto",
                                 choices = unique(tabela_dnevni$merilno_mesto)
                             ),
                             tags$style(HTML("
                  .btn {
                    display:block;
                    height: 40px;
                    width: 100px;
                    border: 1px solid red;

                    }

                    ")),
                             actionButton("selectall","Izberi vse/nobenega", style='height: 40px;
  width: 180px'),
                              
                             hr(),
                    
                             actionButton(inputId = "action", label = "Posodobi", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                         ),
                         mainPanel(tabsetPanel(
                             tabPanel(title = "Tabela",
                                      textOutput("datum"),
                                      shiny::dataTableOutput("tabela"),
                                      hr(),
                                      textOutput("textPar"),
                                      hr(),
                                      img(src= 'https://www.celje.info/wp-content/uploads/2020/05/tabela.jpg'),
                                      style = "background-color: black"),
                            
                             
                             tabPanel(
                                 title = "Grafični prikaz",
                                 textOutput("datum2"),
                                 radioButtons(inputId = "parameterGraf", label="Merjeni parameter", choices=c("PM10", "PM2.5", "NO2", "O3"), selected="PM10", inline = TRUE),
                                 plotOutput("grafPar"),
                                 style = "background-color: black"
                             ),
                             
                             
                             tabPanel(
                                 title= "Opis merjenih parametrov",
                                 style="background-color: black",
                                 radioButtons(inputId = "parameter", label="Merjeni parameter", choices=c("PM10", "PM2.5", "NO2", "O3", "SO2", "Indeks onesnazenosti zraka"), selected="PM10", inline = TRUE),
                                 span(textOutput("opisPar"),
                                 hr(),
                                 img(src= 'https://www.celje.info/wp-content/uploads/2020/05/tabela.jpg')
                                 )
                             )
                         ))
                     ))
            
        )
        
    )
)
    
