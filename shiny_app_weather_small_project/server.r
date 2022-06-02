
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(XML)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT) #za tabelo
library(textreadr)
library(rvest)
library(stringr)

Sys.setlocale("LC_CTYPE", "Slovenian") # da pravilno izpisuje šumnike


shinyServer(function(input, output, session){
    
    observe({
        if (input$selectall%%2 != 0)
        {
            updateCheckboxGroupInput(session,inputId = "merilno_mesto",
                                     label = "Izberi merilno mesto",
                                     choices = unique(tabela_dnevni$merilno_mesto)
                                     )
        }
        else
        {
            updateCheckboxGroupInput(session,inputId = "merilno_mesto",
                                     label = "Izberi merilno mesto",
                                     choices = unique(tabela_dnevni$merilno_mesto),
                                     selected = unique(tabela_dnevni$merilno_mesto)
                                     )
        }
    })

    #PRIKAZ PODATKOV DNEVNI IN URNI
    # podatki za dnevno onesnaženost zraka (vceraj)
    pod_dnevni <- "http://www.arso.gov.si/xml/zrak/ones_zrak_dnevni_podatki_zadnji.xml"
    
    #vsake 1000 milisekund preveri če so podatki posodobljeni in jih posodobi če niso + naredi dataframe
    df_dnevni <- reactivePoll(1000, session, 
                         #ta funkcija vrne čas ko je bil pod_dnevni zadnjič spremenjen
                         checkFunc = function() {
                             if (file.exists(pod_dnevni))
                                 file.info(pod_dnevni)$mtime[1]
                             else
                                 ""
                         },
                         #ta funkcija vrne vsebino df_dnevni
                         valueFunc = function() {
                             pod <- xmlParse(pod_dnevni)
                             df <- xmlToDataFrame(pod, nodes=getNodeSet(pod, "//postaja"))
                             df[,c(1:5, 7:8)]
                             #colnames(df) <- c("merilno mesto", "datum", "pm10 dnevno", "pm2.5 dnevno", "o3 urni max","o3 8-urni max", "no2 urni max")
                         }
                         )
    
    
    #podatki za urno onesnazenost zraka (zadnja ura)
    pod_urni <- "http://www.arso.gov.si/xml/zrak/ones_zrak_urni_podatki_zadnji.xml"
    
    #vsake 1000 milisekund preveri če so podatki posodobljeni in jih posodobi če niso + naredi dataframe
    df_urni <- reactivePoll(1000, session,
                            #ta funkcija vrne čas ko je bil pod_urni zadnjič spremenjen
                              checkFunc = function() {
                                  if (file.exists(pod_urni))
                                      file.info(pod_urni)$mtime[1]
                                  else
                                      ""
                              },
                            #ta funkcija vrne vsebino df_urni
                              valueFunc = function() {
                                  pod <- xmlParse(pod_urni)
                                  df <- xmlToDataFrame(pod, nodes=getNodeSet(pod, "//postaja"))
                                  df[,c(1:3, 6, 7, 4, 5, 8)]
                              }
                            )
    
    #tabela podatkov dnevni ali urni
    output$tabela <- shiny::renderDataTable(options = list(searching = FALSE,paging = FALSE,language = list(
        zeroRecords = "Izberite vsaj eno merilno mesto")),{
        
        
        input$action
        isolate({
        
        if(input$cas=="povprecje vceraj") {
            
            df_dnevni() %>%
                filter(merilno_mesto %in% input$merilno_mesto) %>%
                select(-datum) %>% 
                rename("Merilno mesto" = merilno_mesto, "PM10 [ug/m3]" = pm10_dnevna, "PM2.5 [ug/m3]" = pm2.5_dnevna, "O3 [ug/m3]" = o3_max_urna, "NO2 [ug/m3]" = no2_max_urna, "SO2 [ug/m3]" = so2_dnevna)
                }
        
        else{
            df_urni() %>%
                filter(merilno_mesto %in% input$merilno_mesto) %>%
                select(-datum_od, -datum_do) %>%
                rename("Merilno mesto" = merilno_mesto, "PM10 [ug/m3]" = pm10, "PM2.5 [ug/m3]" = pm2.5, "O3 [ug/m3]" = o3, "NO2 [ug/m3]" = no2, "SO2 [ug/m3]" = so2)
            
        }
            
        })
    })
    
    
    
    #NAPOVED
    
    pm10_napoved <- "http://www.arso.gov.si/zrak/kakovost%20zraka/podatki/PM10_napoved.html"
    
    #vsake 1000 milisekund preveri če so podatki posodobljeni in jih posodobi če niso 
    napoved_danes <- reactivePoll(1000, session,
                            #ta funkcija vrne čas ko je bil pm10_napoved zadnjič spremenjen
                            checkFunc = function() {
                                if (file.exists(pm10_napoved))
                                    file.info(pm10_napoved)$mtime[1]
                                else
                                    ""
                            },
                            #ta funkcija vrne vsebino napoved_danes
                            valueFunc = function() {
                                pod <- read_html(pm10_napoved)
                                pod %>% html_nodes("p:nth-child(7)")%>% html_text(trim = TRUE) #s tem izberemo tisti del na spletni strani ki ga želimo
                            }
                          )
    napoved_jutri <- reactivePoll(1000, session,
                                  #ta funkcija vrne čas ko je bil pm10_napoved zadnjič spremenjen
                                  checkFunc = function() {
                                      if (file.exists(pm10_napoved))
                                          file.info(pm10_napoved)$mtime[1]
                                      else
                                          ""
                                  },
                                  #ta funkcija vrne vsebino napoved_jutri
                                  valueFunc = function() {
                                      pod <- read_html(pm10_napoved)
                                      pod %>% html_nodes("p:nth-child(9)")%>% html_text(trim = TRUE) #s tem izberemo tisti del na spletni strani ki ga želimo
                                  }
    )
    
    pm10_napoved_danes <- "http://www.arso.gov.si/zrak/kakovost%20zraka/podatki/pm10_danes.jpg"
    #vsake 1000 milisekund preveri če so podatki posodobljeni in jih posodobi če niso 
    napoved_danes_slika <- reactivePoll(1000, session,
                                  #ta funkcija vrne čas ko je bil pm10_napoved_danes zadnjič spremenjen
                                  checkFunc = function() {
                                      if (file.exists(pm10_napoved_danes))
                                          file.info(pm10_napoved_danes)$mtime[1]
                                      else
                                          ""
                                  },
                                  #ta funkcija vrne vsebino napoved_danes_slika
                                  valueFunc = function() {
                                     pm10_napoved_danes
                                  }
    )
    
    pm10_napoved_jutri <- "http://www.arso.gov.si/zrak/kakovost%20zraka/podatki/pm10_jutri.jpg"
    #vsake 1000 milisekund preveri če so podatki posodobljeni in jih posodobi če niso 
    napoved_jutri_slika <- reactivePoll(1000, session,
                                        #ta funkcija vrne čas ko je bil pm10_napoved_jutri zadnjič spremenjen
                                        checkFunc = function() {
                                            if (file.exists(pm10_napoved_jutri))
                                                file.info(pm10_napoved_jutri)$mtime[1]
                                            else
                                                ""
                                        },
                                        #ta funkcija vrne vsebino napoved_jutri_slika
                                        valueFunc = function() {
                                            pm10_napoved_jutri
                                        }
    )
    
    #napoved danes text
    output$napoved <- renderText({
        if(input$napoved_cas=="danes") {
            napoved_danes() 
        }
        else{
            napoved_jutri() 
                 
        } 
    })
    
    #napoved danes slika
    output$napoved_slika <- renderUI({
        if(input$napoved_cas=="danes") {
            tags$img(src = napoved_danes_slika())
        }
        else{
            tags$img(src = napoved_jutri_slika()) 
            
        } 
    })
    
    # UPORABA APLIKACIJE

    output$uporaba <- renderText({
       "V tej aplikaciji lahko spremljate trenutno stanje onesnaženosti zraka za prejšno uro in prejšnji dan ter napoved onesnaženosti zraka za danes in jutri.
        Pri trenutnem stanju lahko izberete katero merilno mesto v Sloveniji vas zanima in ali želite grafični ali tabelarni prikaz. Pri napovedi pa lahko prebrete še priporočila, ki veljajo ob povišanih vrednostih onesnaženosti zraka."
    })
        
    # PRIPOROČILA
    
    output$priporocila <- renderText({
        pod <- read_html(pm10_napoved)
        pod %>% html_nodes("ul~ p , .vsebina li , p:nth-child(11)")%>% html_text(trim = TRUE)
    })
    
    
    # TEXT pod tabelo
    
    output$textPar <- renderText({"V Sloveniji za meritve kakovosti zunanjega zraka opravlja Agencija Republike Slovenije za okolje (ARSO). V tabeli so prikazani 4 glavni parametri kakovosti zunanjega zraka, in sicer podatki o povprecnih koncentracijah PM10, PM2.5, O3, NO2 in SO2."})
    
    # DATUMI (čas meritev)
    

    pod <- xmlParse(pod_dnevni)
    tabela_dnevni <- xmlToDataFrame(pod, nodes = getNodeSet(pod, "//postaja"))
    

    pod2 <- xmlParse(pod_urni)
    tabela_urni <- xmlToDataFrame(pod2, nodes=getNodeSet(pod2, "//postaja"))
    
    output$datum  <- renderText({
        input$action
        isolate({
        
        if
        (input$cas=="povprecje vceraj"){
            paste("Datum: ", format(as.Date(tabela_dnevni$datum[1], format="%Y-%m-%d"),"%d.%m.%Y"))
        }
        else{
            paste("Casovno obdobje od ", tabela_urni$datum_od[1], "do", tabela_urni$datum_do[1])
        }
        })
    })
    
    
    output$datum2  <- renderText({

            
            if
            (input$cas=="povprecje vceraj"){
                paste("Datum: ", format(as.Date(tabela_dnevni$datum[1], format="%Y-%m-%d"),"%d.%m.%Y"))
            }
            else{
                paste("Casovno obdobje od ", tabela_urni$datum_od[1], "do", tabela_urni$datum_do[1])
            }

    })

    
    # GRAFICNI PRIKAZ PODATKOV
    
    


    
    output$grafPar <- renderPlot({
        
        
        
        if(input$cas=="povprecje vceraj") {
        tabela_dnevni %>% 
            rename(PM10 = pm10_dnevna, "PM2.5" = pm2.5_dnevna, "O3" = o3_max_urna, "NO2" = no2_max_urna) %>%     
            select(merilno_mesto, input$parameterGraf) %>%
            filter(merilno_mesto %in% input$merilno_mesto) %>%
            rename (Parameter = 2) %>% 
            mutate(Parameter = as.numeric(gsub("[^0-9.-]", "", Parameter))) %>%  
            drop_na %>% 
            mutate(Parameter = as.numeric(Parameter))%>% 
            ggplot(aes(merilno_mesto, Parameter)) +
            geom_point(size=2.5, aes(colour = cut(Parameter,                 
                                        if(input$parameterGraf=="PM10"){
                                            c(0,20,40,75,100,Inf)
                                        }
                                        else if(input$parameterGraf=="PM2.5"){
                                            c(0,10,20,40,50,Inf)
                                        }
                                        else if(input$parameterGraf=="NO2"){
                                            c(0, 40,100,200,400,Inf)
                                        }
                                  
                                        else if(input$parameterGraf=="O3"){
                                            c(0,60,120,180,140,Inf)
                                        }
                                        
                                        
                                        ))) +
            xlab("Merilno mesto") + ylab("Koncentracija onesnazila [ug/m3]") +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            scale_color_manual(name = "Kakovost zraka",
                                  values = c("blue", "green", "yellow", "orange", "red"), labels=c("Zelo dobra", "Dobra", "Srednja", "Slaba", "Zelo slaba"), drop=F)+
            theme_light() +
            scale_x_discrete(guide = guide_axis(angle = 90)) 
        }
        
   
        
        else{
            tabela_urni %>% 
                rename("PM10" = pm10, "PM2.5" = pm2.5, "O3" = o3, "NO2" = no2) %>% 
                select(merilno_mesto, input$parameterGraf) %>%
                filter(merilno_mesto %in% input$merilno_mesto) %>%
                rename (Parameter = 2) %>%     
                mutate(Parameter = as.numeric(gsub("[^0-9.-]", "", Parameter))) %>%  
                drop_na %>% 
                ggplot(aes(merilno_mesto, Parameter)) +
                geom_point(aes(colour = cut(Parameter,                 
                                            if(input$parameterGraf=="PM10"){
                                                c(0,20,40,75,100, Inf)
                                            }
                                            else if(input$parameterGraf=="PM2.5"){
                                                c(0,10,20,40,50, Inf)
                                            }
                                            else if(input$parameterGraf=="NO2"){
                                                c(0, 40,100,200,400, Inf)
                                            }
                                            else if(input$parameterGraf=="SO2"){
                                                c(0,100,200,350,500, Inf)
                                            }
                                            else if(input$parameterGraf=="O3"){
                                                c(0,60,120,180,140, Inf)
                                            }
                                            
                                            
                ))) +
                xlab("Merilno mesto") + ylab("Koncentracija onesnazila [ug/m3]") +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                scale_color_manual(name = "Kakovost zraka",
                                   values = c("blue", "green", "yellow", "orange", "red"), labels=c("Zelo dobra", "Dobra", "Srednja", "Slaba", "Zelo slaba"), drop=F)+
                theme_light() +
                scale_x_discrete(guide = guide_axis(angle = 90)) 
                
        }
        
        
    })


    # OPIS MERJENIH PARAMETROV
    output$opisPar <- renderText({

        if(input$parameter %in% c("PM10", "PM2.5")) {
            paste("PM (Particulate Matter) predstavlja zmes trdnih delcev in kapljic organskega ali anorganskega izvora, ki lebdijo v zraku. Delci nastajajo zaradi naravnih procesov (zemlja, soli morja, prah zaradi požarov v naravi, erozija kamenin, vulkanski prah, cvetni prah) ali človekovih aktivnosti (delci iz motorjev z notranjim izgorevanjem, promet, kmetijstvo, gradbišča, proizvodnja cementa, sežigalnice odpadkov, elektrarne, tobačni dim, prah iz kurilnih naprav). Delce PM delimo po velikosti na v 3 skupine: grobi oz. PM10 (premer 10-2.5 um), fin oz. PM2.5 (2.5-0.1 um) in ultra fini oz PM0.1 (<0.1 um). PM10 imajo sposobnost, da se usedejo relativno hitro, medtem ko PM2.5 in PM0.1 ostanejo v zraku precej dlje. Večja onesnaženost s PM delci je navadno v urbanih naseljih. Predvsem pozimi v kotlinah zaradi temperaturne inverzije in povečane uporabe kurilnih naprav koncentracija PM delci večkrat presegajo mejno dnevno koncentracijo. Prašni delci škodljivo vplivajo na dihala, kožo, oči ter prebavni sistem.")
        }

        else if (input$parameter == "NO2") {
            paste("Dušikov dioksid je rdečkastorjav plin, ki nastaja z oksidacijo dušikovega oksida NO s kisikom iz zraka. Koncentracija NO2 v zraku je visoka predvsem v urbanih predelih in je odvisna od vremenskih razmer, prisotnosti ozona ter količine prometa, ki je eden izmed glavnih virov NO2. Posebej visoke koncentracije lahko pričakujemo v hladnih zimskih dneh z malo vetra. Urna mejna koncentracija za varovanje zdrava ljudi je 200 ug/m3. Višje koncentracije NO2 prizadenejo predvsem kronične astmatike in bronhitike in vodijo do pojava kašlja, bronhitisa, oslabitve imunskega sistema in povečanja alergijskih reakcij.")
        }
        
        else if (input$parameter == "O3"){
            paste("Ozon je plin, sestavljen iz treh atomov kisika. Ozon v stratosferi predstavlja ščit pred sončnim UV sevanjem. Povišane koncentracije ozona v troposferi so škodljive za zdravje ljudi. Med človeške vire ozona oziroma njegovih predhodnikov sodijo izpuhi motornih vozil, industrijske emisije ter hlapi goriv in topil. Povišane koncentracije ozona so predvsem poleti.")
        }

        else if (input$parameter == "SO2"){
            paste("žveplov dioksid je najpogostejši izmed žveplovih oksidov. Je brezbarven plin z dražečim vonjem. Povzroča kisel dež. Glavni naravni viri SO2 so izbruhi požarov in gozdni požari. Med človeške vire onesnaženja sodijo izgorevanje goriv (nafta, premog) in industrijski procesi (predelava rud, elektrarne, rafinerije nafte). Urna mejna koncentracija za varovanje zdrava ljudi je 350 ug/m3. Onesnaženje škodi tako okolju kot zdravju ljudi, saj povzroča kisli dež in težave z dihali.")
        }
        
        else if (input$parameter == "Indeks onesnazenosti zraka"){
            paste("Indeks kakovosti zunanjega zraka se izračunava na podlagi ravni onesnaženosti s štirimi onesnaževali: delci PM10, NO2, SO2 in O3. Za vsako onesnaževalo se po določenem algoritmu vsako uro izračuna vrednost indeksa, pri čemer skupni indeks določa onesnaževalo z najvišjo vrednostjo indeksa. Za O3, NO2 in SO2 se pri izračunu upoštevajo zadnje urne ravni onesnaževal, v primeru delcev PM10 pa uteženo 12 urno drseče povprečje. Na podlagi izračunane vrednosti indeksa se stanje kakovosti zraka uvrsti v enega od štirih razredov: dobra, mejna, slaba in zelo slaba kakovost zraka. Z razredi so povezane tudi barve, dobra kakovost zraka se prikazuje z zeleno barvo, mejna z rumeno, slaba z oranžno in zelo slaba z rdečo barvo. Pričakujemo, da bo v zimskem obdobju indeks kakovosti zunanjega zraka določala raven delcev PM10, poleti pa raven ozona. Ker na vseh merilnih mestih ne izvajamo meritev vseh onesnaževal, se kakovost zraka pozimi prikazuje samo za merilna mesta, kjer so na voljo meritve delcev PM10, poleti pa za merilna mesta, kjer potekajo meritve ozona. Vir: ARSO.")
        }

    })

        
    
    
    
    
    
    
    
})
