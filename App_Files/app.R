library(stringr)
library(plotly)
library(tidyverse)

library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(shinyjs)
library(shinyBS)

ui <- dashboardPage(

# Dashboard Header
  dashboardHeader(
    title = "E.PathDash",
    titleWidth = 300
  ),
  skin = "green",
  

# Dashboard Sidebar
  dashboardSidebar(
    width = 300,
    br(),
    div(class="welcome-button",
        actionButton("welcome_button","Welcome!")
    ),
    sidebarMenu(
      
      #Drop-down menu: let user pick study by associated species and strain
      #Also provide button that allows user to reset filters
      menuItem("Study Filters", tabName = "Study Filters",
               
               br(),
               div(style="text-align:center", tags$h6("Filtering by any of the categories below will update")),
               div(style="text-align:center", tags$h6("the list of studies in 'Study Selection'")),
               
               selectizeInput(
                 inputId = "Filter_Species_Type",
                 label = "Filter by species type:",
                 choices = c("Data loading..."),
                 selected="",
                 multiple = FALSE
               ),
               
               br(),
               
               selectizeInput(
                 inputId = "Filter_Strain",
                 label = "Filter by species strain:",
                 choices = c("Select species..."),
                 multiple = TRUE
               ),
               
               br(),
               actionButton("Reset_Species_Filters", "Reset Filters"),
               br(),
               
               radioButtons(
                 inputId = "Select_Study",
                 label = "Select study by GEO Accession ID",
                 choices = c("Data loading..."),
                 selected = character(0)
               ),
               br()
      ), #Study Filters menuItem
      
      menuItem("Additional Information", tabName = "Additional Information",
               
               br(),
               downloadButton(
                 outputId = "Download_Species_Info",
                 label = "Download table of supported species",
                 icon = icon("download"),
                 style = "color: black; margin-left: 15px; margin-bottom: 5px;"
               ),
               br(),
               br(),
               downloadButton(
                 outputId = "Download_Study_Info",
                 label = "Download table of NCBI Studies",
                 icon = icon("download"),
                 style = "color: black; margin-left: 15px; margin-bottom: 5px;"
               ),
               
               br()
      ),# Additional Information
      menuItem("Other DartCF Applications",
               br(),
               uiOutput(outputId = "ESKAPE_Link"),
               br(),
               uiOutput(outputId = "CFSeq_Link"),
               br()
      ), # Other Apps
      menuItem("User Guide",
              br(),
              downloadButton(
                outputId="User_Guide",
                label="Download User Guide",
                icon=icon("download"),
                style = "color: black; margin-left: 15px; margin-bottom: 5px;"
              ),
              br()
      )
    ) #Sidebar menu
 
  ), # dashboardSidebar
  

# Dashboard Body
  dashboardBody(
    useShinyjs(),
    
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style-sheet.css")
    ),
    
    # Page Info Buttons
    div(class="info-buttons",
        actionButton("studypg_info","Study Explorer Info")
    ),
    div(class="info-buttons",
        actionButton("pathpg_info","Pathway Explorer Info")
    ),
    div(class="info-buttons",
        actionButton("termpg_info","Term Explorer Info")
    ),
    div(class="info-buttons",
        actionButton("comppg_info","Study Comparison Info")
    ),
    
    fluidRow(id="Tab_Buttons",
      column(width=12,
        bsButton("Study_View",
                 type="action",
                 label="Study Explorer",
                 size="large",
                 icon=icon("file"),
                 style="success"
                 ),
        bsButton("Path_View",
                 type="action",
                 label="KEGG Pathway Explorer",
                 size="large",
                 style="basic",
                 icon=icon("random")
                 ),
        bsButton("Term_View",
                 type="action",
                 label="GO Term Explorer",
                 size="large",
                 style="basic",
                 icon=icon("random")
                 ),
        bsButton("Comparison_View",
                 type="action",
                 label="Study Comparison",
                 size="large",
                 style="basic",
                 icon=icon("right-left")
                 )
      )
    ), # Tab_Buttons
    

# Study View --------------------------------------------------------------
    fluidRow(id="Summary_Row_1",
      br(),
      br(),
      br(),
      box( width=6,
        title="Study Information",
        id="Study_Info",
        uiOutput(outputId = "Study_Info_Text"),
        br(),
        downloadButton("Raw_Data_Download","Download Count Table and Design Matrix"),
        br(),
      ),
      box(width=6,
        title="Selected Comparison",
        id="Summary_Selection_1",
        selectizeInput(inputId = "Comparison_Selection_1",
                       label= "Select Treatment Comparison",
                       choices=c("Select a comparison")
                       ),
        br(),
        downloadButton("Comp_LogFC_Data","Download Differential Gene Expression Data")
        )
    ), #Summary_Row_1
    
    fluidRow(id ="Summary_Row_2",
      box(collapsible = T, 
        background = "green",
        id="Summary_Boxplot_1",
        title="Significant KEGG Pathways",
        uiOutput(outputId = "Boxplot_1_Text"),
        shinycssloaders::withSpinner(
          plotlyOutput("KEGGBoxPlot"), type=5, color="#154224"
        )
      ),
      box(collapsible = T,
        background="green",
        id="Summary_BoxPlot_2",
        title="Significant GO Terms",
        uiOutput(outputId = "Boxplot_2_Text"),
        shinycssloaders::withSpinner(
          plotlyOutput("GOBoxPlot"), type=5, color="#154224"
        )
      )
    ), # Summary_Row_2
    
    fluidRow(id= "Summary_Row_3",
      box( width=12,collapsible = T,
        height = "100%",
        id="Summary_Table_1",
        title="KEGG Pathways",
        shinycssloaders::withSpinner(
          DT::dataTableOutput("KEGG_Pathways", width="100%"), type=5, color="#154224"
        ),
        br(),
        downloadButton("Table_1","Download Table")
      )
    ),# Summary_Row_3
    
    fluidRow(id="Summary_Row_4",
             box( width=12,collapsible = T,
                  height = "100%",
                  id="Summary_Table_2",
                  title="GO Terms",
                  shinycssloaders::withSpinner(
                    DT::dataTableOutput("GO_Terms",width="100%"), type=5, color="#154224"
                  ),
                  br(),
                  downloadButton("Table_2","Download Table")
             )
    ), #Summary_Row_4

# Pathway View ------------------------------------------------------------
    fluidRow(id="Pathway_Row_1",
             br(),
             br(),
             br(),
             box(width=12,
                 id="Pathway_Selection",
                 title="Selected KEGG Pathway",
                 selectizeInput(inputId = "Pathway_Selection_1",
                                label="Select a pathway of interest",
                                choices=c("Select a pathway"),
                                selected="Select a pathway"
                )
             )
    ), # Pathway_Row_1
    fluidRow(id="Pathway_Row_2",
              box(width=12,collapsible = T,
                  height="100%",
                  id="Pathway_Study_Table",
                  title="Studies Showing Activated/Repressed Pathway",
                  shinycssloaders::withSpinner(
                    DT::dataTableOutput("Pathway_Studies",width = "100%"),type=5, color="#154224"
                    ),
                  br(),
                  downloadButton("Path_Table","Download table")
                  )
    ), # Pathway_Row_2
    fluidRow(id="Pathway_Row_3",
             box(width=12, collapsible = T,
                 height = "100%",
                 id="Pathway_Volcano_Plot",
                 title="Pathway Volcano Plot",
                 uiOutput(outputId = "Pathway_Selection_Text"),
                 shinycssloaders::withSpinner(
                   plotlyOutput("PathVolcanoPlot"),type=5,color="#154224"
                 ),
                 br(),
                 downloadButton("Path_LogFC_Download","Download Data")
              )
    ),# Pathway_Row_3

# Term View ---------------------------------------------------------------
    fluidRow(id="Term_Row_1",
             br(),
             br(),
             br(),
             box(width=12,
                 id="Term_Selection",
                 title="Selected GO Term",
                 selectizeInput(inputId = "Term_Selection_1",
                                label="Select a GO term of interest",
                                choices=c("Select a term"),
                                selected="Select a term"
                                )
             )
    ), # Term_Row_1
    fluidRow(id="Term_Row_2",
             box(width = 12, collapsible = T,
                 height="100%",
                 id="Term_Study_Table",
                 title="Studies Showing Activated/Repressed GO Term",
                 shinycssloaders::withSpinner(
                   DT::dataTableOutput("Term_Studies", width="100%"), type=5, color="#154224"
                 ),
                 br(),
                 downloadButton("Term_Table","Download table")
             )
     ), # Term_Row_2
    fluidRow(id="Term_Row_3",
             box(width=12, collapsible = T,
                 height = "100%",
                 id="Term_Volcano_Plot",
                 title="Term Volcano Plot",
                 uiOutput(outputId = "Term_Selection_Text"),
                 shinycssloaders::withSpinner(
                   plotlyOutput("TermVolcanoPlot"),type=5,color="#154224"
                 ),
                 br(),
                 downloadButton("Term_LogFC_Download","Download Data")
             )
    ), # Term_Row_3

# Comparison View ---------------------------------------------------------
    fluidRow(id="Comparison_Row_1",
             br(),
             br(),
             br(),
             box(width=6,
                 height="100%",
                 id="Study_Selection_1",
                 title="Select Study",
                 selectizeInput(inputId = "Study_1",
                                label = "Select a study for comparison",
                                choices = c("Select a study")
                                )
                 ),
             box(width=6,
                 height = "100%",
                 id="Study_Selection_2",
                 title="Select Study",
                 selectizeInput(inputId = "Study_2",
                                label = "Select a study for comparison",
                                choices = c("Select a study")
                                )
                 )
    ),
    fluidRow(id="Comparison_Row_2",
             br(),
             box(width=12,collapsible = T,
                 height="100%",
                 id="Comparison_Intersect_Pathways",
                 title="Common Significant Processes",
                 shinycssloaders::withSpinner(
                   DT::dataTableOutput("Comparison_Pathways",width = "100%"),type=5, color="#154224"
                 ),
                 br(),
             )
    ),
    fluidRow(id="Comparison_Row_3",
             box(collapsible = T,
                 background="green",
                 id="Comparison_Barchart1",
                 title= "Process Activation/Repression",
                 shinycssloaders::withSpinner(
                   plotlyOutput("Study1_Barchart"), type=5, color="#154224"
                 )
              ),
             box(collapsible = T,
                 background="green",
                 id="Comparison_Barchart2",
                 title= "Process Activation/Repression",
                 shinycssloaders::withSpinner(
                   plotlyOutput("Study2_Barchart"), type=5, color="#154224"
                 )
             )
    ) # comparison view
  ) # dashboardBody
) # dashboardPage

server <- function(input, output, session) {
  
  rv <- shiny::reactiveValues(
    study_pathway_df=NULL,
    study_term_df=NULL,
    comp_pathways_df=NULL
  )
  
  #Hide windows until study is selected
  {
    shinyjs::hide(id = "Summary_Row_2")
    shinyjs::hide(id = "Summary_Row_3")
    shinyjs::hide(id = "Summary_Row_4")
    shinyjs::hide(id = "Pathway_Row_1")
    shinyjs::hide(id = "Pathway_Row_2")
    shinyjs::hide(id = "Pathway_Row_3")
    shinyjs::hide(id = "Term_Row_1")
    shinyjs::hide(id = "Term_Row_2")
    shinyjs::hide(id = "Term_Row_3")
    shinyjs::hide(id = "Select_Study")
    shinyjs::hide(id="Select_Study")
    shinyjs::hide(id = "Comparison_Row_1")
    shinyjs::hide(id = "Comparison_Row_2")
    shinyjs::hide(id = "Comparison_Row_3")
  }

  # Load data
  {
    load("Data.Rdata")
  }
  

# Initialize Variables ----------------------------------------------------
  {
    ##
    Species_List<-c()
    for (species in names(ESKAPE_Results)){
      Species_List<-append(Species_List, unlist(strsplit(species,"_"))[1])
    }
    names(ESKAPE_Results)<-Species_List
    names(ESKAPE_Sig_Lists)<-Species_List
    ##
    
    ##
    Studies_List<-c()
    for (i in 1:length(ESKAPE_Results)){
      Studies_List<-append(Studies_List,names(ESKAPE_Results[[i]]))
    }
    ##
    
    ##
    Strain_List<-c()
    names(Study_Metadata)<-sapply(names(Study_Metadata),function(x){unlist(str_split(x,"_"))[1]},USE.NAMES = F)
    for (species in names(Study_Metadata)){
      strains<-lapply(Study_Metadata[[species]],function(x){
        return(unlist(strsplit(x$Strain,","))[1])
      })
      index_name<-unlist(strsplit(species,"_"))[1]
      Strain_List[[index_name]]<-unique(unlist(strains,use.names=F))
    }
    ##
    
    ##
    KEGGPath_List<-c()
    GOTerm_List<-c()
    for (species in names(ESKAPE_Sig_Lists)){
      paths<-c()
      terms<-c()
      for (study in names(ESKAPE_Sig_Lists[[species]])){
        for (comp in names(ESKAPE_Sig_Lists[[species]][[study]])){
          paths<-append(paths,unique(names(ESKAPE_Sig_Lists[[species]][[study]][[comp]][["KEGG"]])))
          terms<-append(terms,unique(names(ESKAPE_Sig_Lists[[species]][[study]][[comp]][["GO"]])))
        }
      }
      KEGGPath_List[[species]]<-unique(paths)
      GOTerm_List[[species]]<-unique(terms)
    }
    ##
    
    ##
    Study_Metadata_DF<-data.frame(
      "Date"=c("x"),
      "GEO.Accession"=c("x"),
      "Title"=c("x"),
      "Description"=c("x"),
      "Strain"=c("x"),
      "Medium"=c("x"),
      "Treatment"=c("x"),
      "Genotype"=c("x"),
      "Link"=c("x")
    )
    iter<-0
    for (species in names(Study_Metadata)){
      for (study in names(Study_Metadata[[species]])){
        iter<-iter+1
        Study_Metadata_DF[iter,]<-append(Study_Metadata[[species]][[study]][,1:8],c(Study_Metadata[[species]][[study]]$LinkNCBI))
      }
    }
    ##
    
    ##
    All_Sig_Paths<-list()
    for (species in names(ESKAPE_Sig_Lists)){
      for (study in names(ESKAPE_Sig_Lists[[species]])){
        kegg_paths<-c()
        go_terms<-c()
        for (comp in names(ESKAPE_Sig_Lists[[species]][[study]])){
          kegg_paths<-append(kegg_paths,names(ESKAPE_Sig_Lists[[species]][[study]][[comp]][["KEGG"]]))
          go_terms<-append(go_terms,names(ESKAPE_Sig_Lists[[species]][[study]][[comp]][["GO"]]))
        }
        All_Sig_Paths[[species]][[study]][["KEGG"]]<-unique(kegg_paths)
        All_Sig_Paths[[species]][[study]][["GO"]]<-unique(go_terms)
      }
    }
    ##
    
    ##
    pal<-c("green","blue")
    pal<-setNames(pal,c("green","blue"))
    ##
    
    ## KEGG Mapper Code Variables
    URL_query_sep <- "?"
    URL_base <- "https://www.kegg.jp/kegg-bin/show_pathway"
    URL_query_args_chars <- "&multi_query="
    URL_multi_query_sep <- "%0d%0a" # strangely this is carriage return, line feed
    UpColor <- "green"
    DownColor <- "blue"#"%2354F338"
    DefaultColor <- "white"
    No_Color ="&nocolor=1"
    
    KeggMapperURL <- function(Path_name, strain, Results){
      
      PathGeneTable <- KEGG[[strain]][KEGG[[strain]]$KEGG_paths == Path_name,]
      mergedResults <- merge(PathGeneTable, Results, by.x = "Entry", by.y = "Uniprot" )
      KEGGgenesUp <- mergedResults[mergedResults$logFC > 0, "Cross.reference..KEGG."]
      KEGGgenesDown <- mergedResults[mergedResults$logFC <= 0, "Cross.reference..KEGG."]
      
      URL_multi_queryUp <- paste(paste(KEGGgenesUp, UpColor, sep ='+'),
                                 collapse = URL_multi_query_sep)
      URL_multi_queryDown <- paste(paste(KEGGgenesDown, DownColor, sep ='+'),
                                   collapse = URL_multi_query_sep)
      
      strain_dict<-list("PA14"="pau","PAO1"="pae","VPI-5482"="bth","USA300"="saa","Newman"="sae","SK36"="ssa")
      myId = paste0(strain_dict[[strain]], Path_List[[strain]]$ID[which(Path_List[[strain]]$Name == Path_name)])
      
      URL_map_arg <- paste("map", myId, sep = '=')
      paste(URL_base,
            URL_query_sep,
            URL_map_arg,
            URL_query_args_chars,
            paste(URL_multi_queryUp,URL_multi_queryDown, collapse = URL_multi_query_sep ), 
            No_Color, sep = '')
    }
    ##
    
    ## Boxplot Pre-processing 
    KEGGProcessing<-function(species, study, comp){
      
      KEGGSigList<-ESKAPE_Sig_Lists[[species]][[study]][[comp]][["KEGG"]]
      if(length(KEGGSigList) > 0){
        KEGGSigFC <- lapply(KEGGSigList, function(x){
          Study_LogFC_Genes[[study]][[comp]][which(unique(Study_LogFC_Genes[[study]][[comp]]$Uniprot) %in% x), "logFC"]
        })
        
        KEGGSigMedians <- lapply(KEGGSigFC, function(x){median(x[[1]])})
        Cols <- c(1:length(KEGGSigMedians))
        Meds <- unlist(KEGGSigMedians) < 0
        
        Cols[Meds == TRUE] <- "blue"
        Cols[Meds == FALSE] <- "green"
        
        box_colors<-c()
        paths<-c()
        for (i in 1:length(KEGGSigFC)){
          box_colors<-append(box_colors, rep(Cols[i],length(unlist(KEGGSigFC[[i]]))))
          paths<-append(paths,rep(names(KEGGSigFC[i]),length(unlist(KEGGSigFC[[i]]))))
        }
        
        KEGGSigFC_DF<-data.frame("Path"=paths,
                                 "logFC"=unlist(KEGGSigFC,use.names = F),
                                 "color"=box_colors)
        return(KEGGSigFC_DF)
      }# END IF (significant KEGG pathways)
      return(NULL)
    }
    
    GOProcessing<-function(species, study, comp){
      
      GOSigList<-ESKAPE_Sig_Lists[[species]][[study]][[comp]][["GO"]]
      if(length(GOSigList) > 0){
        GOSigFC <- lapply(GOSigList, function(x){
          Study_LogFC_Genes[[study]][[comp]][which(unique(Study_LogFC_Genes[[study]][[comp]]$Uniprot) %in% x), "logFC"]
        })
        
        GOSigMedians <- lapply(GOSigFC, function(x){median(x[[1]])})
        Cols <- c(1:length(GOSigMedians))
        Meds <- unlist(GOSigMedians) < 0
        
        Cols[Meds == TRUE] <- "blue"
        Cols[Meds == FALSE] <- "green"
        
        box_colors<-c()
        paths<-c()
        for (i in 1:length(GOSigFC)){
          box_colors<-append(box_colors, rep(Cols[i],length(unlist(GOSigFC[[i]]))))
          paths<-append(paths,rep(names(GOSigFC[i]),length(unlist(GOSigFC[[i]]))))
        }
        
        GOSigFC_DF<-data.frame("Path"=paths,
                               "logFC"=unlist(GOSigFC,use.names = F),
                               "color"=box_colors)
        return(GOSigFC_DF)
      }# END IF (Significant GO terms)
      return(NULL)
    }
    ##
    
    ## Extract unique pairs of treatment comparisons (w/in study)
    GetUnique<-function(species,study){
      study_comps<-c()
      for (comparison in names(ESKAPE_Results[[species]][[study]])){
        trts <- str_split_1(comparison, " - ")
        reverse<-paste(trts[2],trts[1],sep=" - ")
        if (reverse %in% study_comps){
          next
        }else{
          study_comps<-append(study_comps,comparison)
        }
      }
      return(study_comps)
    }
    ##
    
    ## Update button styles
    UpdateButtons <- function(status){
      
      updateButton(session, "Study_View",
                   label="Study Explorer",
                   size="large",
                   icon=icon("file"),
                   style=if (status[1]) "success" else "basic"
      )
      updateButton(session, "Path_View",
                   label="KEGG Pathway Explorer",
                   size="large",
                   icon=icon("random"),
                   style=if (status[2]) "success" else "basic"
      )
      updateButton(session,"Term_View",
                   label="GO Term Explorer",
                   size="large",
                   icon=icon("random"),
                   style=if (status[3]) "success" else "basic"
      )
      updateButton(session,"Comparison_View",
                   label="Study Comparison",
                   size="large",
                   icon=icon("right-left"),
                   style=if (status[4]) "success" else "basic"
      )
    }
    ##
    
  } # Initialize Variables and Functions
  

# Initial Filter Updates --------------------------------------------------
  observe({
    updateSelectizeInput(session, inputId = "Filter_Species_Type", choices = Species_List,selected="")
  })
  
  # Hide study filter until species is selected
  shinyjs::hide(id="Select_Study")
  
  # Filter species type
  observeEvent(input$Filter_Species_Type, {
    vec_strain<-c()
    vec_study<-c()
    for (species in input$Filter_Species_Type){
      vec_strain<-append(vec_strain,Strain_List[[species]])
      vec_study<-append(vec_study,names(Study_Metadata[[species]]))
    }
    
    if ((input$Filter_Species_Type != "") & (input$Filter_Species_Type != "Data loading...")){
      shinyjs::show(id="Select_Study")
    }
    observe({
      updateSelectizeInput(session, inputId = "Filter_Strain",choices=vec_strain, selected="")
      updateRadioButtons(session, inputId = "Select_Study", choices = vec_study, selected = character(0))
      
      updateSelectizeInput(session, inputId = "Pathway_Selection_1",choices = KEGGPath_List[[input$Filter_Species_Type]],selected="Select a pathway")
      updateSelectizeInput(session, inputId = "Term_Selection_1", choices=GOTerm_List[[input$Filter_Species_Type]],selected="Select a term")
      
      updateSelectizeInput(session, inputId = "Study_1", choices = vec_study, selected = "")
      updateSelectizeInput(session, inputId = "Study_2", choices = vec_study, selected = "")
    })
  }) # Filter species type
  
  # Filter strain
  observeEvent(input$Filter_Strain, {
    vec_study<-c()
    for (species in input$Filter_Species_Type){
      for (study in Study_Metadata[[species]]){
        if (unlist(strsplit(study$Strain,","))[1] %in% input$Filter_Strain){
          vec_study<-append(vec_study,study$GEO.Accession)
        }
      }
    }
    observe({
      updateRadioButtons(session, inputId = "Select_Study", choices = vec_study, selected = character(0))
      updateSelectizeInput(session, inputId = "Study_1", choices = vec_study, selected = "")
      updateSelectizeInput(session, inputId = "Study_2", choices = vec_study, selected = "")
    })
  }) # Filter strain
  
  # Reset species filter
  observeEvent(input$Reset_Species_Filters, {
    
    shinyjs::hide(id = "Select_Study")
    
    observe({
      updateSelectizeInput(session,inputId = "Filter_Species_Type", choices = Species_List,selected = "")
    })
    observe({
      updateSelectizeInput(session, inputId = "Filter_Strain",choices=c("Select species..."),selected="")
    })
    observe({
      updateSelectizeInput(session, inputId = "Comparison_Selection_1",choices = c("Select a comparison"),selected = "Select a comparison")
    })
    observe({
      updateSelectizeInput(session,inputId="Pathway_Selection_1",choices=c("Select a pathway"),selected="Select a pathway")
    })
    observe({
      updateSelectizeInput(session,inputId="Term_Selection_1",choices=c("Select a term"),selected = "Select a term")
    })
    
    {
      shinyjs::show(id = "Summary_Row_2")
      shinyjs::hide(id = "Summary_Row_3")
      shinyjs::hide(id = "Summary_Row_4")
      shinyjs::hide(id="Study_Info")
      shinyjs::hide(id="Pathway_Row_1")
      shinyjs::hide(id="Pathway_Row_2")
      shinyjs::hide(id = "Pathway_Row_3")
      shinyjs::hide(id="Term_Row_1")
      shinyjs::hide(id = "Term_Row_2")
      shinyjs::hide(id = "Term_Row_3")
      shinyjs::hide(id = "Comparison_Row_1")
      shinyjs::hide(id = "Comparison_Row_2")
      shinyjs::hide(id = "Comparison_Row_3")
    }
  }) # Reset species filter
  
  # Download supported species info
  {
    output$Download_Species_Info<- downloadHandler(
      filename="Supported_Species.csv",
      content = function(file){ 
        file.copy("Supported_Species.csv",file)
        },
      contentType = "text/csv"
    )
  }
  
  # Select study 
  observeEvent(input$Select_Study, {

    Duplicates_Removed<-GetUnique(input$Filter_Species_Type,input$Select_Study)
    
    observe({
      updateSelectizeInput(session, inputId = "Comparison_Selection_1", choices = Duplicates_Removed, selected = Duplicates_Removed[1])
    })
    
    {
      shinyjs::show(id = "Summary_Row_2")
      shinyjs::show(id = "Summary_Row_3")
      shinyjs::show(id = "Summary_Row_4")
      shinyjs::show(id="Study_Info")
      shinyjs::hide(id="Pathway_Row_1")
      shinyjs::hide(id="Pathway_Row_2")
      shinyjs::hide(id = "Pathway_Row_3")
      shinyjs::hide(id="Term_Row_1")
      shinyjs::hide(id = "Term_Row_2")
      shinyjs::hide(id = "Term_Row_3")
      shinyjs::hide( id = "Comparison_Row_1")
      shinyjs::hide(id = "Comparison_Row_2")
      shinyjs::hide(id = "Comparison_Row_3")
    }
    
    UpdateButtons(c(1,0,0,0))
    
    link_text<-paste('<a href="',Study_Metadata[[input$Filter_Species_Type[[1]]]][[input$Select_Study[1]]][,'LinkNCBI'],'" target="_blank">',input$Select_Study[1],"</a>",sep="")
    title_text<-paste(link_text,"<br><b>",Study_Metadata[[input$Filter_Species_Type[[1]]]][[input$Select_Study[1]]][,'Title'],"</b><br>",sep="")
    all_text<-paste(title_text,Study_Metadata[[input$Filter_Species_Type[[1]]]][[input$Select_Study[1]]][,'Description'],sep=" ")
    output$Study_Info_Text<-renderUI(HTML(all_text))
    
  }) # Select study 
  
  {
    output$Download_Study_Info<-downloadHandler(
      filename="Studies.csv",
      content = function(file){
        write.csv(Study_Metadata_DF,file,row.names = F)
      },
      contentType = "text/csv"
    )
  }
  

# Links -------------------------------------------------------------------
  {
    eskape_url<-"<a href=http://scangeo.dartmouth.edu/ESKAPE/ class='link-text' target='_blank'>ESKAPE Act PLUS</a>"
    output$ESKAPE_Link<-renderUI(HTML(eskape_url))
    cfseq_url<-"<a href=http://scangeo.dartmouth.edu/CFSeq/ class='link-text' target='_blank'>CF-Seq</a>"
    output$CFSeq_Link<-renderUI(HTML(cfseq_url))
  }

# Modals -------------------------------------------------------------
  {
    observeEvent(input$welcome_button,{
      showModal(modalDialog(
        title="Welcome to E.PathDash",
        "E.PathDash, ESKAPE Pathogen Pathway Dashboard, enables researchers to analyze  cellular pathways using high-throughput gene expression data for ESKAPE pathogens from the Gene Expression Omnibus. The data has been processed for cellular pathway analysis of KEGG pathways and GO terms, common annotation databases.",
        br(),
        br(),
        'This tool facilitates comparison of pathway activation across experimental conditions and studies, as well as targeted data exploring using a pathway of interest. For more information on each dashboard page, please refer to the "Study Explorer Info", "Pathway Explorer Info", "Term Explorer Info" and "Study Comparison Info" buttons.',
        br(),
        br(),
        strong('Start your search by selecting a bacterial species and strain(s) of interest from "Study Filters".'),
        br(),
        easyClose = TRUE
      ))
    })
    observeEvent(input$studypg_info, {
      showModal(modalDialog(
        title = "Study Explorer Page Use",
        strong("How to Use"),
        br(),
        "Select a bacterial species, strain and GEO study in the filter sidebar. The page will populate with metadata and KEGG pathway/GO Term activation/repression information corresponding to the GEO study.",
        br(),
        br(),
        strong("Page Components"),
        br(),
        em("1. Study Information:"),
        "Link to GEO study page, title and short study description",
        br(),
        br(),
        em("2. Selected Comparison:"),
        "Dropdown to select the study comparisons for which to view pathway/term analysis information. The pathways and terms are activated or repressed in the first sample group compared to the second sample group.",
        br(),
        br(),
        em("3. Significant KEGG Pathways:"),
        "A boxplot showing the distribution of the logFC values for genes in that significantly activated or repressed KEGG pathway. Significance is determined by a FDR value of <0.05 after conducting pathway analysis.",
        br(),
        br(),
        em("4. Significant GO Terms:"),
        "A boxplot showing the distribution of the logFC values for genes in that significantly activated or repressed GO term. Significance is determined by a FDR value of <0.05 after conducting pathway analysis.",
        br(),
        br(),
        em("5. KEGG Pathways:"),
        "A table of all KEGG pathways analyzed with corresponding binomial test statistic, median logFC value for genes in the pathway, P-value and FDR value for significance level of activation/repression. Pathway names link to functional pathway maps, with differentially expressed genes colored according to over (green) or under (blue) expression.",
        br(),
        br(),
        em("6. GO Terms:"),
        "A table of all GO terms analyzed with corresponding binomial test statistic, median logFC value for genes in the pathway, P-value and FDR value for significance level of activation/repression.",
        easyClose = TRUE
      ))
    })
    observeEvent(input$pathpg_info, {
      showModal(modalDialog(
        title = "KEGG Pathway Explorer Page Use",
        strong("How to Use"),
        br(),
        "Explore the compendium of study data by a KEGG pathway of interest. Select a bacterial species for which to view KEGG pathways. After selecting a particular pathway, view studies and treatment comparisons in which that pathway was significantly activated or repressed.",
        br(),
        br(),
        strong("Page Components"),
        br(),
        em("1. Pathway selection:"),
        "Dropdown of KEGG pathways for bacterial species specified in filter.",
        br(),
        br(),
        em("2. Study Table:"),
        "List of studies for which the selected KEGG pathway was significantly activated or repressed. Information includes corresponding treatment comparison within the study for which the pathway was activated/repressed and median logFC of genes within the pathway. Study ID links to GEO page describing study.",
        br(),
        br(),
        em("3. Pathway Volcano Plot:"),
        "Selecting a row from the study table shows a volcano plot for the genes within that pathway. Hovering over a point shows the Uniprot ID of the gene along with its logFC between comparison groups and transformed P-value.",
        easyClose = TRUE
      ))
    })
    observeEvent(input$termpg_info, {
      showModal(modalDialog(
        title="GO Term Explorer Page Use",
        strong("How to Use"),
        br(),
        "Explore the compendium of study data by a GO term of interest. Select a bacterial species for which to view GO terms. After selecting a particular term, view studies and treatment comparisons in which that term was significantly activated or repressed.",
        br(),
        br(),
        strong("Page Components"),
        br(),
        em("1. Term selection:"),
        "Dropdown of GO terms for bacterial species specified in filter.",
        br(),
        br(),
        em("2. Study Table:"),
        "List of studies for which the selected GO term was significantly activated or repressed. Information includes corresponding treatment comparison within the study for which the term was activated/repressed and median logFC of genes within the term. Study ID links to GEO page describing study.",
        br(),
        br(),
        em("3. Term Volcano Plot:"),
        "Selecting a row from the study table shows a volcano plot for the genes within that term. Hovering over a point shows the Uniprot ID of the gene along with its logFC between comparison groups and transformed P-value.",
        easyClose = TRUE
      ))
    })
    observeEvent(input$comppg_info, {
      showModal(modalDialog(
        title="Study Comparison Page Use",
        strong("How to Use"),
        br(),
        "Explore common activated or repressed KEGG pathways and GO terms between studies. After filtering GEO studies by bacterial species and strain, you can select two studies for which to compare pathway analyses.",
        br(),
        br(),
        strong("Page Components"),
        br(),
        em("1. Select Studies:"),
        "Dropdowns of GEO studies to compare.",
        br(),
        br(),
        em("2. Pathway/Term Table:"),
        "List of KEGG pathways and GO terms that are significantly activated or repressed in both selected studies.",
        br(),
        br(),
        em("3. Process Activation/Repression Bar Plots:"),
        "Selecting a pathway/term from the table populates bar plots representing the pathway or term’s activation/repression level across the comparisons in each study. Activation/repression level is symbolized by the median logFC value for genes in the pathway/term under the comparison. Hovering over the bars shows the corresponding sample comparison.",
        easyClose = TRUE
      ))
    })
  }
  

# User Guide --------------------------------------------------------------
  {
      output$User_Guide<-downloadHandler(
        filename="EPathDashUserGuide.pdf",
        content = function(file){ 
          file.copy("EPathDashUserGuide.pdf",file)
        },
        contentType="application/pdf"
      )
  }
# Study View --------------------------------------------------------------
  
  # Download Raw Data
  {
    output$Raw_Data_Download<-downloadHandler(
      filename = function(){
        paste(input$Select_Study,"Raw_Data.zip",sep="_")
      },
      content = function(file){
        count_file<-"Count_Table.csv"
        design_file<-"Design_Matrix.csv"
        write.csv(Raw_Count_Data[[input$Select_Study]][["Count_Table"]],count_file,row.names = F)
        write.csv(Raw_Count_Data[[input$Select_Study]][["Design_Matrix"]],design_file,row.names = F)
        zip(file,files=c(count_file,design_file))
      },
      contentType = "application/zip"
    )
  }

  observeEvent(input$Comparison_Selection_1,{
    
    if(input$Comparison_Selection_1 != "Select a comparison"){

      species<-input$Filter_Species_Type[1]
      study<-input$Select_Study[1]
      comp<-input$Comparison_Selection_1[1]
      strain<-unlist(strsplit(Study_Metadata[[species]][[study]][["Strain"]],","))[1]
      
      KEGGSigFC_DF<-KEGGProcessing(species,study,comp)
      
      Study_Results<-ESKAPE_Results[[species]][[study]]
      KEGG_Results<-Study_Results[[comp]][["KEGG"]]
      
      # KEGG box plot
      if (!is.null(KEGGSigFC_DF)){
        output$Boxplot_1_Text<-renderUI(HTML(""))
        output$KEGGBoxPlot<- renderPlotly({
          plot_ly(data=KEGGSigFC_DF,
                  y=~Path,
                  x=~logFC,
                  color = ~color,
                  colors=pal,
                  type = "box",
                  showlegend=F
          ) %>% layout(
            yaxis=list(tickangle=-45),
            plot_bgcolor="#D8E0E5",
            paper_bgcolor="#D8E0E5"
          ) %>% style(
            hoverinfo = "none"
          )
        })
      }else{
        output$Boxplot_1_Text<-renderUI(HTML("<b>No significant KEGG pathways</b>"))
        output$KEGGBoxPlot <- renderPlotly({})
      }
      
      # KEGG Table
      colnames(KEGG_Results)<-c("Path", "Binomial Test Est.", "Median gene logFC", "P-value", "FDR")
      
      myLinks <- lapply(KEGG_Results$Path, function(x){KeggMapperURL(x,strain,Study_LogFC_Genes[[study]][[comp]])})
      Links_df <- data.frame("Link" = unlist(myLinks),
                             "Name" = KEGG_Results$Path)
      KEGG_Results$Path <- apply(Links_df, 1, function(x){
        paste('<a href="', x[1],'" target="_blank">', 
              x[2], '</a>', collapse="</br>")
      })
      
      output$KEGG_Pathways <- DT::renderDataTable(
        KEGG_Results,
        rownames=F,
        selection='none',
        escape=F
      )
      
      # Download KEGG Table
      output$Table_1 <- downloadHandler(
        filename = function(){
          paste(input$Select_Study, comp, "KEGG.csv",sep="_")
        },
        content = function(file){
          write.csv(Study_Results[[comp]][["KEGG"]], file, row.names = F)
        },
        contentType = "text/csv"
      )

      # Processing GO
      GO_Results<-Study_Results[[comp]][["GO"]]
      GOSigFC_DF<-GOProcessing(species, study, comp)

      # GO Boxplot 
      if (!is.null(GOSigFC_DF)){
        output$Boxplot_2_Text<-renderUI(HTML(""))
        output$GOBoxPlot<- renderPlotly({
          plot_ly(data=GOSigFC_DF,
                  y=~Path,
                  x=~logFC,
                  color = ~color,
                  colors=pal,
                  type = "box",
                  showlegend=F
          ) %>% layout(
            yaxis=list(tickangle=-45),
            plot_bgcolor="#D8E0E5",
            paper_bgcolor="#D8E0E5"
          ) %>% style(
            hoverinfo="none"
          )
        })
      }else{
        output$Boxplot_2_Text<-renderUI(HTML("<b>No significant GO terms</b>"))
        output$GOBoxPlot <- renderPlotly({})
      } 
      
      # GO Table
      colnames(GO_Results)<-c("GO Term", "Binomial Test Est.", "Median gene logFC", "P-value", "FDR")
      output$GO_Terms <- DT::renderDataTable(
        GO_Results,
        rownames=F
      )

      # Download GO Table
      output$Table_2 <- downloadHandler(
        filename = function(){
          paste(input$Select_Study, input$Comparison_Selection_1, "GO.csv",sep="_")
        },
        content = function(file){
          write.csv(Study_Results[[comp]][["GO"]], file, row.names = F)
        },
        contentType = "text/csv"
      )
      
      
      # Download LogFC data
      output$Comp_LogFC_Data <- downloadHandler(
        filename=function(){
          paste(study, comp, "LogFC.csv",sep="_")
        },
        content=function(file){
          write.csv(Study_LogFC_Genes[[study]][[comp]],file,row.names = F)
        },
        contentType="text/csv"
      )
    } # END OUTERMOST IF (A comparison was selected)
    else{ # Reset Species Button Actions
      
      output$KEGGBoxPlot<-NULL
      output$GOBoxPlot<-NULL
      output$KEGG_Pathways<-NULL
      output$GO_Terms<-NULL
      output$Study_Info_Text<-NULL
      
      shinyjs::hide(id="Study_Info")
      shinyjs::hide(id = "Summary_Row_2")
      shinyjs::hide(id = "Summary_Row_3")
      shinyjs::hide(id = "Summary_Row_4")
    }
  }) # Populate single study window
  
  # User navigates via Study View button
  {
    observeEvent(input$Study_View, {
      # Update visible panes
      {
        shinyjs::show(id="Study_Info")
        shinyjs::show(id = "Summary_Row_1")
        shinyjs::show(id = "Summary_Row_2")
        shinyjs::show(id = "Summary_Row_3")
        shinyjs::show(id = "Summary_Row_4")
        shinyjs::hide(id = "Pathway_Row_1")
        shinyjs::hide(id = "Pathway_Row_2")
        shinyjs::hide(id = "Pathway_Row_3")
        shinyjs::hide(id = "Term_Row_1")
        shinyjs::hide(id = "Term_Row_2")
        shinyjs::hide(id = "Term_Row_3")
        shinyjs::hide( id = "Comparison_Row_1")
        shinyjs::hide(id = "Comparison_Row_2")
        shinyjs::hide(id = "Comparison_Row_3")
      }
      UpdateButtons(c(1,0,0,0))
    }) # observeEvent: Study_View
  }
  

# KEGG Pathway View -------------------------------------------------------
  {
    observeEvent(input$Path_View, {
      {
        shinyjs::hide(id="Study_Info")
        shinyjs::hide(id = "Summary_Row_1")
        shinyjs::hide(id = "Summary_Row_2")
        shinyjs::hide(id = "Summary_Row_3")
        shinyjs::hide(id = "Summary_Row_4")
        shinyjs::hide(id = "Term_Row_1")
        shinyjs::hide(id = "Term_Row_2")
        shinyjs::hide(id = "Term_Row_3")
        shinyjs::hide( id = "Comparison_Row_1")
        shinyjs::hide(id = "Comparison_Row_2")
        shinyjs::hide(id = "Comparison_Row_3")
        shinyjs::show(id = "Pathway_Row_1")
        shinyjs::show(id = "Pathway_Row_2")
        shinyjs::show(id = "Pathway_Row_3")
      }
      
      UpdateButtons(c(0,1,0,0))
    }) # observeEvent Path_View
    
    observeEvent(input$Pathway_Selection_1, {
      
      if ((input$Pathway_Selection_1 != "Select a pathway") & (input$Pathway_Selection_1 != "")){
        
        shinyjs::show(id = "Pathway_Row_2")
        output$PathVolcanoPlot<-NULL
        output$Pathway_Selection_Text<-NULL
        
        pathway<-input$Pathway_Selection_1
        Species_Sig_List<-ESKAPE_Sig_Lists[[input$Filter_Species_Type]]
        Pathway_Results<-list()
        for (study in names(Species_Sig_List)){
          sig_comps<-c()
          for (comp in names(Species_Sig_List[[study]])){
            if (pathway %in% names(Species_Sig_List[[study]][[comp]][["KEGG"]])){
              sig_comps<-append(sig_comps, comp)
            }
          }
          if (length(sig_comps)>0){
            Pathway_Results[[study]]<-sig_comps
          }
        }
        
        Pathway_Results_DF<-data.frame(
          "Study"=c('x'),
          "Treatment Comparison"=c('x'),
          "Median logFC"=c('x')
        )
        Pathway_Results_Download<-Pathway_Results_DF
        duplicates<-c()
        for (study in names(Pathway_Results)){
          for (comp in unlist(Pathway_Results[[study]])){
            x<-ESKAPE_Results[[input$Filter_Species_Type]][[study]][[comp]][["KEGG"]][ESKAPE_Results[[input$Filter_Species_Type]][[study]][[comp]][["KEGG"]]$Path==pathway,]$Median_fold_change
            link_text<-paste('<a href="',Study_Metadata[[input$Filter_Species_Type[[1]]]][[study]][,'LinkNCBI'],'" target="_blank">',study,"</a>",sep="")
            trts <- str_split_1(comp," - ")
            reverse<-paste(trts[2],trts[1],sep=" - ")
            if (reverse %in% duplicates){
              next
            }else{
              Pathway_Results_DF[nrow(Pathway_Results_DF)+1,]<-c(link_text,comp,x)
              Pathway_Results_Download[nrow(Pathway_Results_Download)+1,]<-c(study,comp,x)
              duplicates<-append(duplicates,comp)
            }
          }
        }
        output$Pathway_Studies <- DT::renderDataTable(
          Pathway_Results_DF[2:nrow(Pathway_Results_DF),],
          rownames=F,
          escape=F,
          selection="single"
        )
        
        # Download table
        output$Path_Table <-downloadHandler(
          filename=function(){
            paste(pathway,"KEGG.csv",sep="_")
          },
          content=function(file){
            write.csv(Pathway_Results_Download[2:nrow(Pathway_Results_Download),],file,row.names = F)
          },
          contentType = "text/csv"
        )
        
        rv$study_pathway_df <- Pathway_Results_DF[2:nrow(Pathway_Results_DF),]
      
      } #END IF
      else{ # Reset Button hit
        output$Pathway_Studies<-NULL
      }

    }) # observeEvent Pathway_Selection_1
    
    observeEvent(input$Pathway_Studies_rows_selected, {
      
      study <- str_extract(rv$study_pathway_df[input$Pathway_Studies_rows_selected,]$Study,"GSE[0-9]+")
      comp <- rv$study_pathway_df[input$Pathway_Studies_rows_selected,2]
      strain <- str_split_1(Study_Metadata[[input$Filter_Species_Type]][[study]]$Strain,",")[1]
      path_genes<-KEGG[[strain]][KEGG[[strain]]$KEGG_paths == input$Pathway_Selection_1,]$Entry
      logfc_data<-Study_LogFC_Genes[[study]][[comp]][Study_LogFC_Genes[[study]][[comp]]$Uniprot %in% path_genes,]
      output$PathVolcanoPlot <- renderPlotly({
        plot_ly(
          data=logfc_data,
          x=~logFC,
          y=~-log10(PValue),
          type="scatter",
          mode="markers",
          name="Genes in\nKEGG Pathway",
          text=~paste("Uniprot ID:",Uniprot,"<br>-log10(Pvalue):",-log10(PValue),"<br>log2FC:",logFC),
          hoverinfo="text",
          marker=list(color="blue")
        ) %>% 
          layout(
            plot_bgcolor="#D8E0E5",
            paper_bgcolor="#D8E0E5",
            xaxis=list(title="Log2 Fold Change"),
            showlegend=TRUE
          )
      })
      output$Pathway_Selection_Text<-renderUI(paste(input$Pathway_Selection_1,study,comp,sep=" // "))
      
      output$Path_LogFC_Download<-downloadHandler(
        filename = function(){
          paste(input$Pathway_Selection_1,study,comp,"gene_logfc.csv",sep="_")
        },
        content = function(file){
          write.csv(logfc_data,file,row.names = F)
        },
        contentType = "text/csv"
      )
    })
  } #KEGG Pathway View
  

# GO Term View ------------------------------------------------------------

  {
    observeEvent(input$Term_View, {
      {
        shinyjs::hide(id="Study_Info")
        shinyjs::hide(id = "Summary_Row_1")
        shinyjs::hide(id = "Summary_Row_2")
        shinyjs::hide(id = "Summary_Row_3")
        shinyjs::hide(id = "Summary_Row_4")
        shinyjs::hide(id = "Pathway_Row_1")
        shinyjs::hide(id = "Pathway_Row_2")
        shinyjs::hide(id = "Pathway_Row_3")
        shinyjs::hide( id = "Comparison_Row_1")
        shinyjs::hide(id = "Comparison_Row_2")
        shinyjs::hide(id = "Comparison_Row_3")
        shinyjs::show(id = "Term_Row_1")
        shinyjs::show(id = "Term_Row_2")
        shinyjs::show(id = "Term_Row_3")
      }
      
      UpdateButtons(c(0,0,1,0))
    })
    
    observeEvent(input$Term_Selection_1, {
      if ((input$Term_Selection_1 != "Select a term") & (input$Term_Selection_1 != "")){
        
        shinyjs::show(id = "Term_Row_2")
        output$TermVolcanoPlot<-NULL
        output$Term_Selection_Text<-NULL
        
        term<-input$Term_Selection_1
        Species_Sig_List<-ESKAPE_Sig_Lists[[input$Filter_Species_Type]]
        Term_Results<-list()
        for (study in names(Species_Sig_List)){
          sig_comps<-c()
          for (comp in names(Species_Sig_List[[study]])){
            if (term %in% names(Species_Sig_List[[study]][[comp]][["GO"]])){
              sig_comps<-append(sig_comps, comp)
            }
          }
          if (length(sig_comps)>0){
            Term_Results[[study]]<-sig_comps
          }
        }
        
        Term_Results_DF<-data.frame(
          "Study"=c('x'),
          "Treatment Comparison"=c('x'),
          "Median logFC"=c('x')
        )
        Term_Results_Download<-Term_Results_DF
        duplicates<-c()
        for (study in names(Term_Results)){
          for (comp in unlist(Term_Results[[study]])){
            x<-ESKAPE_Results[[input$Filter_Species_Type]][[study]][[comp]][["GO"]][ESKAPE_Results[[input$Filter_Species_Type]][[study]][[comp]][["GO"]]$GO_term==term,]$Median_fold_change
            link_text<-paste('<a href="',Study_Metadata[[input$Filter_Species_Type[[1]]]][[study]][,'LinkNCBI'],'" target="_blank">',study,"</a>",sep="")
            trts <- str_split_1(comp," - ")
            reverse<-paste(trts[2],trts[1],sep=" - ")
            if (reverse %in% duplicates){
              next
            }else{
              Term_Results_DF[nrow(Term_Results_DF)+1,]<-c(link_text,comp,x)
              Term_Results_Download[nrow(Term_Results_Download)+1,]<-c(study,comp,x)
              duplicates<-append(duplicates,comp)
            }
          }
        }
        
        rv$study_term_df<-Term_Results_DF[2:nrow(Term_Results_DF),]
        
        output$Term_Studies <- DT::renderDataTable(
          Term_Results_DF[2:nrow(Term_Results_DF),],
          rownames=F,
          escape=F,
          selection="single"
        )
        # Download table
        {
          output$Term_Table<-downloadHandler(
            filename = function(){
              paste(term,"GO.csv",sep="_")
            },
            content = function(file){
              write.csv(Term_Results_Download[2:nrow(Term_Results_Download),],file,row.names = F)
            },
            contentType = "text/csv"
          )
        }
        
      } #END IF
      else{ # Reset Button hit
        output$Term_Studies<-NULL
      }
    })#observeEvent: Term_Selection_1
    
    observeEvent(input$Term_Studies_rows_selected, {
      
      study <- str_extract(rv$study_term_df[input$Term_Studies_rows_selected,]$Study,"GSE[0-9]+")
      comp <- rv$study_term_df[input$Term_Studies_rows_selected,2]
      strain <- str_split_1(Study_Metadata[[input$Filter_Species_Type]][[study]]$Strain,",")[1]
      term_genes<-GO[[strain]][GO[[strain]]$Gene.ontology..GO. == input$Term_Selection_1,]$Entry
      logfc_data<-Study_LogFC_Genes[[study]][[comp]][Study_LogFC_Genes[[study]][[comp]]$Uniprot %in% term_genes,]
      
      output$TermVolcanoPlot <- renderPlotly({
        plot_ly(
          data=logfc_data,
          x=~logFC,
          y=~-log10(PValue),
          type="scatter",
          mode="markers",
          name="Genes in\nGO Term",
          text=~paste("Uniprot ID:",Uniprot,"<br>-log10(Pvalue):",-log10(PValue),"<br>logFC:",logFC),
          hoverinfo="text",
          marker=list(color="blue")
        ) %>% 
          layout(
            plot_bgcolor="#D8E0E5",
            paper_bgcolor="#D8E0E5",
            xaxis=list(title="Log2 Fold Change"),
            showlegend=TRUE
          )
      })
      output$Term_Selection_Text<-renderUI(paste(input$Term_Selection_1,study,comp,sep=" // "))
      
      output$Term_LogFC_Download<-downloadHandler(
        filename = function(){
          paste(input$Term_Selection_1,study,comp,"gene_logfc.csv",sep="_")
        },
        content = function(file){
          write.csv(logfc_data,file,row.names = F)
        },
        contentType = "text/csv"
      )
    }) #observeEvent Term Studies rows selected
  } 
  

# Comparison View ---------------------------------------------------------

  {
    observeEvent(input$Comparison_View, {
      {
        shinyjs::hide(id="Study_Info")
        shinyjs::hide(id = "Summary_Row_1")
        shinyjs::hide(id = "Summary_Row_2")
        shinyjs::hide(id = "Summary_Row_3")
        shinyjs::hide(id = "Summary_Row_4")
        shinyjs::hide(id = "Pathway_Row_1")
        shinyjs::hide(id = "Pathway_Row_2")
        shinyjs::hide(id = "Pathway_Row_3")
        shinyjs::hide(id = "Term_Row_1")
        shinyjs::hide(id = "Term_Row_2")
        shinyjs::hide(id = "Term_Row_3")
        shinyjs::show( id = "Comparison_Row_1")
        shinyjs::show(id = "Comparison_Row_2")
        shinyjs::show(id = "Comparison_Row_3")
      }
    
      UpdateButtons(c(0,0,0,1))
    })
    
    observeEvent(input$Study_1, {
      species<-input$Filter_Species_Type
      study1<-input$Study_1
      study2<-input$Study_2
      if (study1 != ""){
        if (study2 == ""){
          Intersection_Paths_DF<-data.frame("Process"=All_Sig_Paths[[species]][[study1]][["KEGG"]])
          Intersection_Terms_DF<-data.frame("Process"=All_Sig_Paths[[species]][[study1]][["GO"]])
        }else{
          Intersection_Paths_DF<-data.frame("Process"=intersect(All_Sig_Paths[[species]][[study1]][["KEGG"]],All_Sig_Paths[[species]][[study2]][["KEGG"]]))
          Intersection_Terms_DF<-data.frame("Process"=intersect(All_Sig_Paths[[species]][[study1]][["GO"]],All_Sig_Paths[[species]][[study2]][["GO"]]))
        }
        annotation<-append(rep("KEGG",nrow(Intersection_Paths_DF)),rep("GO",nrow(Intersection_Terms_DF)))
        Intersection_DF<-bind_rows(Intersection_Paths_DF,Intersection_Terms_DF) %>% 
          mutate("Annotation"=annotation)
        
        output$Comparison_Pathways <- DT::renderDataTable(
          Intersection_DF,
          rownames=F,
          escape=F,
          selection="single"
        )
        
        rv$comp_pathways_df <- Intersection_DF
        output$Study2_Barchart<-NULL
        output$Study1_Barchart<-NULL
      }
    }) #observeEvent Study_1
    
    observeEvent(input$Study_2,{
      species<-input$Filter_Species_Type
      study1<-input$Study_1
      study2<-input$Study_2
      if (study2 != ""){
        if (study1 == ""){
          Intersection_Paths_DF<-data.frame("Process"=All_Sig_Paths[[species]][[study2]][["KEGG"]])
          Intersection_Terms_DF<-data.frame("Process"=All_Sig_Paths[[species]][[study2]][["GO"]])
        }else{
          Intersection_Paths_DF<-data.frame("Process"=intersect(All_Sig_Paths[[species]][[study1]][["KEGG"]],All_Sig_Paths[[species]][[study2]][["KEGG"]]))
          Intersection_Terms_DF<-data.frame("Process"=intersect(All_Sig_Paths[[species]][[study1]][["GO"]],All_Sig_Paths[[species]][[study2]][["GO"]]))
        }
        annotation<-append(rep("KEGG",nrow(Intersection_Paths_DF)),rep("GO",nrow(Intersection_Terms_DF)))
        Intersection_DF<-bind_rows(Intersection_Paths_DF,Intersection_Terms_DF) %>% 
          mutate("Annotation"=annotation)
        output$Comparison_Pathways <- DT::renderDataTable(
          Intersection_DF,
          rownames=F,
          escape=F,
          selection="single"
        )
        
        rv$comp_pathways_df <- Intersection_DF
        output$Study2_Barchart<-NULL
        output$Study1_Barchart<-NULL
      }
    }) #observeEvent Study_2
    
    observeEvent(input$Comparison_Pathways_rows_selected, {
      
      species<-input$Filter_Species_Type
      study1<-input$Study_1
      study2<-input$Study_2
      study_1_comps<-GetUnique(input$Filter_Species_Type,study1)
      study_2_comps<-GetUnique(input$Filter_Species_Type,study2)

      Intersection_DF <- rv$comp_pathways_df
      path <- Intersection_DF[input$Comparison_Pathways_rows_selected,]$Process
      annotation <- Intersection_DF[input$Comparison_Pathways_rows_selected,]$Annotation
      col_name<- if(annotation=="KEGG") "Path" else "GO_term"
      fc1 <- c()
      fc2 <- c()
      for (comp in study_1_comps){
        fc1<-append(fc1,ESKAPE_Results[[species]][[study1]][[comp]][[annotation]][ESKAPE_Results[[species]][[study1]][[comp]][[annotation]][,col_name]==path,]$Median_fold_change)
      }
      for (comp in study_2_comps){
        fc2<-append(fc2,ESKAPE_Results[[species]][[study2]][[comp]][[annotation]][ESKAPE_Results[[species]][[study2]][[comp]][[annotation]][,col_name]==path,]$Median_fold_change)
      }
      
      bar_chart_data1<-data.frame("comparison"=study_1_comps,"Median_logFC"=fc1)
      colors<-rep(NA,nrow(bar_chart_data1))
      colors[bar_chart_data1$Median_logFC>0]<-"limegreen"
      colors[bar_chart_data1$Median_logFC<0]<-"blue"
      
      bar_chart_data2<-data.frame("comparison"=study_2_comps,"Median_logFC"=fc2)
      colors2<-rep(NA,nrow(bar_chart_data2))
      colors2[bar_chart_data2$Median_logFC>0]<-"limegreen"
      colors2[bar_chart_data2$Median_logFC<0] <- "blue"
      
      output$Study1_Barchart<-renderPlotly({
        plot_ly(
          data=bar_chart_data1,
          x=~comparison,
          y=~Median_logFC,
          type="bar",
          opacity=0.85,
          marker = list(color = colors)
        ) %>% layout(xaxis=list(showticklabels=FALSE,title="Treatment Comparison"),
                     yaxis=list(title="Median Gene logFC"),
                     margin=3,
                     title=list(text=paste(path,input$Study_1,sep=", "),x=0.05),
                     plot_bgcolor="#D8E0E5",
                     paper_bgcolor="#D8E0E5")
      })
      output$Study2_Barchart<-renderPlotly({
        plot_ly(
          data=bar_chart_data2,
          x=~comparison,
          y=~Median_logFC,
          type="bar",
          opacity=0.85,
          marker = list(color = colors2)
        ) %>% layout(xaxis=list(showticklabels=FALSE,title="Treatment Comparison"),
                     yaxis=list(title="Median Gene logFC"),
                     margin=3,
                     title=list(text=paste(path,input$Study_2,sep=", "),x=0.05),
                     plot_bgcolor="#D8E0E5",
                     paper_bgcolor="#D8E0E5")
      })
   
    }) #observeEvent row selection
  }
}

# Run the application 
shinyApp(ui = ui, server = server)