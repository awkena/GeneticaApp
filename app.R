#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyFeedback)
library(shinyjs)
library(shinythemes)
library(MASS)
library(reactable)

ui <- fluidPage(theme = shinytheme("united"),
                shinyjs::useShinyjs(),
      
      tabsetPanel(
        
        tabPanel("Genetic Crosses",
       # Application description
       h4(tags$b(div('This app generates the genotypes of progenies of genetic 
    crosses for any number of gene loci in DIPLOID species.')), 
          style = 'color: black;',align = 'center'),
       
       hr(),
       tags$img(src = "AGRILOGO.jpg", height = 55, width = 50, 
                align = 'right'),
       
       tags$img(src = "KNUST_logo.jpg", height = 60, width = 50, 
                align = 'right'),
       
       h3(tags$strong("GENETICA"), style = 'color: green', align = 'center'),
       
       fluidRow(
         
         div(
           h4(tags$strong("Input panel"), style = 'text-indent: 50px;'),
           
           column(4,
                  
            helpText("Enter LETTERS of the English alphabets ONLY.", 
                     style = 'color: red;'),
            helpText("Use a SINGLE letter for each allele at a locus.", 
                     style = 'color: red;'),
            helpText("Genotypes MUST be in multiples of two letters.", 
                     style = 'color: red;'),
            
            textInput("Geno_m", "Enter genotype of maternal parent:",
                      value = "AaBb"),
            
            textInput("Geno_p", "Enter genotype of paternal parent:",
                      value = "AaBb"),
            
            tags$head(
              tags$style(HTML('#Generate{background-color: green}'))),
            actionButton("Generate", "Generate gametes", icon("sync")),
            
            tags$head(
              tags$style(HTML('#Cross{background-color: green}'))),
            actionButton("Cross", "Cross parents", icon("times")),
            
           ),
        
           div(
             column(1,
                )
           ),
           
           div(
             column(6,
                    
              tags$head(tags$style("#gamete_m{color: blue;
                font-size: 15px;
                font-style: bold;
                font-family: Verdana;
                word-wrap: break-word;
                }"
                                         
                    )),
                    
            h4(tags$strong("Maternal gamete(s):"),
               br(), br(),
            htmlOutput("gamete_m")), 
            
            tags$br(),
            tags$head(tags$style("#gamete_p{color: blue;
                font-size: 15px;
                font-style: bold;
                font-family: Verdana;
                word-wrap: break-word;
                }"
                                         
                    )),
                    
              h4(tags$strong("Paternal gamete(s):"),
              br(), br(),
              htmlOutput("gamete_p")),  
             )
           )
         
       
        )
        ),
       hr(),
       
       fluidRow(
         
         div(tags$strong("Punnett square output"), style = 'text-indent: 40px;
             font-size: 150%'),
           
           column(10,
                  
                  reactableOutput("punnett"),
                  
          br(),br(),
          
          div(tags$strong("Summary of progeny genotypes in the Punnett square:"), 
              style = "text-indent: 20px; color: blue; font-size: 120%"),
          reactableOutput("summary"),   
          
          br(), br(),
          
          tags$head(tags$style("#progeny{color: blue;
                font-size: 18px;
                font-style: bold;
                font-family: Verdana;
                word-wrap: break-word;
                }"
                               
          )),
          htmlOutput("progeny")
       
           )
         
       ),
       
      
       hr(),
      
       
tags$h5("Copyright ", HTML("&copy"),"2021 | Alexander Wireko Kena, PhD. | Benjamin Annor, PhD. |
        Stephen Amoah, PhD. | Richard Akromah, PhD.",  
               align = 'center', style = 'color: blue')
       
       ),
                 
 tabPanel("Help",
          
    tags$strong(h3("INSTRUCTIONS TO USERS:", style = 'color: black;')),
    
    hr(),
    
    tags$ol(type = '1',
    tags$li("Enter the genotypes of the maternal and paternal parents."),
    tags$li("Enter letters of the English alphabets only."),
    tags$li("Use a single letter to represent each allele at a locus."),
    tags$li("DO NOT mix letters with numbers to represent an allele at any locus."),
    tags$li("Diploid species must have TWO (2) alleles at each locus."),
    tags$li("Thus, inputted parental genotypes must be in multiples of two letters."),
    tags$li("Click on the Generate gametes button to see the gametes of the parents."),
    tags$li("Cross the two parents by clicking on the Cross parents button to 
            view the punnett square.")
    ),
          
          )
 
      )
)


# Define server logic 
server <- function(input, output){
  
  # Check whether locus will segregate using this function
  Seg <- function(x){
    
    if(x[[1]] == x[[2]]){paste("Locus", toupper(x[[1]]), "is homozygous, hence will not segregate")
    } else(paste("Locus", toupper(x[[1]]), "is heterozygous, hence will segregate"))}
  
# Function to order allele combinations -- Dominant allele first
  JJ <- function(x){
    
    if(x[[1]] == x[[2]]){paste0(x[[1]], x[[1]])
    } else if(x[[1]] != x[[2]] & toupper(x[[1]]) == toupper(x[[2]])){
      paste0(toupper(x[[1]]), tolower(x[[1]]))} else if(x[[1]] != x[[2]] & 
                    toupper(x[[1]]) != toupper(x[[2]])){paste0(x[[1]], x[[2]])}
  }
  # Generate Maternal gametes
  Loci_m <- reactive({
    
    input$Generate
    
    req(input$Geno_m, input$Generate, cancelOutput = TRUE)
    
    GenoM <- isolate({gsub(" ", "", input$Geno_m)})
    
    # Obtain number of loci from inputted maternal genotype 
    NLoci_m <- isolate({nchar(GenoM)/2 }) # Diploid species
    
    # Split inputted maternal genotype into individual loci
    LociM <- isolate({strsplit(GenoM, "(?<=.{2})", perl = TRUE)})
    
    # Split genotypes at individual maternal loci into alleles
    Loci.lsm <- isolate({strsplit(unlist(LociM), "")})
   
    # Rename each locus using LETTERS
    names(Loci.lsm) <- isolate({paste0("Locus_", LETTERS[1:NLoci_m])})
    isolate({Loci.lsm})
  })
    
  gamet_m <- reactive({
      
    req(input$Generate, Loci_m(), input$Geno_m, cancelOutput = TRUE)
    
      input$Generate
    
      Gamete1 <- isolate({expand.grid(Loci_m())})
    
    gam_m <- isolate({as.character(unique(interaction(Gamete1[1:nrow(Gamete1),], 
                                                      sep = "")))})
    isolate({gam_m})
  
  })
  
  # Output maternal gametes
  observeEvent(input$Generate,{
    
    req(Loci_m(), gamet_m(), input$Geno_m, cancelOutput = TRUE)
    
    output$gamete_m <- renderUI({
      str3 <- isolate({unname(unlist(sapply(Loci_m(), Seg)))})
      
      str1 <- isolate({paste(gamet_m(), collapse = ", ")})
      
      str2 <- isolate({if(length(gamet_m()) == 1){paste("There is only one unique maternal gamete:")
      } else if(length(gamet_m()) > 1){paste("There are", length(gamet_m()), 
                                               "unique maternal gametes:")}})
      
      isolate({HTML(paste(paste(str3,  collapse = "; "), str2, str1,  sep = "<br/> <br/>"))
      })
    })
  })
    
  # Generate Paternal gametes
  Loci_p <- reactive({
    
    input$Generate
    
    req(input$Geno_p, input$Generate, cancelOutput = TRUE)
    
    GenoP <- isolate({gsub(" ", "", input$Geno_p)})
    
    # Obtain number of loci from inputted genotype 
    NLoci_p <- isolate({nchar(GenoP)/2}) # Diploid species
    
    # Split inputted genotype into individual loci
    LociP <- isolate({strsplit(GenoP, "(?<=.{2})", perl = TRUE)})
    
    # Split genotypes at individual loci into alleles
    Loci.lsp <- isolate({strsplit(unlist(LociP), "")})
    
    # Rename each locus using LETTERS
    names(Loci.lsp) <- isolate({paste0("Locus_", LETTERS[1:NLoci_p])})
    isolate({Loci.lsp})
  })
    
  gamet_p <- reactive({
    req(input$Geno_p, input$Generate, Loci_p(), cancelOutput = TRUE)
    
    input$Generate
    
    Gamete2 <- isolate({expand.grid(Loci_p())})
    
    gam_p <- isolate({as.character(unique(interaction(Gamete2[1:nrow(Gamete2),], 
                                                      sep = "")))})
    isolate({gam_p})
  
  })
  
  # Output Paternal gametes
  observeEvent(input$Generate,{
    req(Loci_p(), gamet_p(),input$Geno_p, cancelOutput = TRUE)
    
    output$gamete_p <- renderUI({
      str3 <- isolate({unname(unlist(sapply(Loci_p(), Seg)))})
      
      str1 <- isolate({paste(gamet_p(), collapse = ", ")})
      
      str2 <- isolate({if(length(gamet_p()) == 1){paste("There is only one unique paternal gamete:")
      } else if(length(gamet_p()) > 1){paste("There are", length(gamet_p()), 
                                             "unique paternal gametes:")}
      })
      
      isolate({HTML(paste(paste(str3,  collapse = "; "), str2, str1,  sep = "<br/> <br/>"))
      }) 
    })
  })
  
  # Generate punnett square
  
  pun1 <- reactive({
    req(input$Generate, gamet_m(), gamet_p(), input$Geno_m, input$Geno_p,
        cancelOutput = TRUE)
    
    input$Cross
    input$Generate
    
    # Create an empty matrix to save results
    
    pun <- isolate({matrix(NA, nrow = length(gamet_m()), ncol = length(gamet_p()))
      })
    rownames(pun) <- isolate({gamet_m()})
    colnames(pun) <- isolate({gamet_p()})
   
     isolate({
    for(i in rownames(pun)){
      for(j in colnames(pun)){
        
        Geno <- paste0(sort(unlist(strsplit(c(i, j),""))),
                       collapse = "")
        
        EE <- strsplit(Geno, "(?<=.{2})", perl = TRUE)
        FF <- strsplit(unlist(EE), "")
        KK <- paste0(sapply(FF, JJ), collapse = "")
        
        pun[i,j] <- KK
        
      }
      
      
    }
     })
      isolate({pun})
  
  })
  
  observeEvent(input$Cross,{
    
    req(input$Generate, pun1(), input$Geno_m, input$Geno_p, cancelOutput = TRUE)
    
    output$punnett <- renderReactable({
    isolate({reactable(data.frame(pun1()), fullWidth = F, defaultPageSize = 10,
      columns = list(.rownames = colDef(name = "Gametes", sortable = TRUE)),
      
              bordered = T, highlight = T, resizable = T)
    })
    })
      
  })
  
  observeEvent(input$Cross,{
    req(input$Generate, pun1(), input$Geno_m, input$Geno_p, cancelOutput = TRUE)
    
    output$summary <- renderReactable({
      
      # Progeny genotypic frequencies in fractions
      Geno.freq <- isolate({fractions(rev(table(pun1()))/sum(table(pun1())))})
      tab1 <- isolate({data.frame(Geno.freq)})
      
      if(dim(tab1)[1] > 1 & dim(tab1)[2] > 1){
        tab1$freq2 <- isolate({as.character(fractions(tab1$Freq))})
      } else(tab1$freq2 <- "1/1")
      
      if(dim(tab1)[1] > 1 & dim(tab1)[2] > 1){
        colnames(tab1) <- isolate({c("Genotype", "Freq", "Frequency")})
      } else(colnames(tab1) <- isolate({c("Freq", "Frequency")}))
      
      tab1 <- isolate({data.frame(tab1)})
      
      if(dim(tab1)[1] > 1 & dim(tab1)[2] > 1){
      isolate({reactable(tab1, fullWidth = F, defaultPageSize = 9, bordered = T,
                         highlight = T, rownames = T, resizable = T)
      })
      } else(isolate({reactable(tab1, fullWidth = F, defaultPageSize = 9, bordered = T,
          columns = list(.rownames = colDef(name = "Genotype", sortable = TRUE)),                  
                                highlight = T, rownames = T, resizable = T)}))
    })
  })
  
  # Punnett square summaries
  observeEvent(input$Cross,{
    
    req(Loci_p(), gamet_p(), input$Generate, cancelOutput = TRUE)
    
    output$progeny <- renderUI({
      
      str4 <- isolate({paste("The number of unique gentoypes in the punnett square
                             is", length(unique(noquote(pun1()))))})
      str5 <- isolate({paste("The minimum population size is", 
                             sum(table(pun1())))})
      
      isolate({HTML(paste(str4, str5,  sep = "<br/> <br/>"))
    })
    
    
    }) 
  })
  
}


# Run the application 
shinyApp(ui = ui, server = server)

