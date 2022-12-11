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
library(reshape2)
library(ggplot2)

ui <- fluidPage(theme = shinytheme("united"),
      shinyjs::useShinyjs(),
            
     tabsetPanel(
              
      tabPanel("Cross",
      # Application description
      h4(tags$b(div('This app generates the genotypes of progenies of genetic 
                    crosses for any number of gene loci in DIPLOID species.')), 
                          style = 'color: black;',align = 'center'),
                       
   hr(),
   
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
  
 tabPanel("Phenotype",
   fluidRow(
    
    div(id = "dom",
            
        tags$style( "#dom{
         color: orange;
         font-size: 14px;
         font-style: bold;
         font-family: Verdana;
        "),
            
      column(5,
      
      checkboxInput("dom", "Assume Complete Dominance", value = TRUE),
      
      helpText("For the best output, the number of loci should not exceed 4.", 
               style = 'color: red;'),
        
                 tags$head(
               tags$style(HTML('#pheno{background-color: green}'))),
             actionButton("pheno", "View phenotypes"),
      
      helpText("Generate gametes and crosses first in the Cross panel.", 
               style = 'color: red;'),
      )
  ),
  
  div(
    column(6,offset = 1,
           
           tags$head(tags$style("#pheno_ratio{color: blue;
  font-size: 15px;
  font-style: bold;
  font-family: Verdana;
  word-wrap: break-word;
  }"
                                
           )),
           
           h4(tags$strong("Phenotypic ratio:"),
              br(), br(),
              htmlOutput("pheno_ratio")), 
           
  
  ))),
  
  
  hr(),
  
  fluidRow(
    div(tags$strong("Punnett square output"), style = 'text-indent: 40px;
             font-size: 150%'),
    column(12,
    plotOutput("comp_plot")
    
  )
  )
  
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
  
  # Obtain number of loci from inputted maternal genotype
  NLoci <- eventReactive (input$Generate,{
     #req(input$Geno_m, input$Generate, cancelOutput = TRUE)
    
    nchar(gsub(" ", "", input$Geno_m))/2 }) # Diploid species
  
  # Generate Maternal gametes
  Loci_m <- eventReactive (input$Generate,{
 
    # req(input$Geno_m, input$Generate, cancelOutput = TRUE)
    
    GenoM <- gsub(" ", "", input$Geno_m)
    
    # Obtain number of loci from inputted maternal genotype 
    # NLoci_m <- nchar(GenoM)/2 # Diploid species
    
    # Split inputted maternal genotype into individual loci
    LociM <- strsplit(GenoM, "(?<=.{2})", perl = TRUE)
    
    # Split genotypes at individual maternal loci into alleles
    Loci.lsm <- strsplit(unlist(LociM), "")
    
    # Rename each locus using LETTERS
    names(Loci.lsm) <- paste0("Locus_", LETTERS[1:NLoci()])
    
    Loci.lsm
  })
  
  gamet_m <- eventReactive (input$Generate,{
    
    #req(input$Generate, Loci_m(), input$Geno_m, cancelOutput = TRUE)
    
    Gamete1 <- expand.grid(Loci_m())
    
    gam_m <- as.character(unique(interaction(Gamete1[1:nrow(Gamete1),], 
                                                      sep = "")))
    gam_m
    
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
  Loci_p <- eventReactive (input$Generate,{
    
    #req(input$Geno_p, input$Generate, cancelOutput = TRUE)
    
    GenoP <- gsub(" ", "", input$Geno_p)
    
    # Obtain number of loci from inputted genotype 
    # NLoci_p <- nchar(GenoP)/2 # Diploid species
    
    # Split inputted genotype into individual loci
    LociP <- strsplit(GenoP, "(?<=.{2})", perl = TRUE)
    
    # Split genotypes at individual loci into alleles
    Loci.lsp <- strsplit(unlist(LociP), "")
    
    # Rename each locus using LETTERS
    names(Loci.lsp) <- paste0("Locus_", LETTERS[1:NLoci()])
    Loci.lsp
  })
  
  gamet_p <- eventReactive (input$Generate,{
    #req(input$Geno_p, input$Generate, Loci_p(), cancelOutput = TRUE)
    
      Gamete2 <- expand.grid(Loci_p())
    
    gam_p <- as.character(unique(interaction(Gamete2[1:nrow(Gamete2),], 
                                                      sep = "")))
    gam_p
    
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
  
  pun1 <- eventReactive(input$Cross,{
    #req(input$Generate, gamet_m(), gamet_p(), input$Geno_m, input$Geno_p,
        #cancelOutput = TRUE)
    
    
    # input$Generate
    
    # Create an empty matrix to save results
    
    pun <- matrix(NA, nrow = length(gamet_m()), ncol = length(gamet_p()))
    
    rownames(pun) <- gamet_m()
    colnames(pun) <- gamet_p()
    
    
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
    
    pun
    
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
  
  
<<<<<<< HEAD
  #' Function to convert genotypes in punnett table to phenotypic
  #' groups
  
=======
  # Function to convert genotypes in punnett table to phenotypic
  # groups
>>>>>>> e196f6efcb03d97f9cf212681d86bc2db668e25f
  nn <- function(x, y){
    for(i in 1:length(x)){
      
      for(j in 1:y){
        if(x[[i]][j] == toupper(x[[i]][j])){ substr(x[[i]][j],2,2) <- '_'
        
        }else if(x[[i]][j] != toupper(x[[i]][j]) & x[[i]][j] != tolower(x[[i]][j]))
        { substr(x[[i]][j],2,2) <- '_'
        } else(x[[i]][j] <- x[[i]][j])
      }
      x[[i]] <- paste0(x[[i]], collapse = "")
    }
    return(x)
  }
  
<<<<<<< HEAD
  
  #' Convert genotypes in Punnett square to phenotypic groups assuming complete dominance
  #' Split genotypes in Punnett square into individual loci
=======
  #' Convert genotypes in Punnett square to phnotypic groups assuming complete dominance
  #' Split genotypes in punnett square into individual loci
>>>>>>> e196f6efcb03d97f9cf212681d86bc2db668e25f
  #' Output is a data frame object
  testdf <- eventReactive(input$pheno,{
    
    input$pheno
    
    if(input$dom == TRUE){
    
    x <- strsplit(as.vector(pun1()), "(?<=.{2})", perl = TRUE)
      
    

     # Use nn function to convert genotypes to phenotypic groups

     aa <- nn(x, NLoci())
     
     #' Melt punnett square into a data frame to plot genotypes and 
     #' color based on phenotypic groups
     
     bb <-  reshape2::melt(pun1(), value.name = "geno")
     
     #' Add phenotypic group for each genotype
     
     dd <- as.vector(unlist(aa))
     
     bb$phenotype <- dd
     
       bb
       }
  
  })
  
  #' Color genotypes in punnett square based on phenotype
  #' Assuming complete dominance at each locus
  #' Easiest way is to use the ggplot2 package
  
  observeEvent(input$pheno,{
    
    req(input$pheno, pun1(), cancelOutput = TRUE)
    
    
    output$comp_plot <- renderPlot({
  
  p <- ggplot2::ggplot(testdf(), aes(x = Var1, y = Var2, 
                          label = geno, fill=phenotype)) +
    isolate({theme(axis.text=element_text(size=15*2/NLoci()),
          axis.title=element_text(size=14,face="bold"))}) +
    isolate({labs(x = 'Female gamete', y = 'Male gamete',
         title = 'Phenotypic classes for genotypes, assuming complete dominance at all loci.',
         subtitle = paste('This cross involved', NLoci(), 'gene loci.'))}) +
    theme(plot.title = element_text(family = "serif",
                                    color = "blue",
                                    size = 20,
                                    hjust = 0.5),
          plot.subtitle = element_text(hjust = 0,
                                       family = "serif",
                                       face = "bold",
                                       color = "darkviolet",
                                       size = 18)) +
    isolate({geom_text(colour = "black", size = 6*2/NLoci(), fontface = 'bold',
              hjust = 0.5)}) +
    geom_tile(alpha = 0.5)
  p
  
    })
    
  })
  
  observeEvent(input$pheno,{
    
    req(input$pheno, pun1(), cancelOutput = TRUE)
  
 
  output$pheno_ratio <- renderUI({
    
    str7 <- isolate({paste("Assuming complete dominance, there would be", 
                  length(unique(testdf()$phenotype)), "phenotypic classes.")})
    pr <- paste(sort(table(testdf()$phenotype), decreasing = T), collapse = ':')
    
  str8 <- isolate({paste("The phenotypic ratio is", pr)})
    
    isolate({HTML(paste(str7, str8,  sep = "<br/> <br/>"))
    })
  
  }) 
  
  })
}


# Run the application 
shinyApp(ui = ui, server = server)

