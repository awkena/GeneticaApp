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
              
        helpText("Enter LETTERS of the English alphabet ONLY.", 
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
                           
                           
  tags$h5("Copyright ", HTML("&copy"),"2021 | Alexander Wireko Kena, PhD.",  
                                   align = 'center', style = 'color: blue')
                           
),
  
 tabPanel("Phenotype",
   fluidRow(
    br(),
    div(id = "dom",
            
        tags$style( "#dom{
         color: orange;
         font-size: 14px;
         font-style: bold;
         font-family: Verdana;
        "),
            
      column(5,
      
      selectInput(inputId = "type", 
                  label = "Select gene interaction",
                  choices = c("Independent assortment" = "IndepA",
                               "Dominant epistasis" = "DE",
                              "Recessive epistasis" = "RE",
                              "Duplicate dominant epistasis" = "DDE",
                              "Duplicate recessive epistasis" = "DRE",
                              "Dominant and recessive epistasis" = "DnRE"),
                  selected = "Independent assortment"),
      
      tags$strong(helpText("For the best display for independent assortment, the number of loci should not exceed 4.", 
               style = 'color: red;')),
      
      tags$strong(helpText("For all digenic interactions, the number of loci = 2.", 
               style = 'color: red;')),
      
      tags$head(
               tags$style(HTML('#pheno{background-color: green}'))),
             actionButton("pheno", "View phenotypes"),
      
      tags$strong(helpText("Generate gametes and crosses first in the Cross panel.", 
               style = 'color: red;')),
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
  
  #' Function to split inputted parental genotypes into alleles
  #' Input data is a string object of parental genotypes with no white spaces
  #' Output is a lsit object of split alleles at each locus
  split_geno <- function(x){
    
    geno <- gsub(" ", "", x) # Remove white spaces
    
    NLoci <- nchar(geno)/2 # Diploid species
    
    # Split inputted genotype into individual loci
    Loci <- strsplit(geno, "(?<=.{2})", perl = TRUE)  
    
    # Split genotypes at individual loci into alleles
    Loci.ls <- strsplit(unlist(Loci), "")
    
    # Rename each locus using LETTERS
    names(Loci.ls) <- paste0("Locus_", LETTERS[1:NLoci])
    
    out <- list(Loci = Loci.ls, NLoci = NLoci)
    
    return(out)
  }
  
  
  #' Check whether locus will segregate using this function
  #' Input data is list object of split alleles for each locus
  #' Output is a string object
  
  Seg <- function(x){
    
    if(x[[1]] == x[[2]]){paste("Locus", toupper(x[[1]]), "is homozygous, hence will not segregate")
    } else(paste("Locus", toupper(x[[1]]), "is heterozygous, hence will segregate"))}
  
  #' Function to generate parental gametes
  #' input data is a list object of split alleles for each locus
  #' Output is a factor object of parental gametes
  gamete <- function(x){
    
    gg <- expand.grid(x) # Generates allele combinations
    
    # Find unique gametes; output is a factor object of unique gametes
    gam <- unique(interaction(gg[1:nrow(gg),], sep = ""))
    
    return(gam)
  }
  
  #' Function to order allele combinations -- Dominant allele first
  #' Input data is a list object; output is a vector of strings
  JJ <- function(x){
    
    if(x[[1]] == x[[2]]){paste0(x[[1]], x[[1]])
    } else if(x[[1]] != x[[2]] & toupper(x[[1]]) == toupper(x[[2]])){
      paste0(toupper(x[[1]]), tolower(x[[1]]))} else if(x[[1]] != x[[2]] & 
            toupper(x[[1]]) != toupper(x[[2]])){paste0(x[[1]], x[[2]])}
    }
  
  #' Function to check if locus is heterozygous
  #' Input data should be a list object; output is a logical object
  is.het <- function(x){
    if(x[[1]] != x[[2]]){
      het <- TRUE
    } else if(x[[1]] == x[[2]]){
      het <- FALSE
    }
    return(het)
    
  }
  
  #' Function to convert genotypes in punnett table to phenotypic
  #' groups
  #' Input data is a list object; output is a list object
  
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
  
  #' Function to convert phenotypic data to digenic epistatic types
  #' Input data is a data frame; output is a data frame
  
  epist <- function(x, loci.m = NULL, loci.p = NULL, NLoci = 2,
                    type = c('DDE', 'DRE', 'DE', 'RE', 'DnRE')){
     
    x <- as.data.frame(x)
    
    mm <- all(sapply(loci.m, is.het))
    pp <- all(sapply(loci.p, is.het))
    
    stopifnot('Both parents must be dihybrids, 
                           and number of loci = 2' = c(NLoci == 2, mm == pp))
    
    # Extract unique phenotypes under independent assortment
    uu <- sort(unique(x$phenotype), decreasing = F)
    
    # Re-code phenotypes based on the type of epistasis
    if(type == 'DE'){
      epi <- ifelse(x$phenotype == uu[1] | x$phenotype == uu[2], 
                    paste(uu[1], uu[2], sep = '||'), ifelse(x$phenotype == uu[3], 
                                                            paste(uu[3]), paste(uu[4])))
      x$epistasis <- epi
      
    } else if(type == 'RE'){
      epi <- ifelse(x$phenotype == uu[1], paste(uu[1]),
                    ifelse(x$phenotype == uu[3] | x$phenotype == uu[4], 
                           paste(uu[3], uu[4], sep = '||'), paste(uu[2])))
      x$epistasis <- epi
      
    } else if(type == 'DDE'){
      epi <- ifelse(x$phenotype == uu[1] | x$phenotype == uu[2] |
                      x$phenotype == uu[3], paste(uu[1], uu[2], uu[3], sep = '||'),
                    paste(uu[4]))
      x$epistasis <- epi
      
    } else if(type == 'DRE'){
      epi <- ifelse(x$phenotype == uu[2] | x$phenotype == uu[3] |
                      x$phenotype == uu[4], paste(uu[2], uu[3], uu[4], sep = '||'),
                    paste(uu[1]))
      x$epistasis <- epi
      
    } else if(type == 'DnRE'){
      epi <- ifelse(x$phenotype == uu[1] | x$phenotype == uu[2] | 
                      x$phenotype == uu[4], paste(uu[1], uu[2], uu[4], sep = '||'),
                    paste(uu[3]))
      x$epistasis <- epi
    }
    
    return(as.data.frame(x))
  }
  
  #' Function to show phenotypes in Punnett square
  #' Input data is a data frame; output is a plot
  gg_plot <- function(data = NULL, fill = NULL, type = NULL,
                      NLoci = NULL, ratio = NULL) {
    
    if(type == 'IndepA'){
      title <- paste('Independent assortment at all loci',
                     paste0('(', ratio, ')'))
    }else if(type == 'DE'){
      title <- paste('Dominant epistasis', paste0('(', ratio, ')'))
    }else if(type == 'RE'){
      title <- paste('Recessive epistasis', paste0('(', ratio, ')'))
    }else if(type == 'DDE'){
      title <- paste('Duplicate dominant epistasis', paste0('(', ratio, ')'))
    }else if(type == 'DRE'){
      title <- paste('Duplicate recessive epistasis', paste0('(', ratio, ')'))
    }else if(type == 'DnRE'){
      title <- paste('Dominant and recessive epistasis', paste0('(', ratio, ')'))
    }
    p <- ggplot2::ggplot(data, aes(x = Var1, y = Var2, 
                                 label = geno, fill = fill)) +
      theme(axis.text=element_text(size  =15*2/NLoci),
            axis.title=element_text(size = 14, face = "bold")) +
      labs(x = 'Female gamete', y = 'Male gamete', fill = 'Phenotype',
           title = title, 
           subtitle = paste('This cross involved', NLoci, 'gene loci.')) +
      theme(plot.title = element_text(family = "serif",
                                      color = "blue",
                                      size = 20,
                                      hjust = 0.5),
            plot.subtitle = element_text(hjust = 0,
                                         family = "serif",
                                         face = "bold",
                                         color = "darkviolet",
                                         size = 18),
            legend.title = element_text(color = "black", size = 18),
            legend.text = element_text(color = "blue", size = 16)) +
      geom_text(colour = "black", size = 6*2/NLoci, fontface = 'bold',
                hjust = 0.5) +
      geom_tile(alpha = 0.5)
    p
    
  }
  
  # Obtain number of loci from inputted maternal genotype
  NLoci <- eventReactive (input$Generate,{
    
    split_geno(input$Geno_m)$NLoci
    
    }) 
  
  # Generate Maternal gametes
  Loci_m <- eventReactive (input$Generate,{
 
    split_geno(input$Geno_m)$Loci
  })
  
  
  # Maternal gametes
  gamet_m <- eventReactive (input$Generate,{
    
    gamete(Loci_m()) 
    
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
    
    split_geno(input$Geno_p)$Loci
  })
  
  gamet_p <- eventReactive (input$Generate,{
    
    gamete(Loci_p())
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
  
  #' Convert genotypes in Punnett square to phenotypic groups assuming complete dominance
  #' Split genotypes in Punnett square into individual loci
  #' Output is a data frame object
  testdf <- eventReactive(input$pheno,{
    
    input$pheno
    
    if(input$type == 'IndepA' | input$type == 'DE'|input$type == 'RE'| input$type =='DDE'|
       input$type == 'DRE'| input$type == 'DnRE'){
    
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
  
  tt <- eventReactive(input$pheno,{
    
    input$pheno
    
     isolate({epist(testdf(), loci.m = Loci_m(), loci.p = Loci_p(),
              NLoci = NLoci(), type = input$type)
      })
    })
  
  # Get phenotypic ratio when gene interaction is independent assortment
  pr_ida <- eventReactive(input$pheno,{
    
    input$pheno
    paste(sort(table(testdf()$phenotype), decreasing = T), collapse = ':')
    
  })
  
  # Get phenotypic ratio when gene interaction is epistasis
  pr_epis <- eventReactive(input$pheno,{
    
    input$pheno
    paste(table(tt()$epistasis), collapse = ':')
    
  })
  
  #' Color genotypes in punnett square based on phenotype
  #' Assuming complete dominance at each locus
  #' Easiest way is to use the ggplot2 package
  
  observeEvent(input$pheno,{
    
    req(testdf(), cancelOutput = TRUE)
    
    output$comp_plot <- renderPlot({
      isolate({if(input$type == 'IndepA'){
          isolate({gg_plot(data = testdf(), fill = testdf()$phenotype, 
                           NLoci = NLoci(), type = input$type, ratio = pr_ida())  
            })
      }else if (input$type == 'DE'|input$type == 'RE'| input$type =='DDE'|
                input$type == 'DRE'| input$type == 'DnRE'){
        isolate({gg_plot(data = tt(), fill = tt()$epistasis, 
                         NLoci = NLoci(), type = input$type, ratio = pr_epis()) 
          })
      }
      })
   
    })
    
  })
  
  
  observeEvent(input$pheno,{
    
    req(testdf(), cancelOutput = TRUE)
  
 
  output$pheno_ratio <- renderUI({
    
    isolate({ if(input$type == 'IndepA'){
    
      str7 <- isolate({paste("Assuming complete dominance, there would be", 
                  length(unique(testdf()$phenotype)), "phenotypic classes.")})
        
    str8 <- isolate({paste("The phenotypic ratio is", pr_ida())})
    
    isolate({HTML(paste(str7, str8,  sep = "<br/> <br/>"))})
    
    } else if(input$type == 'DE'|input$type == 'RE'| input$type =='DDE' |
               input$type == 'DRE'| input$type == 'DnRE'){
      
      str7 <- isolate({paste("There would be", 
              length(unique(tt()$epistasis)), "phenotypic classes.")})
      
      str8 <- isolate({paste("The digenic epistatic ratio is", pr_epis())})
      
      isolate({HTML(paste(str7, str8,  sep = "<br/> <br/>"))})
      
    } 
    })
  
  }) 
  
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

