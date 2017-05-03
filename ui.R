

## NOTEs
## Maybe don't output excel files, output csv so that we don't have to use rJava
## Align to actual zebrafish genome instead of custom genome
## May have to find R versions of bowtie2 and samtools
library(devtools)

library(ggplot2)
library(CrispRVariants)
library(Rsamtools)
library(ShortRead)
library(rtracklayer)
library(zlibbioc)
library(seqTools)
library(Biostrings)
library(BSgenome)
library(GenomicFeatures)
library(gdata)
library(dplyr)
library(fishplot)
library(flowCore)
library(Rtsne)

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Old Faithful Geyser Data"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
       sliderInput("bins",
                   "Number of bins:",
                   min = 1,
                   max = 50,
                   value = 30)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
       plotOutput("distPlot")
    )
  )
))
