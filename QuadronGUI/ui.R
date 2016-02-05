################################################################################
#//    Aleksandr B. Sahakyan (aleksahak [at] cantab.net), Cambridge 2016     \\#
################################################################################
library(shiny)

shinyUI(pageWithSidebar(

  headerPanel("Quadron - sequence-based predictor of DNA G-quadruplex formation.",
              windowTitle = "QuadronGUI"),

  #***#######################################################################
  sidebarPanel(
      imageOutput("logo", width = "100%", height = "300px", inline = FALSE),

      #--#########################################################
      wellPanel(
        selectInput("nCPU", "Number of CPU (warp) cores:", list("1 core"  = 1,
                                                                "2 cores" = 2,
                                                                "4 cores" = 4,
                                                                "8 cores" = 8,
                                                               "12 cores" = 12,
                                                               "16 cores" = 16,
                                                               "24 cores" = 24))
      )
      #--#########################################################
  ),
  #***#######################################################################

  #***#######################################################################
  mainPanel(

      #--#########################################################
      conditionalPanel(
        condition = "output.text == 'Idle'",
        h5("A program to predict sequence-based and context dependent strength (G4-seq mismatch levels in quadruplex induced polymersae stalling assay) of quadruplex structures in any DNA sequence. The on-line version is only suitable for small sequences and cannot utilise more than one CPU core. To use the full power of the program, please download the stand-alone version, usable with or without (terminal access) this graphical user interface. For the best performance, utilize as many computing processors as possible. The full instruction sets along with the source data are available from ", a("HERE", href="http://quadron.atgcdynamics.org/"),".")
      ),
      #--#########################################################

      #--#########################################################
      wellPanel(
        h4("Input"),
        fileInput('FastaFile', 'Upload the DNA fasta file:',accept=NULL),
        checkboxInput(inputId = "pasteseq", label = "Paste the sequence instead.", FALSE),
        conditionalPanel(
          condition = "input.pasteseq == true",
          textInput("sequence", label="Place the cursor on the box below and paste the sequence without blank spaces and new lines:", value = "")
        )
      ),
      #--#########################################################

      #--#########################################################
      verbatimTextOutput("text"),
      #--#########################################################
      conditionalPanel(
        condition = "output.text != 'Idle'",
        downloadButton('downloadQuadronOUT', 'Download')
      ),
      #--#########################################################
      tableOutput("out")
  )
  #***#######################################################################

))
################################################################################
