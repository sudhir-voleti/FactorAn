#################################################
#               Factor Analysis                 #
#################################################

library("shiny")
#library("foreign")

shinyUI(fluidPage(
  tags$head(includeScript("google_analytics.js")),
  title = "Factor analysis",
  # titlePanel("Regression Tree"),
  titlePanel(title=div(img(src="logo.png",align='right'),"Factor analysis")),
  
  # Header:
  #titlePanel("Factor analysis"),

  # Input in sidepanel:
  sidebarPanel(
    # Upload data:
    fileInput("file", "Upload input data (csv file with header)"),  
    uiOutput("colList"),
    htmlOutput("fselect"),
    textInput('fname',label = "Enter Factor Name (seperated by comma)"),
    sliderInput("cutoff", "Cut-off for factor loadings(for Plotting only)", min = 0,  max = 1, value = 0.25),
    htmlOutput("xaxis"),
    htmlOutput("yaxis"),
    
    sliderInput("cutoffcorr", "Cut-off for Correlation in Factors vs Variable", min = 0,  max = 1, value = 0.25),
    br()
     # submitButton(text = "Apply Changes", icon("refresh"))
  ),
  # Main:
  mainPanel( 
    
    tabsetPanel(type = "tabs",
                #

tabPanel("Overview",
# h5(p("Factor Analysis")), 
# p("Factor analysis is a statistical method used to describe variability among observed, correlated variables in terms of a potentially lower 
# number of unobserved variables called factors. For example, it is possible that variations in four observed variables mainly reflect the variations 
# in two unobserved variables. Factor analysis searches for such joint variations in response to unobserved latent variables. The observed variables 
# are modelled as linear combinations of the potential factors, plus error terms. The information gained about the interdependencies between 
# observed variables can be used later to reduce the set of variables in a dataset. Computationally this technique is equivalent to low rank 
# approximation of the matrix of observed variables. Factor analysis originated in psychometrics, and is used in behavioral sciences, social sciences,
# marketing, product management, operations research, and other applied sciences that deal with large quantities of data.",align="justify"),
# a(href="http://en.wikipedia.org/wiki/Factor_analysis","- Wikipedia"),
h4(p("How to use this shiny application")),
p("This shiny application require one data input from the user. To do so, click on the Browse (in left side-bar panel) and upload the csv data input file.
  Note that this application can read only csv file(comma delimited file), so if you don't have csv input data file, first convert your data in csv format 
  and then proceed. Make sure you have top row as variable names and first column as respondent id/name in csv file",align="justify"),
p("Once csv file is uploaded successfully, application will fit a factor model with optimal factors from parallel Analysis and various 
results will be showed in the above tabs. In left-side bar panel you can change the parameters value and correspondingly new results 
  will be showed",align="justify"),
br(),
h4(p("Download Sample Input File")),
downloadButton('downloadData', 'Download Example file'),
p("Please note that download will not work with RStudio interface. Download will work only in web-browsers. So open this app in a web-browser and then download the example file. For opening this app in web-browser click on \"Open in Browser\" as shown below -"),
img(src = "example1.png") #, height = 280, width = 400

),
    
                tabPanel("Summary",
                         (h4(p("Correlation"))),
                         (plotOutput("corplot",height = 850, width = 850)),
                         (h4(p("Test Summary"))),(textOutput("text1")),(textOutput("text2")),(textOutput("text3")),
                         (h4(p("Factors Loadings Summary"))),
                          (verbatimTextOutput("mat")),
                         
                         (h4(p("Uniqueness table - "))),
                         (dataTableOutput("uni")),
#                          (textOutput("text4")),
                         plotOutput("plot1",height = 600, width = 850)),
                tabPanel("Loadings",dataTableOutput("loadings")),
                
#                tabPanel("Scores",tableOutput("scores")),   # origi code
                # my edits 16-9-2017 below:
                tabPanel("Scores", 	# tab name
	                      br(),
                        downloadButton('downloadDataX', 
		                        'Download Factor Scores File (Works only in browser)'), 
	                      br(),br(),
	                      dataTableOutput("scores")),
                tabPanel("Parameter Tests",
                         
                         h2("KMO Test for Factor Adequacy"), br(),
                         h4("The Kaiser-Meyer-Olkin (KMO) is a test conducted to examine the strength of the partial correlations
                            (how the factors explain each other) between the variables. KMO values closer to 1.0 are consider 
                            ideal while values less than 0.5 are unacceptable. According to Kaiser’s (1974) guidelines, 
                            a suggested cutoff for determining the factorability of the sample data is KMO ≥ 0.6."),
                         br(),
                         verbatimTextOutput("dummy"),
                         br(),
                         h2("Bartlett’s test of Sphericity"), br(),
                         h4("Bartlett’s Test of Sphericity compares an observed correlation matrix to the
                            identity matrix. An identity correlation matrix means your variables are 
                            unrelated and not ideal for factor analysis. The null hypothesis of the test
                            is that the variables are orthogonal, i.e. not correlated. 
                            
                            Look for a significant test statistic (p.value < 0.05) that 
                            rejects this null (i.e. proves that the correlation matrix is NOT an identity matrix)."),
                         verbatimTextOutput("dummy2"), br(),
                         h2("Determinant for Factorizability"), br(),
                         h4("This is a simple test that looks for a positive determinant value, which implies the matrix in question is invertible and that the factorization can proceed smoothly."),
                         p("The determinant of the correlation matrix for this dataset is :"), br(),
                         verbatimTextOutput("dummy3"),
                         h2("Visualization of Principal Axis Factoring"),
                         plotOutput("plot_PA")
                         ),
                
                tabPanel("Factor vs Variables",plotOutput("plot20",height = 600, width = 850)),
                tabPanel("Factor vs Variables 2",plotOutput("plot2",height = 600, width = 850)),
                tabPanel("Factor vs Users",plotOutput("plot3",height = 600, width = 850)),
                tabPanel("Data",dataTableOutput("table")) 
    )
  ) 
) 
)
