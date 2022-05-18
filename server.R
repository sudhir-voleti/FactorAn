#################################################
#               Factor Analysis                 #
#################################################

library("shiny")
library("nFactors")
library("qgraph")
library("corrplot")
library("magrittr")
library('psych')
library('car')

shinyServer(function(input, output) {

  Dataset <- reactive({
    if (is.null(input$file)) { return(NULL) }
    else{
    Dataset <- as.data.frame(read.csv(input$file$datapath ,header=TRUE, sep = ","))
    rownames(Dataset) = Dataset[,1]
    Dataset1 = Dataset[,2:ncol(Dataset)]
    #Dataset = t(Dataset)
    return(Dataset1)
    }
  })

fname <- reactive({
  if(length(strsplit(input$fname,',')[[1]])==0){return(NULL)}
  else{
    return(strsplit(input$fname,',')[[1]])
  }
})


filtered_dataset <- reactive({if (is.null(input$file)) { return(NULL) }
  else{
    Dataset <- Dataset() %>% dplyr::select(!!!input$selVar)
    return(Dataset)
  }})


output$colList <- renderUI({
  varSelectInput("selVar",label = "Select Variables",data = Dataset(),multiple = TRUE,selectize = TRUE,selected = colnames(Dataset()))
})

# output$table22 <- renderTable ({ 
#   round(cor(Dataset()),2) 
#                                       })
output$corplot = renderPlot({
  if (is.null(input$file)) { return(NULL) }
    else{
      my_data = filtered_dataset()
      cor.mat <- round(cor(my_data),4)
      corrplot(cor.mat, 
               type = "upper",    # upper triangular form
               order = "hclust",  # ordered by hclust groups
               tl.col = "black",  # text label color
               tl.srt = 45,
               method = "number")  
      
    }
  })



output$table <- renderDataTable({ Dataset() },options = list(pageLength=25))
  
nS = reactive ({    
    if (is.null(input$file)) { return(NULL) }
    else{
      ev = eigen(cor(filtered_dataset(), use = 'pairwise.complete.obs'))  # get eigenvalues
      ap = parallel(subject=nrow((filtered_dataset())),var=ncol((filtered_dataset())),rep=100,cent=.05)
      nS = nScree(ev$values, aparallel= ap$eigen$qevpea)
}
})

output$fselect <- renderUI({ 
  if (is.null(input$file)) { return(NULL) }
  else{
  numericInput("fselect", "Number of Factors:", unlist((nS())[1])[3])
  }
  })

fselect = reactive({
  if (is.null(input$file)) { return(NULL) }
  else{
    
  fselect=input$fselect
  return(fselect)
  }
})

fit = reactive ({ 
  if (is.null(input$file)) { return(NULL) }
  else{
    
  fit = factanal(na.omit(filtered_dataset()), fselect() , scores="Bartlett", rotation="varimax");
  return (fit)
  }
  }) 

output$dummy <- renderPrint(
  
  if (is.null(input$file)) { return(NULL) }
  else {
    KMO(r=cor(filtered_dataset()))
  }
  
)

output$dummy2 <- renderPrint(
  
  if (is.null(input$file)) { return(NULL) }
  else {
    cortest.bartlett(filtered_dataset())
  }
  
)

output$dummy3 <- renderPrint(
  
  if (is.null(input$file)) { return(NULL) }
  else {
    det(cor(filtered_dataset()))
  }
  
)

output$plot_PA <- renderPlot(
  if(is.null(input$file)) {return(NULL)}
  else{
    
    fa.none <- fa(r=filtered_dataset(), 
                  nfactors = fselect(), # comes from UI
                  # covar = FALSE, SMC = TRUE,
                  fm="pa", # type of factorAn to use (“pa” == principal axis factoring)
                  max.iter=100, # (50 is the default, but we have changed it to 100
                  rotate="varimax")
    
    fa.diagram(fa.none)
  }
)

output$xaxis <- renderUI({
  
  if (is.null(input$file)) { return(NULL) }
  else {
  if(is.null(fname())){
    n =(fselect())
    list = character(0)
    for (i in 1:n) { 
      temp = paste("Factor",i)
      list = c(list, temp)
    }
    
    selectInput("xaxis", "Choose Factor for plotting on X axis",
                list, selected = "Factor 1")
  }else{
      temp = fname()
       selectInput("xaxis", "Choose Factor for plotting on X axis",
                temp, selected = temp[1])
    
  }
    
    
  }
  
  
})

output$yaxis <- renderUI({
  
  if (is.null(input$file)) { return(NULL) }
  else {
    if(is.null(fname())){
      n =(fselect())
      list = character(0)
      for (i in 1:n) { 
        temp = paste("Factor",i)
        list = c(list, temp)
      }
      list2 = setdiff(list,input$xaxis)
      selectInput("yaxis", "Choose Factor for plotting on Y axis",
                  list2, selected = "Factor 2")
    }else{
      temp = fname()
      temp2 = setdiff(temp,input$xaxis)
      selectInput("yaxis", "Choose Factor for plotting on Y axis",
                  temp2, selected = temp2[1])
      
    }
      }
  
})

f1 = reactive({
  if(is.null(fname())){
    f = input$xaxis
    s <- strsplit(f, "[^[:digit:]]")
    solution <- as.numeric(unlist(s))
    solution <- unique(solution[!is.na(solution)])
    return(solution)
  }else{
index <- match(input$xaxis,fname())
return(index)
}
  
})

f2 = reactive({
  if(is.null(fname())){
    f = input$yaxis
    s <- strsplit(f, "[^[:digit:]]")
    solution <- as.numeric(unlist(s))
    solution <- unique(solution[!is.na(solution)])
    return(solution)
  }else{
    index<-match(input$yaxis,fname()) 
    return(index)
    
   }
  
  
})



output$text1 <- renderText({    
  if (is.null(input$file)) { return(NULL) }
else {
      return(paste("Test of the hypothesis that",(fit())$factors,"factors are sufficient."))}
})

output$text2 <- renderText({
  if (is.null(input$file)) { return(NULL) }
 else{
     return(paste("The chi square statistic is",round((fit())$STATISTIC,3),"on",(fit())$dof," degrees of freedom.")) }
                                   })

output$text3 <- renderText({
  if (is.null(input$file)) { return(NULL) }
 else{
   return(paste("The p-value is",round((fit())$PVAL,3)))
 }
  })
#output$text4 <- renderText({ return(paste("Note - Optimal factors from are:",unlist((nS())[1])[3])) })

output$plot1 = renderPlot({
  if (is.null(input$file)) { return(NULL) }
  else{
  
    fafitfree <- fa(filtered_dataset(), nfactors = fselect(), rotate = "none")
    n_factors <- length(fafitfree$e.values)
    scree     <- data.frame(
      Factor_n =  as.factor(1:n_factors), 
      Eigenvalue = fafitfree$e.values)
    
    scree_plot1 = ggplot(scree, aes(x = Factor_n, y = Eigenvalue, group = 1)) + 
      geom_point() + geom_line() + geom_hline(yintercept = 1) +
      xlab("Number of factors") +
      ylab("Initial eigenvalue") +
      labs( title = "Scree Plot", 
            subtitle = "(Based on the unreduced correlation matrix)")
    
    scree_plot1  
  
  }
  })


output$plot20 = renderPlot({
  if (is.null(input$file)) { return(NULL) }
  else{
    
a = unclass((fit())$loadings)
grp = NULL
for (i in 1:nrow(a)){
  max = max(abs(a[i,]))
  temp0 =  which(abs(a[i,]) == max)
  temp = temp0[1]
  grp = c(grp,temp)
}
grp = matrix(c(grp,seq(1:length(grp))),,2)
rownames(grp) = colnames(filtered_dataset())

gr = vector("list", length = length(table(grp[,1])))
for (i in 1:length(table(grp[,1]))) {
  l1  = grp[(grp[,1] == as.numeric(names(table(grp[,1])[i]))),2]
  gr[[i]][1:length(l1)] = c(l1)   
}

qgraph(cor(filtered_dataset(), use= 'complete.obs'),layout="spring", groups = gr, labels=names(filtered_dataset()), label.scale=F, label.cex = 1, minimum=input$cutoffcorr)
}
})

output$plot2 = renderPlot({
  if (is.null(input$file)) { return(NULL) }
  else{
    
  a0 = (fit())$loadings
  a1 = a0
  for (i in 1:ncol(a1)){  a1[,i] = a0[,i]*(abs(a0[,i]) > input$cutoff)}
  k2 = f1()
  k3 = f2()
  
  factor.plot = function(a0, a1, k2, k3){
    
    load = a0[((a1[, k2] != 0)|(a1[, k3] != 0)), c(k2, k3)]
    
    par(col="black") #black lines in plots
    
    plot(load,type="p",pch=19,col="red", xlim=c(-1, 1), ylim=c(-1, 1),xlab=input$xaxis,ylab = input$yaxis) # set up plot
    
    abline(h=0);abline(v=0)#draw axes
    
    arrows(0,0, x1=load[,1], y1=load[,2], col="blaCK", lwd=1.5);
    
    text(load,labels = rownames(load),cex=1,pos=1)
    
  } # factor.plot() func ends
  
  factor.plot(a0, a1, k2, k3)
  }
})






output$plot3 = renderPlot({
  if (is.null(input$file)) { return(NULL) }
  else{
    
plot(x=(fit())$scores[,(f1())], y=(fit())$scores[,(f2())], type="p", pch=19, col="red",
    xlab = paste0(input$xaxis), ylab = paste0(input$yaxis))   # added this line in edit

text(x=(fit())$scores[,(f1())],y=(fit())$scores[,(f2())],labels=rownames(Dataset()), pos = 2, col="blue", cex=0.8)

abline(h=0); abline(v=0)
}
})

output$loadings <- renderDataTable({ 
  if (is.null(input$file)) { return(NULL) } else{
  # rownames((fit())$loadings) = colnames(Dataset())  # edit 2
  b2 <- unclass((fit())$loadings); rownames(b2) <- NULL;  
  b1 <- data.frame(colnames(filtered_dataset()), round(b2,2));
  names(b1)[1] <- "Variable"# [2:ncol(Dataset())];rownames(b1) <- colnames(Dataset())  # edit 2  
  #-------#
  if(is.null(fname())){return(b1)}
  else{
    names(b1)[c(-1)]<-fname()
    return(b1)
  }
  
  #-------#
  return(b1) # unclass((fit())$loadings)
  }
  })

mat = reactive({
  fact = (fit())
# SS.loadings= colSums(fact$loading*fact$loading)
# Proportion.Var = colSums(fact$loading*fact$loading)/dim(fact$loading)[1]
# Cumulative.Var= cumsum(colSums(fact$loading*fact$loading)/dim(fact$loading)[1])
# mat = rbind(SS.loadings,Proportion.Var,Cumulative.Var)
laodings = print(fact$loadings, digits=2, cutoff=.25, sort = TRUE)
# out = list(Stat = mat, loadings=laodings)
# return(laodings)

})

output$mat <- renderPrint({ 
  if (is.null(input$file)) { return(NULL) }
  else{ mat() }
  
  })

# uni = reactive({ 
# a = matrix(fit()$uniqueness,1,)
# colnames(a) = rownames(as.matrix(fit()$uniqueness))
# rownames(a) = "Uniqueness"
# return(a)
#  })

output$uni <- renderDataTable({ 
  if (is.null(input$file)) { return(NULL) }
  else{ 
    # n = ceiling(length(uni())/3)
    # matrix(uni(), ncol = 1)
    data.frame(Variable = rownames(as.matrix(fit()$uniqueness)), Uniqueness = round(fit()$uniqueness,2))
    }
},options = list(pageLength=10))


output$scores <- renderDataTable({     
  if (is.null(input$file)) { return(NULL) } else{
      # rownames((fit())$scores) = rownames(Dataset()) # edit 3 i made.
      # b0 <- (fit())$scores;   rownames(b0) <- rownames(Dataset()); 
      b2 <- unclass((fit())$scores); rownames(b2) <- NULL; 
      b0 <- data.frame(rownames(filtered_dataset()), round(b2,2)) # else ends    
      names(b0)[1] <- "Variable"
  if(is.null(fname())){return(b0)}
  else{
    names(b0)[c(-1)]<-fname()
    return(b0)
  }
  }
      # unclass((fit())$scores)
      #                             }
})  
  
# my edits 16-sept-2017 below
  output$downloadDataX <- downloadHandler(
    filename = function() { "Fac_scores.csv" },
    content = function(file) {
      write.csv(data.frame(rownames(filtered_dataset()), (fit())$scores), file, row.names = F)
   				 }
	)

output$downloadData <- downloadHandler(
  filename = function() { "Big_five_Survey_data.csv" },
  content = function(file) {
    write.csv(read.csv("data/Big_five_Survey_data.csv"), file, row.names=F, col.names=F)
  }
)
  
})

