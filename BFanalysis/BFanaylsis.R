# live dangerously, get rid of pesky warnings
oldw <- getOption("warn")
options(warn = -1)

shhh <- suppressPackageStartupMessages # stops annoying warnings when loading libraries
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(MASS)
library(Matrix)
library(reshape2)
library(ape) # stats
library(vegan) # stats
library(RColorBrewer)
library(cocor)
library(DescTools)
library(reshape2)
library(grid)
library(ggplotify)



# read the csv data files into a dataframe
files = list.files(pattern="*.csv")
data = sapply(files, read.csv, simplify=FALSE) %>% bind_rows(.id = "id")

colnames(data)

# Select variables we need for analysis 
trial_vars<- c( "participant",
                "Colour_1", "Colour_2", "Colour1", "Colour2", 
                "similarity", "response_time",  "Ecc","screen_size_x",
                "screen_size_y","answer5",
                "viewerdistancecm", 'viewer_distance',"trialnumber","Ecc",'expName')

data <- subset(data, select = trial_vars)


# Create data frame for trials 
dftrials <- subset(data, !is.na(Colour1))

# Label participant number
dftrials$ID <- NA
subjectlist <- unique(dftrials$participant)
k= 0
for (participant in subjectlist){
  k = k + 1
  dftrials$ID[dftrials$participant == participant] <- k
}


# changing color values from RGB to hex for graphing purpose
dftrials$Colour1 <- as.character(dftrials$Colour1)
dftrials$Colour1 <- revalue(dftrials$Colour1, 
                            c(  "1" = '#FF0000',
                                "2" = '#FFAA00',
                                "3" = '#AAFF00',
                                "4" = '#00FF00',
                                "5" = '#00FFA9',
                                "6" = '#00A9FF',
                                "7" = '#0000FF',
                                "8" = '#AA00FF',
                                "9" = '#FF00AA'))
dftrials$Colour2 <- as.character(dftrials$Colour2)
dftrials$Colour2 <- revalue(dftrials$Colour2, 
                            c(  "1" = '#FF0000',
                                "2" = '#FFAA00',
                                "3" = '#AAFF00',
                                "4" = '#00FF00',
                                "5" = '#00FFA9',
                                "6" = '#00A9FF',
                                "7" = '#0000FF',
                                "8" = '#AA00FF',
                                "9" = '#FF00AA'))

# colors for the labels
# red, orange, yellow, green, cyan, cyan-blue, blue, purple, pink
colors <- c('#FF0000','#FFAA00','#AAFF00','#00FF00','#00FFA9','#00A9FF','#0000FF','#AA00FF','#FF00AA')
# can change the way the plot line up
# red, pink, orange, purple, yellow, blue, green, cyan-blue, cyan
#colors <- c('#FF0000','#FF00AA','#FFAA00','#AA00FF','#AAFF00','#0000FF','#00FF00','#00A9FF','#00FFA9')
abcolors <- sort(colors) # this was messing up the asymmetry plot, maybe useful for some other stuff


# factor the dataframes for the plot function
dissimdata2 <- function(dftrials, colors){
  
  # refactor the levels so they can be plotted properly later if need be
  dftrials$Colour1 <- with(dftrials, factor(Colour1, levels = colors))
  dftrials$Colour2 <- with(dftrials, factor(Colour2, levels = colors))
  
  return(dftrials)
}


df2mat_asymmetry_temporal <- function(dftrials){
  
  datatemp <- dissimdata2(dftrials, colors)
  
  # aggregate over the remaining columns of interest
  nmdsdata <- aggregate(datatemp, by = list(datatemp$Colour1, datatemp$Colour2),FUN=mean)
  nmdsdata$Colour1 <- nmdsdata$Group.1
  nmdsdata$Colour2 <- nmdsdata$Group.2
  
  nmdsdata = subset(nmdsdata, select = c("Colour1","Colour2","similarity"))  # get rid of unnecessary columns
  nmdsmatrix <- spread(nmdsdata, Colour1, similarity) # convert the dataframe to a matrix
  nmdsmatrix <- data.matrix(nmdsmatrix) # change first column from colour to number(just some label stuff) 
  nmdsmatrix <- nmdsmatrix[,-1] # get rid of the labels in the first column, it messes up the code
  
  matdf<-  as.data.frame(nmdsmatrix - t(nmdsmatrix)) # calculate the asymmetry
  matdf$colorset <- c(colors) # adding additional column "colorset"
  num_colors <- length(colors)
  matdf <- matdf %>% gather(othercolor,asymmetry ,1:num_colors) # convert the matrix back to the data frame which has the 
  # column "colortset", "othercolor", "asymmetry"
  return(matdf)
}


# Convert dataframe into a list of asymmtery matrix for each subject
dissimgroup2matlist <- function(dftrials){
  subjectlist <- sort(unique(dftrials$ID)) # obtain a list of all the subjects
  mat.list <- list()
  k = 1
  for(ID in subjectlist){
    subjectdf <- dftrials[which(dftrials$ID==ID),]
    subject.mat <- df2mat_asymmetry_temporal(subjectdf)
    mat.list[[k]] <- subject.mat
    k = k + 1
  }
  return(mat.list)
}

mat.list <- dissimgroup2matlist(dftrials)

# Return a list of asymmtery values for each cell 
vals.list.fun <- function(mat.list){
  n.cells <- 81
  n.subjects <- length(mat.list)
  
  k=1
  vals.list <- list()
  
  for(cell in 1:n.cells){
    v <- vector()
    for(subject in 1:n.subjects){
      v <- c(v,mat.list[[subject]][[3]][cell])
    }
    vals.list[[cell]] <- v
  }
  return(vals.list)   
}

vals.list <-vals.list.fun(mat.list)

# T-test 
t.test.mu.zero <- function(cell, out='p-value'){
  test <- t.test(cell,
                 alternative = c("two.sided"),
                 mu = 0,
                 conf.level = 0.95)
  t <- test$statistic
  if(out=='t'){
    return(t)
  } 
  if(out=='p-value'){
    return(p-value)
  } 
}

# p values 
t.test.mu.zero.p <- function(cell){
  test <- t.test(cell,
                 alternative = c("two.sided"),
                 mu = 0,
                 conf.level = 0.95)
  p <- test$p.value
  return(p)
}

# Create a matrix of p-values
p.mat <- function(mat.list){
  mat.vals <- vals.list
  out.mat <- matrix(nrow=9,ncol=9)
  
  for(cell in 1:length(mat.vals)){
    p <- t.test.mu.zero.p(mat.vals[[cell]])
    out.mat[cell] <- p
  }
  print(out.mat)
  return(out.mat)
}

sigmat <- p.mat(mat.list)

# Create a matrix of t-scores
t.mat <- function(mat.list){
  mat.vals <- vals.list
  out.mat <- matrix(nrow=9,ncol=9)
  
  for(cell in 1:length(mat.vals)){
    t <- t.test.mu.zero(mat.vals[[cell]],'t')
    out.mat[cell] <- t
  }
  print(out.mat)
  return(out.mat)
}

out.mat <- t.mat(mat.list)

# Create dataframe of t-scores 
# Select half of the matrix
out.mat[upper.tri(out.mat)] <- NA
# Convert matrix of t scores into a dataframe
dfmat <- as.data.frame((out.mat))
colnames(dfmat) <- c(abcolors) 
dfmat$colorset <- c(abcolors) # adding additional column "colorset"
num_colors <- length(colors)
dfmat <- dfmat %>% gather(othercolor,t ,1:num_colors) # convert the matrix back to the data frame which has the 
# column "colortset", "othercolor", "asymmetry"
# Remove all NA values
dfmat <- dfmat[complete.cases(dfmat), ] #dataframe of t-scores

dfmat$t <- abs(dfmat$t) # Take the absolute asymmetry 


## HYPOTHESIS 1 -  Participants will provide a different similarity 
#judgment for a colour depending on the order of the two colours


# Models I have tried

model.test1 <- stan_glm(x=dfmat$t, family=Gamma()) # Does not work, incorrect formula 

# Gamma distribution 
posterior <- distribution_gamma(x = dfmat$t, n=36, shape=1)
bayesfactor_parameters(posterior)

# Using Ariel function 


max.DSR = 7 

offset.beta <- function(vector, max.DSR){
  offset <- (vector * (length(vector)-1) + 0.5) / length(vector)
  return(offset/max.DSR)
}


h1a <- function(data,summary=FALSE){
  data <- offset.beta(data,max.DSR)
  df <- as.data.frame(data)
  
  model.test <- stan_betareg(
    formula = data~1, # is this the correct model? I saw that beta regressions like the fit we want would be
    # data ~ 2? I think
    data = df,
    iter=10000, # helps to stabilise the estimates
    refresh=0
  )
  
  if(summary){
    print(summary(model.test))
  }else{
    print(model.test)
  }
  HPD <- posterior_interval(model.test, prob = 0.95)
  print(HPD)
  BF <- bayesfactor_parameters(model.test)
  print(BF)
  print(plot(BF, show_intercept=TRUE))
  
  # Converted back to normalised trace
  print("")
  print("Converted back to Normalised Trace")
  intercept <- model.test[[1]][1]
  print(paste('Intercept:',round(undo.logit(intercept,max.DSR),2)))
  print(paste('HPD:',round(undo.logit(HPD[1,],max.DSR),2)))
  return(model.test)
}


h1a(dfmat$t,summary=TRUE) # Does not print the phi graph or calcualte the phi BF, unsure why 

