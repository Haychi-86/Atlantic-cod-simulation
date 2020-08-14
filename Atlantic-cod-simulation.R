# Example code of simulations presented in Wootton et al. (2020) where population growth models are adapted from Wilson et. al. (2018). 
# Life history parameters are taken from published estimates from the Greenland Atlantic cod (Gadus morhua) stock to give population growth curves maximum realism. 
# Population growth is first modeled under baseline parameters and then replicate populations are modeled where the coefficient of variation around size at age is allowed 
# to vary. We assume this variation is constant throughout ontogeny at each level of CV. Modifications to this code can be made to simulate variation in other growth
# parameters or to simulate fishery selection as described in the manuscript.

# First, clear global environment and load required packages

rm(list=ls())
library(ggplot2)
library(dplyr)
#library(Ryacas) # this rearranges equations if you need to back calculate parameters

# Greenland cod stock step 1. Generate growth length-at-age data with variation built in
# *modify the below code if you wish to model the Skagerrak cod stock or any other stock with differing parameters

# Generate a curve that we will use as an average size for each age class 

ages <- 1:25   # Create an integer sequence of ages

Tmat <- 8      # Age of maturity #rounded down from 8.5 from fishbase (Jonsson 1959)

linf <- 1540   # fishbase (Ratz and Stein 1999)
vbk <-  0.06   # fishbase (Ratz and Stein 1999)
t0 <-  -2.88   # fishbase (Ratz and Stein 1999)

# To back-calculate g, h and t1 from equations from the original Wilson (2018) script, substitute parameters in the below equation and solve for g, h and t1
# This allows us to generate population growth curves from data derived from real populations
# linf <- 3 * h/g                                        # Wilson (2018) conversion to VBGF L-infinity
# vbk <- log(1 + g/3)                                    # Wilson (2018) conversion to VBGF kappa
# t0 <- Tmat + log(1 - g * (Tmat - t1)/3)/log(1 + g/3)   # Wilson (2018) conversion to VBGF t0 

# vbk <- log(1 + g/3) 
g <- (exp(vbk)-1)*3   # Proportion of energy in adult phase allocated to reproduction per year

# linf <- 3 * h/g 
h <- (linf/3)*g       # Somatic growth in millimeters per year

# t0 <- Tmat + log(1 - g * (Tmat - t1)/3)/log(1 + g/3)
t1 <- Tmat - (3/g) * (1 - (1 + g/3) ^ (t0 - Tmat))     # Age when size=0 for the juvenile phase

lena_phase1 <- h * (ages - t1)  # Length-at-age for phase 1
lena_phase2 <- linf * (1 - exp(-vbk * (ages - t0)))  # Length-at-age for phase 2
biphasic <- ifelse(ages < Tmat, lena_phase1, lena_phase2)  # if-else statement for which phase a fish is allocating surplus energy

# Plot the population average growth-curve

plot(ages, lena_phase1, ylab = "Size", xlab = "Age")
lines(ages, lena_phase2)
lines(ages, biphasic)
abline(v = Tmat, col = "red", lty = 2)

Cv <- 0.155672 # Cv for Greenland cod - calculated as a proxy from the ICES data on Atlantic cod from the North Sea. Data are analysed by age class and area and then averaged. Cv is assumed to be constant through ages in the model.

Size <- NULL # Generate a vector of sizes

No.inAgeClasses <- 1000 # How many samples of each age class do we want to create? (for survival to be applied to)

for (i in 1:length(ages)) Size <- c(Size, rnorm(No.inAgeClasses, mean=biphasic[i], sd=Cv*biphasic[i])) # A for loop that generates size-distribution with variation

Data <- data.frame(cbind(rep(ages, each = No.inAgeClasses), Size)) # Create a data-frame
colnames(Data) <- c("Age", "Size")

# Plot the results

plot(Data$Age, Data$Size, pch = 21, xlab = "Age (yrs)", 
     ylab = "Length (mm)")

# Combine the vector of sizes with discreet probabilities of survival at a given age to simulate realisitic size-age strutures

surv <- rep(NA, length(ages))  # Create an empty vector
surv[1] <- 1
for (j in 2:max(ages)) {
  surv[j] <- surv[j - 1] * exp(log(((g/1.18) - 1)/-1))
}

# Plot the survival data

plot(Data$Age, Data$Size, pch = 21, xlab = "Age (yrs)", 
     ylab = "Length (mm)")
lines(ages, surv*2000)

SurvivalData <- as.data.frame(cbind(Age=ages, surv)) 

SurvivalChanceData <- merge(Data, SurvivalData, by= "Age")

# Apply survival

SurvivalChanceData$SurvSample <- NA
for ( i in seq_len(nrow(SurvivalChanceData))){
  SurvivalChanceData$SurvSample[i] <- rbinom(1, size=1, prob =  SurvivalChanceData$surv[i])
}

# Remove 'dead' individuals

SurvivalChanceData$SurvSample[SurvivalChanceData$SurvSample < 1] <- NA

SurvivalChanceData <- SurvivalChanceData[!is.na(SurvivalChanceData$SurvSample), ]

# Sample from survived 

No.Samples <- 250 # How many datapoints do we want in our dataset that incorporates survival? - lets sample 250 to try to ensure that it is Cv affecting model performance not sample size (Honsey et al. 2017)

SurvivalChanceData <- SurvivalChanceData[sample(nrow(SurvivalChanceData), No.Samples), ]

plot(SurvivalChanceData$Age, SurvivalChanceData$Size, pch = 21, xlab = "Age (yrs)", 
     ylab = "Length (mm)")
lines(ages, surv*1400)

# Now run code that allows us to vary CV across a number of populations and stores the resulting populations in a dataframe

nPop <- 5 # How many populations do we want to create and analyse later?

Cv <- seq(0.05, 0.45, length = nPop) # What distribution of variation around size at age do we want to create?

# A for-loop that generates size-at-age distrubutions for the number of populaitons specified above

PopData <- NULL
for (j in 1:nPop) {
  Size <- NULL
  # Create size distribution
  
  for (i in 1:length(ages)) Size <- c(Size, rnorm(No.inAgeClasses, mean=biphasic[i], sd=Cv[j]*biphasic[i]))
  
  # Merge sizes with ages
  
  Data <- data.frame(cbind(rep(ages, each = No.inAgeClasses), Size, rep(j, length(ages)), rep(Cv[j], length(ages))), row.names = NULL)
  colnames(Data) <- c("Age", "Size", "No.Pop", "PopCv")
  
  SurvivalData <- as.data.frame(cbind(Age=ages, surv)) 
  
  SurvivalChanceData <- merge(Data, SurvivalData, by= "Age")
  
  # Apply survival function to generated size-at-age data
  
  SurvivalChanceData$SurvSample <- NA
  for ( i in seq_len(nrow(SurvivalChanceData))){
    SurvivalChanceData$SurvSample[i] <- rbinom(1, size=1, prob =  SurvivalChanceData$surv[i])
  }
  
  # Remove 'dead' individuals
  
  SurvivalChanceData$SurvSample[SurvivalChanceData$SurvSample < 1] <- NA
  
  SurvivalChanceData <- SurvivalChanceData[!is.na(SurvivalChanceData$SurvSample), ]
  
  SurvivalChanceData <- select(SurvivalChanceData,-c(SurvSample)) #clean out dataset with this line of code
  
  # Remove sizes less than zero  
  
  SurvivalChanceData$Size[SurvivalChanceData$Size < 0] <- NA 
  
  SurvivalChanceData <- SurvivalChanceData[!is.na(SurvivalChanceData$Size), ]
  
  # Take a random sample from size-at-age data (with survival built in) and repeat 10 times to create 10 random populations for each level of Cv drawn from the same pool of survival data
  
  RepeatSampleData <-SurvivalChanceData
  
  for(i in 1:10){RepeatSampleDataIteration <- RepeatSampleData[sample(nrow(RepeatSampleData), No.Samples), ] 
  RepeatSampleDataIteration$RepNo. <- rep(i, length = nrow(RepeatSampleDataIteration))
  PopData <- rbind(PopData, RepeatSampleDataIteration)
  }
}

tapply(X=round(PopData$Size), INDEX=list(round(PopData$Size)), FUN=length) # Check if the loop worked

# Lets plot our populations

SimulatedCvLengthAgePopPlot <- ggplot(data= PopData, aes(x= Age, y= Size)) + 
  geom_point() +  facet_wrap(No.Pop~RepNo.)  
SimulatedCvLengthAgePopPlot

# Lets select data from the first mat-age input and plot it, which allows us to look more closely at our repeated sampling of survival data at the same level of CV

SimulatedCvLengthAgePopPlotPop1 <- ggplot(data= PopData[with(PopData, No.Pop == "1"),], aes(x= Age, y= Size)) + 
  geom_point() +  facet_wrap(RepNo.~No.Pop)  
SimulatedCvLengthAgePopPlotPop1