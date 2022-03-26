
###############
# Script Info #
###############
# PURPOSE: Try to identify the ecosystem thresholds
# AUTHOR: Masoud A. Rostami
#


############
# PACKAGES #
############

# install.packages("gradientForest", repos="http://R-Forge.R-project.org")
# install.packages("extendedForest", repos="http://R-Forge.R-project.org")
# install.packages("reshape")
# install.packages("dplyr")
library(randomForest)
require(gradientForest)
library(reshape)
library(dplyr)
library(magrittr)
#############
# Load data #
#############


################ read the Harpacticoida families: Respond variables

Sp_foram=read.csv("foram_sps.csv")
dim(Sp_foram)

View(Sp_foram)
Sp_foram=Sp_foram[,2:27]



################ read the Environemnts variables: Predictor variables

abiotic_var=read.csv("foram_phys.csv")

dim(abiotic_var)
abiotic_var=abiotic_var[,2:13]
View(abiotic_var)


##############  Saving on object in RData format
save(Sp_foram, file = "Sp_foram.RData")

save(abiotic_var, file = "abiotic_var.RData")


# To load the data again
load("Sp_foram.RData")
load("abiotic_var.RData")





############################################### 2.2 gradientForest analysis

nSites <- dim(Sp_foram)[1]
nSites
nSpecs <- dim(Sp_foram)[2]
nSpecs
# set depth of conditional permutation
lev <- floor(log2(nSites * 0.368/2))
lev



######################

set.seed(12345)

gf <- gradientForest(cbind(abiotic_var, Sp_foram),predictor.vars = colnames(abiotic_var), response.vars = colnames(Sp_foram),ntree = 2000, transform = NULL, compact = T,
                     nbin = 201, maxLevel = lev, corr.threshold = 0.5, trace=T)

gf


names(importance(gf))

rank(-gf$overall.imp)

gf$result %>% data.frame()

gf$overall.imp2 %>% data.frame()
gf$overall.imp %>% data.frame()

Overall_performance=gf$result %>% data.frame()

gf$res



write.csv(Overall_performance,file = "Overall_performance.csv")

#########
# PLOTS #
#########



############################################ gradient Forest plots

plot(gf, plot.type = "O", col="darkgray", family = 'serif', cex.axis=0.75,cex.lab=2, bty="l")




##### The predictor gradient plots are best presented in order of importance; in this example the
## top 9 predictors are presented in 5 by 5 panels.


most_important <- names(importance(gf))[1:4]
par(mgp = c(2, 0.75, 0))



##The second plot is the splits density plot (plot.type="S"), which shows binned split importance and location on each gradient (spikes), kernel density of splits (black lines), of observations
#(red lines) and of splits standardised by observations density (blue lines). Each distribution integrates to predictor importance. These show where important changes in the abundance of
#multiple species are occurring along the gradient; they indicate a composition change rate.
#Many of the usual plot options can be set in the call



plot(gf, plot.type = "S", imp.vars = most_important,leg.posn = "topright", cex.legend = 1.2, cex.axis = 1,
     cex.lab = 1, line.ylab = 1, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))




####The third plot is the species cumulative plot (plot.type="C", show.overall=F), which
#for each species shows cumulative importance distributions of splits improvement scaled by
#R2 weighted importance, and standardised by density of observations. These show cumulative
#change in abundance of individual species, where changes occur on the gradient, and the species
#changing most on each gradient. Again many of the usual plot options can be set in the call; in
#this example the legend identifies the top 5 most responsive species for each predictor



plot(gf, plot.type = "C", imp.vars = most_important, show.overall = F, legend = T, leg.posn = "topleft",
     leg.nspecies = 5, cex.lab = 1.1, cex.legend = 1,cex.axis = 1.1, line.ylab = 1, 
     par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 1.5, 0.3, 0.5), omi = c(0,0.3, 0, 0)))





## The fourth plot is the predictor cumulative plot (plot.type="C", show.species=F), which
#for each predictor shows cumulative importance distributions of splits improvement scaled by
#R2 weighted importance, and standardised by density of observations, averaged over all species.
#These show cumulative change in overall composition of the community, where changes occur on
#the gradient. Again many of the usual plot options can be set in the call;
#in this example common.scale=T ensures that plots for all predictors have the same y-scale 
#as the most important predictor.



plot(gf, plot.type = "C", imp.vars = most_important,
     show.species = F, common.scale = T, cex.axis = 1,col="gray27",
     cex.lab = 1.2, line.ylab = 0.9, par.args = list(mgp = c(1.5,
                                                             0.5, 0), mar = c(2.5, 1, 1, 0.5), omi = c(0, 0.2, 0, 0)))









# https://www.r-graph-gallery.com/73-box-style-with-the-bty-function.html   # bty = 'L'show different boundary in graph



###### The fifth plot shows the R2 measure of the fit of the random forest model for each family, ordered in various ways

par(oma=c(3,3,3,3),mar=c(10,10,7,7),mfrow=c(1,1),bty = 'L')

plot(gf, plot.type = "P", show.names = T, horizontal = F,cex.axis = 0.8, cex.labels = 0.8,
     line = 2, col="darkgray",cex.main = 0.9, cex=1.9)


















################################## |Removing species lesss than 20





################ read the Harpacticoida families: Respond variables

Sp_foram=read.csv("foram_sps_few_species.csv")
dim(Sp_foram)

#View(Sp_foram)
Sp_foram=Sp_foram[,2:13]



################ read the Environemnts variables: Predictor variables

abiotic_var=read.csv("foram_phys.csv")

dim(abiotic_var)
abiotic_var=abiotic_var[,2:14]
#View(abiotic_var)


##############  Saving on object in RData format
save(Sp_foram, file = "Sp_foram.RData")

save(abiotic_var, file = "abiotic_var.RData")


# To load the data again
load("Sp_foram.RData")
load("abiotic_var.RData")



############################################### 2.2 gradientForest analysis

nSites <- dim(Sp_foram)[1]
nSites
nSpecs <- dim(Sp_foram)[2]
nSpecs
# set depth of conditional permutation
lev <- floor(log2(nSites * 0.368/2))
lev


######################

set.seed(12345)

gf <- gradientForest(cbind(abiotic_var, Sp_foram),predictor.vars = colnames(abiotic_var), response.vars = colnames(Sp_foram),ntree = 2000, transform = NULL, compact = T,
                     nbin = 201, maxLevel = lev, corr.threshold = 0.5, trace=F)

gf



############################################ gradient Forest plots

plot(gf, plot.type = "O", col="darkgray", family = 'serif', cex.axis=0.75,cex.lab=2, bty="l")



most_important <- names(importance(gf))[1:6]
par(mgp = c(2, 0.75, 0))



plot(gf, plot.type = "S", imp.vars = most_important,leg.posn = "topright", cex.legend = 1.2, cex.axis = 1,
     cex.lab = 1, line.ylab = 1, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))




plot(gf, plot.type = "C", imp.vars = most_important, show.overall = F, legend = T, leg.posn = "topleft",
     leg.nspecies = 6, cex.lab = 1.1, cex.legend = 1,cex.axis = 1.1, line.ylab = 1, 
     par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 1.5, 0.3, 0.5), omi = c(0,0.3, 0, 0)))




plot(gf, plot.type = "C", imp.vars = most_important,
     show.species = F, common.scale = T, cex.axis = 1.1,col="darkgray",
     cex.lab = 1.1, line.ylab = 1, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 1.5, 0.3, 0.5), omi = c(0, 0.3, 0, 0)))




###### The fifth plot shows the R2 measure of the fit of the random forest model for each family, ordered in various ways

par(oma=c(3,3,3,3),mar=c(10,10,7,7),mfrow=c(1,1),bty = 'L')

plot(gf, plot.type = "P", show.names = T, horizontal = F,cex.axis = 0.8, cex.labels = 0.8,
     line = 2, col="darkgray",cex.main = 0.9, cex=1.9)


plot(gf, plot.type = "P", show.names = T, horizontal = T,cex.axis = 0.8, cex.labels = 0.8,
     line = 2, col="darkgray",cex.main = 0.9, cex=1.9)





