working_dir = directory = 'C://directory/data/'
setwd(working_dir)
      
install.packages("qtl")
##Remember to install libraries first

library(qtl)

#A tutorial for the R/qtl program is available here: http://www.rqtl.org/tutorials/rqtltour.pdf
#In our lab, we will not go through all steps of checking for genotyping errors that the tutorial does. You may find the tutorial helpful in your own work.

# R/qtl analysis of 500 RILs mean freezing tolerance (survival after freezing) from Oakley et al. (2014)

# The researchers first use quantile-normalized data then re-fit the model with the raw data to get effect size estimates on the raw scale
# Read in and specify the format of the data files. Here, meangeno.csv is the genotypic data and meanpheno.csv is phenotypic data. In the phenotypic data "QN_Mean_survival" represents the quantile normalized data, whereas Mean_survival is the raw data.
# a = Italian genotype, b = Swedish genotype

setwd("~/Documents/GENE 8150/2018/lab10/")
sweden_italy_cross<-read.cross("csvs",dir="",
  genotypes=c("a","b"), 
  genfile="meangeno.csv", 
  phefile= "meanpheno.csv", 
  na.strings=c("NA","-"))
  
# Set the cross type to RIL (recombinant inbred line, generated from selfing)
class(sweden_italy_cross)[1] <-"riself"

# Calculate genotype probabilities, conditional on the available marker data, and calcualte the linkage map from genotypic data. Genotypic probabilities are needed for QTL mapping functions. The argument step increases the step size (in cM) at which probabilities are calculated and determines the step size at which later LOD scores are calculated.
sweden_italy.calc<-calc.genoprob(sweden_italy_cross, 
  step=2, error.prob=.0001, map.function="kosambi",
  stepwidth="max", off.end=0)
  
 
 #The plot function gives you some interesting plots, including the extent of missing data, the genetic map ("Chromosome" should really be labeled "Linkage Groups"), and hisograms of the two phenotypes (the first is quantile normalized survival and the second is the raw data - you'll see that zero inflation rears its ugly head here, too)
  plot(sweden_italy.calc)
  
  ##If you want to plot the map with the marker names listed, you can use this
  plotMap(sweden_italy.calc, show.marker.names=TRUE)
  
  ##If you just want to plot the phenotypic histogram
  
  plotPheno(sweden_italy.calc, pheno.col="Mean_survival")
  
  
 ##This step is not necessary, but it is another way of getting the map
 newmap<-est.map(sweden_italy_cross, verbose=TRUE)
 plotMap(newmap)

  
  
  
##Here, we will calculate recombintion frequencies across markers and plot them. Note that red and yellow corersponds to a small recombination fraction, while blue is the reverse. Do thse patterns make sense to you? 

sweden_italy.calc <- est.rf(sweden_italy.calc)
plotRF(sweden_italy.calc)
#If you were just interested in linkage groupd 1
plotRF(sweden_italy.calc, chr=1)

#Now that we have the genetic map, we can begin with QTL mapping. Please look at the distribution patterns of freezing tolerance (the phenotype) - you can evaluate the histograms you plotted with plotPheno(sweden_italy.calc, pheno.col="Mean_survival") and compare them to Figure 1. It departs from standard assumptions of normality of interval mapping. How do the authors approach this departure in their paper?
  
# Do a one way scan, which performs a single QTL genome scan via interval mapping (IM) using maximum likelihood via the Haley-Knott regression model. You could change "hk" to "em" to do standard interval mapping (Lander and Botstein 1989) using a normal model and the EM algorithm

scan1_Mean_survival <- scanone(sweden_italy.calc, method="hk", pheno.col="Mean_survival")

#You could do the same IM scan with quantile normalized data

scan2_QN_Mean_survival <- scanone(sweden_italy.calc, method="hk", pheno.col="QN_Mean_survival")


#The permutations here can be used to determine the signifcance threshold. This may take 30 seconds or so first for mean survival and then for QN_Mean_Survival
perm_scan1_Mean_survival <- scanone(sweden_italy.calc, method="hk",n.perm=10000, pheno.col="Mean_survival")
perm_scan2_QN_Mean_survival <- scanone(sweden_italy.calc, method="hk",n.perm=10000, pheno.col="QN_Mean_survival")

#These statements give you the LOD significance thresholds from the permutations
np_Mean_survival <- summary(perm_scan1_Mean_survival)
np_Mean_survival

np_QN_survival<-summary(perm_scan2_QN_Mean_survival)
np_QN_survival

# Summarize one way scan (IM) with permutations summary. This summary statement displays the maximum LOD score on each linkage group for which LOD exceeds the threshold found via permutation (which was 2.6 in my run above)

summary(scan1_Mean_survival, perms=perm_scan1_Mean_survival, alpha=0.05)

#And now for the quantile normalized data. Do you notice any differences?

summary(scan2_QN_Mean_survival, perms= perm_scan2_QN_Mean_survival, alpha=0.05)

# Plot one way scan (interval mapping) with LOD threshold from permutations with raw data

plot(scan1_Mean_survival, main="scan1_Mean_survival")
abline(h=np_Mean_survival[1], col="red")

# Let's add on the quantile normalized QTLs for comparison (raw data = black, quantile normalized = blue; permutation threshold for raw=red and for normalized = orange). Do the results seem robust to data transformation?

plot(scan1_Mean_survival, scan2_QN_Mean_survival, col=c("black","blue"), lwd=3, lty=1:2)
abline(h=np_Mean_survival[1], col="red", lty=3)
abline(h= np_QN_survival[1], col="orange", lty=4)

# You might want to zoom in on a chromosome. Let's try chromosome 2 because of the small effect QTL that may be there in the normalized data
plot(scan1_Mean_survival, scan2_QN_Mean_survival, col=c("black","blue"), lwd=3, lty=1:2, chr=c("2"))
abline(h=np_Mean_survival[1], col="red", lty=3)
abline(h= np_QN_survival[1], col="orange", lty=4)



#Let's try composite interval mapping. Here, we are using 4 markers as covariates. What happens if you change that number? The perm_scan3_cim step takes a couple of minutes.

est.cim<-cim(sweden_italy.calc, n.marcovar=4, window=10, pheno.col="Mean_survival")
perm_scan3_cim <- cim(sweden_italy.calc, n.marcovar=4, window=10, pheno.col="Mean_survival",n.perm=1000)



#This code will generate a figure to compare interval mapping with composite interval mapping. In this case, interval mapping results are black solid lines, CIM results are blue dashed lines and the signficiance threshold is the black dashed line (interval mapping) and blue dashed line (composite interval mapping). How do the two results compare?

plot(scan1_Mean_survival, est.cim, col=c("black","blue"), lwd=3, lty=1:2)
add.threshold(scan1_Mean_survival,perm= perm_scan1_Mean_survival, col="black", lty=3)
add.threshold(est.cim,perm= perm_scan3_cim, col="blue", lty=3)


# Do a two way scan. This is a two-dimensional genome scan to identify interacting QTL. For every pair of positions, it calculates a LOD score for the full model (two QTL plus interactions) and a LOD score or the additive model (two QTL but no interaction). 

scan2_QN_Mean_survival <- scantwo(sweden_italy.calc, 
 pheno.col="QN_Mean_survival", model="normal",  method="hk", addcovar=NULL, intcovar=NULL, weights=NULL,
  clean.output=FALSE, clean.nmar=FALSE, clean.distance=FALSE,
  incl.markers=TRUE,  assumeCondIndep=FALSE,  verbose=FALSE)


# Perform permutations for two way scan. The researchers used 10,000 permutations distributed across 20 computational clusters (the relevant text in the manuscript is : "We ran 10,000 permutations to determine LOD thresholds (experimentwise alpha=0.05)" on p. 4307 under QTL mapping). However, that would take a *really* long time on our desktops, so we're going to be rather simplistic and use 500 permutations instead. On my computer, this took 8 mintues. This would be a good time to review the MacKay et al. paper, reread Oakley, and/or start working on the homework questions

QN_Mean_survival_perm<-scantwo(sweden_italy.calc, 
 pheno.col="QN_Mean_survival", model="normal",  method="hk", n.perm=500, 
  addcovar=NULL, intcovar=NULL, weights=NULL,
  clean.output=FALSE, clean.nmar=FALSE, clean.distance=FALSE,
  incl.markers=TRUE,  assumeCondIndep=FALSE,  verbose=FALSE)




###starthere
# Calculate LOD threshold penalties at alpha=0.05 for additive QTL
QN_Mean_survival_pen1.05<-calc.penalties(QN_Mean_survival_perm, alpha=.05)
summary(QN_Mean_survival_pen1.05)

# Perform the stepwise forward-backward regression procedure on quantile normalized data at alpha=0.05 for additive and epistatic QTL, max.qtl is the maximum number of QTL allowed in the model, which should be set to be the actual number of QTL. The relevant text in the paper is under the QTL mapping subheading (p 4307): "... and employed automated stepwise model selection, scanning for additive and epistatic QTL at each step"

EPI_HEAVY_stepwise_QN_Mean_survival<- stepwiseqtl(sweden_italy.calc, method="hk", model="normal", 
pheno.col="QN_Mean_survival", penalties=QN_Mean_survival_pen1.05[1:2], max.qtl=15, covar=NULL, 
  scan.pairs=FALSE, additive.only=F,keeplodprofile=TRUE, keeptrace=TRUE, 
  refine.locations=TRUE, verbose=TRUE)

# Summarize the stepwise model and plot the LOD profiles of significant QTL. How do these results differ from interval mapping (scan one) and composite interval mapping?

summary(EPI_HEAVY_stepwise_QN_Mean_survival)
plotLodProfile(EPI_HEAVY_stepwise_QN_Mean_survival, main="EPI_HEAVY_stepwise_QN_Mean_survival_LodProfile")

# Fit and refine the stepwise model on the quantile normalized data
EPI_HEAVY_fitqtl_QN_Mean_survival<-fitqtl(sweden_italy.calc, pheno.col="QN_Mean_survival", 
  method="hk", model="normal", covar=NULL, 
  qtl=EPI_HEAVY_stepwise_QN_Mean_survival, formula=formula(EPI_HEAVY_stepwise_QN_Mean_survival), 
  get.ests=TRUE, dropone=TRUE)
  summary(EPI_HEAVY_fitqtl_QN_Mean_survival)

# Re-fit this model using the raw data
EPI_HEAVY_fitqtl_Mean_survival<-fitqtl(sweden_italy.calc, pheno.col="Mean_survival", 
  method="hk", model="normal", covar=NULL, 
  qtl=EPI_HEAVY_stepwise_QN_Mean_survival, formula=formula(EPI_HEAVY_stepwise_QN_Mean_survival), 
  get.ests=TRUE, dropone=TRUE)
    summary(EPI_HEAVY_fitqtl_Mean_survival)


# Are there major differences in QTL mapped with quantile normalized data vs. raw data? How do your estimates of effect sizes compare with what Oakley et al present in their paper?




###############################################################
# Get ANOVA tables from fitqtl models
# Use the LOD scores and PVE (% variance explained) from the quantile-normalized data
# Use the allelic effect sizes from the raw data (multiply these by 2 to get genotypic effects (RILs))
# Layout changes in R/qtl for n QTL >1 vs. n QTL = 1
# if statements extract the bits you want regardless of number of qtl found

if(length(EPI_HEAVY_stepwise_QN_Mean_survival$name)>1) {summary(EPI_HEAVY_fitqtl_QN_Mean_survival)[1:2]} else {summary(EPI_HEAVY_fitqtl_QN_Mean_survival)[1]}

if(length(EPI_HEAVY_stepwise_QN_Mean_survival$name)>1) {summary(EPI_HEAVY_fitqtl_Mean_survival)[3]} else {summary(EPI_HEAVY_fitqtl_Mean_survival)[2]}


# For loop to print Bayesian 95% credible intervals for each QTL in model

for (i in 1:length(EPI_HEAVY_stepwise_QN_Mean_survival$name)) {
	
	print(bayesint(EPI_HEAVY_stepwise_QN_Mean_survival, prob=0.95, qtl.index=i, expandtomarkers=F))
	
}



# Code to extract LOD profiles
# This example extracts the LOD profile for the second QTL in the model
lp.EPI_HEAVY_stepwise_QN_Mean_survival<-attr(EPI_HEAVY_stepwise_QN_Mean_survival,"lodprofile") 
lp.EPI_HEAVY_stepwise_QN_Mean_survival[2]
