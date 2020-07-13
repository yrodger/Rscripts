####################################
###       AMOVA and Fst          ###
####################################

setwd("//ad.monash.edu/home/User006/yrod0001/Desktop/PhD/Chapter_3_Tetraploids")

############ Analysis of Molecular Variance using poppr #################
library('adegenet')
library('dartR')
library('poppr')
library('ade4')

## Load genind object

load("tetraploidsONLY_genind.rdata")


#Make sure your genind object has strata defined - this will be useful in specifying the hierarchical structure of the AMOVA

pops <- read.csv("tet_strata.csv")

addStrata(tet.gi, pops)
strata(tet.gi) <- pops$pop

## Perform AMOVA 

tet.amova <- poppr.amova(tet.gi, ~pop)

tet.amova

amova.test <- randtest(tet.amova)

amova.test


############ Pairwise Fst using hierfstat #################

library(dartR)
library(adegenet)
library(hierfstat)

## Load genind object

load("tetraploidsONLY_genind.rdata")

#Convert genind object to hierfstat object
tet.hier <- genind2hierfstat(tet.gi)


#Have to edit tet.hier to change population names to numbers 
tet.hier$pop <- as.character(tet.hier$pop)
tet.hier$pop[tet.hier$pop == "Ararat"] <- 1
tet.hier$pop[tet.hier$pop == "Bannockburn"] <- 2
tet.hier$pop[tet.hier$pop == "Dobie_Bridge"] <- 3
tet.hier$pop[tet.hier$pop == "Glenelg"] <- 4
tet.hier$pop[tet.hier$pop == "Middle_Creek_subpop1"] <- 5
tet.hier$pop[tet.hier$pop == "Middle_Creek_subpop3"] <- 6
tet.hier$pop[tet.hier$pop == "Rokewood"] <- 7
tet.hier$pop[tet.hier$pop == "Wickcliffe"] <- 8
tet.hier$pop[tet.hier$pop == "Woodnaggerak"] <- 9
tet.hier$pop <- as.numeric(tet.hier$pop)
save(tet.hier, file = "tethierfstat.rdata")

#Compute the fst matrix
fst <- genet.dist(tet.hier, method = "WC84")
fst.m <- as.matrix(fst)
write.csv(fst.m, file = "tetraploidfst.csv")

## Bootstrapping
boot <- boot.ppfst(dat=rut.hier, nboot = 100, quant = c(0.025,0.975)) # p 0.05
boot99 <- boot.ppfst(dat=rut.hier, nboot = 100, quant = c(0.001,0.999)) # p 0.001
save(boot, file = "fstbootstrap.rdata")

# Save the Upper and Lower confidence limits
write.csv(boot99$ul, "fstbootstrapUI99.csv")
write.csv(boot99$ll, "fstbootstrapLI99.csv")

write.csv(boot$ul, "fstbootstrapUI.csv")
write.csv(boot$ll, "fstbootstrapLI.csv")

#If the range of LL - UL encompasses 0, the pairwise Fst is not significant.
