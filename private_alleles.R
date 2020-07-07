###############################
### Private allele analysis ###
###############################

### Code written by Yael Rodger ###

### Option 1: using dartR 

## If you have DArT data from Diversity Arrays Technology you will be able to import it as a genlight object
## Alternatively, if you have a genind object (adegenet) you can convert it to a genlight object using gi2gl()

library(dartR)
load("gl12filtered6.rdata) # This is my dataset: the genlight object gl12filtered6

#To calculate the pairwise private alleles
pairwise.pa <- gl.report.pa.pop(gl12filtered6)

#To calculate the total number of private alleles possessed by one or other population, as a matrix
pa.dist <- gl.dist.pop(gl12filtered6, method = "pa")
write.csv(pa.dist, "private allele population matrix.csv")

#Private allele analysis between regions

#First merge your populations into regions
glparegion <- gl.merge.pop(gl12filtered6, old=c("Bredbo", "Michelago", "Captains Flat", "Capital Hill",
                                                "Red Hill", "Stirling Ridge", "Crace", "Queanbeyan",
                                                "Majura", "Gundary"), new="ACTNSW")
glparegion <- gl.merge.pop(glparegion, old=c("Truganina", "St Albans"), new="Victoria")


### Option 2: using poppr

## Data format here is a genind object (adegenet)

library(poppr)

pral <- private_alleles(gi_private, form= alleles ~ . , report="table", count=TRUE) #The table format is more helpful

write.csv(pral, file = "private allele count regions.csv")
#With this you can identify the loci that have private alleles and how many. This allows you to calculate the frequency of alleles in each population/region (I did this in Excel).
