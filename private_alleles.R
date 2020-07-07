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






library(poppr)

gi_private <- gl2gi(glparegion)

pral <- private_alleles(gi_private, form = alleles ~ ., report = "data.frame", level = "population", count.alleles = TRUE, drop = FALSE)

write.csv(pral, file = "private allele count REGIONS.csv")

Privates <- private_alleles(gi_private, form= alleles ~ . , report="table", count=TRUE)

write.csv(Privates, file = "private allele Joe table.csv")

###>20% allele threshold in both populations
##Start with your imported dart data as genlight object. Calcualte the frequency of each allele at each locus:

frequencies <- gl.alf(glparegion)

##Create a column in this new object with the name of each locus
frequencies$row.names <- glparegion$loc.names

## Subset this frequency list, leaving only loci with with <0.2 frequency for either allele
to_remove <- subset(frequencies, alf1 < 0.2 | alf2 <0.2)

##Remove these < 0.2 frequency loci from your genlight object

filtered_glparegion <- gl.drop.loc(glparegion, to_remove$row.names, v = 2)

nLoc(filtered_glparegion) #390 loci

paregion20 <- gl.report.pa.pop(filtered_glparegion)

gipraregion20 <- gl2gi(filtered_glparegion)

pa_filtered <- private_alleles(gipraregion20, form= alleles ~ . , report="table", count=TRUE)

write.csv(pa_filtered, file = "private allele filetered 20 Joe table.csv")
