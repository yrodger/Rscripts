
#calculating pop stats of genetic diversity in adegenet (using full dataset)
```{r pressure, echo=FALSE}
library(adegenet)
library(pegas)
library(hierfstat)
dataPop <-read.fstat.data("Rutidosis_fulldataset12pops_FSTAT.dat")
stats <- basic.stats(data = dataPop)
save(stats, file = "stats_hierfstat.rdata")
load("stats_hierfstat.rdata")
N <- stats$n.ind.samp #calculating 
write.csv(N, file = "nindsamp_hierfstat.csv")
Ho <-stats$Ho
Hs <-stats$Hs
Fis <-stats$Fis
perloc <- stats$perloc

colMeans(Ho, na.rm=TRUE) #prints a table of observed Heterozygosity
Hosd <- apply(Ho, 2, sd)
colMeans(Hs, na.rm=TRUE) #prints a table of gene diversity
colMeans(Fis, na.rm=TRUE) #value derived from the H0 and Hs

```
#Calculates per locus allelic richness
```{r pressure, echo=FALSE}
AR_full_dataset <- allelic.richness(dataPop, diploid=TRUE)  #scales AR to the smallest pop size in the dataset
write.csv(AR_full_dataset, "AR_full_dataset_YR.csv")

##Analysis of genetic diversity in diveRsity (using reduced dataset - no missing data or else no Ho calculated)
```{r pressure, echo=FALSE}
#requires conversion into GENEPOP file format; done via PGDSpider from Structure format
#PGDSpider had issues with pop names that included spaces. GENEpop file was given a name in the first lane to comply (see below)
#NOTE: Genepop input file must have a title for diveRsity to run, will output error message otherwise:
#"Error in all_alleles[[i]] : subscript out of bounds".

library(diveRsity)

divBasic(infile = "Rutidosis_fulldataset12popsGENEPOP.txt", outfile = "Rutidosis12_diveRsity", gp = 3, bootstraps = 1, HWEexact = FALSE)
