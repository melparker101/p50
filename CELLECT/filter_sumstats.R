# Use a high number of cores (e.g. 10)

# Load in libraries
library(data.table)
library(stringr)

# Read in infertility sumstats
filename <- 'female_infertility_analysis1_UKBB_Finngen_EstBB_GandH_noMACfilter_March20231.out'
Infertility1_F_Finngen <- fread(filename)

# Read in hormone sumstats
hormone_sumstats <- list()
for(i in seq_along(filelist)) {
  hormone_sumstats[] <- name
}


filelist <- list.files(pattern='.txt')
# hormone_sumstats <- list()
for(i in seq_along(filelist)){
  name <- str_replace(filelist[i],"_filtered.txt","")
  sumstats <- fread(filelist[i])
  assign(name, sumstats)
  hormone_sumstats <- append(hormone_sumstats,name)
}

# Clean up
rm(sumstats)

# Make a MarkerName column for hormone sumstats
# Columns in hormome sumstats include: CHROM, GENPOS, Allele1, Allele2
# touppper() is probably not necessary but leave in for robustness
for(file in seq_along(filelist)) {
  sumstats <- str_replace(filelist[file],"_filtered.txt","")
  n <- nrow(sumstats)
  for (i in 1:n){
    	a1 <- toupper(sumstats$Allele1[i])
	    a2 <- toupper(sumstats$Allele2[i])
      chrom <- sumstats$CHROM[i]
      position <- sumstats$GENPOS[i]
      alleles <- paste(sort(c(a1,a2)), collapse="_")
      sumstats[i,"key"]<- paste(chrom,":",position,"_",alleles,sep="")
    }
}

library(dplyr)
# Try to do this with dplyr()
df_before <- df %>% 
   mutate(a1=sum(a+1), .before = )


FSH_F_EUR <- FSH_F_EUR %>% 
              mutate(MarkerName=paste(CHROM,":",position,"_",alleles,sep=""))

# Make a function for sorting alleles in alphabetical order and joining with an underscore

sort_alleles <- function(a1,a2) {
  a1 <- toupper(a1)
  a2 <- toupper(a2)
  sorted_alleles <- paste(sort(c(a1,a2)), collapse="_")
  return(sorted_alleles)
}

sort_alleles <- function(a1, a2) {
  a1 <- toupper(a1)
  a2 <- toupper(a2)
  sorted_alleles <- paste(sort(c(a1, a2))[1:2], collapse = "_")
  return(sorted_alleles)
}

test %>% mutate(MarkerName=paste(CHROM,GENPOS,sort_alleles(Allele1,Allele2), sep=":"), .before = ID)

test <- head(FSH_F_EUR)

FSH_F_EUR %>% mutate(MarkerName=paste(CHROM,GENPOS,sep=":",sort_alleles(Allele1,Allele2)), .before = ID)

# debug

