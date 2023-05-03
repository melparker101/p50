#############################################################
## Add marker name column to hormone sumstats
## melodyjparker14@gmail.com - April 2023
## This script adds a MarkerName column to sumstats data
## Recommended cores > 10
## This code is actually quite slow - use the add_MarkerName_col_hormones.sh instead
## Not fully tested, but a template for adding a marker name column using R
#############################################################

##############################
# 1 - Load librairies
##############################
library(data.table)
library(stringr)
library(dplyr)
# library(doParallel)
# library(foreach)

############################## 
# 2 - Source files
##############################
# Infertility summary stats
filename <- 'female_infertility_analysis1_UKBB_Finngen_EstBB_GandH_noMACfilter_March20231.out'

# Hormones summary stats
# filelist <- list.files(pattern='.txt')

##############################
# 3 - Define functions
##############################
# Function to sort alleles in alphabetical order and concatenate with an underscore
sort_alleles_vec <- function(a1, a2) {
  a1 <- toupper(a1)
  a2 <- toupper(a2)
  alleles <- matrix(c(a1, a2), ncol = 2)
  sorted_alleles <- apply(alleles, 1, function(x) paste(sort(x), collapse = "_"))
  return(sorted_alleles)
}
			  
##############################
# 4 - Start main code
##############################
# Read in infertility sumstats
Infertility1_F_Finngen <- fread(filename)

# Read in hormone sumstats
sumstats_list = c()
for(i in seq_along(filelist)){
  name <- str_replace(filelist[i],"_filtered.txt","")
  sumstats <- fread(filelist[i])
  assign(name, sumstats)
  name <- data.table(name)
  sumstats_list <- append(sumstats_list,name)
}

# Clean up
rm(sumstats)
		  
# This is the fastest method that I've tested in R, however, awk is way way faster so use awk beforehand
# It still took ~10 mins per data table - awk took ~15 seconds for each
# Even if this is run in parallel rather than a loop it is still a lot slower
start <- Sys.time()
for(sumstats in sumstats_list){
		sumstats[, MarkerName := paste(CHROM, GENPOS, sort_alleles_vec(Allele1, Allele2), sep = ":")]
    # Outfile name = markernameR_sumstats.txt
    out <- paste0("MNR_",sumstats,".txt")
    fwrite(sumstats,out)
}
print( Sys.time() - start )	

# Debugging required for running in parallel
if(FALSE) {
  start <- Sys.time()
  n=length(filelist)
  foreach(i=1:n) %dopar% {
	  sumstats <- sumstats_list[i]
	  sumstats[, MarkerName := paste(CHROM, GENPOS, sort_alleles_vec(Allele1, Allele2), sep = ":")]
	}
  print(Sys.time() - start)
}
