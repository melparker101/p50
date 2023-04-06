#############################################################
## Filter Sumstats
## melodyjparker14@gmail.com - April 2023
## This script adds a MarkerName column to sumstats data and then
## Use a high number of cores (e.g. 10)
#############################################################

##############################
# 1 - Load librairies
##############################
library(data.table)
library(stringi)
library(dplyr)

############################## 
# 2 - Source files
##############################
# Infertility summary stats
filename <- 'female_infertility_analysis1_UKBB_Finngen_EstBB_GandH_noMACfilter_March20231.out'

# Hormones summary stats
# filelist <- list.files(pattern='.txt)
filelist <- list.files(pattern='MN_')
			  
##############################
# 3 - Start main code
##############################
# Read in infertility sumstats
Infertility1_F_Finngen <- fread(filename)

# Read in hormone sumstats
sumstats_list = c()
for(i in seq_along(filelist)){
  sumstats <- fread(filelist[i])
  name <- stri_replace_all_regex(filelist[i], pattern=c('MN_', '_filtered.txt'), replacement=c('',''), vectorize=FALSE)
  assign(name, sumstats)
  sumstats_list <- append(sumstats_list,name) 
}

# Clean up
rm(sumstats)

# Check all of the data tables we have loaded in
for(i in 1:length(sumstats_list)) {print(head(get(sumstats_list[i])))}

# Loop through hormone sumstats to filter
for (name in sumstats_list) {
  sumstats <- get(name)
  print("Before filtering for MAF")
  print(c(name,dim(sumstats)[1]))
  
  # Filter for MAF > 0.01
  sumstats_filtered <- sumstats %>% filter(MAF > 0.01)
  
  print("After filtering for MAF")
  print(c(name,dim(sumstats_filtered)[1]))
  
  # Overwrite original object with filtered object
  assign(name, sumstats_filtered, envir = .GlobalEnv)
}

# Now filter infertility sumstats 
# Add MAF column
Infertility1_F_Finngen[, MAF := pmin(Freq1, 1-Freq1)]

# Filter for MAF > 0.01
Infertility1_F_Finngen <- Infertility1_F_Finngen %>% filter(MAF>0.01)

# Add infertility to sumstats list
sumstats_list <- append(sumstats_list,"Infertility1_F_Finngen")

for(i in 1:length(sumstats_list)) {print(get(sumstats_list[i]))}

# Write sumstats to file
for(name in c(sumstats_list)){
	sumstats <- get(name)
	out <- paste0("MAF_filtered/MAFfiltered_",name,".txt")
	fwrite(sumstats,out,sep='\t')
}
