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
path <- "data/sumstats/other/"

# Infertility summary stats
infert_filename <- "Ncol_Infertility1_F_EUR.txt"

# Hormones summary stats
filelist <- list.files(path,pattern='MN_')
			  
##############################
# 3 - Start main code
##############################
# Read in infertility sumstats
Infertility1_F_EUR <- fread(paste0(path,infert_filename))

# Read in hormone sumstats
sumstats_list = c()
for(i in seq_along(filelist)){
  sumstats <- fread(paste0(path,filelist[i]))
  name <- stri_replace_all_regex(filelist[i], pattern=c('MN_Ncol_', '.txt'), replacement=c('',''), vectorize=FALSE)
  assign(name, sumstats)
  sumstats_list <- append(sumstats_list,name) 
}

# Clean up
rm(sumstats)

# Make a sample size table manually for hormones
# col1 <- paste0(sumstats_list)
# col2 <- c(19117, 16610, 54522, 5169, 197234, 381985)
# df <- data.frame(col1, col2)
# colnames(df) <- c("sumstats","sample_size")
# fwrite(df,"sample_size.txt",sep="\t",col.names=FALSE)

# Check all of the data tables we have loaded in
for(i in 1:length(sumstats_list)) {print(head(get(sumstats_list[i])))}
head(Infertility1_F_EUR)

# Now filter infertility sumstats 
# Add MAF column
Infertility1_F_EUR[, MAF := pmin(Freq1, 1-Freq1)]

# Add infertility to sumstats list
sumstats_list <- append(sumstats_list,"Infertility1_F_EUR")

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

# Change column names of infertility sumstats in preparation for munging
setnames(Infertility1_F_EUR, "P-value", "PVALUE")
setnames(Infertility1_F_EUR, "StdErr", "SE")

# Check data before writing to file
for(i in 1:length(sumstats_list)) {print(get(sumstats_list[i]))}

# Write sumstats to file
for(name in c(sumstats_list)){
	sumstats <- get(name)
	out <- paste0(path,"MAFfiltered_MN_Ncol_",name,".txt")
	fwrite(sumstats,out,sep='\t')
}

# End
