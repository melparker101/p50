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
library(stringr)
library(dplyr)
library(purrr)

############################## 
# 2 - Source files
##############################
# Infertility summary stats
filename <- 'female_infertility_analysis1_UKBB_Finngen_EstBB_GandH_noMACfilter_March20231.out'
# Hormones summary stats
filelist <- list.files(pattern='.txt')

##############################
# 3 - Define functions
##############################
# Make a function for sorting alleles in alphabetical order and joining with an underscore
sort_alleles <- function(a1,a2) {
  a1 <- toupper(a1)
  a2 <- toupper(a2)
  sorted_alleles <- paste(sort(c(a1,a2)), collapse="_")
  return(sorted_alleles)
}

##############################
# 4 - Start main code
##############################
# Read in infertility sumstats
Infertility1_F_Finngen <- fread(filename)

# Read in hormone sumstats
for(i in seq_along(filelist)){
  name <- str_replace(filelist[i],"_filtered.txt","")
  sumstats <- fread(filelist[i])
  assign(name, sumstats)
}

# Clean up
rm(sumstats)

#####################################

FSH_F_EUR <- FSH_F_EUR %>% 
              mutate(MarkerName=paste(CHROM,":",position,"_",alleles,sep=""))


FSH_F_EUR

test %>% mutate(MarkerName=paste(CHROM,GENPOS,sort_alleles(Allele1,Allele2), sep=":"), .before = ID)

# it will also have the side effect of grouping and outputting a tibble, instead of a data.frame
# https://stackoverflow.com/questions/21818181/applying-a-function-to-every-row-of-a-table-using-dplyr
test <- as.data.frame(test %>%
  rowwise() %>%
  mutate(MarkerName = paste(CHROM, GENPOS, sort_alleles(Allele1, Allele2), sep = ":"), .before = ID))



#########################

test <- head(FSH_F_EUR)

test %>% 
  mutate(MarkerName = paste(CHROM,GENPOS,map2_chr(Allele1, Allele2, sort_alleles), sep = ":"), .before = ID)
#########
# Do a comparison test of both of these

start <- Sys.time()
FSH_F_EUR_edited <- as.data.frame(FSH_F_EUR %>%
  rowwise() %>%
  mutate(MarkerName = paste(CHROM, GENPOS, sort_alleles(Allele1, Allele2), sep = ":"), .before = ID))
print( Sys.time() - start )

# Time difference of 18.23314 mins

start <- Sys.time()
FSH_F_EUR_edited2 <- FSH_F_EUR %>% 
  mutate(MarkerName = paste(CHROM,GENPOS,map2_chr(Allele1, Allele2, sort_alleles), sep = ":"), .before = ID)
print( Sys.time() - start )

# Time difference of 14.10276 mins

# Still a bit slow but we can just send them off
# An awk command will be too complicated and hard to read/understand for anyone else



# Try it with data.table
FSH_F_EUR_edited_dt <- FSH_F_EUR 

start <- Sys.time()
FSH_F_EUR_edited_dt[, MarkerName := paste(CHROM, GENPOS, unlist(Map(sort_alleles, Allele1, Allele2)), sep = ":")]
print( Sys.time() - start )

test[, MarkerName := paste(CHROM,GENPOS,map2_chr(Allele1, Allele2, sort_alleles), sep = ":"), .before = ID]
test[, MarkerName := paste(CHROM,GENPOS, sort_alleles(Allele1, Allele2), sep = ":"), .before = ID]

 DT <- data.table(New = c(1:10), DT)
DT[New = rpois(20,2)]


# add a new column as the first column
test[, MarkerName:= paste(CHROM,GENPOS, do.call(sort_alleles(Allele1, Allele2)), sep = ":")]
test <- test[, c(MarkerName,colnames(test)[-1])]

# view the updated data table
dt

pmin <- function (...) UseMethod("pmin")
pmin.default <- function (...) {
  res <- .Primitive("pmin")(...)
  if (anyNA(res)) res[is.na(res)] <- -Inf
  res
}
#######################
#try this
FSH_F_EUR_edited_dt <- FSH_F_EUR 

sort_alleles_vec <- function(a1, a2) {
  a1 <- toupper(a1)
  a2 <- toupper(a2)
  alleles <- matrix(c(a1, a2), ncol = 2)
  sorted_alleles <- apply(alleles, 1, function(x) paste(sort(x), collapse = "_"))
  return(sorted_alleles)
}
			  
start <- Sys.time()
FSH_F_EUR_edited_dt[, MarkerName := paste(CHROM, GENPOS, sort_alleles_vec(Allele1, Allele2), sep = ":")]
print( Sys.time() - start )
# Time difference of 10.81209 mins

start <- Sys.time()
FSH_F_EUR_edited_dt[, MarkerName := paste(CHROM, GENPOS, sort_alleles_vec(Allele1, Allele2), sep = ":")]
print( Sys.time() - start )			  

# Other suggestions
			  
# 1.
set(test, j = "MarkerName", value = paste(test$CHROM, test$GENPOS, sort_alleles(test$Allele1, test$Allele2), sep = ":"))
			  
# 2. Use stringr's str_c() function: The str_c() function from the stringr package is a faster alternative to base R's paste() function. 
set(test, j = "MarkerName", value = str_c(test$CHROM, ":", test$GENPOS, ":", sort_alleles(test$Allele1, test$Allele2)))
			 
#3. Subset the data
subset_test <- test[, .(CHROM, GENPOS, Allele1, Allele2)]
set(subset_test, j = "MarkerName", value = paste(subset_test$CHROM, subset_test$GENPOS, sort_alleles(subset_test$Allele1, subset_test$Allele2), sep = ":"))
set(test, i = which(test$MarkerName %in% subset_test$MarkerName), j = "MarkerName", value = subset_test$MarkerName)

			  
			  
# Try 2
start <- Sys.time()			  
set(FSH_F_EUR_set, j = "MarkerName", value = str_c(FSH_F_EUR_set$CHROM, ":", FSH_F_EUR_set$GENPOS, ":", sort_alleles(FSH_F_EUR_set$Allele1, FSH_F_EUR_set$Allele2)))
print( Sys.time() - start )	
# took ages so I shopped it			  
			  
############################
			  
start <- Sys.time()
FSH_F_EUR_edited_dt[, MarkerName := paste(CHROM, GENPOS, unlist(Map(sort_alleles, Allele1, Allele2)), sep = ":")]
print( Sys.time() - start )
# Time difference of 16.78605 mins
			  

############
# MAF > 1%
dt <- data.table(dt)
dt[, MAF := pmin(Freq1, 1-Freq1)]
dt %>% filter(MAF>0.01)
