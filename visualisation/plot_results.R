#############################################################
## Creating boxplots from CELLECT results
## melodyjparker14@gmail.com - Apr 2023
## Run from p50 dir
## Plot size needs adjusting for different number of phenotypes.
#############################################################

##############################
# 0 - Load librairies
##############################
library(dplyr)  # for mutate 
library(data.table)

############################## 
# 1 - Choose data
##############################
# edit this based on chosen data/phenotypes

# phenotypes need to be consistent with CELLECT output directory names
# for the alternative annotations, add "_ALT" on the end

phenotype_table <- fread("sumstats/sample_size.txt", select = 1)

phenotypes <- c()
for(i in phenotype_table){phenotypes <- append(phenotypes,i)}

############################## 
# 2 - Source file
##############################
# paste general results data file
data_type <- paste("CELLECT", toupper(data_source), toupper(tissue_type), toupper(RNAseq_type), sep="-")
# "CELLECT-ATLAS-ADIPO-SN"

# input files need to be inside loop in main code
output_file <- paste("CELLECT/plots/", data_source, "_", tissue_type, "_", RNAseq_type, sep="")

# CELLECT output directory
library(stringr)

datasets <- list.files("counts")

##############################
# 3 - Define functions
##############################
# function for adding -log(p value) column
resultsCol <- function(x) {
  mutate(prioritization, "{x}" := -log10(pvalue))
}

##############################
# 4 - Start code
##############################

for (dataset in datasets){

# dataset <- datasets[1]

  # Extract last three characters of dataset
  ds <- str_sub(dataset,-3,-1)

  cellect_out_dir <- paste0('CELLECT_OUT_',ds)

  # Create a joint results table for a dataset
  for (phenotype in phenotypes){
    input_file <- paste(cellect_out_dir, phenotype, "CELLECT-LDSC/results/prioritization.csv", sep="/")
    # Read in results
    prioritization <- read.csv(input_file)
    # Use function to add -log(pvalue) column
    prioritization <- resultsCol(phenotype)
    # Add phenotype to joint results table
    if (exists('joint_results')){
      joint_results <- full_join(df1, prioritization[,c(3,7)], by = "annotation")
    } else {
      joint_results <- prioritization[,c(3,7)]    
    }
  }

  joint_results <- joint_results[order(tolower(joint_results$annotation)),]
  rownames(joint_results) <- joint_results[,1]
  joint_results <- joint_results[,-1]
  matrix <- as.matrix(joint_results)
  rm(joint_results)

  # Choose bonferroni-adjusted significance threshold
  n <- nrow(matrix)
  m <- ncol(matrix)
  sig <-0.05
  b_sig = -log10((sig/n))
  y_lim <- ceiling(max(matrix))

  output_file <- paste0('plots/CELLECT_plot_hormones_',dataset,".pdf")

    # Generate barplot
    title <- paste0("CELLECT results with hormones sumstats and dataset ",dataset)

    pdf(file=output_file, width=30)
    {par(lwd = 2)
    barplot(matrix, space=c(0,8), width = 1, border=T, beside = TRUE, col = terrain.colors(n), xlim = c(6, y_lim +150), ylim = c(0, y_lim), legend = TRUE, args.legend = list(ncol = 4, bty = "n", xjust=1, x = "top"))
    abline(h=b_sig,lty = 2)
    title(ylab = '-log(p value)', line=2.2, font.lab=1)
    title(main = title, adj = 0.5)
    abline(h=0)
    }
 
    dev.off()
  rm(matrix)
    
}
