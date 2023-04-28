#######################################################
# Adding an N column to summary stats data using the direction column and n case/control per study
# melodyjparker14@gmail.com - Apr 23
#######################################################

##############################
# Load libraries
##############################
library(dplyr)
library(data.table)

##############################
# Define file names
##############################
in_file <- "data/sumstats/other/Infertility1_F_EUR_directions.txt"
out_file <- "data/sumstats/other/Infertility1_F_EUR.txt"

##############################
# Functions
##############################
calcNeff <- function(n_cases, n_controls){ 
  N_eff <- 2 / (1 / n_cases + 1 / n_controls)
  return(N_eff)
}

##############################
# Start code
##############################
# Read in the direction file - this contains 0s and 1s
direc <- fread(in_file)

# N Case/controls/Neff
FinnGen_N_Cases <- 14759
FinnGen_N_Controls <- 111583
FinnGen_N_eff <- calcNeff(FinnGen_N_Cases, FinnGen_N_Controls)  # 26069.77

UKBB_N_Cases <- 3447
UKBB_N_Controls <- 245591
UKBB_N_eff <- calcNeff(UKBB_N_Cases, UKBB_N_Controls)  # 6798.578

EstBB_N_Cases <- 12027
EstBB_N_Controls <- 111395
EstBB_N_eff <- calcNeff(EstBB_N_Cases, EstBB_N_Controls)  # 21710.03
  
infert2 <- infert %>% 
  mutate(N = (direc$FinnGen * FinnGen_N_eff) + (direc$UKBB * UKBB_N_eff) + (direc$EstBB * EstBB_N_eff))

infert3 <- infert2 %>% 
  mutate(N2 = calcNeff( (direc$FinnGen*FinnGen_N_Cases + direc$UKBB*UKBB_N_Cases + direc$EstBB*EstBB_N_Cases), (direc$FinnGen*FinnGen_N_Controls + direc$UKBB*UKBB_N_Controls + direc$EstBB*EstBB_N_Controls)))

# It can be shown that by splitting the case-control data into subsets, the sum of the effective sample sizes over the subsets is always ≤ the effective sample size of the whole data
# https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS9.html

# Make a scatter plot of the two N-eff columns for the different combinations of studies
df <- data.frame (FinnGen  = c(1,0,1,1),
                  UKBB = c(1,1,0,1),
                  EstBB = c(0,1,1,1)
                  )

df <- df %>% 
  mutate(N = (FinnGen * FinnGen_N_eff) + (UKBB * UKBB_N_eff) + (EstBB * EstBB_N_eff)) %>%
  mutate(N2 = calcNeff( (FinnGen*FinnGen_N_Cases + UKBB*UKBB_N_Cases + EstBB*EstBB_N_Cases), (FinnGen*FinnGen_N_Controls + UKBB*UKBB_N_Controls + EstBB*EstBB_N_Controls)))
  
plot(df$N,df$N2, main="Comparison of N-eff when calculated in two different ways",
   xlab="N ", ylab="N2 ", pch=19)

fwrite(infert3,outfile,sep=" ")

##############################
# End
##############################
