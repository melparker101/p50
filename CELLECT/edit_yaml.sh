#!/bin/bash

########################################################################
# This script uses sed to create a CELLECT config yaml file per phenotype
# It copies and edits the template config file 'config_template.yml'
########################################################################

# Out directory
CONFIG=config_127

mkdir -p "$CONFIG"

for phenotype in $(cut -f1 sumstats/sample_size.txt)
do
  sumstats_file=munged_"$phenotype".sumstats.gz
  config_out="$CONFIG"/config_"$phenotype".yml
  
  # Add backslashes before special characters ready for double quoted sed
  # sumstats_file=${sumstats_file////\\/}
  
  # Copy the template
  cp config_template.yml "$config_out"
  
  # Edit yaml
  # Use double quotes to expand out the shell variables
  sed -i "s/<phenotype>/$phenotype/" $config_out
  sed -i "s/<sumstats_id>/$phenotype/" $config_out
  sed -i "s/<sumstats_file>/$sumstats_file/" $config_out
done
