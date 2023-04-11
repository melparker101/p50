###################################################
# This script is set up to edit the config yaml file to run CELLECT per phenotype
# melodyjparker14@gmail.com
#
# Arg 1 = Phenotype
# Arg 2 = Dataset
#
# There are three key-value pairs that need editing:
# - BASE_OUTPUT_DIR
# - SPECIFICITY_INPUT
# - GWAS_SUMSTATS
###################################################

import yaml
import sys
import pandas as pd
import os

# Extract phenotypes from the index file
index_file = 'sumstats/sample_size.txt'
df = pd.read_csv(index_file,sep='\t',usecols=[0],header=False)
phenotypes = df.iloc[:,0].values.tolist()
# phenotypes = ['FSH_F_EUR', 'LH_F_EUR', 'Oestradiol_F_EUR', 'Progesterone_F_EUR', 'Testosterone_F_EUR', 'Testosterone_sex_comb_EUR']

# List all specificity files in the esmu dir
specificity_files = os.listdir("esmu")
# Only keep files with ensembl id (ens) and mean (esmu); filter out files with gene symbols (sym) or standard deviation (essd)
specificity_files = list(filter(lambda x:'ens.esmu' in x, specificity_files))
# specificity_files = ['GSE202601_ens.esmu.csv.gz', 'GSE118127_ens.esmu.csv.gz', 'GSE213216_ens.esmu.csv.gz']

input_file = config.yml

# Make a config directory to store new files in
isExist = os.path.exists("config")
if isExist == 'False':
    os.mkdir("config")

# Read in yaml
with open(input_file,'r') as f:
    y = yaml.safe_load(f)

# Create new yaml files per phenotype
for phenotype in phenotypes:
    output_file = 'config/config_' +  phenotype + '.yml'
    # Edit 1
    key = 'BASE_OUTPUT_DIR'
    value = 'CELLECT_OUT/' + phenotype
    y[key] = value
    # Edit 2
    key = 'SPECIFICITY_INPUT'
    subkey = 'id' 
    value = dataset
    y[key]['id'] = 'new_test_id'
    subkey = 'path'
    value = specificity_file 

    y['key'] = value
    # Edit 3
    y['key'] = value
    with open(output_file,'w') as f:
        yaml.safe_dump(y,f)



############################################################
# Check for correct number of arguments
if len(sys.argv) != 3:
    print("Usage: python edit_yaml.py <phenotype> <dataset>")
    exit()
    
phenotype = sys.argv[1]  # e.g. LH
dataset = sys.argv[2]  # e.g. GSE118127

# Yaml files
input_file = config.yml
output_file = 'config_' +  phenotype + '.yml'

# Only use F EUR hormone sumstats for now
specificity_path = 'esmu/'
specificity_file = specificity_path + '/' + dataset + '_ens.esmu.csv.gz'  # e.g. esmu/GSE118127_ens.esmu.csv.gz
sumstats_path = 'sumstats/munged/'
sumstats_file = sumstats_path + 'munged_' phenotype '_F_EUR.sumstats.gz'  # e.g. sumstats/munged/munged_LH_F_EUR.sumstats.gz

# Load YAML file
with open(input_file, 'r') as file:
    data = yaml.safe_load(file)

# Modify YAML data
key, value = key_value_pair.split('=')
data[key] = value

# Write the updated YAML data back to file
with open(output_file, 'w') as file:
    yaml.safe_dump(data, file)
    
#######################
# Things we want to edit

# BASE_OUTPUT_DIR: CELLECT-EXAMPLE
key = 'BASE_OUTPUT_DIR'
value = 'CELLECT_OUT/' + phenotype
data[key] = value

#####

# Start with one specificity csv (dataset) for now
# Use the ENS

# SPECIFICITY_INPUT:
#   - id: mousebrain-test
#     path: example/mousebrain-test.csv # this test file contains 2 annotations

key = 'SPECIFICITY_INPUT'
subkey = 'id' 
value = dataset

key = 'SPECIFICITY_INPUT'
subkey = 'path'
value = specificity_file 

#####

# GWAS_SUMSTATS:
#   - id: EA3_Lee2018
#     path: example/EA3_Lee2018.sumstats.gz 
#   - id: BMI_Yengo2018
#     path: example/BMI_Yengo2018.sumstats.gz

key = 'GWAS_SUMSTATS'
subkey = 'id' 
value = phenotype

key = 'GWAS_SUMSTATS'
subkey = 'path'
value = sumstats_file 

key, value = key_value_pair.split('=')

data
    y['SPECIFICITY_INPUT']['id'] = 'new_test_id'
    print(yaml.dump(y, default_flow_style=False, sort_keys=False))


#################################
import yaml
import sys

# Check for correct number of arguments
if len(sys.argv) != 4:
    print("Usage: python edit_yaml.py <input_file> <output_file> <key_value_pair>")
    exit()

input_file = sys.argv[1]
output_file = sys.argv[2]
key_value_pair = sys.argv[3]

# Load YAML file
with open(input_file, 'r') as file:
    data = yaml.safe_load(file)

# Modify YAML data
key, value = key_value_pair.split('=')
data[key] = value

# Write the updated YAML data back to file
with open(output_file, 'w') as file:
    yaml.safe_dump(data, file)







#################################



with open("config2.yml") as f:
    y = yaml.safe_load(f)
    y['db']['admin']['password'] = 'new_admin_pass'
    print(yaml.dump(y, default_flow_style=False, sort_keys=False))

print(sys.argv[1])

with open("config.yaml") as f:
     list_doc = yaml.safe_load(f)

for sense in list_doc:
    if sense["name"] == "sense2":
         sense["value"] = 1234

with open("config2.yaml", "w") as f:
    yaml.dump(list_doc, f)
    
with open("config2.yml") as f:
    y = yaml.safe_load(f)
    y['SPECIFICITY_INPUT']['id'] = 'new_test_id'
    print(yaml.dump(y, default_flow_style=False, sort_keys=False))
    
    
BASE_OUTPUT_DIR:
SPECIFICITY_INPUT:
  
https://stackoverflow.com/questions/29518833/editing-yaml-file-by-python
  https://stackoverflow.com/questions/53351840/how-to-edit-a-yaml-file-in-python-repeatedly-using-a-variable-without-defining-a
    https://stackoverflow.com/questions/63581308/edit-yaml-file-with-bash
    
