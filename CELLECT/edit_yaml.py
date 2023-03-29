#################
# Arg 1 = Phenotype
# Arg 2 =
# Arg 3 = 
#################

with open("config2.yml") as f:
    y = yaml.safe_load(f)
    y['db']['admin']['password'] = 'new_admin_pass'
    print(yaml.dump(y, default_flow_style=False, sort_keys=False))

import yaml
import sys

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
    
