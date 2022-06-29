# interva
Python implementation of the InterVA algorithm for assigning causes of death to verbal autopsy data

## Importing package and installing dependencies
To import this package's functions:  
`from interva.interva5 import InterVA5`  
  
To install all package dependencies, run:  
`python setup.py install`  
  
## Example data
To access example data from the package:  
`va_data_csv = pkgutil.get_data("interva", "data/randomva5.csv")`  
`va_data = read_csv(BytesIO(va_data_csv))`  
  
## Creating and running an InterVA5 object
To initialize an InterVA5 object:  
`iv5out = InterVA5(va_data, hiv="h", malaria="l", write=False, directory="VA test", filename="VA5_result", output="extended", append=False, return_checked_data=True)`  
  
To run the InterVA5 algorithm on the InterVA5 object:  
`run_output = iv5out.run()`  
  
To access the algorithm output:  
`id_output = run_output["ID"]`  
`va5_output = run_output["VA5"]`  
`malaria_output = run_output["Malaria"]`  
`hiv_output = run_output["HIV"]`  
`checked_data_output = run_output["checked_data"]`  
  
To get likelihoods of HIV or malaria as a cause of death from the InterVA5 object:  
`hiv = iv5out.get_hiv()`  
`malaria = iv5out.get_malaria()`  
  
To set likelihoods ("h", "l", "v") of HIV or malaria as a cause of death for the InterVA5 object:  
`iv5out.set_hiv("l")`  
`iv5out.set_malaria("v")`  
  
To get ids from the InterVA5 object:  
`ids = iv5out.get_ids()`  
  