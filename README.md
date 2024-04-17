# interva

[![image](https://img.shields.io/pypi/pyversions/interva)](https://pypi.org/project/interva/)
[![pytest](https://github.com/verbal-autopsy-software/interva/actions/workflows/python-package.yml/badge.svg)](https://github.com/verbal-autopsy-software/interva/actions)

Python implementation of the InterVA (version 5) algorithm for assigning causes of death to verbal autopsy data.  This package replicates the R
version [InterVA5](https://github.com/verbal-autopsy-software/InterVA5).


## Importing package and installing dependencies

To install all package dependencies, run:  

```python
pip install interva
```

To import this package's functions:  

```python
from interva.interva5 import InterVA5
```
## Example data

To access example data from the package:  

```python
from interva.interva5 import get_example_input
va_data = get_example_input()
```

```python
va_data
       ID i004a i004b i019a i019b i022a  ... i454o i455o i456o i457o i458o i459o
0      d1     .     .     y     .     y  ...     n     n     n     n     n     n
1      d2     .     .     .     y     y  ...     n     n     n     n     n     n
2      d3     .     .     y     .     .  ...     n     n     n     n     n     n
3      d4     .     .     .     y     .  ...     n     n     n     n     n     n
4      d5     .     .     y     .     .  ...     n     n     n     n     n     n
..    ...   ...   ...   ...   ...   ...  ...   ...   ...   ...   ...   ...   ...
195  d196     .     .     .     y     .  ...     n     n     n     n     n     n
196  d197     .     .     y     .     y  ...     n     n     n     n     n     n
197  d198     .     .     y     .     y  ...     n     n     n     n     n     n
198  d199     .     .     .     y     y  ...     n     n     n     n     n     n
199  d200     .     .     .     y     y  ...     n     n     n     n     n     n

[200 rows x 354 columns]
```
  
## Creating and running an InterVA5 object

To initialize an InterVA5 object (results will be written to `VA_outpt/VA5_results.csv`):

```python
iv5out = InterVA5(va_data, hiv="h", malaria="l", write=True, directory="VA_output", filename="VA5_result", output="extended", append=False, return_checked_data=True)
```
  
To run the InterVA5 algorithm on the InterVA5 object:  

```python
iv5out.run()
```

```python
Using Probbase version: probbase v19 20210720
..........10% completed
..........20% completed
..........30% completed
..........40% completed
..........50% completed
..........60% completed
..........70% completed
..........80% completed
..........90% completed
..........100% completed
100% completed
```

Get cause-specific mortality fractions (CSMFs) for the top 8 causes:

```python
iv5out.get_csmf(top=8)
```

```python
HIV/AIDS related death               0.194408
Undetermined                         0.136583
Digestive neoplasms                  0.083285
Other and unspecified infect dis     0.063096
Renal failure                        0.061253
Reproductive neoplasms MF            0.053655
Other and unspecified cardiac dis    0.047557
Stroke                               0.045583
dtype: float64
```

getters for InterVA5 parameters:  

```python
iv5out.get_hiv()
iv5out.get_malaria()
iv5out.get_ids()
```

```python
HIV parameter is h
Malaria parameter is l
ids
0        d1
1        d2
2        d3
3        d4
4        d5
       ...
195    d196
196    d197
197    d198
198    d199
199    d200
Name: ID, Length: 200, dtype: object
```

To change the parameters, use setters and re-run.   

```python
iv5out.set_hiv("v")
iv5out.set_malaria("v")
iv5out.run()
```
