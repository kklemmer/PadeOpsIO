# PadeOpsIO Overview

PadeOpsIO is a tool for loading, visualizing, and analyzing data from [PadeOps](https://github.com/FPAL-Stanford-University/PadeOps), an LES and DNS box code. 

## Dependencies

Aside from standard python library packages, PadeOpsIO requries: <br>
* `numpy` <br>
* `scipy` <br>
* `matplotlib` <br>
* `f90nml` (This package isn't required unless loading data from Fortran; conveinently parses input namelist files)

## Installation

Download the repo from github. The module can be imported into a Python script with: 
```
import padeopsIO
```
If the module needs to be added to the system path, use `sys.path.append`. 

## Usage

PadeOpsIO is used for visualization of and analysis for output data from PadeOps. For more, see the [quick start](https://github.com/kirbyh/PadeOpsIO/blob/main/padeopsIO/padeopsIO_quickstart.ipynb). 

Data can be instanteous data: 

![image](https://user-images.githubusercontent.com/8905274/197601106-86fd32e4-52dc-4cf5-bcc3-d1bd664cdc08.png)

or time-averaged: 

![image](https://user-images.githubusercontent.com/8905274/197600994-47325c6d-89f3-4d09-9a44-1a0822fe81b5.png)
