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
If the module needs to be added to the system path, us `sys.path.append`. 

## Usage

PadeOpsIO is used for visualization of and analysis for output data from PadeOps. 

Data can be instanteous data: 
![u_504](https://user-images.githubusercontent.com/8905274/197599663-3869848d-0b21-4759-896f-fecfdfe0cfdf.png)

or time-averaged data: 
![image](https://user-images.githubusercontent.com/8905274/197600823-d5d019ce-923b-4f60-af80-c2111efbcc90.png)
