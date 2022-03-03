# Binding energy app

This repository holds an app written in python to plot and export binding energy values of atomic elements from Z = 1 - 92. This app considers data from four sources contained in the ```data/``` folder, and they are described below.

### ```experimental```
Experimental measurements from solids, compiled by [Williams (1995)](https://xdb.lbl.gov/Section1/Sec_1-1.html) from references 
1. [Bearden and Burr (1967)](https://doi.org/10.1103/RevModPhys.39.125), 
2. [Photoemission in Solids I (1978)](https://doi.org/10.1007/3-540-08685-4), and 
3. [Fuggle and Mårtensson (1980)](https://doi.org/10.1016/0368-2048(80)85056-0) .

### ```dirac-fock``` 
Results by [Desclaux (1973)](https://doi.org/10.1016/0092-640X(73)90020-X) computed with the Dirac-Fock method.

### ```perturbative```
Relativistic values computed with the parametric potential method by [M. Klapisch (1971)](https://doi.org/10.1016/0010-4655(71)90001-4) implemented in the [```HULLAC```](https://doi.org/10.1016/S0022-4073(01)00066-8) code. Some of these results have been published in [Mendez et al. (2019)](https://doi.org/10.1016/j.nimb.2019.02.002).

### ```hartree-fock```
Theoretical calculations with the ```HF86``` code by [C. Froese Fischer (1987)](https://doi.org/10.1016/0010-4655(87)90053-1).



**IMPORTANT**: The app has not been packaged (yet) and it requires ```dash``` ```dash-bootstrap-components```, ```jupyter-dash```, ```periodictable```, ```pandas```, ```numpy```, ```scipy```, ```os```, ```re``` and ```nbformat``` to work. At this point, the user should install them manually (for example, using pip or conda).

