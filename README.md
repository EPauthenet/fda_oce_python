# fda_oce_python
Functional Data Analysis of hydrographic profiles for python (fda_oce_python)

**Functional Data Analysis** is a set of tools to study curves or functions. Here we see vertical hydrographic profiles of several variables (temperature, salinity, oxygen,...) as curves and apply a functional principal component analysis (FPCA) in the multivaraite case to reduce the dimensionality of the system. The classical case is done with couples of temperature and salinity. It can be used for front detection, water mass identification, unsupervised or supervised classification, model comparison, data calibration ... This toolbox is also available in [Matlab](https://github.com/EPauthenet/fda.oceM) and [R](https://github.com/EPauthenet/fda.oce) with a doi [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4073123.svg)](https://doi.org/10.5281/zenodo.4073123)

**Installation of fda.oce using Anaconda:

If you have installed Anaconda, you can install fda_oce_python in a new environment *fda_env* by executing in your terminal the commands:
``` r
conda create -n fda_env --yes
conda activate fda_env
conda install python=3.6 --yes
conda install -c r r-devtools rpy2 r-fda --yes
conda install -c conda-forge jupyterlab --yes
conda install -c conda-forge numpy --yes
conda install -c conda-forge scipy --yes
conda install -c conda-forge matplotlib --yes
```
Note that we use python 3.6 because at this date, the R module fda is not compatible with later versions of python.

Before opening jupyter lab (or jupyter notebook), make sure you select the right conda environment. This can be done from the Anaconda graphical interface, or using the command line:
``` r
conda activate fda_env
```

Download the github folder fda_oce_python, move there in your terminal and open jupyter lab:
``` r
conda activate fda_env
jupyter lab
```

You should now be able to execute the code in the fda_oce_python.ipynb file. Note that you may have to modify the variable *os.environ['R_HOME']* defined below in order to access R from python.

See the Jupyter notebook : https://github.com/EPauthenet/fda_oce_python/blob/master/fda_oce_python.ipynb

**Installation of fda.oce using an existing R install** If you have an installation of R obtained differently, you may install devtools by executing in R the command:
``` r
install.packages("devtools")
install.packages("fda")
```
You should then set the python variable *os.environ['R_HOME']* using the path to the home of R, in order to access R from python.


*Author*:
Fabien Roquet (fabien.roquet@gu.se)
Etienne Pauthenet (etienne.pauthenet@locean-ipsl.upmc.fr)

*References*: 
- Pauthenet et al. (2019) The thermohaline modes of the global ocean. Journal of Physical Oceanography, [10.1175/JPO-D-19-0120.1](https://doi.org/10.1175/JPO-D-19-0120.1)
- Pauthenet et al. (2018) Seasonal meandering of the Polar Front upstream of the Kerguelen Plateau. Geophysical Research Letters, [10.1029/2018GL079614](https://doi.org/10.1029/2018GL079614)
- Pauthenet et al. (2017) A linear decomposition of the Southern Ocean thermohaline structure. Journal of Physical Oceanography, [10.1175/JPO-D-16-0083.1](http://dx.doi.org/10.1175/JPO-D-16-0083.1)
- Ramsay, J. O., and B. W. Silverman, 2005: Functional Data Analysis. 2nd Edition Springer, 426 pp., Isbn : 038740080X.

