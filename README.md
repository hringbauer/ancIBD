# ancIBD
This GitHub repository contains development code for the Python package `ancIBD` that is [available via PyPi](https://pypi.org/project/ancIBD/). Please find the official documentation on [readthedocs](https://ancibd.readthedocs.io). 

This [github repository](https://github.com/hringbauer/ancIBD) contains the development code (in `./package/`), and also the code underlying the `ancIBD` publication (in particular in `./notebook/`).

The notes below are only for developers of this package who run development versions.

### Produce C extension
To run the code, you need to have the C extension that implements the forward/backward pass. When installing the Python package, this is done automatically. However, if you are a developer and want to test the latest changes to the packages, you will have to do that manually.

### On Leipzig Cluster
Simply cd into the python folder, and then cythonize:

cd package/ancIBD/
cythonize -a -i cfunc.pyx

### On o2 cluster 

### Python 3.13  (via Jupyter server)
Run a manual pip install:

cd package
python3 -m pip install --user .

### Python3.7 (outdated since update)
Do this from my Python Environment (loaded with `envpython37`). 
Then switch to the correct folder and build the extensions using the following commands:

envpython37  #load the python environment
cd package/ancIBD/  # or path where the cfunc.pyx file is located
module load gcc/6.2.0
module load python/3.7.4
cythonize -a -i cfunc.pyx

Warning: cythonizing can be killed on the head node (probably due to resource constraints). In that case, delete intermediate files and try again. Mid-term to-do: Need to build a batchable script.

### Python 3.8 (outdated since update):
cd package/ancIBD
module load gcc/9.2.0
module load python/3.8.12
cythonize -a -i cfunc.pyx








