# hapBLOCK
This github repository contains development code for the Python package `ancIBD` that is available via PyPi [here](https://pypi.org/project/ancIBD/). Please find the official instructions on Readthedocs [here](https://ancibd.readthedocs.io). The notes below are for the developers of this package who run development versions.

### Depracation warning
The code in `./python3` is all legacy code and will be removed soon. All Python code has now been moved to `./package/ancIBD` - and continues to be updated there.

### Produce C extension
In order to run the code, you need to have the C extension that implements the forward/backward pass.

When installing the Python package, this is done automatically. However, if you are a developer and want to test the latest changes to the packages you will have to do that manually.

### On Leipzig Cluster
Simply cd into the python folder, and then cythonize:

cd package/ancIBD/
cythonize -a -i cfunc.pyx

### On o2 cluster: Do this from my Python Environment (loaded with `envpython37`). 
Then switch to th correct folder and build the extensions via:

envpython37  #load the python environment
cd python3/  # or path where the cfunc.pyx file is located
module load gcc/6.2.0
module load python/3.7.4
cythonize -a -i cfunc.pyx



