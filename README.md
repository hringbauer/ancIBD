# hapBLOCK
This package calls IBD with an HMM from phased/imputed data.

I currently work with the following versions (but should work on a much broader set):
module load gcc/6.2.0
module load python/3.7.4

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



