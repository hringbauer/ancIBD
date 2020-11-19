# hapBLOCK
This package calls IBD with an HMM from phased/imputed data.

I currently work with the following versions (but should work on a much broader set):
module load gcc/6.2.0
module load python/3.7.4

### Produce C extension
In order to run the code, you need to build the C extension that implements the forward/backward pass.

On o2 cluster: Do this from my Python Environment (loaded with `envpython37`). 
Then switch to th correct folder and build the extensions via:

cd python3/
module load gcc/6.2.0
module load python/3.7.4
cythonize -a -i cfunc.pyx

