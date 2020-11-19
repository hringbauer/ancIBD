# hapBLOCK
This package calls IBD with an HMM from phased/imputed data.

It is tested with the follwowing versions (but should work on a much broader set):
module load gcc/6.2.0
module load python/3.7.4

### Produce C extension
On cluster: Do this from my Python Environment (loaded with envpython37)

cd python3/
module load gcc/6.2.0
module load python/3.7.4
cythonize -a -i cfunc.pyx

