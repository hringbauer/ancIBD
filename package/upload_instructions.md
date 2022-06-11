# Instructions how to upload package to PYPI
### June 2022, Harald Ringbauer


This instruction note is a brief summary of update and upload instructions from https://packaging.python.org/tutorials/packaging-projects/

### Go to base folder
cd /home/hr97/hringbauer/git/hapBLOCK/package

On Harvard O2 cluster:  
envpython37

### Run Tests of Expected Behavior
Use `/notebooks/tests/unit_tests.ipynb` to run tests of expected behavior of hapROH

### Create the Source Package 
Update version in setup.py to next version number

### Update setuptools. 
Delete previous ./dist/* (alternatively be specific below what to upload):  

rm ./dist/*

Run the setup file:
python3 setup.py sdist

### Upload to the Sources (copy into shell, to interactively do it!)
### For full PyPi server
python3 -m twine upload dist/* 
### Alternatively: Upload on test server (for testing)
python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/* 

# Bonus Material and Checks

## To test whether extensions builds
python3 setup.py build_ext --inplace

## Test the Python test server package
python3 -m pip install --index-url https://test.pypi.org/simple/ ancIBD

# Further Documentation 
### To install via pip:
Follow instructions on pypi site of `hapROH`.

### for packaging: 
https://packaging.python.org/tutorials/packaging-projects/

### For C extension (here cython):

### for version numbers:
https://www.python.org/dev/peps/pep-0440/
