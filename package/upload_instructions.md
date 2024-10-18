# Instructions how to upload package to PYPI
### June 2022, Harald Ringbauer


This instruction note gives a summary of update and upload instructions described in detail in https://packaging.python.org/tutorials/packaging-projects/. It also contains specific instructions for the Cluster environments of the developers of ancIBD.

### Go to base folder
On Leipzig cluster:
cd /mnt/archgen/users/hringbauer/git/hapBLOCK/package

On Harvard O2 cluster:
cd /n/groups/reich/hringbauer/git/hapBLOCK/package/

## [with Python 3.7 on O2]
envpython37

## [with Python 3.8 on O2]
module load gcc/9.2.0
module load python/3.8.12

### Create the Source Package 
Update version in setup.py to next version number

### Update setuptools. 
Delete previous ./dist/* (alternatively be specific below what to upload):  

rm ./dist/*

### Run the setup file:
python3 setup.py sdist

### [Optional] Run local install to finish tests (e.g. bash commands):
python3 -m pip install ./

### Run Tests of Expected Behavior
Use `/notebooks/tests/unit_tests.ipynb` to run tests of `ancIBD`.

### Upload to the Sources to official PyPi server:
python3 -m twine upload dist/* 


## Check Installation in Leipzig
python3 -m pip install --user --upgrade --no-deps --force-reinstall ancIBD

# Bonus Material and Checks
## Alternatively: Upload on test server (for testing)
python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/* 

## Manually test whether extensions builds
python3 setup.py build_ext --inplace

## Test the Python test server package
python3 -m pip install --index-url https://test.pypi.org/simple/ ancIBD

# Further Documentation 
### To install via pip:
Follow instructions on PyPi site of `ancIBD`.

### Further Reading:
- For packaging: 
https://packaging.python.org/tutorials/packaging-projects/

- For Versioning:
https://www.python.org/dev/peps/pep-0440/

- To set up Token for Twine Upload to PyPi
https://pypi.org/manage/account/token/
