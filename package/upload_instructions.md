# Instructions how to upload package to PYPI
### June 2022, Harald Ringbauer


This instruction note is a brief summary of update and upload instructions from https://packaging.python.org/tutorials/packaging-projects/

### Go to base folder
On Leipzig cluster:
cd /mnt/archgen/users/hringbauer/git/hapBLOCK/package

On Harvard O2 cluster:
cd /home/hr97/hringbauer/git/hapBLOCK/package
envpython37

### Run Tests of Expected Behavior
Use `/notebooks/tests/unit_tests.ipynb` to run tests of expected behavior of `ancIBD`.

### Create the Source Package 
Update version in setup.py to next version number

### Update setuptools. 
Delete previous ./dist/* (alternatively be specific below what to upload):  

rm ./dist/*

Run the setup file:
python3 setup.py sdist

### Upload to the Sources to official PyPi server:
python3 -m twine upload dist/* 


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
