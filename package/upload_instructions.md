# Instructions on how to upload the ancIBD package to PYPI
### 2026, Harald Ringbauer

This instruction note summarizes the update and upload instructions described in detail at https://packaging.python.org/tutorials/packaging-projects/. It also contains specific instructions for the developers' Cluster environments in ancIBD.

I distribute the source - and no wheels. Because the Cython extension makes it non-Pythonic, building the wheel is nontrivial.

## First steps
- On Leipzig HPC Cluster (primary source to create package since 2025):
Activate Python environment and go to the hapROH package folder:

pyenvhpc312
cd /mnt/archgen/users/hringbauer/git/hapBLOCK/package/

- On Harvard O2 cluster:
cd /n/groups/reich/hringbauer/git/hapBLOCK/package/

-- [with Python 3.7 on O2]
envpython37

-- [with Python 3.8 on O2]
module load gcc/9.2.0
module load python/3.8.12

### Create the Source Package 
Update version in `./pyproject.toml` to next version number and update `./change_log.md`

### Remove prior version: 
Delete previous ./dist/* (alternatively, be specific below what to upload):  

rm ./dist/*


### Run the setup file:

### For local test
pip3 install ./

### Build package (only source, no wheel)
%%% [Legacy <=2025 v0.7] python3 setup.py sdist
python3 -m build --sdist


### [Optional] Run local install to finish tests (e.g., bash commands):
python3 -m pip install ./

### Run Tests of Expected Behavior
Use `/notebooks/tests/unit_tests.ipynb` to run tests of `ancIBD` (on O2 Harvard cluster).
- Leipzig Unit test in progress in `/notebooks/tests/unit_tests_leipzig.ipynb`

### Upload to the Sources to the official PyPi server:
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
Follow instructions on the PyPi site of `ancIBD`.

### Further Reading:
- For packaging: 
https://packaging.python.org/tutorials/packaging-projects/

- For Versioning:
https://www.python.org/dev/peps/pep-0440/

- To set up Token for Twine Upload to PyPi
https://pypi.org/manage/account/token/
