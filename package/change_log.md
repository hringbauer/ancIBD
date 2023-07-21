## VERSION NUMBER, DATE, AUTHOR
List of updates to ancIBD.

## 0.5, July 21st, Harald Ringbauer, Yilei Huang
- Added command line support for ancIBD (as described on readthedocs)
- Fixed bug on some setups of function vcf_to_1240K_hdf when intermediate .vcf was compressed but not flagged in path_vcf as .gz. We now choose whether to compress based on the path automatically. 

## 0.3a2, April 12th, Harald Ringbauer
- Added Filter to only existing GP data when loading (use case: HDF5 with missing data). Disregards all SNPs where at least some data is missing.

## 0.3a1, Feb 12th 2023, Harald Ringbauer
- Fix Bug of np.float and np.int being not supported in newer version of numpy.

## 0.2a2
Minor bug fixes:
- Only import ancIBD functions from package.
- Removed several unneeded import
- Added Dependency for scikit-allel

## 0.2a1
After several bug fixes and making sure that the package structure works (in particular imports), this is the first version on the official PyPi Server.

## 0.1a10, May 27th 2022, Harald Ringbauer
This is the first ever ancIBD version, on the Test PyPi Server. Happy birthday. May your IBD be long and prosper.
- Includes the first three Vignette Notebooks
- Introduces the ./package/ancIBD structure