Installation
===============

The software package ``ancIBD`` is on the official Python package repository (pypi). Therefore you can install it using ``pip``,
::
    python3 -m pip install ancIBD

This software package distributes source code that is compiled during installation - and pip installs ``ancIBD`` automatically. 

Note for experts: The ``setup.py`` contains the relevant information used by ``pip`` for the installation, and one can also compile the relevant C code using ``Cython`` manually.



Upgrading    
************
If you already have a ``ancIBD`` release installed via pip and wish to upgrade to the latest stable release, you can do so by adding ``--upgrade``:
::
    pip install --upgrade ancIBD
    
c Extension
************
For performance reasons, the heavy lifting of the algorithms is coded into c methods (``cfunc.c``). This "extension" is built automatically during installation from ``cfunc.pyx`` via the package cython (when ``CYTHON=True`` in setup.py, the default setting). If you set ``CYTHON=False``, then the extension is directly compiled from ``cfunc.c`` (experimental, not tested on all platforms).


Dependencies
************
The basic Python package dependencies are sufficient for the core functions of  ``ancIBD``. We kept the strictly required dependencies low to avoid creating dependency conflicts. When ``ancIBD`` is installed, the following dependent Python packages should be automatically installed without any action on your part:

* ``numpy`` for calculations with numerical arrays at C speed 
* ``pandas`` for handling databases and tables at C speed 
* ``h5py`` for handling hdf5, a file format with partial I/O
* ``psutil`` for process monitoring

Some downstream and advanced functionalities require additional packages, e.g. ``matplotlib` for plotting. You can install them manually via pip, in case you are missing those you will be alerted by import errors.
