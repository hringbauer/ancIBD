Installation
===============

You can find the software package ``ancIBD`` on the official pypi Python package repository. 

Therefore you can install it using ``pip``,
::
    python3 -m pip install ancIBD

The package distributes source code that is compiled during installation - and pip installs ``ancIBD`` automatically. 

For experts: The ``setup.py`` contains the relevant information used by ``pip`` for the installation, and one can also compile the relevant C code using ``Cython`` manually.



Upgrading    
************
If you already have a ``ancIBD`` release installed via pip and wish to upgrade to the latest stable release, you can do so by adding ``--upgrade``:
::
    pip install --upgrade ancIBD
    
c Extension
************
For performance reasons, the heavy lifting of the algorithms is coded into c methods (``cfunc.c``). The package is set up so that this "extension" is built during installation. This is done automatically from ``cfunc.pyx`` via the package cython (when ``CYTHON=True`` in setup.py, the default setting). You can also set ``CYTHON=False``, then the extension is directly compiled from ``cfunc.c`` (experimental, not tested on all platforms).


Dependencies
************
The basic dependencies are sufficient for the core functions of the algorithms and are kept low to avoid package dependency conflicts. When ``ancIBD`` is installed, the following dependent Python packages should be automatically installed without any action on your part:

* ``numpy`` for calculaions with numerical arrays at C speed 
* ``pandas`` for handling databases and tables at C speed 
* ``h5py`` for handling hdf5, a file format with partial I/O
* ``psutil`` for process monitoring
