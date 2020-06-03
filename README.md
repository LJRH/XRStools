XRStools (with APS LERIX support)
===========================================

Main development website: https://gitlab.esrf.fr/mirone/XRStools
This version of XRStools is maintained by Luke Higgins at the University of Leeds (when I get chance) for use with X-ray Raman data from both the European Synchrotron (ESRF) and the Advanced Photon Source (APS). The original code is developed at the ESRF by Alessandro Mirone and Christoph Sahle.

The LERIX-working branch of this code has an extra class, `xrs_read.read_lerix()`, which can be edited at line  3250 in `xrs_read.py`. Any helpful edits/suggestions are welcome.

References
----------

* Ch.J. Sahle, A. Mirone, J. Niskanen, J. Inkinen, M. Krisch, and S. Huotari:Journal of Synchrotron Radiation 22, No. 2 (2015): 400-409.

Installation
------------

With Conda (recommended)
........

Conda installation is buggy but works but definitely works on Ubuntu (my main driver) and has been tested on Windows/MacOS, but is still development code (BEWARE). Tested with conda version: 4.7.12

Step 0 - Install conda -- [miniconda comes with less bulk](https://conda.io/en/latest/miniconda.html)

Step 1 - Install git and clone the XRStools repository:

```shell
    conda install git
    git clone git@github.com:LJRH/XRStools.git
    git checkout LERIX-working #this step for APS LERIX support
```

Step 2 - Change directory to the cloned XRStools repo and create a new conda virtual environment called "xrstools":

```shell
    cd <...>/XRStools/
    conda create -n xrstools python=2.7 --file requirements.txt
```

Step 3 - Install other required libraries including praxes's version of pymca (newer versions don't work yet):

```shell
    conda activate xrstools
    conda install -yc praxes pymca
    conda install -yc conda-forge silx
```

Step 4 - Intitiate XRStools (this might generate errors, ignore these.)

```shell
    python setup.py install [--user]
```

Step 5 - Open iPython and check that you are able to load important modules:

```shell
    ipython
```

```python
    from XRStools import xrs_read xrs_utilities xrs_extraction
```

With PIP
........

As most Python packages, pyFAI is available via [PIP](https://pip.pypa.io/en/stable/)::
```shell
   pip install xrstools [--user]
```
Provide the *--user* to perform an installation local to your user.
Under UNIX, you may have to run the command via *sudo* to gain root access an
perform a system wide installation.



Documentation
-------------

Documentation & an example for [using XRStools with APS LERIX can be found here](https://github.com/LJRH/XRStools/tree/LERIX-working/LERIX-Guide).
Documentation can be build using this command and Sphinx (installed on your computer)::

    python setup.py build build_doc
