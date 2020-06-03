Author: L.Higgins
Date: 16 Jan 20

-------------------------------------------------------------------------------
## Introduction

LERIX (Low Energy Resolution Inelastic X-ray) spectrometer is the instrument in operation at [beamline 20-ID](https://www.aps.anl.gov/Spectroscopy/Beamlines/20-ID) of the Advanced Photon Source (APS). The data collected from X-ray Raman scattering spectroscopy (XRSS) experiments at the APS has a very different format (Labview prodcing ASCII text docs) to the more data-efficient ESRF "spec" format. Therefore, in order to use the feature-rich XRStools software, a small add-on (read_lerix) was made. This guide should help those who want to analyze their data using this tool. It is is development software, made by a very amateur coder in python - it will be buggy. Any help from those with more skill is welcome.

-------------------------------------------------------------------------------
## Getting Started - an example

Included in the LERIX-guide directory is some test data of a graphite pellet on a spinner from 20-ID. The experiment was performed using the Si(311) mono at 9890 eV with XRS data collected for the C K edge.

### Loading relevant modules
If XRStools has been installed properly then on opening an iPython command line the following modules should be able to be called:

```python
    from XRStools import xrs_read xrs_utilities xrs_extraction
    import numpy as np
    import matplotlib.pyplot as plt
    import os
```

If at this stage there are errors such as 'conda: module not found', try and install the missing modules. The next step would be to rebuild XRStools `python setup.py install`.

If successful you should see the following output:

```
>>>>>>>>  use_PyMca  True
>>>>>>>>  use_PyMca  True
/opt/conda/envs/christoph6/lib/python2.7/site-packages/PyMca5/PyMcaGui/plotting/MaskImageWidget.pyc
Could not load PyTango
>>>>>>>>  use_PyMca  True

############################# Welcome to XRStools #############################
# If you are using this software, please cite the following work:             #
# Ch.J. Sahle, A. Mirone, J. Niskanen, J. Inkinen, M. Krisch, and S. Huotari: #
# "Planning, performing, and analyzing X-ray Raman scattering experiments."   #
# Journal of Synchrotron Radiation 22, No. 2 (2015): 400-409.                 #
###############################################################################
```
### Setting the Experiment Directory
At this point we set the directory of the test data and read the experiment:

```python
path = '../XRStools/LERIX-guide/graphite/' #change this to the test data
graphite = xrs_read.read_lerix(exp_dir=path)
```

Once the data is loaded we can check to see which scans have been loaded into the class:

```python
graphite.elastic_scans
graphite.nixs_scans
graphite.wide_scans
```
### Loading Elastic Scans and Checking FWHM & E<sub>0</sub>
Once we know the scans have been loaded, we can load the elastic scans:

```python
graphite.load_elastics(analyzers='all')
```

Now that the elastic scans are loaded, we can access the cenom (centre-of-mass) dictionary which contains the FWHM and centre (E0) of the elastic peaks. This is performed using the `xrs_utilities.fwhm()` function and fits a gaussian function to the elastic scans. The cenom dictionary can be accessed in the following ways:

```python
graphite.cenom_dict
# # or
graphite.cenom_dict['Analyzer03']
# # or
graphite.cenom_dict['Analyzer03']['elastic0001']
# # or
graphite.cenom_dict['Analyzer03']['elastic0001']['fwhm'] # or ['e0']
```

### Loading NIXS and Wide scans, then Joining them Together
At this point we can begin to load the NIXS and WIDE scans. N.B. wide scans are important as they allow proper subtraction of any Compton background present in the XRSS scan.

```python
graphite.load_nixs(exp_dir=path)
graphite.load_wides(exp_dir=path)
```

At this stage, the NIXS and wide data should be merged together. Sometimes though, I only perform this step before background subtraction, since scans are easier to compare without the Compton background.

```python
graphite.join_nixs_wide(scaling='auto')
```
### Plotting the Analyzers and selecting best q's


-------------------------------------------------------------------------------
# Dictionary of 'read_lerix' functions
## read_lerix()

Variable      | Type   |  Description
------------- | ------ | -------------
exp_dir       | String | Unix path to location of directory containing nrixs, elastic and wide scans.
energy_column | string | `default='energy'`, a part of the header name of the ASCII column holding the energy data.
monitorcolumn | string | `default= "__ enc"`, a part of the header name of the ASCII column holding the monitor (I0) data. (Note M.B reccomends using the encoded channel.)
NIXS_name     | string | `default= "nixs"`, name of the files containing nixs data.
elastic_name  | string | `default= "elastic"`, name of the files containing elastic data.
wide_name  | string |`default= "wide"`, name of the files containing wide data.

## read_lerix.load_elastics(self,exp_dir=None,scans='all',analyzers='all')
Module to load the `elastic_scans` and analyse them to produce the average E0 and FWHM for each analyzer/crystal pair. Rejects scans that do not meet 500 counts per seconda or have more than 5000 cps, as this is the range of linearity of the detectors. Loads all the elastic scans in the experiment directory by default.

Variable      | Type   |  Description
------------- | ------ | -------------
exp_dir       | String | `default=None`; typically the elastic scans are in the main directory, can be changed to load elastic scans from another directory.
scans         | String/List | default='all'; It is possible to load a list of scans e.g. [1,2,5] where 1 is scan .0001.
analyzers     | String/List | `default='all'`; Usually imports all the elastic scans, this is important because they are typically collected with different filters to protect the detectors. However a list of elastic scans to use can be given e.g. [1,2,5]

## read_lerix.update_cenom(self, analyzers="all")
Updates the `cenom_dict` for a set of selected analyzer crystals e.g. `[0,2,18]`. Before performing Compton subtraction, this should be run with "all" otherwise the subtraction will produce an error.

Variable      | Type   |  Description
------------- | ------ | -------------
analyzers     | String/list | `default='all'`; Chosen analyzers could be a list e.g. `[0,2,18]`.


## read_lerix.load_nixs(self,exp_dir=None,scans='all',analyzers='all')
Module to load the `nixs_scans` into XRStools.

Variable      | Type   |  Description
------------- | ------ | -------------
exp_dir       | String | `default=None`; typically the NIXS scans are in the main directory, can be changed to load NIXS scans from another directory.
scans         | String/List | `default='all'`; It is possible to load a list of scans e.g. [1,2,5] where 1 is scan .0001.
analyzers     | String/List | `default='all'`; It is possible to load a list of scans e.g. [1,2,5] where 1 is scan .0001.

## read_lerix.load_wides(self,exp_dir=None,scans='all',analyzers='all',join=False)
Module to load the `wide_scans` into XRStools.

Variable      | Type   |  Description
------------- | ------ | -------------
exp_dir       | String | `default=None`; typically the NIXS scans are in the main directory, can be changed to load NIXS scans from another directory.
scans         | String/List | `default='all'`; It is possible to load a list of scans e.g. [1,2,5] where 1 is scan .0001.
analyzers     | String/List | `default='all'`; It is possible to load a list of scans e.g. [1,2,5] where 1 is scan .0001.

## read_lerix.join_nixs_wide(self,scaling='auto')
Joins together the loaded nixs_scans and wide_scans so that the compton (below edge) AND the nixs scans become one. This is important for compton subtraction.  TO DO: works but is ugly and can be shortened with a sensible loop.

Variable      | Type   |  Description
------------- | ------ | -------------
scaling       | String/Float | `default='auto'`; scales wide to nixs data can be `auto`, `none` or a `Float`. If a float is given the Wide part of the data is multiplied by the float in order to avoid 'gaps' in intensity at the join between the wide and nixs data.

## read_lerix.scan_info(self, f):
Gets the scan number, name, type and file extention for a filename or path. The scan assumes typical format e.g. elastic.0001, nixs.0001
returns:
[0] -> scan number (e.g. 0001)
[1] -> scan name (e.g. nixs0001)
[2] -> scan_type (e.g. nixs)
[3] -> file name (e.g. nixs.0001)"""

Variable      | Type   |  Description
------------- | ------ | -------------
f             | String | Filename or path to a scan, returns useful information.

## read_lerix.save_H5(self,H5name='20ID_APS_data.H5')
Saves all the loaded data as a H5 file. Works but could do with some refinement, e.g. saving the background subtracted data as well.

Variable      | Type   |  Description
------------- | ------ | -------------
H5name        | String | `default='20ID_APS_data.H5'`; Automatically saves the H5 file in the `self.path` location and appends the H5name as the filename.

## read_lerix.plot_data(self,analyzer=False)
Plots a GUI with radio-buttons able to show the data for each analyzer, can be a quick tool for choosing the best scans/analysers to use - but is not very well developed and uses pyqt=4.

Variable      | Type   |  Description
------------- | ------ | -------------
Analyzer      | String | `default='False'`; Tool to plot each analyzer individually, choose analyzer as 'Analyzer01' etc.
