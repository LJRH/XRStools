Code Invocation
===============
  * Some of the XRStools capabilities can be accessed by invocation of the *XRS_swissknife* script, providing as input a file in the *yaml* format.
  * To use the wizard the suggested instruction is ::

         XRS_wizard  --wroot ~/software/XRStoolsSuperResolution/XRStools/WIZARD/methods/

    the wroot argument tells where extra workflow can be found. In the above instruction we give workflows in the home source directory. This is practical because the wizard allows to edit them online and the modification will remain in the sources. or to access extra workflows that are not coming with the main disribution.
  

  * Depending on the details of your installation, you have the *XRS_swissknife* script sitting somewhere in a directory. Check the *Installation* page to see how to set PYTHONPATH and PATH in case of a local installation.

    *The following documentation has been generated automatically from the comments found in the code*.


GENERALITIES about XRS_swissknife
---------------------------------
.. automodule:: XRS_swissknife
    :members: generality_doc 


Super Resolution
----------------

to fit optical responses of all the  analysers (you selected a ROI for) and the pixel response based on a foil scan
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

embedded doc :

.. automodule:: XRS_swissknife
   :members:  superR_fit_responses 


to extrapolate to a larger extent the ROIS and the foils scan, thus to cover a larger sample 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

embedded doc :

.. automodule:: XRS_swissknife
   :members:   superR_recreate_rois


to calculate the scalar product between a foil scan and a sample, for futher use in the inversion problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

embedded doc :	       

.. automodule:: XRS_swissknife
   :members:  superR_scal_deltaXimages 


Other features
--------------

.. automodule:: XRS_swissknife
    :members: help, create_rois,load_scans,Extraction, HFspectrum
e_rois   

