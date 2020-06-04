from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from .Wizard import *
from . import Wizard 

# from .methods.RIXS_spectra import winfo_RIXS_spectra_extraction_deconvolve 
# from .methods.RIXS_spectra import winfo_RIXS_spectra_extraction_preparation 
# from .methods.RIXS_spectra import winfo_RIXS_spectra_extraction 
# from .methods.predictions import winfo_xrsprediction

## from .methods.predictions   winfo_response_denoiser import metodo as metodo_reponse_denoiser
installation_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),"..")
animation_filename =  os.path.join(installation_dir,'data/shadok_eb-003.gif')


def main():
    import sys, getopt
    
    usage = "USAGE : XRS_wizard --shift_the_reference <0/1>(defaults 1) --do_deconvolution <0/1>(defaults 0) --wroot  <extra_w_path>"

    shift_the_reference = 1
    do_deconvolution = 0

    try:
        opts, args = getopt.getopt(sys.argv[1:],"h",["shift_the_reference=","do_deconvolution=","wroot="])
    except getopt.GetoptError:
        print( usage)
        sys.exit(2)
    print( opts, args)
    extrawpaths=[]
    for opt, arg in opts:
        if opt == '-h':
            print( usage)
            sys.exit()
        elif opt in ("--shift_the_reference"):
            shift_the_reference = int(arg)
        elif opt in ("--do_deconvolution"):
            do_deconvolution = int(arg)
        elif opt in ("--wroot"):
            extrawpaths.append(arg)



    options = { "shift_the_reference" :shift_the_reference, "do_deconvolution":do_deconvolution         } 

#     metodo_RIXS_spectra_extraction_preparation = winfo_RIXS_spectra_extraction_preparation.getMethod(options)
#     metodo_RIXS_spectra_extraction = winfo_RIXS_spectra_extraction.getMethod(options)
#     metodo_xrsprediction = winfo_xrsprediction.getMethod(options)


    import collections

    metodi = collections.OrderedDict()
    
#     metodi["RIXS_spectra_extraction_preparation"] = metodo_RIXS_spectra_extraction_preparation
#     metodi["RIXS_spectra_extraction"] = metodo_RIXS_spectra_extraction
#     metodi["XRS Prediction"] = metodo_xrsprediction

    
    print( animation_filename)
    wizardMain(animation_filename, metodi, extrawpaths, options)    
        

if __name__=="__main__":
    main()







