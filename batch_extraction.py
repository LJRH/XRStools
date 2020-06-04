# Hi Alessandro,
# here's a file with 1. column ROI# and 2. column E0 (in keV).
# There were in total two samples we tried to measure at many energies (I also attached the logbook for clearity):
# 1) this is the little insect trapped in amber
# Scans: #440 - #589, the scheme was as follows:
# - scans 440 - 464 were around the elastic line, 2 scans per energy (i.e. 2 sty scans for two different values of stz at each energy point), the energy grid was 0.5 eV
# - scans 464 - 589 were around the C K-edge, 2 scans per energy (same as for elastic)
# 2) this was a little piece of burnt prehistoric wood
# Scans: #1047 - #1141, the scheme is the same as above, but 0.25 eV energy steps:
# - scans 1047 - 1195: scans around the elastic line (2 sty scans per energy at 2 different stz positions)
# - scans 1095 - 1141, scans around the C K-edge (2 sty scans per energy at 2 different stz positions)
# the reference scans (25 micron Kapton foil) are:
# # 591 E = 9.686
# # 592 E = 9.9712
# # 593 E = 9.9745
# # 594 E = 9.9776
# # 595 E = 9.986
# # 596 E = 9.956
# # 597 E = 9.946
# # 598 E = 10.006
# # 599 E = 10.026
# let me know if all this makes sense (or not). Thanks a lot!
# Christoph
import numpy as np
import h5py
import glob

cenom = np.array(
    [ 0,   9.68715366 ,
      1,   9.68646557 ,
      2,   9.68662118 ,
      3,   9.68626895 ,
      4,   9.68712335 ,
      5,   9.68650054 ,
      6,   9.68603224 ,
      7,   9.68585970 ,
      8,   9.68648548 ,
      9,   9.68587371 ,
     10,   9.68528659 ,
     11,   9.68522053 ,
     12,   9.68575133 ,
     13,   9.68453831 ,
     14,   9.68495099 ,
     15,   9.68573933 ,
     16,   9.68612397 ,
     17,   9.68537876 ,
     18,   9.68550170 ,
     19,   9.68629621 ,
     20,   9.68686062 ,
     21,   9.68620125 ,
     22,   9.68611576 ,
     23,   9.68710568 ,
     24,   9.68663625 ,
     25,   9.68615702 ,
     26,   9.68634820 ,
     27,   9.68703000 ,
     28,   9.68593704 ,
     29,   9.68544516 ,
     30,   9.68526687 ,
     31,   9.68596309 ,
     32,   9.68507243 ,
     33,   9.68492916 ,
     34,   9.68469521 ,
     35,   9.68589879 ,
     36,   9.68620847 ,
     37,   9.68638990 ,
     38,   9.68647063 ,
     39,   9.68743229 ,
     40,   9.68610182 ,
     41,   9.68520195 ,
     42,   9.68558591 ,
     43,   9.68649882 ,
     44,   9.68551274 ,
     45,   9.68487196 ,
     46,   9.68511143 ,
     47,   9.68583333 ,
     48,   9.68539694 ,
     49,   9.68503707 ,
     50,   9.68508553 ,
     51,   9.68637998 ,
     52,   9.68646747 ,
     53,   9.68580797 ,
     54,   9.68568058 ,
     55,   9.68629252 ,
     56,   9.68729703 ,
     57,   9.68647016 ,
     58,   9.68634066 ,
     59,   9.68693394 ,
     60,   9.68626225 ,
     61,   9.68617870 ,
     62,   9.68654891 ,
     63,   9.68712219 ,
     64,   9.68598538 ,
     65,   9.68545723 ,
     66,   9.68553384 ,
     67,   9.68648591 ,
     68,   9.68543486 ,
     69,   9.68492219 ,
     70,   9.68508252 ,
     71,   9.68554590  ])
cenom=np.reshape(cenom,[-1,2])
Enominal = np.median(cenom[:,1])
cenom[:,1] -= Enominal






import os

def process_input(s, go=0):
    open("input_tmp_%d.par"%go, "w").write(s)
    if not(go%12):
        os.system("XRS_swissknife  input_tmp_%d.par "%go)
    else:
        os.system("XRS_swissknife  input_tmp_%d.par &"%go)
        

def select_rois(roi_scan_num=592):
    
    inputstring = """
    create_rois:
        expdata : /data/id20/inhouse/data/run5_17/run7_ihr/hydra 
        scans : [{roi_scan_num}] 
        roiaddress : roi_{roi_scan_num}.h5:/extracted/ROI_AS_SELECTED 
        filter_path : mask.h5:/FILTER_MASK/filter 
"""
    s=inputstring.format(roi_scan_num = roi_scan_num )
    process_input(s)

def extract_sample_givenrois(roi_scan_num=592, nick_name="insect_ck", target_file ="signals.h5"  , Start = 464, End = 589, Thickness = 2  ):

    inputstring = """
    loadscan_2Dimages :
      expdata : /data/id20/inhouse/data/run5_17/run7_ihr/hydra 
      roiaddress : roi_{roi_scan_num}.h5:/extracted/ROI_AS_SELECTED 
      monitor_column : izero/0.000001 
      scan_interval : [{start}, {end}] 
      energy_column : sty 
      signaladdress : {target_file}:/{where}/_{start}_{end} 

      sumto1D  : 0
      monitorcolumn : izero/0.000001
    """
  
    for start in range(Start,End, Thickness):    
        s=inputstring.format(start=str(start), end=str(start+2) , where= nick_name ,roi_scan_num = roi_scan_num, target_file =target_file  )
        process_input(s)

  

class InterpInfo:
    def __init__(self, cenom, interp_file, source, target):
        volum_list = list(interp_file[source].keys())
        scan_num_list = np.array([ int( t.split("_") [1]) for t in volum_list])
        ene_list      = np.array([    interp_file[source][vn]["scans"]["Scan%d"%sn ]["motorDict"]["energy"].value                   for vn,sn in     zip(volum_list, scan_num_list   )   ])

        
        print ( " ecco la scannumlist " , scan_num_list)
        print (" ecco ene_list", ene_list)
        
        
        self.volum_list    =  volum_list
        self.scan_num_list =  scan_num_list
        self.ene_list      =  ene_list

        order = np.argsort(  self.ene_list    )
        
        self.ene_list  = self.ene_list [order]
        self.scan_num_list  = self.scan_num_list [order]
        self.volum_list  = [ self.volum_list [ii]  for ii in order  ] 
        
        self.interp_file=interp_file
        self.source= source
        self.target = target 
        self.cenom=cenom
        
    def interpola(self):
        # print ( " ECCO I DATI ")
        # print (  self.ene_list  ) 
        # print (  self.cenom   )
        # raise
        for t_vn, t_sn, t_ene in zip(self.volum_list,  self.scan_num_list, self.ene_list    ):
            rois_coeffs={}
            for roi_num, de in enumerate(     self.cenom   ):
                if  t_ene+de < self.ene_list .min() or t_ene+de > self.ene_list .max():
                    continue

                print ( t_ene+de, self.ene_list .min() ,self.ene_list .max() )
                
                i0 = np.searchsorted(   self.ene_list    , t_ene+de )-1
                assert(i0>=0)
                i1=i0+1
                print (i0, i1, len(self.ene_list))
                print (self.ene_list) 
                assert(i1<len( self.ene_list ))

                DE = (  self.ene_list[i1] -  self.ene_list[i0]   )
                df = (  t_ene+de  -  self.ene_list[i0]   )/ DE
                
                rois_coeffs[ roi_num  ] =  [   i0,(1-df)   , i1,df         ]
            print ( " for reinterpolation of ", t_vn ," interpolation scheme is the following ",  rois_coeffs)


            self.interp_file.require_group(self.target)
            
            if (self.target+"/"+t_vn) in self.interp_file:
                pass
                ## del  self.interp_file[ self.target+"/"+t_vn   ]
            else:
                self.interp_file.copy(self.source+"/"+t_vn, self.interp_file, name =  self.target+"/"+t_vn  )


            fscans = self.interp_file[ self.target+"/"+t_vn   ]["scans"]
            keys_list = list(  fscans.keys() )
            print ( " keylist ", keys_list)
            for k in keys_list:
                if k[:3]=="ROI":
                    if int(k[3:]) not in rois_coeffs:
                        print (" rimuovo ", k)
                        del fscans[k]
            for sn in range(t_sn, t_sn+2):
                fScan = fscans["Scan%d"% sn]
                keys_list = list(  fScan.keys() )
                for k in keys_list:
                    if k!="motorDict":
                        if int(k) not in rois_coeffs:
                            print (" rimuovo da scans", k)
                            del fScan[k]


            # TESTED till here
                            
            for sn in range(t_sn, t_sn+2):
                fScan = fscans["Scan%d"% sn]
                keys_list = list(  fScan.keys() )
                for k in keys_list:
                    if k!="motorDict":
                        assert( int(k)  in rois_coeffs)
                        k = int(k)
                        i0,f0,i1,f1 = rois_coeffs[k]

                        matrix0 = self.interp_file[self.source][self.volum_list[i0]  ]["scans"]["Scan%d"%( self.scan_num_list[i0]+sn-t_sn)  ][str(k)]["matrix"][:]
                        matrix1 = self.interp_file[self.source][self.volum_list[i1]  ]["scans"]["Scan%d"%( self.scan_num_list[i1]+sn-t_sn)  ][str(k)]["matrix"][:]
                        monitor = np.ones( matrix0.shape[0],"f" )
                        newmatrix = f0* matrix0+f1*matrix1
                        
                        if "matrix" in fScan[str(k)] :
                            del fScan[str(k)]["matrix"]
                        if "monitor" in fScan[str(k)] :
                            del fScan[str(k)]["monitor"]
                        if "monitor_divider" in fScan[str(k)] :
                            del fScan[str(k)]["monitor_divider"]
                            
                        fScan[str(k)]["matrix"] = newmatrix
                        fScan[str(k)]["monitor"] = monitor
                        fScan[str(k)]["monitor_divider"] = 1.0
                        

def get_reference(roi_scan_num=592):
    inputstring = """
    loadscan_2Dimages :
       expdata : /data/id20/inhouse/data/run5_17/run7_ihr/hydra 
       roiaddress : roi_{roi_scan_num}.h5:/extracted/ROI_AS_SELECTED 
       monitor_column : izero/0.000001 
       scan_interval : [{roi_scan_num},{roi_scan_num_plus1} ] 
       signaladdress : calibration_scan 
       isolateSpot : 7 
  
       sumto1D  : 0
       energycolumn : 'stx'
       monitorcolumn : izero/0.000001
    """

    s=inputstring.format(  roi_scan_num= roi_scan_num, roi_scan_num_plus1=  roi_scan_num+1  )
    process_input( s ) 



def get_scalars( signals_file = "signals.h5"  , Start = 464,  Thickness = 2 , roi_scan_num=592 , nick = "reinterp_insect_ck"):
    inputstring = """
    superR_scal_deltaXimages :
       sample_address : {signals_file}:/{nick}/_{start}_{end}/scans
       delta_address : roi_{roi_scan_num}.h5:/extracted/ROI_AS_SELECTED/calibration_scan/scans/Scan{roi_scan_num}
       nbin : 5
       optional_solution : 
       target_address : volumes.h5:/{nick}/_{start}_{end}/scal_prods
    """
    
    s=inputstring.format(start=Start, end=Start+Thickness  , roi_scan_num = roi_scan_num,   nick=nick , signals_file = signals_file )
    process_input(s)
    
    # inputstring = """
    # superR_scal_deltaXimages :
    #    sample_address : {signals_file}:/{nick}/_{start}_{end}/scans
    #    delta_address : roi_{roi_scan_num}.h5:/extracted/ROI_AS_SELECTED/calibration_scan/scans/Scan{roi_scan_num}
    #    nbin : 1
    #    optional_solution : 
    #    roi_keys  : [14,33]
    #    target_address : volumes.h5:/{nick}_14_33/_{start}_{end}/scal_prods
 
    # """
    # s=inputstring.format(start=Start, end=Start+Thickness  , roi_scan_num = roi_scan_num,   nick=nick , signals_file = signals_file )
    # process_input(s)
    
    # inputstring = """
    # superR_scal_deltaXimages :
    #    sample_address : {signals_file}:/{nick}/_{start}_{end}/scans
    #    delta_address : roi_{roi_scan_num}.h5:/extracted/ROI_AS_SELECTED/calibration_scan/scans/Scan{roi_scan_num}
    #    nbin : 1
    #    optional_solution : 
    #    roi_keys  : [67, 58]
    #    target_address : volumes.h5:/{nick}_67_58/_{start}_{end}/scal_prods
 
    # """
    # s=inputstring.format(start=Start, end=Start+Thickness  , roi_scan_num = roi_scan_num,   nick=nick , signals_file = signals_file )
    # process_input(s)

    # inputstring = """
    # superR_scal_deltaXimages :
    #    sample_address : {signals_file}:/{nick}/_{start}_{end}/scans
    #    delta_address : roi_{roi_scan_num}.h5:/extracted/ROI_AS_SELECTED/calibration_scan/scans/Scan{roi_scan_num}
    #    nbin : 1
    #    optional_solution : 
    #    roi_keys  : [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    #    target_address : volumes.h5:/{nick}_0_12/_{start}_{end}/scal_prods
 
    # """
    # s=inputstring.format(start=Start, end=Start+Thickness  , roi_scan_num = roi_scan_num,   nick=nick , signals_file = signals_file )
    # process_input(s)
   
def get_volume(nick = "reinterp_insect_ck",  Start = 464,  Thickness = 2, go=[1] ):
    inputstring = """    
     superR_getVolume :
        scalprods_address : volumes.h5:/{nick}/_{start}_{end}/scal_prods
        target_address : volumes_{go}.h5:/{nick}/_{start}_{end}/volume
        niter : 400
        beta : 6e-06
        eps : 2e-06
        debin : [1, 1]
    """
    s=inputstring.format(start=Start, end=Start+Thickness,   nick=nick, go=go[0] )
    print ( " INPUT ", s)
    process_input(s, go[0])
    go[0] = (go[0]+1)
    
    
def reshuffle(   volumefile  = "volumes.h5",   nick = "reinterp_insect_ck"     ):

    h5file_root = h5py.File( volumefile ,"r+" )
    h5file = h5file_root[nick]
    scankeys = list( h5file.keys())
    scankeys.sort()
    volumes = []
    for k in scankeys:
        if k[:1]!="_":
            continue
        print( k)
        if "volume" in h5file[k]:
            volumes.append( h5file[k]["volume"]  )
    # volume = np.concatenate(volumes,axis=0)
    volume = np.array(volumes)
    if "concatenated_volume" in h5file:
        del h5file["concatenated_volume"]
    h5file["concatenated_volume"]=volume
    h5file_root.close()
    

## THE FOLLOWING PART IS THE RELEVANT ONE

    
if(0):   # ROI selection and reference scan
    select_rois(roi_scan_num=592)



if(0): # SAMPLE extraction    
    extract_sample_givenrois(roi_scan_num=592, nick_name="insect_ck", target_file ="signals.h5"  , Start = 464, End = 588, Thickness = 2  )
    # extract_sample_givenrois(roi_scan_num=592, nick_name="wood_ck", target_file ="signals.h5"  , Start =1095 , End =1141 , Thickness = 2  )
    # extract_sample_givenrois(roi_scan_num=592, nick_name="insect_elastic", target_file ="signals.h5"  , Start = 442, End = 464, Thickness = 2  )
    # extract_sample_givenrois(roi_scan_num=592, nick_name="wood_elastic", target_file ="signals.h5"  , Start = 1047, End = 1095, Thickness = 2  )



if(0):    # INTERPOLATION 
    # interp_file = h5py.File("signals_interpolated.h5","r+")
    interp_file = h5py.File("signals.h5","r+")
    i_info = InterpInfo(  cenom[:,1] , interp_file, "insect_ck" , "reinterp_insect_ck" )
    i_info.interpola()
    interp_file.close()

if(0):  # of course we need the REFERENCE SCAN
    get_reference(roi_scan_num=592)

    
if(0):    ## The scala products, which define the equation to invert
    
    for start in range(464,588,2):
        get_scalars( signals_file = "signals.h5"  , Start = start,  Thickness = 2 , roi_scan_num=592 ,nick="reinterp_insect_ck"  )
    # for start in range(1095,1141,2):
    #     get_scalars( signals_file = "signals.h5"  , Start = start,  Thickness = 2 , roi_scan_num=592 ,nick="reinterp_wood_ck"  )

        
if(0):  # inversion of the equations  
    # for start in range(464,589,2):
    #     get_volume(nick = "reinterp_insect_ck",  Start = start,  Thickness = 2)
    #     get_volume(nick = "reinterp_insect_ck_67_58",  Start = start,  Thickness = 2)
    #     get_volume(nick = "reinterp_insect_ck_14_33",  Start = start,  Thickness = 2)
    # for start in range(1095,1141,2):
    for start in range(464,588,2):
        #get_volume(nick = "reinterp_insect_ck_67_58",  Start = start,  Thickness = 2)
        #get_volume(nick = "reinterp_insect_ck_14_33",  Start = start,  Thickness = 2)
        #get_volume(nick = "reinterp_insect_ck_0_12",  Start = start,  Thickness = 2)
        get_volume(nick = "reinterp_insect_ck",  Start = start,  Thickness = 2)

if(1):
    fl = glob.glob("volumes_*.h5")
    target = h5py.File("volumes.h5","r+"     )
    for  fn in fl:
        source =  h5py.File( fn ,"r"     )
        keylist = list(  source.keys() )
        
        for k in keylist:
            keylist2 = list(  source[k].keys() )
            for k2 in keylist2:
                print(" copiando ", k,k2, " da ", fn)
                if k2 +"/volume"   in target[k]:
                    del target[k][k2 +"/volume" ]
                source[k].copy( k2+"/volume" , target[k], name =  k2 +"/volume"  )
 

    
        
        
if(1): # putting everything in a 4D volume  Dimensions : ispectra,z,y,x
    
    reshuffle( volumefile  = "volumes.h5",   nick = "reinterp_insect_ck"    )
    #reshuffle( volumefile  = "volumes.h5",   nick = "reinterp_insect_ck_67_58"    )
    #reshuffle( volumefile  = "volumes.h5",   nick = "reinterp_insect_ck_14_33"    )


        
