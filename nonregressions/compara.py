from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import h5py
import fabio
import os
from matplotlib.colors import LogNorm
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
from matplotlib import patches
from six.moves import range
from six.moves import zip
SHOW = 0

# Temperature as defined in spslut
cdict = {'red': ((0.0, 0.0, 0.0),
                 (0.5, 0.0, 0.0),
                 (0.75, 1.0, 1.0),
                 (1.0, 1.0, 1.0)),
         'green': ((0.0, 0.0, 0.0),
                   (0.25, 1.0, 1.0),
                   (0.75, 1.0, 1.0),
                   (1.0, 0.0, 0.0)),
         'blue': ((0.0, 1.0, 1.0),
                  (0.25, 1.0, 1.0),
                  (0.5, 0.0, 0.0),
                  (1.0, 0.0, 0.0))}

#Do I really need as many colors?
temperature = LinearSegmentedColormap('temperature',cdict, 65536)

filea = "demo_foilroi.h5"
fileb = "demo_newrois_scan.h5"
filer = "demo_responses_0.0_0.2_Tfreezed.h5"
filer = "demo_responses_0.0_2.0_Tfreezed.h5"


totf =  fabio.open("FIGS/cumulata.edf","r").data

h5resp = h5py.File(filer, "r")
pixelr  = h5resp["response"][:]

group_scan = "/ROI_FOIL/foil_scan/scans/Scan273/"
group_rois = "/ROI_FOIL/rois_definition/rois_dict/"
group_scanb = "/ROI_FOIL/scan_foil/scans/Scan273/"

from pylab import *

def  get_data(filea, group_scan,chiave  ):
    h5f = h5py.File(filea, "r")
    h5  = h5f[group_scan]
    data = h5[str(chiave)]["matrix"][:]
    h5f.close()
    return data

def get_rois_list(filea,group_rois ):
    h5f = h5py.File(filea, "r")
    h5  = h5f[group_rois]
    chiavi = list(h5.keys())
    chiavi.sort()
    
    origins = []
    sezioni = []
    for c in chiavi:
        h5b = h5[c]
        origins.append(h5b["origin"][:])
        if "sezioneold" in h5b:
            sezioni.append(h5b["sezioneold"][:])
        else:
            sezioni.append(None) 
    chiavi = [ int( c [-2:]) for c in chiavi    ]
    h5f.close()
    return chiavi, origins, sezioni

roiskeys_list_a,origins_list_a, tmp = get_rois_list(filea,group_rois )
roiskeys_list_b,origins_list_b,sezione_list  = get_rois_list(fileb,group_rois )

assert(roiskeys_list_a==roiskeys_list_b)

imshow(pixelr, interpolation='nearest')


if SHOW:
    show()
else:
    savefig('FIGS/PR.png', dpi=200, bbox_inches='tight')

    
def mostra(chiave,data_a, data_b , data_b_ori, sez, pointer, totf, shor, inset_b):   
    statsa =  (data_a.sum(axis=-1)).sum(axis=-1)
    statsb =  (data_b.sum(axis=-1)).sum(axis=-1)
    statsbb =  (data_b_ori.sum(axis=-1)).sum(axis=-1)

    d = "FIGS/"+str(chiave)+"/"
    if not os.path.exists(d):
        os.makedirs(d)

    fig = plt.figure()
    # ax = fig.add_subplot(111, autoscale_on=False, xlim=(-1, 5), ylim=(-4, 3))
    ax = fig.add_subplot(111)
    ax.imshow(totf, interpolation='nearest',  cmap=temperature, norm=LogNorm(vmin=1 ) ) #  vmin=0.01, vmax=1

    dim1,dim2 = totf.shape
    print( dim1,dim2)
    
    ax.annotate('Analyser '+str(chiave), xy=(ob[1], ob[0]) , xycoords='data', color = "red",
                xytext=(0,0*dim1//2), textcoords ='data',  #  textcoords='offset points',
                size=20, 
                # bbox=dict(boxstyle="round", fc="0.8"),
                arrowprops=dict(arrowstyle="->",# arrowstyle="fancy",
                                # fc="0.6", ec="none",
                                # patchB=el,
                                connectionstyle="angle,angleA=-120,angleB=0,rad=10", color="red")
                                # connectionstyle="angle3,angleA=180,angleB=-90", color = "red")
                )

    
    for thick in range(4) :
        ax.add_patch(
            patches.Rectangle(
                (ob[1]-thick, ob[0]-thick),
                shor[2]+2*thick,
                shor[1]+2*thick,
                fill=False    ,
                color = "red" # remove background
            )
        )
        
    if SHOW :
        show()
    else:
        fig.savefig(d+'roi.png', dpi=200, bbox_inches='tight')

    if (chiave == 10000) :
        for chiave_, oa_, ob_, sez_ in zip( roiskeys_list_a, origins_list_a, origins_list_b , sezione_list )  :
            
            shor_ = get_data(fileb, group_scanb,chiave_ ).shape 
         
            for thick in range(4) :
                ax.add_patch(
                    patches.Rectangle(
                        (ob_[1]-thick, ob_[0]-thick),
                        shor_[2]+2*thick,
                        shor_[1]+2*thick,
                        fill=False    ,
                        color = "red" # remove background
                    )
                )           
        fig.savefig('FIGS/tutteroi.png', dpi=200, bbox_inches='tight')

        
    fig = plt.figure()
    inset = totf[ ob[0]-3*0: ob[0]+shor[1]+3*0,  ob[1]-3*0: ob[1]+shor[2]+3*0       ]
    ax = fig.add_subplot(111)
    ax.imshow( inset, interpolation='nearest',  cmap=temperature, norm=LogNorm(vmin=1 ) ) #  vmin=0.01, vmax=1
    # for thick in range(4) :        
    #     ax.add_patch( patches.Rectangle((3-thick, 3-thick),shor[2]+2*thick,shor[1]+2*thick,fill=False,color = "red" ) )
                      
    if SHOW :
        show()
    else:
        fig.savefig(d+'insetbig.png', dpi=200, bbox_inches='tight')


    fig = plt.figure()
    inset = zeros( [ shor[1]+6 , shor[2]+6    ], "f")
    inset[ 3*0: inset_b.shape[0]+3*0, 3*0: inset_b.shape[1]+3*0   ] = inset_b
    ax = fig.add_subplot(111)
    ax.imshow(inset, interpolation='nearest',  cmap=temperature, norm=LogNorm(vmin=1 ) ) #  vmin=0.01, vmax=1
    # for thick in range(4) :
    #     ax.add_patch( patches.Rectangle((3-thick, 3-thick),shor[2]+2*thick,shor[1]+2*thick,fill=False,color = "red" ) )
    if SHOW :
        show()
    else:
        fig.savefig(d+'insetsmall.png', dpi=200, bbox_inches='tight')

    fig = plt.figure()
    inset = totf[ ob[0]-3*0: ob[0]+shor[1]+3*0,  ob[1]-3*0: ob[1]+shor[2]+3*0       ]

    spia = inset_b.sum(axis=0)
    dove = (spia>0).astype("f")
    # dove =  roll(dove,1) +  roll(dove,-1)+roll(dove,2) +  roll(dove,-2)+roll(dove,3) +  roll(dove,-3)
    
    #inset[:, dove>0  ] = inset_b[:,dove>0]
    inset[ inset_b>0  ] = inset_b[inset_b>0]
    ax = fig.add_subplot(111)
    ax.imshow(inset, interpolation='nearest',  cmap=temperature, norm=LogNorm(vmin=1 ) ) #  vmin=0.01, vmax=1
    # for thick in range(4) :
    #     ax.add_patch( patches.Rectangle((3-thick, 3-thick),shor[2]+2*thick,shor[1]+2*thick,fill=False,color = "red" ) )
    if SHOW :
        show()
    else:
        fig.savefig(d+'insetcollage.png', dpi=200, bbox_inches='tight')


    fig = plt.figure()
    
    data  = h5resp[str(chiave)]["data"][:]+1
    data[0,0]=0
    imshow(data+1.0, interpolation='nearest', norm=LogNorm(vmin=1.1 ))
    if SHOW :
        show()
    else:
        fig.savefig(d+'optical.png', dpi=200, bbox_inches='tight')
        
    x = arange(len(statsa))+sez
    xx = arange(len(statsbb))
    
    fig = plt.figure()

    plot(x,statsa )
    plot(xx,statsbb)
    # plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')
    if SHOW :    
        show()
    else:
        fig.savefig(d+'largeplot.png', dpi=200, bbox_inches='tight')
        
    fig = plt.figure()

    plot(x,statsa )
    plot(x,statsb )
    # plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')
    if SHOW :    
        show()
    else:
        fig.savefig(d+'smallplot.png', dpi=200, bbox_inches='tight')

    def mostra_massimo(bl, names):
        b=bl[0]
        k = argmax(b)
        j = k%b.shape[1]
        i = k//b.shape[1]
        
        for b, name in zip(bl, names):
            inset = zeros([7,7])

            t = b[max(0,i-3):i+4, max(0,j-3):j+4]
            
            inset[ -min(0,i-3):-min(0,i-3) +t.shape[0] , -min(0,j-3):-min(0,j-3) +t.shape[1] ] = t

            
            print( "b.shape",inset.shape)
            imshow(inset, interpolation='nearest',  cmap=temperature, norm=LogNorm(vmin=1 ) )

            if SHOW :
                show()
            else:
                fig.savefig(  name , dpi=200, bbox_inches='tight')
        
    mostra_massimo(  [data_b_ori[0]], [d+"estremo_inizio.png"] )
    mostra_massimo( [data_b[0] , data_a[0]]  , [d+"segmento_inizio.png"  , d+"segmentoExp_inizio.png"    ] )
    mostra_massimo( [data_b[-1]  , data_a[-1] ], [d+"segmento_fine.png"  , d+"segmentoExp_fine.png"    ] )
    mostra_massimo( [data_b_ori[-1] ] , [d+"estremo_fine.png"] )

    print( " faccio chiave " , chiave)

for chiave, oa, ob, sez in zip( roiskeys_list_a, origins_list_a, origins_list_b , sezione_list )  :

    # chiave = 3
    
    data_a = get_data(filea, group_scan,chiave  ) 
    data_b = get_data(fileb, group_scanb,chiave  ) 

    data_b_ori = array(data_b)

    data_b  =  data_b[sez[0]:sez[1]]


    inset_b = data_b.sum(axis=0)
                      
    c = oa[0]-ob[0], oa[1]-ob[1]
    sh = data_a.shape[1:]

    print( " CORNER A " , oa)
    print( " CORNER B " , ob)
    
    print( " SHAPE A ",  data_a.shape)
    print( " SHAPE B ",  data_b.shape)
    
    by0 = max(   oa[0],ob[0]  )
    bx0 = max(   oa[1],ob[1]  )
    
    by1 = min(   oa[0]+ data_a.shape[1]   ,ob[0] + data_b.shape[1] )
    bx1 = min(   oa[1] + data_a.shape[2]    ,ob[1] + data_b.shape[2] )

    print( " DATAB ",  max(0,by0-ob[0]) ,  by1-ob[0]      ,  max(0,bx0-ob[1] ),bx1-ob[1] )
    print( " DATA1 ",  max(0,by0-oa[0]) ,  by1-oa[0]      ,  max(0,bx0-oa[1] ),bx1-oa[1] )
                          
    data_b  =  data_b[  :  ,  max(0,by0-ob[0] ):  by1-ob[0]      ,  max(0,bx0-ob[1] ):bx1-ob[1]       ]
    data_a  =  data_a[  :  ,  max(0,by0-oa[0] ):  by1-oa[0]      ,  max(0,bx0-oa[1] ):bx1-oa[1]       ]
    
    print( sez)
    print( oa)
    print( ob)

    mostra( chiave ,data_a, data_b , data_b_ori    ,sez[0], ob , totf, data_b_ori.shape, inset_b)







