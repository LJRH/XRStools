from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from .xrstools import calctools
from pylab import *

tmao = calctools.erkale('/home/christoph/programs/erkale-runs/tmao-6M/clusters_4A/results/cluster_mol','/trans_four-1.50.dat',0,980,20,stepformat=3)
tmao.cut_rawspecs(500,600)
tmao.broaden_lin()
tmao.sum_specs()
tmao.plot_spec()

urea = calctools.erkale('/home/christoph/programs/erkale-runs/urea-6M/clusters_4A/results/cluster_mol','/trans_four-1.50.dat',0,980,20,stepformat=3)
urea.cut_rawspecs(500,600)
urea.sum_specs()
urea.broaden_lin()


water = calctools.erkale('/home/csahle/programs/erkale-runs/tip4p-water/clusters_6A/cluster_mol','/xrs/trans_four-1.50.dat',0,3000,50,stepformat=3)
water.cut_rawspecs(500,600)
water.broaden_lin() #params=[0.8, 8, 537.5, 550]
water.sum_specs()
water.norm_area(530,560)
water.plot_spec()

ion()
plot(urea.energy,urea.sqw,tmao.energy,tmao.sqw)

# TMAO 4M, only H2O oxygens
tmao4M_OW = calctools.erkale('/home/csahle/programs/erkale-runs/tmao-4M/clusters_5A/cluster_mol','/xrs/trans_four-1.50.dat',80,980,20,stepformat=3)
tmao4M_OW.cut_rawspecs(500,600)
tmao4M_OW.broaden_lin()
tmao4M_OW.sum_specs()
tmao4M_OW.norm_area(530,560)
tmao4M_OW.plot_spec()

# TMAO 8M, only H2O oxygens
tmao8M_OW = calctools.erkale('/home/csahle/programs/erkale-runs/tmao-8M/clusters_5A/cluster_mol','/xrs/trans_four-1.50.dat',160,1120,20,stepformat=3)
tmao8M_OW.cut_rawspecs(500,600)
tmao8M_OW.broaden_lin()
tmao8M_OW.sum_specs()
tmao8M_OW.norm_area(530,560)
tmao8M_OW.plot_spec()

# UREA 4M, only H2O oxygens
urea4M_OW = calctools.erkale('/home/csahle/programs/erkale-runs/urea-4M/clusters_5A/cluster_mol','/xrs/trans_four-1.50.dat',80,1080,20,stepformat=3)
urea4M_OW.cut_rawspecs(500,600)
urea4M_OW.broaden_lin()
urea4M_OW.sum_specs()
urea4M_OW.norm_area(530,560)
urea4M_OW.plot_spec()

# UREA 8M, only H2O oxygens
urea8M_OW = calctools.erkale('/home/csahle/programs/erkale-runs/urea-8M/clusters_4A/cluster_mol','/xrs/trans_four-1.50.dat',160,1080,20,stepformat=3)
urea8M_OW.cut_rawspecs(500,600)
urea8M_OW.broaden_lin()
urea8M_OW.sum_specs()
urea8M_OW.norm_area(530,560)
urea8M_OW.plot_spec()

plot(water.energy,water.sqw,tmao4M_OW.energy,tmao4M_OW.sqw,tmao8M_OW.energy,tmao8M_OW.sqw,urea4M_OW.energy,urea4M_OW.sqw,urea8M_OW.energy,urea8M_OW.sqw)












