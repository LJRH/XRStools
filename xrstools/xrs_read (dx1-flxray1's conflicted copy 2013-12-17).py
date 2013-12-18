#!/usr/bin/python
# Filename: xrs_read.py

from helpers import * 

# should have classes id20_read, p01_read, lerix_read, etc. 

import os
import numpy as np
import array as arr
from itertools import groupby
from scipy.integrate import trapz
from scipy.interpolate import interp1d
from pylab import *
from scipy import signal
from scipy.ndimage import measurements

from matplotlib.widgets import Cursor
import matplotlib.pyplot as plt

__metaclass__ = type # new style classes

class read_id20:
	"""
	a class for loading data form SPEC files and the according edf files from the old id16
	"""
	def __init__(self,absfilename,energycolumn='energy_cc',monitorcolumn='monitor'):
		self.scans         = {} # was a dictionary before
		if not os.path.isfile(absfilename):
			raise Exception('IOError! No such file, please check filename.')
		self.path          = os.path.split(absfilename)[0] + '/'
		self.filename      = os.path.split(absfilename)[1]
		self.scannumbers   = []
		self.EDF_PREFIXh   = 'edf/h_'
		self.EDF_PREFIXv   = 'edf/v_'
		self.EDF_POSTFIX   = '.edf'
		self.DET_PIXEL_NUMx = 768
		self.DET_PIXEL_NUMy = 256

		# which column in the SPEC file to be used for the energy and monitor
		self.encolumn      = energycolumn.lower()
		self.monicolumn    = monitorcolumn.lower()
		# here are the attributes of the old rawdata class
		self.eloss    = []	# common eloss scale for all analyzers
		self.energy   = []	# common energy scale for all analyzers
		self.signals  = []	# signals for all analyzers
		self.errors   = []	# poisson errors
		self.qvalues  = []	# for all analyzers
		self.groups   = {}      # groups (such as 2 'elastic', or 5 'edge1', etc.)
		self.cenom    = []
		self.E0       = []
		self.tth      = []
		self.resolution = []    # list of FWHM of the elastic lines for each analyzer
		# some new attributes, specific to id20
		# raw data signals sorted by analyzer module:
		self.VDsignals = []
		self.VUsignals = []
		self.VBsignals = []
		self.HRsignals = []
		self.HLsignals = []
		self.HBsignals = []
		# raw data errors sorted by analyzer module:
		self.VDerrors = []
		self.VUerrors = []
		self.VBerrors = []
		self.HRerrors = []
		self.HLerrors = []
		self.HBerrors = []
		# q-values sorted by analyzer module:
		self.VDqvalues = []
		self.VUqvalues = []
		self.VBqvalues = []
		self.HRqvalues = []
		self.HLqvalues = []
		self.HBqvalues = []
		# tth-values sorted by analyzer module:
		self.VDtth = []
		self.VUtth = []
		self.VBtth = []
		self.HRtth = []
		self.HLtth = []
		self.HBtth = []
		# TTH offsets from center of V and H modules
		# tth offsets in one direction (horizontal for V-boxes, vertical for H-boxes)
		self.TTH_OFFSETS1   = np.array([5.0, 0.0, -5.0, 5.0, 0.0, -5.0, 5.0, 0.0, -5.0, 5.0, 0.0, -5.0])
		# tth offsets in the other direction (horizontal for H-boxes, vertical for V-boxes)
		self.TTH_OFFSETS2   = np.array([-9.71, -9.75, -9.71, -3.24, -3.25, -3.24, 3.24, 3.25, 3.24, 9.71, 9.75, 9.71]) 
		# input
		self.rois     = []	# can be set once the rois are defined to integrate scans directly

	def readscan(self,scannumber):
		"""
		returns the data, motors, counter-names, and edf-files from the calss's specfile of "scannumber"
		should use PyMca's routines to read the Spec- and edf-files
		"""
		if self.path: # if a path was provided upon class construction, use this one
			path = self.path
		else: # if no path was provided
			print 'please provide a path to the SPEC-file'
			return
		if not self.filename:
			print 'please provide a SPEC-file name!'
			return
		else:
			fn = path + self.filename
		# try loading variables from file
		print 'parsing edf- and SPEC-files of scan No. %s' % scannumber
		data, motors, counters = specread(fn,scannumber)
		edfmatsh = np.array(np.zeros((len(counters['ccdno']),self.DET_PIXEL_NUMy,self.DET_PIXEL_NUMx)))
		edfmatsv = np.array(np.zeros((len(counters['ccdno']),self.DET_PIXEL_NUMy,self.DET_PIXEL_NUMx)))

		edfmats  = np.array(np.zeros((len(counters['ccdno']),self.DET_PIXEL_NUMy*2,self.DET_PIXEL_NUMx)))
		
		for m in range(len(counters['ccdno'])):
			ccdnumber = counters['ccdno'][m]
			edfnameh   = path + self.EDF_PREFIXh + self.filename + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
			edfnamev   = path + self.EDF_PREFIXv + self.filename + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
			edfmatsh[m,:,:] = (edfread_test(edfnameh))
			edfmatsv[m,:,:] = (edfread_test(edfnamev))

			edfmats[m,0:self.DET_PIXEL_NUMy,:] = edfmatsv[m,:,:]
			edfmats[m,self.DET_PIXEL_NUMy:,:]  = edfmatsh[m,:,:]

		# add the scannumber to self.scannumbers, if not already present
		if not scannumber in self.scannumbers:
			self.scannumbers.extend([scannumber])
		return data, motors, counters, edfmats

	def loadscan(self,scannumbers,scantype='generic'):
		"""
		loads the files belonging to scan No "scannumber" and puts it into an instance
		of the container-class scan. the default scantype is 'generic', later the scans
		will be grouped (then added) based on the scantype
		"""
		# make sure scannumbers are iterable (list)
		if not isinstance(scannumbers,list):
			scannums = []
			scannums.append(scannumbers)
		else:
			scannums = scannumbers 
		for number in scannums:
			scanname = 'Scan%03d' % number
			data, motors, counters, edfmats = self.readscan(number)
			# can assign some things here already (even if maybe redundant)
			monitor   = counters[self.monicolumn]
			monoangle = 1 # have to check for this later ... !!! counters['pmonoa']
			energy    = counters[self.encolumn]
			# create an instance of "scan" class for every scan
			onescan = scan(edfmats,number,energy,monoangle,monitor,counters,motors,data,scantype)
			# assign one dictionary entry to each scan 
			self.scans[scanname] = onescan

	def loadloop(self,begnums,numofregions):
		"""
		loads a whole loop of scans based on their starting scannumbers and the number of single
		scans in a single loop
		begnums      = scannumbers of the first scans of each loop (is a list)
		numofregions = number of scans in each loop (integer)
		"""
		typenames = []
		for n in range(numofregions):
			typenames.append('edge'+str(n+1))
		numbers = []
		for n in range(len(begnums)):
			for m in range(numofregions):
				numbers.append((begnums[n]+m))
		typenames = typenames*len(begnums)
		for n in range(len(typenames)):
			thenumber = []
			thenumber.append(numbers[n]) # make the scannumber into an interable list, sorry it's messy
			self.loadscan(thenumber,typenames[n])

	def loadelastic(self,scann):
		"""
		loads a scan using the loadscan function and sets the scantype attribute to 'elastic'
		"""
		self.loadscan(scann,'elastic')
	
	def loadlong(self,scann):
		"""
		loads a scan using the loadscan function and sets the scantype attribute to 'long'
		"""
		self.loadscan(scann,'long')

	def createrois(self,scannumbers):
		"""
		create rois object from this class and scannumbers
		"""
		rois_object = rois(self.scans,scannumbers)
		return rois_object

	def getautorois(self,scannumbers,kernel_size=5,colormap='jet',threshfrac=100):
		"""
		define ROIs automatically using median filtering and a variable threshold
		"""
		rois_object = self.createrois(scannumbers)
		rois_object.getautorois(kernel_size,colormap,thresholdfrac=threshfrac)
		self.rois = rois_object.rois

	def getautorois_eachdet(self,scannumbers,kernel_size=5,colormap='jet',threshfrac=100):
		"""
		define ROIs automatically one detector at a time using median filtering and a variable threshold
		"""
		rois_object = self.createrois(scannumbers)
		rois_object.getautorois_eachdet(kernel_size,colormap,thresholdfrac=threshfrac)
		self.rois = rois_object.rois

	def getlinrois(self,scannumbers,numrois=72,logscaling=True,height=5,colormap='jet'):
		"""
		define ROIs by clicking two points
		"""
		rois_object = self.createrois(scannumbers)
		rois_object.getlinrois(numrois,logscaling,height,colormap)
		self.rois = rois_object.rois

	def getzoomrois(self,scannumbers,numrois=72,logscaling=True,colormap='jet'):
		"""
		define ROIs by zooming in on a region of interest
		"""
		rois_object = self.createrois(scannumbers)
		rois_object.getzoomrois(numrois,logscaling,colormap)
		self.rois = rois_object.rois

	def getlinedgerois(self,scannumbers,energyval,numrois=72,logscaling=True,height=5,colormap='jet'):
		"""
		define ROIs by clicking two points on a matrix, which is a difference of edf-files above and 
		below an energy value
		"""
		if not isinstance(scannumbers,list):
			scannums = []
			scannums.append(scannumbers)
		else:
			scannums = scannumbers
		# search for correct index in the first of the given scans
		scanname = 'Scan%03d' % scannums[0]
		index = self.scans[scanname].energy.flat[np.abs(self.scans[scanname].energy - energyval).argmin()]
		rois_object = self.createrois(scannumbers)
		rois_object.getlinedgerois(index,numrois,logscaling,height,colormap)
		self.rois = rois_object.rois

	def getzoomedgerois(self,scannumbers,energyval,numrois=72,logscaling=True,colormap='jet'):
		"""
		define ROIs by zooming in on a matrix, which is a difference of edf-files above and 
		below an energy value
		"""
		if not isinstance(scannumbers,list):
			scannums = []
			scannums.append(scannumbers)
		else:
			scannums = scannumbers
		# search for correct index in the first of the given scans
		scanname = 'Scan%03d' % scannums[0]
		index = self.scans[scanname].energy.flat[np.abs(self.scans[scanname].energy - energyval).argmin()]
		rois_object = self.createrois(scannumbers)
		rois_object.getlinedgerois(index,numrois,logscaling,colormap)
		self.rois = rois_object.rois

	def saverois(self,filename):
		"""
		saves the rois into an file using pickle
		filename = absolute path to file and filename
		"""
		import pickle
		f = open(filename, 'wb')
		theobject = self.rois
		pickle.dump(theobject, f,protocol=-1)
		f.close()

	def loadrois(self,filename):
		"""
		loads a file written with the saverois-function using pickle
		filename = absolute path to file and filename
		"""
		import pickle
		f = open(filename,'rb')
		self.rois = pickle.load(f)
		f.close()

	def plotrois(self):
		"""
		returns a plot of the ROI shapes
		"""
		figure()
		thematrix = np.zeros((self.DET_PIXEL_NUMy*2.0,self.DET_PIXEL_NUMx))
		for roi in self.rois:
			for index in roi:
				thematrix[index[1],index[0]] = 1
		theimage = imshow(thematrix)
		show()

	def getrawdata(self):
		"""
		goes through all instances of the scan class and calls it's applyrois method
		to sum up over all rois
		"""
		if not self.rois:
			print 'please define some ROIs first.'
			return
		for scan in self.scans:
			if len(self.scans[scan].edfmats):
				print ("integrating "+scan)
				self.scans[scan].applyrois(self.rois)

	def loadscandirect(self,scannumbers,scantype='generic'):
		"""
		loads a scan without saving the edf files in matrices 
		needs ROIs as input (this needs testing)
		"""
		# make sure scannumbers are iterable (list)
		if not isinstance(scannumbers,list):
			scannums = []
			scannums.append(scannumbers)
		else:
			scannums = scannumbers 
		# check if there are ROIs defined		
		if not self.rois:
			print 'please define some ROIs first'
			return
		for number in scannums:
			scanname = 'Scan%03d' % number
			data, motors, counters, edfmats = self.readscan(number)
			# can assign some things here already (even if maybe redundant)
			monitor   = counters[self.monicolumn]
			monoangle = 1 # counters['pmonoa'] # this still needs checking
			energy    = counters[self.encolumn]
			# create an instance of "scan" class for every scan
			onescan = scan(edfmats,number,energy,monoangle,monitor,counters,motors,data,scantype)
			onescan.applyrois(self.rois)
			print 'deleting edf-files of scan No. %03d' % number
			onescan.edfmats = [] # delete the edfmats
			self.scans[scanname] = onescan

	def loadloopdirect(self,begnums,numofregions):
		"""
		loads a scan without saving the edf files in matrices and sets
		the scantype attribute to 'edge1', 'edge2', ...
		"""
		typenames = []
		for n in range(numofregions):
			typenames.append('edge'+str(n+1))
		numbers = []
		for n in range(len(begnums)):
			for m in range(numofregions):
				numbers.append((begnums[n]+m))
		typenames = typenames*len(begnums)
		for n in range(len(typenames)):
			thenumber = []
			thenumber.append(numbers[n]) # make the scannumber into an interable list, sorry it's messy
			self.loadscandirect(thenumber,typenames[n])

	def loadelasticdirect(self,scann):
		"""
		loads a scan without saving the edf files in matrices and sets
		the scantype attribute to 'elastic'
		"""
		self.loadscandirect(scann,'elastic')

	def loadlongdirect(self,scann):
		"""
		loads a scan without saving the edf files in matrices and sets
		the scantype attribute to 'long'
		"""
		self.loadscandirect(scann,'long')

	def deletescan(self,scannumbers):
		"""
		deletes the instance 'scann' of the scans class from the dictionary of scans and the number scann from the scannumbers list
		"""
		numbers = []
		if not isinstance(scannumbers,list):
			numbers.append(scannumbers)
		else:
			numbers = scannumbers
		for number in numbers:
			scanname = 'Scan%03d' % number
			del(self.scans[scanname])
			self.scannumbers.remove(number)

	def getspectrum(self):
		"""
		groups the instances of the scan class by their scantype attribute, 
		adds equal scans (each group of equal scans) and appends them 
		"""
		# find the groups
		allgroups = findgroups(self.scans)
		for group in allgroups:
			# self.makegroup(group)
			onegroup = makegroup_nointerp(group)
			self.groups[onegroup.gettype()] = onegroup
		self.energy,self.signals,self.errors = appendscans(self.groups)

	def getspectrum_test(self):
		"""
		groups the instances of the scan class by their scantype attribute, 
		adds equal scans (each group of equal scans) and appends them
		trying to correct each scan on energy scale to elastic closest to it in scnanumber (but smaller) 
		"""
		# find the groups
		allgroups   = findgroups(self.scans)
		elasticnums = [] # list with scannumbers of all elastic lines taken
		cenoms      = [] # list of lists with all center of masses for all elastic lines and ROIs
		# catch the information of all elastic lines (scannumber and center of mass)
		for group in allgroups:
			if group[0].gettype() == 'elastic':
				for n in range(len(group)):
					elasticnums.append(group[n].number)
					cenoms.append(group[n].cenom)

		for group in allgroups:
			# create groups with the info from the elastic line
			onegroup = makegroup_test(group,elasticnum=elasticnums,cenom=cenoms)
			self.groups[onegroup.gettype()] = onegroup
		self.energy,self.eloss,self.signals,self.errors,self.resolution = appendscans_test(self.groups)

	def geteloss(self):
		"""
		finds the center of mass of each roi and interpolates the signals onto a common grid
		"""	
		if not 'elastic' in self.groups:
			print 'please load/integrate at least one elastic scan first!'
			return
		else:
			for n in range(len(self.rois)):
				self.cenom.append(trapz(self.groups['elastic'].signals[:,n]*self.groups['elastic'].energy,x=self.groups['elastic'].energy)/trapz(self.groups['elastic'].signals[:,n],x=self.groups['elastic'].energy))
				try:
					FWHM,x0 = fwhm((self.groups['elastic'].energy - self.cenom[n])*1e3,self.groups['elastic'].signals[:,n])
					self.resolution.append(FWHM)
				except:
					self.resolution.append(0) # need a more sofisticated way of finding the FWHM
			self.E0 = np.mean(self.cenom)
			self.eloss = (self.energy - self.cenom[0])*1e3 # energy loss in eV
			# append the FWHM to self.resolution for each analyzer
			
			for n in range(len(self.rois)):
				# inserting zeros at beginning and end of the vectors to avoid interpolation errors
				x = (self.energy-self.cenom[n])*1e3
				x = np.insert(x,0,-1e10,)
				x = np.insert(x,-1,1e10)
				y = self.signals[:,n]
				y = np.insert(y,0,0)
				y = np.insert(y,-1,0)
				f = interp1d(x,y, bounds_error=False,fill_value=0.0)
				self.signals[:,n] = f(self.eloss)
			# set the eloss scale for all scans too but do not interpolate (using self.E0 for all)
			for number in self.scannumbers:
				scanname = 'Scan%03d' % number
				self.scans[scanname].eloss = (self.scans[scanname].energy - self.cenom[0])*1e3
				for n in range(len(self.rois)):
					# inserting zeros at beginning and end of the vectors to avoid interpolation errors
					x = (self.scans[scanname].energy-self.cenom[n])*1e3
					x = np.insert(x,0,-1e10,)
					x = np.insert(x,-1,1e10)
					y = self.scans[scanname].signals[:,n]
					y = np.insert(y,0,0)
					y = np.insert(y,-1,0)
					f = interp1d(x,y, bounds_error=False,fill_value=0.0)
					self.scans[scanname].signals[:,n] = f(self.scans[scanname].eloss)

	def gettths(self,rvd=0.0,rvu=0.0,rvb=0.0,rhl=0.0,rhr=0.0,rhb=0.0,order=[0,1,2,3,4,5]):
		"""
		uses the defined TT_OFFSETS of the read_id20 class to set all scattering angles tth
		from the mean angle avtth of the analyzer modules
		rhl = mean tth angle of HL module (default is 0.0)
		rhr = mean tth angle of HR module (default is 0.0)
		rhb = mean tth angle of HB module (default is 0.0)
		rvd = mean tth angle of VD module (default is 0.0)
		rvu = mean tth angle of VU module (default is 0.0)
		rvb = mean tth angle of VB module (default is 0.0)
		"""
		# reset all values, just in case mean angles are redefined (values are otherwise appended to existing values)
		self.VDtth = []
		self.VUtth = []
		self.VBtth = []
		self.HRtth = []
		self.HLtth = []
		self.HBtth = []
		# horizontal modules
		# HL (motor name rhl)
		v_angles = self.TTH_OFFSETS1 
		h_angles = self.TTH_OFFSETS2 + rhl
		for n in range(len(h_angles)):
			self.HLtth.append( np.arccos(np.cos(np.radians(h_angles[n]))*np.cos(np.radians(v_angles[n])))*180.0/np.pi)
		# HR (motor name rhr)
		v_angles = self.TTH_OFFSETS1 
		h_angles = self.TTH_OFFSETS2 + rhr
		for n in range(len(h_angles)):
			self.HRtth.append( np.arccos(np.cos(np.radians(h_angles[n]))*np.cos(np.radians(v_angles[n])))*180.0/np.pi)
		# HB (motor name rhb)
		v_angles = self.TTH_OFFSETS1 
		h_angles = self.TTH_OFFSETS2 + rhb
		for n in range(len(h_angles)):
			self.HBtth.append( np.arccos(np.cos(np.radians(h_angles[n]))*np.cos(np.radians(v_angles[n])))*180.0/np.pi)

		# vertical modules
		# VD
		v_angles = self.TTH_OFFSETS2 + rvd
		h_angles = self.TTH_OFFSETS1 
		for n in range(len(h_angles)):
			self.VDtth.append( np.arccos(np.cos(np.radians(v_angles[n]))*np.cos(np.radians(h_angles[n])))*180.0/np.pi)
		# VU
		v_angles = self.TTH_OFFSETS2 + rvu
		h_angles = self.TTH_OFFSETS1 
		for n in range(len(h_angles)):
			self.VUtth.append( np.arccos(np.cos(np.radians(v_angles[n]))*np.cos(np.radians(h_angles[n])))*180.0/np.pi)
		# VB
		v_angles = self.TTH_OFFSETS2 + rvb
		h_angles = self.TTH_OFFSETS1 
		for n in range(len(h_angles)):
			self.VBtth.append( np.arccos(np.cos(np.radians(v_angles[n]))*np.cos(np.radians(h_angles[n])))*180.0/np.pi)

		# list of TTH values
		tth = [self.VDtth, self.VUtth, self.VBtth, self.HRtth, self.HLtth, self.HBtth]
		# list all TTH values in one long list ordered by the 'order'-keyword
		for n in order:
			self.tth.extend(tth[n])

	def getqvals(self,invangstr=False):
		"""
		calculates q-values from E0 and tth values in either atomic units (defalt) or
		inverse angstroms 
		"""
		theqvals = np.zeros_like(self.signals)
		if invangstr:
			for n in range(len(self.tth)):
				theqvals[:,n] = momtrans_inva(self.E0+self.eloss/1e3,self.E0,self.tth[n]) 
		else:
			for n in range(len(self.tth)):
				theqvals[:,n] = momtrans_au(self.E0+self.eloss/1e3,self.E0,self.tth[n])
		self.qvalues = theqvals

	def getqvals_energy(self,energy):
		"""
		returns all q-values at a certain energy loss
		"""
		ind = np.abs(self.eloss - energy).argmin()
		return self.qvalues[ind,:]

	def copy_edf_files(self,scannumbers,destdir):
		"""
		copies all edf-files from scan with scannumber or scannumbers into directory 'destdir'
		"""
		import shutil
		numbers = []
		if not isinstance(scannumbers,list):
			numbers.append(scannumbers)
		else:
			numbers = scannumbers
		fn = self.path + self.filename
		for n in range(len(numbers)):
			data, motors, counters = specread(fn,numbers[n])
			for m in range(len(counters['ccdno'])):
				ccdnumber = counters['ccdno'][m]
				edfnameh   = self.path + self.EDF_PREFIXh + self.filename + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
				edfnamev   = self.path + self.EDF_PREFIXv + self.filename + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
				shutil.copy2(edfnameh, destdir)
				shutil.copy2(edfnamev, destdir)

	def make_posscan_image(self,scannumber,motorname):
		"""
		loads a scan from a sample position scan (x-scan, y-scan, z-scan), lets you choose a zoomroi and constructs a 2D image from this
		scannumber = number of the scan
		motorname  = string that contains the motor name (must be the same as in the SPEC file)
		"""
		plt.clf()
		# load the scan
		data, motors, counters, edfmats = self.readscan(scannumber)
		# the scan motor
		position = counters[motorname.lower()]
		# create a roi object
		roi_object = self.createrois(scannumber)
		# get one zoomroi 
		roi_object.getzoomrois_frommatrix(sumx(edfmats),1)
		roixinds = []
		roiyinds = []
		for n in range(len(roi_object.rois[0])):
			roixinds.append(roi_object.rois[0][n][0])
			roiyinds.append(roi_object.rois[0][n][1])
		# go through all edf files of the scan, sum over the height of the roi and stack the resulting lines into a matrix
		axesrange = [0,roiyinds[-1],position[0],position[-1]]
		image = (np.sum(edfmats[:,np.amin(roiyinds):np.amax(roiyinds)+1,np.amin(roixinds):np.amax(roixinds)+1],axis=1))
		plt.close()
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.imshow(np.log(image),extent=axesrange)
		ax.set_aspect('auto')
		plt.xlabel('pixel along the beam')
		ylabelstr = motorname.lower() + ' position [mm]'
		plt.ylabel(ylabelstr)
		plt.show()

	def animation(self,scannumber,logscaling=True,timeout=-1,colormap='jet'):
		"""
		shows the edf-files of a scan
		"""
		if isinstance(scannumber,list):
			if len(scannumber)>1:
				print 'this only works for a single scan, sorry'
				return
			else:
				scannumber = scannumber[0]
		scanname = 'Scan%03d' % scannumber
		edfmats = self.scans[scanname].edfmats
		scanlen  = np.shape(edfmats)[0]
		plt.ion()
		plt.clf()
		for n in range(scanlen):
			if logscaling:
				theimage = plt.imshow(np.log(edfmats[n,:,:]))
			else:
				theimage = plt.imshow(edfmats[n,:,:])
			plt.xlabel('detector x-axis [pixel]')
			plt.ylabel('detector y-axis [pixel]')
			if timeout<0:
				titlestring = 'Frame No. %d' % (n+1) + ' of %d' % scanlen + ', press key or mouse botton to continue' 
				plt.title()
			else:
				titlestring = 'Frame No. %d' % (n+1) + ' of %d' % scanlen + ', updating every %2.2f ' % timeout + ' seconds'
				plt.title(titlestring)
			theimage.set_cmap(colormap)
			plt.draw()
			plt.waitforbuttonpress(timeout=timeout)

	def showscans(self,scannumbers,columns=None):
		"""
		plots spectra from a scan or series of scans
		"""
		numbers = []
		if not isinstance(scannumbers,list):
			numbers.append(scannumbers)
		else:
			numbers = scannumbers
		plt.ion()
		plt.clf()
		for number in numbers:
			scanname = 'Scan%03d' % number
			# show all individual ROIs if no columns are specified
			if not columns:
				y = np.array(self.scans[scanname].signals)
				for n in range(len(self.rois)):
					y[:,n] = y[:,n]/self.scans[scanname].monitor
				try:
					titlestring = 'all analyzer signals of scan No. ' + scanname
					plt.title(titlestring)
					plt.plot(self.scans[scanname].eloss,y)
					plt.xlabel('energy loss [eV]')
					plt.ylabel('signal [arb. units]')
				except:
					titlestring = 'all analyzer signals of scan No. ' + scanname
					plt.title(titlestring)
					plt.plot(self.scans[scanname].energy,y)
					plt.xlabel('primary energy [keV]')
					plt.ylabel('signal [arb. units]')
			# sum over the columns if they are provided
			if columns:
				if isinstance(columns,list):
					y = np.zeros_like(self.scans[scanname].eloss)
					for n in columns:
						y += self.scans[scanname].signals[:,n]
					y = y/self.scans[scanname].monitor
				else:
					y = self.scans[scanname].signals[:,columns]
				try:
					titlestring = 'all analyzer signals of scan No. ' + scanname
					plt.title(titlestring)
					plt.plot(self.scans[scanname].eloss,y)
					plt.xlabel('energy loss [eV]')
					plt.ylabel('signal [arb. units]')
				except:
					titlestring = 'all analyzer signals of scan No. ' + scanname
					plt.title(titlestring)
					plt.plot(self.scans[scanname].energy,y)
					plt.xlabel('primary energy [keV]')
					plt.ylabel('signal [arb. units]')

	def showloops(self,scannumbers,looplength):
		"""
		plots the results one loop at a time
		"""
		pass

	def plotspectrum(self,columns=None):
		"""
		make a plot of the raw spectrum of the entire instance, sum over selected analyzers if porvided
		"""
		if len(self.eloss) == 0:
			print 'please provide elastic line scans first and run the geteloss() method'
			return
		if not columns:
			plt.plot(self.eloss,self.signals)
			plt.title('all analyzers vs. energy loss')
			plt.xlabel('energy loss [eV]')
			plt.ylabel('normalized signal [arb. units]')
		else:
			plt.plot(self.eloss,np.sum(self.signals[:,columns],axis=1))
			plt.title('all analyzers vs. energy loss')
			plt.xlabel('energy loss [eV]')
			plt.ylabel('normalized signal [arb. units]')

	def printlength(self,scannumbers):
		numbers = []
		if not isinstance(scannumbers,list):
			numbers.append(scannumbers)
		else:
			numbers = scannumbers

		for i in numbers:
			name = 'Scan%03d' % i
			print 'length of scan %03d ' %i + ' is ' + str(len(self.scans[name].energy))


class read_id16:
	"""
	a class for loading data form SPEC files and the according edf files from the old id16
	"""
	def __init__(self,absfilename,energycolumn='energy_cc',monitorcolumn='monitor'):
		self.scans         = {} # was a dictionary before
		if not os.path.isfile(absfilename):
			raise Exception('IOError! No such file, please check filename.')
		self.path          = os.path.split(absfilename)[0] + '/'
		self.filename      = os.path.split(absfilename)[1]
		self.scannumbers   = []
		self.EDF_PREFIX    = 'edf/'
		self.EDF_POSTFIX   = '.edf'
		self.DET_PIXEL_NUM = 256
		self.TTH_OFFSETS   = np.array([-13.0,-6.5,-6.5,0.0,0.0,6.5,6.5,13.0,13.0])
		# which column in the SPEC file to be used for the energy and monitor
		self.encolumn      = energycolumn.lower()
		self.monicolumn    = monitorcolumn.lower()
		# here are the attributes of the old rawdata class
		self.eloss    = []	# common eloss scale for all analyzers
		self.energy   = []	# common energy scale for all analyzers
		self.signals  = []	# signals for all analyzers
		self.errors   = []	# poisson errors
		self.qvalues  = []	# for all analyzers
		self.groups   = {}  # groups (such as 2 'elastic', or 5 'edge1', etc.)
		self.cenom    = []
		self.E0       = []
		self.tth      = []
		# input
		self.rois     = []	# can be set once the rois are defined to integrate scans directly
		# intermediates
		self.rawsignals = []	# structure of rawdata (energies, raw-signals, monitor, etc.)
		self.rawenergy  = []
		self.rawmonitor = []

	def readscan(self,scannumber):
		"""
		returns the data, motors, counter-names, and edf-files from the calss's specfile of "scannumber"
		should use PyMca's routines to read the Spec- and edf-files
		"""
		if self.path: # if a path was provided upon class construction, use this one
			path = self.path
		else: # if no path was provided
			print 'please provide a path to the SPEC-file'
			return
		if not self.filename:
			print 'please provide a SPEC-file name!'
			return
		else:
			fn = path + self.filename
		#try: #try loading variables from file
		print ("parsing edf- and SPEC-files of scan No. %s" % scannumber)
		data, motors, counters = specread(fn,scannumber)
		edfmats = np.array(np.zeros((len(counters['ccdno']),self.DET_PIXEL_NUM,self.DET_PIXEL_NUM)))
		for m in range(len(counters['ccdno'])):
			ccdnumber = counters['ccdno'][m]+1
			edfname   = path + self.EDF_PREFIX + self.filename + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
			edfmats[m,:,:] = edfread(edfname)
		self.scannumbers.extend([scannumber]) # keep track of the scnanumbers
		return data, motors, counters, edfmats

	def loadscan(self,scannumbers,scantype='generic'):
		"""
		loads the files belonging to scan No "scannumber" and puts it into an instance
		of the container-class scan. the default scantype is 'generic', later the scans
		will be grouped (then added) based on the scantype
		"""
		# make sure scannumbers are iterable (list)
		if not isinstance(scannumbers,list):
			scannums = []
			scannums.append(scannumbers)
		else:
			scannums = scannumbers 
		for number in scannums:
			scanname = 'Scan%03d' % number
			data, motors, counters, edfmats = self.readscan(number)
			# can assign some things here already (even if maybe redundant)
			monitor   = counters[self.monicolumn]
			monoangle = counters['pmonoa']
			energy    = counters[self.encolumn]
			# create an instance of "scan" class for every scan
			onescan = scan(edfmats,number,energy,monoangle,monitor,counters,motors,data,scantype)
			# assign one dictionary entry to each scan 
			self.scans[scanname] = onescan

	def loadloop(self,begnums,numofregions):
		"""
		loads a whole loop of scans based on their starting scannumbers and the number of single
		scans in a single loop
		begnums      = scannumbers of the first scans of each loop (is a list)
		numofregions = number of scans in each loop (integer)
		"""
		typenames = []
		for n in range(numofregions):
			typenames.append('edge'+str(n+1))
		numbers = []
		for n in range(len(begnums)):
			for m in range(numofregions):
				numbers.append((begnums[n]+m))
		typenames = typenames*len(begnums)
		for n in range(len(typenames)):
			thenumber = []
			thenumber.append(numbers[n]) # make the scannumber into an interable list, sorry it's messy
			# self = self.loadscan(thenumber,typenames[n])
			self.loadscan(thenumber,typenames[n])

	def loadelastic(self,scann):
		"""
		loads a scan using the loadscan function and sets the scantype attribute to 'elastic'
		"""
		self.loadscan(scann,'elastic')
	
	def loadlong(self,scann):
		"""
		loads a scan using the loadscan function and sets the scantype attribute to 'long'
		"""
		self.loadscan(scann,'long')

	def createrois(self,scannumbers):
		"""
		create rois object from this class and scannumbers
		"""
		rois_object = rois(self.scans,scannumbers)
		return rois_object

	def getautorois(self,scannumbers,kernel_size=5,colormap='jet'):
		"""
		define ROIs automatically using median filtering and a variable threshold
		"""
		rois_object = self.createrois(scannumbers)
		rois_object.getautorois(kernel_size,colormap)
		self.rois = rois_object.rois

	def getlinrois(self,scannumbers,numrois=9,logscaling=True,height=5,colormap='jet'):
		"""
		define ROIs by clicking two points
		"""
		rois_object = self.createrois(scannumbers)
		rois_object.getlinrois(numrois,logscaling,height,colormap)
		self.rois = rois_object.rois

	def getzoomrois(self,scannumbers,numrois=9,logscaling=True,colormap='jet'):
		"""
		define ROIs by zooming in on a region of interest
		"""
		rois_object = self.createrois(scannumbers)
		rois_object.getzoomrois(numrois,logscaling,colormap)
		self.rois = rois_object.rois

	def getlinedgerois(self,scannumbers,energyval,numrois=9,logscaling=True,height=5,colormap='jet'):
		"""
		define ROIs by clicking two points on a matrix, which is a difference of edf-files above and 
		below an energy value
		"""
		if not isinstance(scannumbers,list):
			scannums = []
			scannums.append(scannumbers)
		else:
			scannums = scannumbers
		# search for correct index in the first of the given scans
		scanname = 'Scan%03d' % scannums[0]
		index = self.scans[scanname].energy.flat[np.abs(self.scans[scanname].energy - energyval).argmin()]
		rois_object = self.createrois(scannumbers)
		rois_object.getlinedgerois(index,numrois,logscaling,height,colormap)
		self.rois = rois_object.rois

	def getzoomedgerois(self,scannumbers,energyval,numrois=9,logscaling=True,colormap='jet'):
		"""
		define ROIs by zooming in on a matrix, which is a difference of edf-files above and 
		below an energy value
		"""
		if not isinstance(scannumbers,list):
			scannums = []
			scannums.append(scannumbers)
		else:
			scannums = scannumbers
		# search for correct index in the first of the given scans
		scanname = 'Scan%03d' % scannums[0]
		index = self.scans[scanname].energy.flat[np.abs(self.scans[scanname].energy - energyval).argmin()]
		rois_object = self.createrois(scannumbers)
		rois_object.getlinedgerois(index,numrois,logscaling,colormap)
		self.rois = rois_object.rois

	def getrawdata(self):
		"""
		goes through all instances of the scan class and calls it's applyrois method
		to sum up over all rois
		"""
		if not self.rois:
			print 'please define some ROIs first.'
			pass
		for scan in self.scans:
			if len(self.scans[scan].edfmats):
				print ("integrating "+scan)
				self.scans[scan].applyrois(self.rois)

	def loadscandirect(self,scannumbers,rois,scantype='generic'):
		"""
		loads a scan without saving the edf files in matrices 
		needs ROIs as input (this needs testing)
		"""
		for number in scannumbers:
			scanname = 'Scan%03d' % number
			data, motors, counters, edfmats = self.readscan(number)
			# can assign some things here already (even if maybe redundant)
			monitor   = counters[self.monicolumn]
			monoangle = counters['pmonoa']
			energy    = counters[self.encolumn]
			# create an instance of "scan" class for every scan
			onescan = scan(edfmats,number,energy,monoangle,monitor,counters,motors,data,scantype)
			onescan.applyrois(rois.rois)
			onescan.edfmats = [] # delete the edfmats
			self.scans[scanname] = onescan
		# if there are no rois yet, set the rois	
		if not self.rois:
			self.rois = rois

	def loadloopdirect(self,begnums,numofregions,rois):
		"""
		loads a scan without saving the edf files in matrices and sets
		the scantype attribute to 'edge1', 'edge2', ...
		"""
		typenames = []
		for n in range(numofregions):
			typenames.append('edge'+str(n+1))
		numbers = []
		for n in range(len(begnums)):
			for m in range(numofregions):
				numbers.append((begnums[n]+m))
		typenames = typenames*len(begnums)
		for n in range(len(typenames)):
			thenumber = []
			thenumber.append(numbers[n]) # make the scannumber into an interable list, sorry it's messy
			self.loadscandirect(thenumber,typenames[n])

	def loadelasticdirect(self,scann,rois):
		"""
		loads a scan without saving the edf files in matrices and sets
		the scantype attribute to 'elastic'
		"""
		self.loadscandirect(scann,rois,'elastic')
	
	def loadlongdirect(self,scann,rois):
		"""
		loads a scan without saving the edf files in matrices and sets
		the scantype attribute to 'long'
		"""
		self.loadscandirect(scann,rois,'long')

	def getspectrum(self):
		"""
		groups the instances of the scan class by their scantype attribute, 
		adds equal scans (each group of equal scans) and appends them 
		"""
		# find the groups 
		allgroups = findgroups(self.scans)
		for group in allgroups:
			# self.makegroup(group)
			onegroup = makegroup(group)
			self.groups[onegroup.gettype()] = onegroup
		self.energy,self.signals,self.errors = appendscans(self.groups)

	def geteloss(self):
		"""
		finds the center of mass of each roi and interpolates the signals onto a common grid
		"""	
		if not 'elastic' in self.groups:
			print 'please load/integrate at least one elastic scan first!'
		else:
			for n in range(len(self.rois)):
				self.cenom.append(trapz(self.groups['elastic'].signals[:,n]*self.groups['elastic'].energy,x=self.groups['elastic'].energy)/trapz(self.groups['elastic'].signals[:,n],x=self.groups['elastic'].energy))
			self.E0 = np.mean(self.cenom)
			self.eloss = (self.energy - self.cenom[0])*1e3 # energy loss in eV
			# interpolate all signals onto the energy scale of the first analyzer
			#if self.signals.any():
			for n in range(len(self.rois)):
				f = interp1d((self.energy-self.cenom[n])*1e3,self.signals[:,n], bounds_error=False,fill_value=0.0)
				self.signals[:,n] = f(self.eloss)

	def gettths(self,avtth):
		"""
		uses the defined TT_OFFSETS of the read_id20 class to set all scattering angles tth
		from the mean angle avtth of the analyzer module
		"""
		self.tth = np.absolute(avtth + self.TTH_OFFSETS)

	def getqvals(self,invangstr=False):
		"""
		calculates q-values from E0 and tth values in either atomic units (defalt) or
		inverse angstroms 
		"""
		theqvals = np.zeros_like(self.signals)
		if invangstr:
			for n in range(len(self.tth)):
				theqvals[:,n] = momtrans_inva(self.E0+self.eloss/1e3,self.E0,self.tth[n]) 
		else:
			for n in range(len(self.tth)):
				theqvals[:,n] = momtrans_au(self.E0+self.eloss/1e3,self.E0,self.tth[n])
		self.qvalues = theqvals

	def getqvals_energy(self,energy):
		"""
		returns all q-values at a certain energy loss
		"""
		ind = np.abs(self.eloss - energy).argmin()
		return self.qvalues[ind,:]

	def copy_edf_files(self,scannumbers,destdir):
		"""
		copies all edf-files from scan with scannumber or scannumbers into directory 'destdir'
		"""
		import shutil
		numbers = []
		if not isinstance(scannumbers,list):
			numbers.append(scannumbers)
		else:
			numbers = scannumbers
		fn = self.path + self.filename
		for n in range(len(numbers)):
			data, motors, counters = specread(fn,numbers[n])
			for m in range(len(counters['ccdno'])):
				ccdnumber = counters['ccdno'][m]+1
				edfname   = self.path + self.EDF_PREFIX + self.filename + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
				shutil.copy2(edfname, destdir)

class rois:
	"""
	a class to define ROIs
	ROIs are saved in self.rois as a list with the same number of entries as ROIs defined. Each entry is a list of (x,y)-coordinate tuples
	"""
	def __init__(self,scans,scannumbers):
		self.scandata = scans # dictionary with instances of the onescan class
		if not isinstance(scannumbers,list):
			thescannumbers = []
			thescannumbers.append(scannumbers)
		else:
			thescannumbers = scannumbers  
		self.scannums = thescannumbers
		self.rois     = []
		self.roinums  = []
		self.roikind  = []
		self.roiheight = []
		self.roiwidth  = []
		self.DET_PIXEL_NUM = 256
		# adding list of lists with x- and y-indices, useful for 2D things (compare Simos Matlab scripts)
		self.roisx = [] 
		self.roisy = []

	def preparemats(self):
		"""
		sums and squeezes all edf-files of a scan into one matrix 
		"""
		# take shape of the first edf-matrices
		dim = np.shape(self.scandata['Scan%03d' % self.scannums[0]].edfmats)
		#edfmatrices = np.zeros((1,self.DET_PIXEL_NUM,self.DET_PIXEL_NUM))
		edfmatrices = np.zeros((1,dim[1],dim[2]))
		for number in self.scannums:
			scanname = 'Scan%03d' % number
			edfmatrices = np.append(edfmatrices,self.scandata[scanname].edfmats,axis=0)
		return sumx(edfmatrices)
	
	def prepareedgemats(self,index):
		"""
		difference between two summed and squeezed edf-files of a scan from below and above energyval 
		"""
		# take shape of the first edf-matrices 
		dim = np.shape(self.scandata['Scan%03d' % self.scannums[0]].edfmats)
		abovemats = np.zeros((1,dim[1],dim[2]))
		belowmats = np.zeros((1,dim[1],dim[2]))
		# abovemats = np.zeros((1,self.DET_PIXEL_NUM,self.DET_PIXEL_NUM))
		# belowmats = np.zeros((1,self.DET_PIXEL_NUM,self.DET_PIXEL_NUM))
		for number in self.scannums:
			scanname = 'Scan%03d' % number
			# abovemats = np.append(abovemats,self.scandata[scanname].edfmats[self.scandata[scanname].energy >  energyval,:,:],axis=0)
			# belowmats = np.append(belowmats,self.scandata[scanname].edfmats[self.scandata[scanname].energy <= energyval,:,:],axis=0)
			abovemats = np.append(abovemats,self.scandata[scanname].edfmats[index:,:,:],axis=0)
			belowmats = np.append(belowmats,self.scandata[scanname].edfmats[:index,:,:],axis=0)
		return np.absolute(np.squeeze(np.sum(abovemats,axis=0))-np.squeeze(np.sum(belowmats,axis=0)))

	def getlinrois(self,numrois,logscaling=True,height=5,colormap='jet'):
		"""
		define ROIs by clicking two points
		"""

		fig = plt.figure(figsize=(8, 6))
		ax = fig.add_subplot(111, axisbg='#FFFFCC')

		thematrix = self.preparemats()
		plt.title('Click start and endpoints of the linear ROIs you whish.',fontsize=14)		
		if logscaling:
			theimage = plt.imshow(np.log(np.absolute(thematrix)))
		else:
			theimage = plt.imshow(thematrix)
		# Choose a color palette
		theimage.set_cmap(colormap)
		plt.colorbar()

		therois   = [] # list of lists of tuples with (x,y) coordinates
		theroisx  = [] # list of lists of x-indices
		theroisy  = [] # list of lists of y-indices

		cursor = Cursor(ax, useblit=True, color='red', linewidth=1 )

		for m in range(numrois):
			plt.clf()
			if m>0:
				if therois[m-1]:
					for index in therois[m-1]:
						thematrix[index[1],index[0]] = np.amax(thematrix)
			print ("chose two points as endpoints for linear ROI No. %s" % (m+1))

			if logscaling:
				theimage = imshow(np.log(np.absolute(thematrix)))
			else:
				theimage = imshow(thematrix)	
			theimage.set_cmap(colormap)
			plt.draw()	

			endpoints = []
			endpts    = ginput(2)
			for point in endpts:
				point = [round(element) for element in point]
				endpoints.append(np.array(point))
			roix = np.arange(endpoints[0][0],endpoints[1][0])
			roiy = [round(num) for num in np.polyval(np.polyfit([endpoints[0][0],endpoints[1][0]],[endpoints[0][1],endpoints[1][1]],1),roix)]

			theheight = np.arange(-height,height)
			
			eachroi  = []
			eachroix = []
			eachroiy = []
			for n in range(len(roix)):
				for ind in theheight:
					eachroi.append((roix[n],roiy[n]+ind))
					eachroix.append(roix[n])
					eachroiy.append(roiy[n]+ind)
			therois.append(eachroi)
			theroisx.append(eachroix)
			theroisy.append(eachroiy)

			self.roiheight.append(len(theheight)) # save the hight of each roi (e.g. for imaging)
			self.roiwidth.append(endpoints[1][0]-endpoints[0][0]) # save the width of each roi
			del endpoints
			del eachroi
			del eachroix
			del eachroiy
	
		plt.show()		
		self.roikind  = 'linear'
		self.rois     = therois
		self.roisx    = theroisx
		self.roisy    = theroisy
		self.roinums  = numrois

	def getzoomrois(self,numrois,logscaling=True,colormap='jet'):
		"""
		define ROIs by zooming in on plot
		"""
		thematrix = self.preparemats()
		limmax    = np.shape(thematrix)[1]
		therois   = [] # list of lists of tuples with (x,y) coordinates
		theroisx  = [] # list of lists of x-indices
		theroisy  = [] # list of lists of y-indices
		title('Zoom around ROI, change to shell, and press Return')
		if logscaling:
			theimage = imshow(np.log(np.absolute(thematrix)))
		else:
			theimage = imshow(thematrix)	
		theimage.set_cmap(colormap)
		ion()	
		for m in range(numrois):
			clf()
			if m>0:
				if therois[m-1]:
					for index in therois[m-1]:
						thematrix[index[1],index[0]] = np.amax(thematrix)
			title('Zoom around ROI, change to shell, and press Return')
			if logscaling:
				theimage = imshow(np.log(np.absolute(thematrix)))
			else:
				theimage = imshow(thematrix)	
			theimage.set_cmap(colormap)
			draw()	
						
			thestring = 'Zoom around ROI No. %s and press enter to continue.' % (m+1)
			wait = raw_input(thestring)
								
			limits = np.floor(axis())
			#axis([0.0,self.DET_PIXEL_NUM,self.DET_PIXEL_NUM,0.0]) # check this for the larger detector size!
			
			inds = limits < 1
			limits[inds] = 1
			inds = limits > limmax
			limits[inds] = limmax
			print 'selected limits are: ', limits

			T = np.zeros(np.shape(thematrix))
			T[limits[3]:limits[2],:] = 1
			t = np.zeros(np.shape(thematrix))
			t[:,limits[0]:limits[1]] = 1
			indsx,indsy = np.where(t*T == 1)
			eachroi  = []
			eachroix = []
			eachroiy = []
			for n in range(len(indsx)):
				eachroi.append((indsy[n],indsx[n]))
				eachroix.append(indsx[n])
				eachroiy.append(indsy[n])
			therois.append(eachroi)
			theroisx.append(eachroix)
			theroisy.append(eachroiy)
		show()
		self.roikind  = 'zoom'
		self.rois     = therois
		self.roisx    = theroisx
		self.roisy    = theroisy
		self.roinums  = numrois

	def getzoomrois_frommatrix(self,matrix,numrois,logscaling=True,colormap='jet'):
		"""
		define ROIs by zooming in on plot
		"""
		thematrix = matrix
		limmax = np.shape(thematrix)[1]
		therois = []
		title('Zoom around ROI, change to shell, and press Return')
		if logscaling:
			theimage = imshow(np.log(np.absolute(thematrix)))
		else:
			theimage = imshow(thematrix)	
		theimage.set_cmap(colormap)
		ion()	
		for m in range(numrois):
			
			title('Zoom around ROI, change to shell, and press Return')
			if logscaling:
				theimage = imshow(np.log(np.absolute(thematrix)))
			else:
				theimage = imshow(thematrix)	
			theimage.set_cmap(colormap)
			ion()	
						
			thestring = 'Zoom around ROI No. %s and press enter to continue.' % (m+1)
			wait = raw_input(thestring)
								
			limits = np.floor(axis())
			axis([0.0,self.DET_PIXEL_NUM,self.DET_PIXEL_NUM,0.0])
			
			inds = limits < 1
			limits[inds] = 1
			inds = limits > limmax
			limits[inds] = limmax
			print 'selected limits are: ', limits

			T = np.zeros(np.shape(thematrix))
			T[limits[3]:limits[2],:] = 1
			t = np.zeros(np.shape(thematrix))
			t[:,limits[0]:limits[1]] = 1
			indsx,indsy = np.where(t*T == 1)
			eachroi = []
			for n in range(len(indsx)):
				eachroi.append((indsy[n],indsx[n]))
			therois.append(eachroi)
			
		show()
		self.roikind  = 'zoom'
		self.rois     = therois
		self.roinums  = numrois

	def getlinedgerois(self,index,numrois=9,logscaling=True,height=5,colormap='jet'):
		"""
		define ROIs by clicking two points in difference picture
		"""
		# thematrix = self.prepareedgemats(energyval) # change this to use an index
		thematrix = self.prepareedgemats(index) # change this to use an index
		if logscaling:
			theimage = imshow(np.log(np.absolute(thematrix)))
		else:
			theimage = imshow(thematrix)
		# Choose a color palette
		theimage.set_cmap(colormap)
		colorbar()

		therois = []
		for m in range(numrois):
			print ("chose two points as endpoints for linear ROI No. %s" % (m+1))
			endpoints = []
			endpts    = ginput(2)
			for point in endpts:
				point = [round(element) for element in point]
				endpoints.append(np.array(point))
			roix = np.arange(endpoints[0][0],endpoints[1][0])
			roiy = [round(num) for num in np.polyval(np.polyfit([endpoints[0][0],endpoints[1][0]],[endpoints[0][1],endpoints[1][1]],1),roix)]
			del endpoints
			theheight = np.arange(-height,height)
			eachroi = []
			for n in range(len(roix)):
				for ind in theheight:
					eachroi.append((roix[n],roiy[n]+ind))
			therois.append(eachroi)
			del eachroi			
		show()		
		self.roikind  = 'linear'
		self.rois     = therois
		self.roinums  = numrois

	def getzoomedgerois(self,index,numrois=9,logscaling=True,colormap='jet'):
		"""
		define ROIs by zooming in on plot
		"""
		# thematrix = self.prepareedgemats(energyval)
		thematrix = self.prepareedgemats(index)
		limmax = np.shape(thematrix)[1]
		therois = []
		title('Zoom around ROI, change to shell, and press Return')
		if logscaling:
			theimage = imshow(np.log(np.absolute(thematrix)))
		else:
			theimage = imshow(thematrix)
		theimage.set_cmap(colormap)
		ion()
		for m in range(numrois):
			
			title('Zoom around ROI, change to shell, and press Return')
			if logscaling:
				theimage = imshow(np.log(np.absolute(thematrix)))
			else:
				theimage = imshow(thematrix)
			theimage.set_cmap(colormap)
			ion()	
						
			thestring = 'Zoom around ROI No. %s and press enter to continue.' % (m+1)
			wait = raw_input(thestring)
								
			limits = np.floor(axis())
			axis([0.0,self.DET_PIXEL_NUM,self.DET_PIXEL_NUM,0.0])
			
			inds = limits < 1
			limits[inds] = 1
			inds = limits > limmax
			limits[inds] = limmax
			print 'selected limits are: ', limits
			
			T = np.zeros(np.shape(thematrix))
			T[limits[3]:limits[2],:] = 1
			t = np.zeros(np.shape(thematrix))
			t[:,limits[0]:limits[1]] = 1
			indsx,indsy = np.where(t*T == 1)
			eachroi = []
			for n in range(len(indsx)):
				eachroi.append((indsy[n],indsx[n]))
			therois.append(eachroi)
		show()
		self.roikind  = 'zoom'
		self.rois     = therois
		self.roinums  = numrois

	def getautorois(self,kernel_size=5,colormap='jet',thresholdfrac=100):
		"""
		define ROIs by choosing a threshold through a slider bar on the plot window
		"""
		thematrix = self.preparemats() # the starting matrix to plot
		title('Crank the threshold and close the plotting window, when you are satisfied.',fontsize=14)		
		ax = subplot(111) 
		subplots_adjust(left=0.05, bottom=0.2)
		thres0 = 0 # initial threshold value
		theimage = imshow(np.log(np.absolute(thematrix)))
		theimage.set_cmap(colormap)
		colorbar()		
		thresxcolor = 'lightgoldenrodyellow'
		thresxamp  = axes([0.2, 0.10, 0.55, 0.03], axisbg=thresxcolor)
		sthres = Slider(thresxamp, 'Threshold', 0.0, np.floor(np.amax(thematrix)/thresholdfrac), valinit=thres0)
		
		def update(val):
			thres = sthres.val
			newmatrix = signal.medfilt2d(thematrix, kernel_size=kernel_size)
			belowthres_indices = newmatrix < thres
			newmatrix[belowthres_indices] = 0			
			labelarray,numfoundrois = measurements.label(newmatrix)
			print str(numfoundrois) + ' ROIs found!'	
			# organize rois
			therois = []
			for n in range(numfoundrois):
				rawindices = np.nonzero(labelarray == n+1)
				eachroi = []				
				for m in range(len(rawindices[0])):
					eachroi.append((rawindices[1][m],rawindices[0][m]))
				therois.append(eachroi)
			self.rois = therois
			self.roinums = numfoundrois
			theimage.set_data(newmatrix)
			draw()
		sthres.on_changed(update)
		show()
		self.roikind = 'auto'

	def getautorois_eachdet(self,kernel_size=5,colormap='jet',thresholdfrac=100):
		"""
		autoroi, detector for detector for ID20 (6 detector setup) so that different 
		thresholds can be choosen for the different detectors
		define ROIs by choosing a threshold through a slider bar on the plot window
		"""
		self.rois = []
		self.roinums = 0
		wholematrix = np.array(self.preparemats()) # full 6-detector image
		imagelims = [[0,256,0,256],[0,256,256,512],[0,256,512,768],[256,512,0,256],[256,512,256,512],[256,512,512,768]] 

		for n in range(len(imagelims)):
			thematrix = wholematrix[imagelims[n][0]:imagelims[n][1],imagelims[n][2]:imagelims[n][3]]
			title('Crank the threshold and close the plotting window, when you are satisfied.',fontsize=14)		
			ax = subplot(111) 
			subplots_adjust(left=0.05, bottom=0.2)
			thres0 = 0 # initial threshold value
			theimage = imshow(np.log(np.absolute(thematrix)))
			theimage.set_cmap(colormap)
			colorbar()		
			thresxcolor = 'lightgoldenrodyellow'
			thresxamp  = axes([0.2, 0.10, 0.55, 0.03], axisbg=thresxcolor)
			sthres = Slider(thresxamp, 'Threshold', 0.0, np.floor(np.amax(thematrix)/thresholdfrac), valinit=thres0)
		
			def update(val):
				thres = sthres.val
				newmatrix = signal.medfilt2d(thematrix, kernel_size=kernel_size)
				belowthres_indices = newmatrix < thres
				newmatrix[belowthres_indices] = 0			
				labelarray,numfoundrois = measurements.label(newmatrix)
				print str(numfoundrois) + ' ROIs found!'	
				# organize rois
				therois = []
				for n in range(numfoundrois):
					rawindices = np.nonzero(labelarray == n+1)
					eachroi = []				
					for m in range(len(rawindices[0])):
						eachroi.append((rawindices[1][m],rawindices[0][m]))
					therois.append(eachroi)
				self.rois.append(therois[-1])
				theimage.set_data(newmatrix)
				draw()

			sthres.on_changed(update)
			show()
		self.roikind = 'auto'

	def saverois(self,filename):
		"""
		save the ROIs in file with name filename
		"""	
		f = open(filename, 'wb')
		theobject = [self.rois, self.roinums, self.roikind]
		pickle.dump(theobject, f,protocol=-1)
		f.close()

	def loadrois(self,filename):
		"""
		load ROIs from file with name filename
		"""
		f = open(filename,'rb')
		theobject = pickle.load(f)
		f.close()
		self.rois    = theobject[0]
		self.roinums = theobject[1]
		self.roikind = theobject[2]
		
	def deleterois(self):
		"""
		delete the existing ROIs
		"""
		self.rois    = []
		self.roinums = []
		self.roikind = []

	def plotrois(self):
		"""
		returns a plot of the ROI shapes
		"""
		figure()
		thematrix = np.zeros((self.DET_PIXEL_NUM,self.DET_PIXEL_NUM))
		for roi in self.rois:
			for index in roi:
				thematrix[index[1],index[0]] = 1
		theimage = imshow(thematrix)
		show()

# horizontal zoom selector
"""
ax = subplot(111)
ax.plot(x,y)

def onselect(vmin, vmax):
    print vmin, vmax
span = SpanSelector(ax, onselect, 'horizontal')
"""

# example for rectangle selector from widgets
"""
from matplotlib.widgets import RectangleSelector

fig, ax = plt.subplots()
x = np.random.normal(size=1000)
y = np.random.normal(size=1000)
c = np.zeros((1000, 3))
c[:, 2] = 1  # set to blue
points = ax.scatter(x, y, s=20, c=c)

def selector_callback(eclick, erelease):
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    global c
    c[(x >= min(x1, x2)) & (x <= max(x1, x2))
      & (y >= min(y1, y2)) & (y <= max(y1, y2))] = [1, 0, 0]
    points.set_facecolors(c)
    fig.canvas.draw()
    
    
selector = RectangleSelector(ax, selector_callback,
                             drawtype='box', useblit=True,
                             button=[1,3], # don't use middle button
                             minspanx=5, minspany=5,
                             spancoords='pixels')
"""
# another one
"""
from matplotlib.widgets import  RectangleSelector
from pylab import *

def onselect(eclick, erelease):
  'eclick and erelease are matplotlib events at press and release'
  print ' startposition : (%f, %f)' % (eclick.xdata, eclick.ydata)
  print ' endposition   : (%f, %f)' % (erelease.xdata, erelease.ydata)
  print ' used button   : ', eclick.button

def toggle_selector(event):
    print ' Key pressed.'
    if event.key in ['Q', 'q'] and toggle_selector.RS.active:
        print ' RectangleSelector deactivated.'
        toggle_selector.RS.set_active(False)
    if event.key in ['A', 'a'] and not toggle_selector.RS.active:
        print ' RectangleSelector activated.'
        toggle_selector.RS.set_active(True)

x = arange(100)/(99.0)
y = sin(x)
fig = figure
ax = subplot(111)
ax.plot(x,y)

toggle_selector.RS = RectangleSelector(ax, onselect, drawtype='line')
connect('key_press_event', toggle_selector)
show()

"""

class read_p01:
	"""
	a class for loading data fio, DSC, and binary files from PetraIII beamline P01 (experimental)
	"""
	def __init__(self,prefix):
		self.scans         = {}
		self.prefix        = prefix
		self.scannumbers   = []
		self.DET_PIXEL_NUM = 256
		self.TTH_OFFSETS   = np.array([])
		self.eloss    = []	# common eloss scale for all analyzers
		self.energy   = []	# common energy scale for all analyzers
		self.signals  = []	# signals for all analyzers
		self.errors   = []	# poisson errors
		self.qvalues  = []	# for all analyzers
		self.groups   = {}  # groups (such as 2 'elastic', or 5 'edge1', etc.)
		self.cenom    = []
		self.E0       = []
		self.tth      = []
		# input
		self.rois     = []		# can be set once the rois are defined to integrate scans
		# intermediates
		self.rawsignals = []	# structure of rawdata (energies, raw-signals, monitor, etc.)
		self.rawenergy  = []
		self.rawmonitor = []
	
	def loadscan(self,scannumbers,reps=None,usemonoa=False,scantype='generic'):
		monicol  = 6
		encol    = 0
		anglecol = 0
		# make sure scannumbers are iterable (list)
		if not isinstance(scannumbers,list):
			scannums = []
			scannums.append(scannumbers)
		else:
			scannums = scannumbers
		if not reps:
			for number in scannums:
				scanname = 'Scan%03d' % number
				fiodata, mats, mats1, mats2 = readp01scan(self.prefix,number)
				# can assign some things here already (even if maybe redundant)
				print np.shape(fiodata)				
				monitor   = fiodata[:,monicol]
				monoangle = fiodata[:,anglecol]				
				if usemonoa:
					# need to flipup energy scale since monoangle is scanned from small to large angles
					inds    = np.argsort(energy_monoangle(fiodata[:,anglecol]))
					fiodata = fiodata[inds,:]						
					energy  = energy_monoangle(fiodata[:,anglecol])
					mats    = mats[inds,:,:]
					mats1   = mats1[inds,:,:]
					mats2   = mats2[inds,:,:]
					monitor = monitor[inds]
					monoangle= monoangle[inds]
				else:
					energy = fiodata[:,encol]/1.0e3
				
				# create an instance of "scan" class for every scan
				onescan = scan(mats,number,energy,monoangle,monitor,[],[],fiodata,scantype)
				# assign one dictionary entry to each scan 
				self.scans[scanname] = onescan
		else:
			for rep in reps:
				scanname = 'Scan%03d' % scannums[0] + 'r%01d' % rep
				fiodata, mats, mats1, mats2 = readp01scan_rep(self.prefix,scannums[0],rep)
				# can assign some things here already (even if maybe redundant)
				monitor   = fiodata[:,monicol]
				monoangle = fiodata[:,anglecol]				
				if usemonoa:
					# need to flipup energy scale since monoangle is scanned from small to large angles
					inds    = np.argsort(energy_monoangle(fiodata[:,anglecol]))
					fiodata = fiodata[inds,:]						
					energy  = energy_monoangle(fiodata[:,anglecol])
					mats    = mats[inds,:,:]
					mats1   = mats1[inds,:,:]
					mats2   = mats2[inds,:,:]
					monitor = monitor[inds]
					monoangle= monoangle[inds]
				else:
					energy = fiodata[:,encol]/1.0e3

				# create an instance of "scan" class for every scan
				onescan = scan(mats,scannums[0],energy,monoangle,monitor,[],[],fiodata,scantype)
				# assign one dictionary entry to each scan 
				self.scans[scanname] = onescan

	def loadelastic(self,scann,reps=None):
		self = self.loadscan(scann,reps,scantype='elastic')
	
	def loadedge(self,scann,reps=None,usemonoa=True):
		self = self.loadscan(scann,reps,usemonoa,scantype='edge')
	
	def loadlong(self,scann,reps=None):
		self = self.loadscan(scann,reps,scantype='long')

	def createrois(self,scannumbers):
		"""
		create rois object from this class and scannumbers
		"""
		rois_object = rois(self.scans,scannumbers)
		return rois_object		

	def getzoomrois(self,scannumbers,numrois=12,logscaling=True,colormap='jet'):
		rois_object = self.createrois(scannumbers)
		rois_object.getzoomrois(numrois,logscaling,colormap)
		self.rois = rois_object.rois

	def getautorois(self,scannumbers,kernel_size=5,colormap='jet'):
		rois_object = self.createrois(scannumbers)
		rois_object.getautorois(kernel_size,colormap)
		self.rois = rois_object.rois

	def getlinrois(self,scannumbers,numrois=9,logscaling=True,height=5,colormap='jet'):
		rois_object = self.createrois(scannumbers)
		rois_object.getlinrois(numrois,logscaling,height,colormap)
		self.rois = rois_object.rois

	def getrawdata(self):
		if not self.rois:
			print 'please define some ROIs first.'
			pass
		for scan in self.scans:
			if len(self.scans[scan].edfmats):
				print ("integrating "+scan)
				self.scans[scan].applyrois(self.rois)

	def getspectrum(self):
		# find the groups 
		allgroups = findgroups(self.scans)
		for group in allgroups:
			# self.makegroup(group)
			onegroup = makegroup(group)
			self.groups[onegroup.gettype()] = onegroup
		self.energy,self.signals,self.errors = appendscans(self.groups)

	def plotoneq(self,scannumbers,whichq,norm=True,reps=False):
		if not rep:
			for number in scannumbers:
				scanname = 'Scan%03d' % number
				if norm:
					plot(self.scans[scanname].energy,self.scans[scanname].signals[:,whichq]/self.scans[scanname].monitor)
				else:
					plot(self.scans[scanname].energy,self.scans[scanname].signals[:,whichq])
			show()
		else:
			for number in scannumbers:
				for rep in reps:
					scanname = 'Scan%03d' % number + 'r%01d' % rep 
					if norm:
						plot(self.scans[scanname].energy,self.scans[scanname].signals[:,whichq]/self.scans[scanname].monitor)
					else:
						plot(self.scans[scanname].energy,self.scans[scanname].signals[:,whichq])
			show()
