from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, sys, re, time, h5py, glob
import numpy as np
from XRStools import xrs_utilities, xrs_scans
from  dateutil.parser import parse as dateparse

TINY = 1.e-7
MAX_FILESIZE = 100*1024*1024  # 100 Mb limit
COMMENTCHARS = '#;%*!$'
NAME_MATCH = re.compile(r"[a-zA-Z_][a-zA-Z0-9_]*(\.[a-zA-Z_][a-zA-Z0-9_]*)*$").match
VALID_SNAME_CHARS = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_'
VALID_NAME_CHARS = '.%s' % VALID_SNAME_CHARS
VALID_CHARS1 = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_'
BAD_FILECHARS = ';~,`!%$@$&^?*#:"/|\'\\\t\r\n (){}[]<>'
GOOD_FILECHARS = '_'*len(BAD_FILECHARS)
RESERVED_WORDS = ('and', 'as', 'assert', 'break', 'class', 'continue',
                'def', 'del', 'elif', 'else', 'eval', 'except', 'exec',
                'execfile', 'finally', 'for', 'from', 'global', 'if',
                'import', 'in', 'is', 'lambda', 'not', 'or', 'pass',
                'print', 'raise', 'return', 'try', 'while', 'with',
                'group', 'end', 'endwhile', 'endif', 'endfor', 'endtry',
                'enddef', 'True', 'False', 'None')


################################################################################
# Functions to get the header attributes
################################################################################
class Lerix:

    def __init__(self):
        self.scans         = {}
        self.keys          = {}

        self.elastic_scans = []
        self.elastic_name  = 'elastic'
        self.nixs_scans    = []
        self.nixs_name     = 'nixs'
        self.wide_scans    = []
        self.wide_name     = 'wide'
        self.scan_name     = []
        self.eloss_avg     = [] #elastic eloss average
        self.signals_avg   = [] #elastic signals average used to plot analyzer resolutions at the end
        self.energy        = []
        self.signals       = []
        self.errors        = []
        #self.groups        = {}
        self.tth           = []
        self.resolution    = {}
        self.E0            = []
        self.cenom         = []
        self.cenom_dict    = {}

    def isValidName(self,filename):
        "input is a valid name"
        if filename in RESERVED_WORDS:
            return False
            tnam = filename[:].lower()
            return NAME_MATCH(tnam) is not None

    def fixName(self,filename, allow_dot=True):
        "try to fix string to be a valid name"
        if self.isValidName(filename):
            return filename
        if self.isValidName('_%s' % filename):
            return '_%s' % filename
            chars = []
            valid_chars = VALID_SNAME_CHARS
            if allow_dot:
                valid_chars = VALID_NAME_CHARS
                for s in filename:
                    if s not in valid_chars:
                        s = '_'
                        chars.append(s)
                        filename = ''.join(chars)
                    # last check (name may begin with a number or .)
                    if not self.isValidName(filename):
                        filename = '_%s' % filename
                        return filename

    def getfloats(self,txt, allow_times=True):
        words = [w.strip() for w in txt.replace(',', ' ').split()]
        mktime = time.mktime
        for i, w in enumerate(words):
            val = None
            try:
                val = float(w)
            except ValueError:
                try:
                    val = mktime(dateparse(w).timetuple())
                except ValueError:
                    pass
                    words[i] = val
                    return words


    def colname(self,txt):
        return fixName(txt.strip().lower()).replace('.', '_')

    def attributes(self,fn):
        fh = open(fn, 'r')
        text = fh.read()
        text = text.replace('\r\n', '\n').replace('\r', '\n').split('\n')
        _labelline = None
        ncol = None
        data, footers, headers = [], [], []
        text.reverse()
        section = 'FOOTER'
        for line in text:
            line = line.strip()
            if len(line) < 1:
                continue

        # look for section transitions (going from bottom to top)
        if section == 'FOOTER' and not None in self.getfloats(line):
            section = 'DATA'
        elif section == 'DATA' and None in self.getfloats(line):
            section = 'HEADER'
            _labelline = line
            if _labelline[0] in COMMENTCHARS:
                _labelline = _labelline[1:].strip()
        # act of current section:
        if section == 'FOOTER':
            footers.append(line)
        elif section == 'HEADER':
            headers.append(line)
        elif section == 'DATA':
            rowdat  = self.getfloats(line)
            if ncol is None:
                ncol = len(rowdat)
            if ncol == len(rowdat):
                data.append(rowdat)

    # reverse header, footer, data, convert to array footers.reverse(headers.reverse(data.reverse(data = np.array(data).transpose()
    # try to parse attributes from header text
        header_attrs = {}
        for hline in headers:
            hline = hline.strip().replace('\t', ' ')
            if len(hline) < 1: continue
            if hline[0] in COMMENTCHARS:
                hline = hline[1:].strip()
            keywds = []
            if ':' in hline: # keywords in  'x: 22'
                words = hline.split(':', 1)
                keywds = words[0].split()
            elif '=' in hline: # keywords in  'x = 22'
                words = hline.split('=', 1)
                keywds = words[0].split()
            if len(keywds) == 1:
                key = self.colname(keywds[0])
                if key.startswith('_'):
                    key = key[1:]
                if len(words) > 1:
                    header_attrs[key] = words[1].strip()
        return header_attrs

    ################################################################################
    # Get Scan Info - returns the number, name, type and f.ext of the scan
    ################################################################################

    def scan_info(self, file):
        """get the scan number, name, type and file extention from the title of
        the scan assuming typical format e.g. elastic.0001, nixs.0001"""
        fn,fext = os.path.splitext(file)
        if str.lower(fn)==str.lower(self.nixs_name):
            scan_type = 'nixs'
        elif str.lower(fn)==str.lower(self.elastic_name):
            scan_type = 'elastic'
        elif str.lower(fn)==str.lower(self.wide_name):
            scan_type = 'wide'
        scan_number = fext.lstrip('.')
        scan_number = int(scan_number)
        scan_name = scan_type + '%04d' %scan_number
        return scan_number, scan_name, scan_type, fext, file

    def sort_dir(self, dir):
        """Returns a list of directory contents after filtering out scans without
        the correct format or size e.g. 'elastic.0001, nixs.0001 '"""
        dir_scans = []
        for file in os.listdir(dir):
            file_lc = str.lower(file)
            fn,fext = os.path.splitext(file_lc)
            if not file_lc.startswith('.'):
                    if fext.lstrip('.').isdigit():
                        if not os.stat(dir + '/' + file).st_size > 8000:
                            print("{} {}".format(">> >> Warning!! skipped empty scan (<8KB): ", file))
                            continue
                        elif not os.stat(dir + '/' + file).st_size < MAX_FILESIZE:
                            print("{} {}".format(">> >> Warning!! skipped huge scan (>100MB): ", file))
                            continue
                        else:
                            if fn==self.nixs_name:
                                dir_scans.append(file)
                            elif fn==self.elastic_name:
                                dir_scans.append(file)
                            elif fn==self.wide_name:
                                dir_scans.append(file)
        sorted_dir = sorted(dir_scans, key=lambda x: os.path.splitext(x)[1])
        return sorted_dir

    def isValidDir(self,dir):
        """Show that the scan directory is valid, that the directory holds a scan
        with the correct elastic name, nixs name and then let the user know if it
        has not found a wide scan. Returns True if valid directory."""
        if not os.path.isdir(dir):
            print('Check the directory you have supplied')
            return False
        elif not os.path.isfile(dir+'/'+self.elastic_name+'.0001'):
            print("The directory you supplied does not have a elastic.0001 file!!! \n If your elastic scan has a different name, please specify as: 'elastic_name'")
            return False
        elif not os.path.isfile(dir+'/'+self.nixs_name+'.0001'):
            print("The directory you supplied does not have a NIXS.0001 file!!! \n If your raman scan has a different name, please specify as: 'NIXS_name'")
            return False
        elif not os.path.isfile(dir+'/'+self.wide_name+'.0001'):
            print("No wide scans found. Continuing...")
            return True
        else:
            return True

    def plot_data(self,analyzer=False):
        """<classObj>.plot_data() Function that can be called to plot the eloss
        data for each channel and build an average by clicking a button.
        Requires matplotlib >2.1"""

        import matplotlib.pyplot as plt
        from matplotlib.widgets import CheckButtons, Cursor

        channels = []
        for analyzer in self.resolution:
            if analyzer.startswith('Analyzer'):
                if self.resolution[analyzer] < 1.0:
                    channels.append(int(analyzer.lstrip('Analyzer'))-1)

        data = np.average(self.signals[:,channels],axis=1)

        fig, ax = plt.subplots()
        ax.plot(self.eloss, data, lw=2)
        ax.set_xlabel('Energy Loss (eV)')
        ax.set_title('Plotting Raman Analysers')
        plt.subplots_adjust(left=0.3)

        rax = plt.axes([0.02, 0.1, 0.2, 0.8])
        check = CheckButtons(rax, ('Analyzer1', 'Analyzer2', 'Analyzer3','Analyzer4','Analyzer5',
        'Analyzer6','Analyzer7','Analyzer8','Analyzer9','Analyzer10','Analyzer11','Analyzer12'
        ,'Analyzer13','Analyzer14','Analyzer15','Analyzer16','Analyzer17','Analyzer18','Analyzer19'),
        (False, False, False, False, False, False, False, False, False, False, False, False, False
        ,False,False,False,False,False,False))


        def func(label):
            is_checked = []
            on = check.get_status()
            for ii in range(19):
                if on[ii]:
                    is_checked.append(ii)
                    #errors.append(np.average(noodle.errors[:,ii],axis=0))
            data = np.average(self.signals[:,is_checked],axis=1)
            ax.clear()
            ax.plot(self.eloss, data, lw=2)
            ax.autoscale(True)
            ax.set_xlabel('Energy Loss (eV)')
            ax.set_title('Plotting Raman Analysers')
            plt.draw()
            cursor = Cursor(ax, useblit=False, color='red', linewidth=2)
            plt.show()
            #print("{} {}".format("Average Error: ", ))

        check.on_clicked(func)
        plt.show()


    def write_H5scanData(self,dir,H5file,averaged='False'):
        """Writes all the scan information into a H5 file named after the sample name. inside
        this H5 directory scans are split into elastic and NIXS and the averaged scans. No support
        yet for wide scans"""
        sample_name = os.path.basename(dir)
        g = H5file.create_group(sample_name) #H5 subgroup with the name of the sample
        H5_ela = g.create_group('elastic') #H5 subgroup for elastics
        H5_xrs = g.create_group('XRS')     #H5 subgroup for NIXS
        all_scans = self.elastic_scans+self.nixs_scans
        for file in all_scans:
            scan_info = self.scan_info(file)
            if scan_info[2] == 'elastic':
                h5group = H5_ela.create_group(scan_info[1])
                h5group.create_dataset("energy",data=self.scans[scan_info[1]].energy)
                h5group.create_dataset("signals",data=self.scans[scan_info[1]].signals)
                h5group.create_dataset("errors",data=self.scans[scan_info[1]].errors)
                h5group.create_dataset("cenoms",data=self.scans[scan_info[1]].cenom)
            elif scan_info[2]=='nixs':
                h5group = H5_xrs.create_group(scan_info[1])
                h5group.create_dataset("energy",data=self.scans[scan_info[1]].energy)
                h5group.create_dataset("signals",data=self.scans[scan_info[1]].signals)
                h5group.create_dataset("eloss",data=self.scans[scan_info[1]].eloss)
                h5group.create_dataset("errors",data=self.scans[scan_info[1]].errors)
                h5group.create_dataset("tth",data=self.scans[scan_info[1]].tth)

        g.create_dataset("energy",data=self.energy)
        g.create_dataset("signals",data=self.signals)
        g.create_dataset("eloss",data=self.eloss)
        g.create_dataset("errors",data=self.errors)
        g.create_dataset("tth",data=self.tth)
        g.create_dataset("Mean Resolutions", data=np.array(self.resolution.items()))

        #Never forget to close an open H5 file!!!
        H5file.close()

    ################################################################################
    # Read Scan
    ################################################################################
    def get_cenoms(self, scan_info):
        """Internal Function to get the centre of mass of the elastic peak and
        the E0 for each elastic scan using XRStools"""
        cenom_list = []
        for analyzer in range(19): #The analyzer channels in the scan ASCII
            self.scans[scan_info[1]].cenom.append(xrs_utilities.find_center_of_mass(self.scans[scan_info[1]].energy,self.scans[scan_info[1]].signals[:,analyzer]))
        cenom_list.append(self.scans[scan_info[1]].cenom)
        self.cenom = [sum(a)/len(a) for a in zip(*cenom_list)]
        self.E0 = np.mean(self.cenom)/1e3

    def get_resolutions(self,scan_numbers):
        """Internal function to get the average resolution of each analyzer and
        average to give a self.resolution over the 19 analyzers. Returns a Dictionary
        of resolutions the mean and each analyser"""
        eloss_running_elastic = []
        signals_running_elastic = []

        if scan_numbers=='all':
            chosen_scans = []
            for number in range(len(self.elastic_scans)):
                scan_info = self.scan_info(self.elastic_scans[number])
                chosen_scans.append(scan_info[4])
        elif type(scan_numbers) is list:
            scan_numbers[:] = [x - 1 for x in scan_numbers] #scan 1 will be the 0th item in the list
            chosen_scans = []
            for number in scan_numbers:
                scan_info = self.scan_info(self.elastic_scans[number])
                chosen_scans.append(scan_info[4])
        else:
            print('scan numbers must be a list of the scans with correct length')
            return

        #populate lists with eloss and signals and then find the average over the whole range
        for file in chosen_scans:
            scan_info = self.scan_info(file)
            eloss_running_elastic.append(self.scans[scan_info[1]].eloss)
            signals_running_elastic.append(self.scans[scan_info[1]].signals)
        self.eloss_avg = np.array([sum(a)/len(a) for a in zip(*eloss_running_elastic)])
        self.signals_avg = np.array([sum(a)/len(a) for a in zip(*signals_running_elastic)])

        #take these average values and find the average FWHM for each analyzer and then find the total average
        for file in chosen_scans:
            resolution = []
            skipped = []
            for analyzer in range(19):
                try:
                    resolution.append(xrs_utilities.fwhm(self.eloss_avg, self.signals_avg[:,analyzer])[0])
                    self.resolution['Analyzer%s'%analyzer] = resolution[analyzer]
                except:
                    skipped.append(analyzer+1)
                    continue
        if len(skipped) > 1:
            print("{} {}".format("Skipped resolution for analyzer/s: ", list(set(skipped))))
        self.resolution['Resolution'] = round(np.mean(resolution),3)


    def read_scans(self,dir,file,valid_elastic='True'):
        """Internal Function that reads the APS data using numpy and finds the cenoms for each elastic
        scan ready to be passed to read_nixs to get eloss"""
        scan_info = self.scan_info(file)
        analyzers = [range(19)]
        try:
            scan_data = np.loadtxt(dir+'/'+file, comments='#')
        except:
            print("{} {}".format("NumPy failed to load scan name: ", file))
            pass
        self.scans[scan_info[1]].energy = np.array(scan_data[:,0]) #this format for np.repeat to work
        self.scans[scan_info[1]].signals = np.array(scan_data[:,5:24])
        self.scans[scan_info[1]].errors  = np.array(np.sqrt(np.absolute(self.scans[scan_info[1]].signals)))
        if scan_info[2]=='elastic':
            self.get_cenoms(scan_info)
            for analyzer in range(19):
                self.scans[scan_info[1]].eloss = np.subtract(self.scans[scan_info[1]].energy,self.scans[scan_info[1]].cenom[analyzer])
        elif scan_info[2]=='nixs' or scan_info[2]=='wide':
            #create empty array with shape energy.v.signals
            eloss = np.zeros(self.scans[scan_info[1]].signals.shape)
            self.scans[scan_info[1]].tth = list(range(9,180,9)) #assign tth to each scan
            self.tth = list(range(9,180,9)) #assign tth to self
            if valid_elastic=='True':
                for analyzer in range(19):
                    self.scans[scan_info[1]].eloss = np.subtract(self.scans[scan_info[1]].energy,self.scans['elastic%04d'%scan_info[0]].cenom[analyzer])
            elif valid_elastic=='False':
                for analyzer in range(19):
                    self.scans[scan_info[1]].eloss = np.subtract(self.scans[scan_info[1]].energy,self.cenom[analyzer])
            else:
                print('valid_elastic is a boolean')
                sys.exit()

    def average_scans(self,scan_numbers='all'):
        """Function to calculate the average eloss, energy, signals and errors over
        all the read scans (default) or over a list of scans e.g. [1,3,5]"""
        energy_running = []
        signals_running = []
        eloss_running = []
        errors_running = []
        if scan_numbers=='all':
            for file in self.nixs_scans:
                scan_info = self.scan_info(file)
                energy_running.append(self.scans[scan_info[1]].energy)
                signals_running.append(self.scans[scan_info[1]].signals)
                eloss_running.append(self.scans[scan_info[1]].eloss)
                errors_running.append(self.scans[scan_info[1]].errors)
            self.energy = np.array([sum(a)/len(a) for a in zip(*energy_running)])
            self.signals = np.array([sum(a)/len(a) for a in zip(*signals_running)])
            self.eloss = np.array([sum(a)/len(a) for a in zip(*eloss_running)])
            self.errors = np.array([sum(a)/len(a) for a in zip(*errors_running)])

        elif type(scan_numbers) is list:
            scan_numbers[:] = [x - 1 for x in scan_numbers] #scan 1 will be the 0th item in the list
            chosen_scans = []
            for number in scan_numbers:
                scan_info = self.scan_info(self.nixs_scans[number])
                chosen_scans.append(scan_info[1])
            print("{} {}".format("Averaging scan numbers: ", chosen_scans))
            for scan in chosen_scans:
                energy_running.append(self.scans[scan_info[1]].energy)
                signals_running.append(self.scans[scan_info[1]].signals)
                eloss_running.append(self.scans[scan_info[1]].eloss)
                errors_running.append(self.scans[scan_info[1]].errors)
            self.energy = np.array([sum(a)/len(a) for a in zip(*energy_running)])
            self.signals = np.array([sum(a)/len(a) for a in zip(*signals_running)])
            self.eloss = np.array([sum(a)/len(a) for a in zip(*eloss_running)])
            self.errors = np.array([sum(a)/len(a) for a in zip(*errors_running)])

        else:
            print("scan_numbers must be blank, 'all' or a list of scan numbers e.g.[1,3,5]")
            sys.exit()

    ################################################################################
    # Begin the reading
    ################################################################################
    def load_scan(self,dir,nixs_name='NIXS',wide_name='wide',elastic_name='elastic',scan_numbers='all',H5=False,H5path=None):
        """Function to load scan data from a typical APS 20ID Non-Resonant inelastic
        X-ray scattering experiment. With data in the form of elastic.0001, allign.0001
        and NIXS.0001. Function reteurns the averaged energy loss, signals, errors, E0
        and 2theta angles for the scans in the chosen directory."""

        self.nixs_name = str.lower(nixs_name)
        self.wide_name = str.lower(wide_name)
        self.elastic_name = str.lower(elastic_name)

        #check dir location
        if not self.isValidDir(dir):
            print('IO Error - sorry about that!')
            sys.exit()
        else:
            pass

        #sort the directory so that scans are in order, determine number of scans
        #open list to be filled with the elastic/nixs scan names
        sorted_dir = self.sort_dir(dir) #returns all the scan names
        number_of_scans = len(glob.glob(dir+'/'+self.nixs_name+'*'))-1 #number of nixs scans
        self.elastic_scans = []
        self.nixs_scans = []
        self.wide_scans = []
        #self.keys = {"eloss":np.array, "energy":np.array, "signals":np.array, "errors":np.array,"E0":np.float, "tth":np.array} #,"resolution":array }

        #split scans into NIXS and elastic and begin instance of XRStools scan class for each scan
        for file in sorted_dir:
                scan_info = self.scan_info(file)
                scan = xrs_scans.Scan()

                if scan_info[2]=='elastic':
                    self.elastic_scans.append(file)
                    self.scans[scan_info[1]] = scan #self.scans {} in class _init_
                    self.scans[scan_info[1]].scan_type = scan_info[2]
                    self.scans[scan_info[1]].scan_number = scan_info[0]

                if scan_info[2]=='nixs':
                    self.nixs_scans.append(file)
                    self.scans[scan_info[1]] = scan #self.scans {} in class _init_
                    self.scans[scan_info[1]].scan_type = scan_info[2]
                    self.scans[scan_info[1]].scan_number = scan_info[0]

                if scan_info[2]=='wide':
                    self.wide_scans.append(file)
                    self.scans[scan_info[1]] = scan #self.scans {} in class _init_
                    self.scans[scan_info[1]].scan_type = scan_info[2]
                    self.scans[scan_info[1]].scan_number = scan_info[0]

                else:
                    continue

        #read elastic scans first to calculate cenom
        # reading all the elastics in the file to improve E0 accuracy and get a
        # good grasp on the scan resolution
        for file in self.elastic_scans:
            scan_info = self.scan_info(file)
            print("{} {}".format("Reading elastic scan named: ", file))
            self.read_scans(dir,file)

        print('I always read all the Elastic scans to improve Resolution and E0 Accuracy\n >> Type <class>.resolution to see the analyzer resolutions.')

        #Read NIXS scans - if there isn't a corresponing elastic scan, subtract the
        #running average cenoms and tell the user.
        for file in self.nixs_scans:
            scan_info = self.scan_info(file)
            corresponding_elastic = dir+'/'+self.elastic_name+scan_info[3]
            if os.path.isfile(corresponding_elastic):
                print("{} {}".format("Reading NIXS scan name: ", file))
                self.read_scans(dir,file)
            elif not os.path.isfile(corresponding_elastic):
                print("{} {} {}".format(">> >> WARNING:", scan_info[1],"has no corresponding elastic - finding eloss by average elastic values!"))
                print("{} {}".format("Reading NIXS scan named: ", file))
                self.read_scans(dir,file,valid_elastic='False')

        for file in self.wide_scans:
            scan_info = self.scan_info(file)
            print("{} {}".format("Reading wide scan named: ", file))
            self.read_scans(dir,file)

        #call function to calculate the average values over the scans - all by default
        self.average_scans(scan_numbers)
        self.get_resolutions(scan_numbers)

        #if the user asks, call function to write all info to H5 file
        if H5:
            if H5path==None:
                H5path = dir
            elif H5path:
                if os.path.isdir(H5path):
                    H5path = H5path
                else:
                    print('H5 path directory does not exist!')

            H5name = '20ID_APS_data.H5'
            saveloc = H5path+'/'+H5name

            if os.path.isfile(saveloc):
                H5file = h5py.File(saveloc, "a")
            else:
                H5file = h5py.File(saveloc, "w")

            self.write_H5scanData(dir,H5file)
            print("{} {}".format("Wrote scan data to H5 file: ", saveloc))

        #let the user know the program has finished
        print('Finished Reading!')

"""To do:
1) Check if the cenoms are close to the E0 specified in the ASCII header, and if not, do not average that scan_name
2) Make the H5 file location more interactive and allow many different samples to be read into the H5file - e.g. check if it exists and if so
write into it. DONE
3) read wide scans - Proving hard, can't do this until I've fixed the header attrbs reading, to pull out the NIXS scan region
4) Header reading and H5 attributes
5) Maybe a file size check to make sure the input file isn't crazy - DONE
6) If file size is less than 5KB, then ignore it from the list - DONE
"""
