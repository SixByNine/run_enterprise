#!/usr/bin/env python
'''Pulsar glitch module for processing glitch data.
    Written by Y.Liu (yang.liu-50@postgrad.manchester.ac.uk).'''

from __future__ import print_function
# import argparse
# import sys
import os
import subprocess
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
# from matplotlib.ticker import FormatStrFormatter
# from scipy.optimize import curve_fit
# from astropy import coordinates as coord
# from astropy import units as u
# from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
# from mpl_toolkits.axes_grid1.inset_locator import mark_inset
# from mpl_toolkits.axes_grid.inset_locator import inset_axes



def mjd2sec(t, epoch):
    '''Convert MJD t into unit of seconds: x = time since epoch in seconds'''
    x = (t-epoch)*86400.0
    return x

def lin_fit(t, mjd, yoffs):
    '''Linear fit for discrete yoffs as a function of mjd at time t'''
    if t < mjd[0]:
        return yoffs[0]
    elif t > mjd[-1]:
        return yoffs[-1]
    else:
        for i in range(len(mjd)):
            if t > mjd[i] and t < mjd[i+1]:
                x1, x2 = mjd[i], mjd[i+1]
                y1, y2 = yoffs[i], yoffs[i+1]
                x=(t-x1)/(x2-x1)
                return (y2-y1)*x + y1

class Pulsar:
    '''A pulsar class consists basic info of pulsar'''
    def __init__(self, parfile, timfile):
        '''Initialize pulsar class with corresponding par file and tim file'''
        self.par = parfile
        self.tim = timfile
        self.psrn = None
        self.ra = None
        self.dec = None
        self.PB = None
        self.F0 = None
        self.F1 = None
        self.F2 = None
        self.start = None
        self.finish = None
        self.pepoch = None
        self.tidypar = None
        self.testpar = None
        self.truth = None
        self.epochfile = None
        self.datafile = None
        self.load_info()

    def load_info(self, tidy=False):
        '''Load basic info of pulsar class from par file and tim file'''
        if tidy==False or self.tidypar==None:
            loadpar = self.par
        else:
            loadpar = self.tidypar
        self.max_glitch = 0
        self.pglep = np.zeros(100)
        self.pglf0 = np.zeros(100)
        self.pglf1 = np.zeros(100)
        self.pglf2 = np.zeros(100)
        self.pglf0d = np.zeros(100)
        self.pgltd = np.ones(100)
        self.pglf0d2 = np.zeros(100)
        self.pgltd2 = np.ones(100)
        self.taug = 200*np.ones(100)
        self.pglf0ins = np.zeros(100)
        self.pglf0tg = np.zeros(100)
        self.redamp = None
        self.redgam = None
        with open(loadpar) as f1:
            for line in f1:
                line = line.strip()
                e = line.split()
                if e[0] == "PSRJ":
                    self.psrn = e[1]
                if e[0] == "RAJ":
                    self.ra = e[1]
                if e[0] == "DECJ":
                    self.dec = e[1]
                if e[0].startswith("GLEP_"):
                    i = int(e[0][5:])
                    self.pglep[i-1] = float(e[1])
                    self.max_glitch = max(i, self.max_glitch)
                if e[0].startswith("GLF0_"):
                    i = int(e[0][5:])
                    self.pglf0[i-1] = float(e[1])
                if e[0].startswith("GLF1_"):
                    i = int(e[0][5:])
                    self.pglf1[i-1] = float(e[1])
                if e[0].startswith("GLF2_"):
                    i = int(e[0][5:])
                    self.pglf2[i-1] = float(e[1])
                if e[0].startswith("GLTD_"):
                    i = int(e[0][5:])
                    self.pgltd[i-1] = float(e[1])
                if e[0].startswith("GLF0D_"):
                    i = int(e[0][6:])
                    self.pglf0d[i-1] = float(e[1])
                if e[0].startswith("GLTD2_"):
                    i = int(e[0][6:])
                    self.pgltd2[i-1] = float(e[1])
                if e[0].startswith("GLF0D2_"):
                    i = int(e[0][7:])
                    self.pglf0d2[i-1] = float(e[1])
                if e[0].startswith("GLF0(T="):
                    i = int(e[0].split('_')[-1])
                    numb = e[0].split('=')[1]
                    self.taug[i-1] = int(numb.split(')')[0])
                if e[0] == "PB":
                    self.PB = float(e[1])
                if e[0] == "F0":
                    self.F0 = float(e[1])
                if e[0] == "F1":
                    self.F1 = float(e[1])
                if e[0] == "F2":
                    self.F2 = float(e[1])
                if e[0] == "START":
                    self.start = float(e[1])
                if e[0] == "FINISH":
                    self.finish = float(e[1])
                if e[0] == "PEPOCH":
                    self.pepoch = float(e[1])
                if e[0] == "TNRedAmp":
                    self.redamp = float(e[1])
                if e[0] == "TNRedGam":
                    self.redgam = float(e[1])
        for i in range(self.max_glitch):
            self.pglf0ins[i] = self.pglf0[i]
            self.pglf0tg[i] = self.pglf1[i] * self.taug[i] + 0.5 * self.pglf2[i] * self.taug[i]**2
            if self.pglf0d[i] != 0:
                if self.taug[i] == 200:
                    self.taug[i] = int(self.pgltd[i])
                self.pglf0ins[i] += self.pglf0d[i]
                self.pglf0tg[i] += self.pglf0d[i]*(np.exp(-self.taug[i]/self.pgltd[i])-1)
            if self.pglf0d2[i] != 0:
                self.pglf0ins[i] += self.pglf0d2[i]
                self.pglf0tg[i] += self.pglf0d2[i]*(np.exp(-self.taug[i]/self.pgltd2[i])-1)
        # Load tim file info
        self.minmjds = 1e6
        self.maxmjds = 0
        self.toanum = 0
        self.toaseries = []
        with open(self.tim) as f2:
            for line in f2:
                e = line.split()
                if len(e) > 2.0 and e[0] != "C" and "-pn" in e:
                    self.toaseries.append(float(e[2]))
                    self.minmjds = min(self.minmjds, float(e[2]))
                    self.maxmjds = max(self.maxmjds, float(e[2]))
                    self.toanum += 1
        self.toaseries.sort()
        self.toaseries = np.array(self.toaseries)
        self.toaspan = self.maxmjds-self.minmjds
        self.cadence = self.toaspan/(self.toanum-1)

    def toa_gap(self, x, gap=0):
        '''Find huge gap in TOAs and create mask array for x in MJD'''
        self.toainterval = []
        for i, toai in enumerate(self.toaseries):
            if i == 0:
                self.toainterval.append(0)
            else:
                self.toainterval.append(toai-self.toaseries[i-1])
        self.toainterval = np.array(self.toainterval)
        if gap<=0:
            gap = 10*self.cadence
        print('Threshold of gap:', gap)
        maski = x<0
        counter = 1
        for i, value in enumerate(self.toainterval):
            if value >= gap:
                print('The No.%i gap (%f) in TOA is from %f to %f'%(counter, value, self.toaseries[i-1], self.toaseries[i]))
                counter += 1
                mask1 = self.toaseries[i-1] < x
                mask2 = x < self.toaseries[i]
                maski += mask1 * mask2
        x_mask = ma.array(x, mask=maski)
        x_inverse = ma.array(x, mask=~maski)
        return x_mask, x_inverse

    def delete_null(self):
        ''' Delete empty entries in glitch parameters'''
        self.pglep = self.pglep[:self.max_glitch]
        self.pglf0 = self.pglf0[:self.max_glitch]
        self.pglf1 = self.pglf1[:self.max_glitch]
        self.pglf2 = self.pglf2[:self.max_glitch]
        self.pglf0d = self.pglf0d[:self.max_glitch]
        self.pgltd = self.pgltd[:self.max_glitch]
        self.pglf0d2 = self.pglf0d2[:self.max_glitch]
        self.pgltd2 = self.pgltd2[:self.max_glitch]
        self.taug = self.taug[:self.max_glitch]
        self.pglf0ins = self.pglf0ins[:self.max_glitch]
        self.pglf0tg = self.pglf0tg[:self.max_glitch]

    def generate_truth(self):
        ''' Generate truth file'''
        self.truth = "trh_"+self.psrn+".txt"
        with open(self.truth, 'w') as f:
            for gi in range(self.max_glitch):
                idx = gi + 1
                f.write('GLEP_%i   %f\n'%(idx, self.pglep[gi]))
                if self.pglf0[gi] != 0:
                    f.write('GLF0_%i   %e\n'%(idx, self.pglf0[gi]))
                if self.pglf1[gi] != 0:
                    f.write('GLF1_%i   %e\n'%(idx, self.pglf1[gi]))
                if self.pglf2[gi] != 0:
                    f.write('GLF2_%i   %e\n'%(idx, self.pglf2[gi]))
                if self.pglf0d[gi] != 0:
                    f.write('GLF0D_%i   %e\n'%(idx, self.pglf0d[gi]))
                if self.pgltd[gi] != 0:
                    f.write('GLTD_%i   %f\n'%(idx, self.pgltd[gi]))
                if self.pglf0d2[gi] != 0:
                    f.write('GLF0D2_%i   %e\n'%(idx, self.pglf0d2[gi]))
                if self.pgltd2[gi] != 0:
                    f.write('GLTD2_%i   %f\n'%(idx, self.pgltd2[gi]))
                if self.pglf0[gi] != 0 and self.pglf1[gi] != 0:
                    #glf0_i = glf0 + glf0d + glf0d2
                    #glf0_T = glf1*t200*86400+glf0d*(np.exp(-t200/gltd)-1)+glf0d2*(np.exp(-t200/gltd2)-1)
                    f.write('GLF0(instant)_%i   %e\n'%(idx, self.pglf0ins[gi]))
                    f.write('GLF0(T=%d)_%i   %e\n'%(self.taug[gi], idx, self.pglf0tg[gi]))
            if all(p is not None for p in [self.redamp, self.redgam]):
                #alpha = redgam
                #P0 = ((redamp**2)/(12*np.pi**2))*(fc**(-alpha))
                f.write('TNRedAmp   %f\n'%self.redamp)
                f.write('TNRedGam   %f\n'%self.redgam)

    def tidy_glitch(self, chop=None):
        ''' Convert double recovery into (old) tempo2 style by creating a new glitch.
            Chop tim file and only keep TOAs from chop days before first glitch to chop days after the last glitch'''
        parlines = []
        glitches = {}
        with open(self.par) as f1:
            for line in f1:
                if line.startswith("GL"):
                    e = line.split()
                    pp = e[0].split("_")
                    i = int(pp[1])
                    param = pp[0]
                    if not i in glitches:
                        glitches[i] = {'turns':0}
                    if param == "GLEP":
                        glitches[i]['epoch'] = float(e[1])
                    #if param == "GLPH":
                        #glitches[i]['turns'] = round(float(e[1]))
                        #e[1] = "{}".format(float(e[1]) - glitches[i]['turns'])
                    glitches[i][param] = " ".join(e[1:])
                else:
                    parlines.append(line)
        #for ig in glitches:
            #print("glitch[{}] epoch {} turns {}".format(ig,glitches[ig]['epoch'],glitches[ig]['turns']))
        gg = sorted(glitches, key=lambda x: glitches[x]['epoch'])
        #print("Time-order is",gg)
        i = 1
        for ig in gg:
            glitches[ig]['newid'] = i
            i += 1
        timlines = []
        with open(self.tim) as f2:
            for line in f2:
                e = line.split()
                if "-pn" in e:
                    epoch = float(e[2])
                    ii = e.index("-pn")
                    pn = int(e[ii+1])
                    for ig in gg:
                        if epoch > glitches[ig]['epoch']:
                            pn -= glitches[ig]['turns']
                    newline = " ".join(e[:ii])+" -pn {} ".format(pn)+(" ".join(e[ii+2:]))
                    if isinstance(chop, int):
                        if self.pglep[0]-chop <= epoch <= self.pglep[-1]+chop:
                            timlines.append(" "+newline+"\n")
                    else:
                        timlines.append(" "+newline+"\n")
                elif len(e)>3:
                    epoch = float(e[3])
                    if isinstance(chop, int):
                        if self.pglep[0]-chop <= epoch <= self.pglep[-1]+chop:
                            timlines.append(line)
                    else:
                        timlines.append(" "+newline+"\n")
                else:
                    timlines.append(line)
        if isinstance(chop, int):
            self.tim = "chp_"+self.psrn+".tim"
        with open(self.tim,"w") as f2:
            f2.writelines(timlines)
        nextid = len(glitches)+1
        for ig in gg:
            if "GLF0D2" in glitches[ig]:
                print("Creating extra glitch {} to account for second recovery term".format(nextid))
                glitches[nextid] = {}
                glitches[nextid]['newid'] = nextid
                glitches[nextid]['turns'] = 0
                glitches[nextid]['epoch'] = glitches[ig]['epoch']
                glitches[nextid]['GLEP'] = glitches[ig]['epoch']
                glitches[nextid]['GLF0D'] = glitches[ig]['GLF0D2']
                glitches[nextid]['GLTD'] = glitches[ig]['GLTD2']
                del glitches[ig]['GLF0D2']
                del glitches[ig]['GLTD2']
                gg.append(nextid)
                nextid += 1
        self.tidypar = "tdy_"+self.par.split('_', 1)[1]
        with open(self.tidypar,"w") as f:
            f.writelines(parlines)
            for ig in gg:
                for param in glitches[ig]:
                    if param in ["epoch","turns","newid"]:
                        continue
                    f.write("{}_{} {}\n".format(param,glitches[ig]['newid'],glitches[ig][param]))


    def noglitch_par(self):
        ''' Make a copy of par file without glitch parameters'''
        parlines = []
        with open(self.par) as f:
            for line in f:
                if line.startswith("GLEP_"):
                    parlines.append(line)
                elif line.startswith("GL") or line.startswith("TN"):
                    continue
                elif any(line.startswith(ls) for ls in ["JUMP", "DM", "PM", "F0", "F1", "F2", "PX"]):
                    parlines.append(line)
                else:
                    e = line.split()
                    if len(e) > 2 and float(e[2]) == 1:
                        # Turn off tempo2 fitting of all other parameters
                        newline = e[0]+'   '+e[1]+'   '+e[3]+'\n'
                        parlines.append(newline)
                    else:
                        parlines.append(line)
        self.testpar = "tst_"+self.par.split('_', 1)[1]
        with open(self.testpar, "w") as newf:
            newf.writelines(parlines)

    def sf_create_global(self, leading, trailing, width, step):
        ''' Set fit timespan using a global par file'''
        try:
            os.remove("global.par")
        except OSError:
            pass
        with open("global.par", 'a') as glob:
            for ge in self.pglep:
                if ge is not None:
                    if leading < ge and trailing >= ge:
                        old_trail = trailing
                        trailing = ge - 0.01*step
                        if trailing - leading < width / 2.0:
                            leading = ge + 0.01*step
                            trailing = old_trail
            glob.write("START {} 1\n".format(leading))
            glob.write("FINISH {} 1".format(trailing))
        return leading, trailing

    def sf_run_fit(self, epoch, F0e=1e-7, F1e=1e-14, F2e=1e-17):
        '''Run tempo2 and fit for parameters'''
        epoch = str(epoch)
        command = ['tempo2', '-f', self.testpar, self.tim, '-nofit', '-global', 'global.par',
                   '-fit', 'F0', '-fit', 'F1', '-fit', 'F2', '-epoch', epoch]
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=None)
        while proc.poll() is None:
            line = proc.stdout.readline().decode('utf-8')
            fields = line.split()
            if len(fields) > 4:
                if fields[0] == "PEPOCH":
                    pepoch = fields[3]
                if fields[0] == "F0":
                    F0 = fields[3]
                    F0_e = fields[4]
                    if not 0 < abs(float(F0_e)) < F0e:
                        return None
                if fields[0] == "F1":
                    F1 = fields[3]
                    F1_e = fields[4]
                    if not 0 < abs(float(F1_e)) < F1e:
                        return None
                if fields[0] == "F2":
                    F2 = fields[3]
                    F2_e = fields[4]
                    if not 0 < abs(float(F2_e)) < F2e:
                        return None
        try:
            return pepoch, F0, F0_e, F1, F1_e, F2, F2_e
        except UnboundLocalError:
            return None

    def sf_main(self, width, step):
        '''Main function for stride fitting'''
        if step <= 0:
            step = 3*int(self.cadence)
        if width < 2*step:
            width = 2*step
        print('Stride fitting with box width %d and step size %d'%(width, step))
        first, last = self.minmjds, self.maxmjds
        leading = first
        trailing = first + width
        counter = 0
        self.epochfile = self.psrn+'_g'+str(self.max_glitch)+'_w'+str(int(width))+'_s'+str(int(step))+'_epoch.txt'
        with open(self.epochfile, 'w') as f1:
            while trailing <= last:
                leading, trailing = self.sf_create_global(leading, trailing, width, step)
                epoch = leading + ((trailing - leading)/2.0)
                print(leading, trailing, epoch, file=f1)
                counter += 1
                leading = first + counter*step
                trailing = first + width + counter*step
        starts, ends, fitepochs = np.loadtxt(self.epochfile, unpack=True)
        self.datafile = self.psrn+'_g'+str(self.max_glitch)+'_w'+str(int(width))+'_s'+str(int(step))+'_data.txt'
        with open(self.datafile, 'w') as f2:
            for i, (start_value, end_value) in enumerate(zip(starts, ends)):
                self.sf_create_global(start_value, end_value, width, step)
                epoch = fitepochs[i]
                out = self.sf_run_fit(epoch)
                os.remove("global.par")
                if out:
                    print(out[0], out[1], out[2], out[3], out[4], out[5], out[6], starts[i], ends[i], file=f2)

    def sf_calculate_data(self, save=True, plot=False):
        '''Load stride fitting results, plot stride fitting results and save to text files'''
        # sft, f0, f0e, f1, f1e, f2, f2e, mjds, mjdf = np.loadtxt(self.datafile, unpack=True)
        sft, f0, f0e, f1, f1e = np.loadtxt(self.datafile, usecols=(0, 1, 2, 3, 4), unpack=True)
        sfx = mjd2sec(sft, self.pepoch)
        sff1, sff2, sfdf2 = self.psr_taylor_terms(sfx)
        sfglf0, sfglf1, sfglf2, sfexp1, sfexp2, sfdglf1, sfdglf2, sfdexp = self.glitch_terms(sft)
        p1y = (f0 - self.F0 - sff1 - sff2)
        p2y = (f0 - self.F0 - sff1 - sff2 - sfglf0 - sfglf1 - sfglf2)
        p3y = (f0 - self.F0 - sff1 - sff2 - sfglf0 - sfglf1 - sfglf2 - sfexp1 - sfexp2)
        p4y = (f1 - self.F1 - sfdf2)
        p5y = (f1 - self.F1 - sfdf2 - sfdglf1 - sfdglf2 + sfdexp)
        if save:
            with open("panel1_{}.txt".format(self.psrn), "w") as file1:
                for i, value in enumerate(sft):
                    file1.write('%f   %e   %e   \n' % (value, 1e6*(p1y[i]), 1e6*f0e[i]))
                file1.close()
            with open("panel2_{}.txt".format(self.psrn), "w") as file2:
                for i, value in enumerate(sft):
                    file2.write('%f   %e   %e   \n' % (value, 1e6*(p2y[i]), 1e6*f0e[i]))
                file2.close()
            with open("panel3_{}.txt".format(self.psrn), "w") as file3:
                for i, value in enumerate(sft):
                    file3.write('%f   %e   %e   \n' % (value, 1e6*(p3y[i]), 1e6*f0e[i]))
                file3.close()
            with open("panel4_{}.txt".format(self.psrn), "w") as file4:
                for i, value in enumerate(sft):
                    file4.write('%f   %e   %e   \n' % (value, 1e10*(p4y[i]), 1e10*f1e[i]))
                file4.close()
            with open("panel5_{}.txt".format(self.psrn), "w") as file5:
                for i, value in enumerate(sft):
                    file5.write('%f   %e   %e   \n' % (value, 1e10*(p5y[i]), 1e10*f1e[i]))
                file5.close()
        if plot:
            plt.errorbar(sft, 1e6*p1y, yerr=1e6*f0e, marker='.', color='k', ecolor='k', linestyle='None')
            plt.show()
            plt.errorbar(sft, 1e6*p2y, yerr=1e6*f0e, marker='.', color='k', ecolor='k', linestyle='None')
            plt.show()
            plt.errorbar(sft, 1e6*p3y, yerr=1e6*f0e, marker='.', color='k', ecolor='k', linestyle='None')
            plt.show()
            plt.errorbar(sft, 1e10*p4y, yerr=1e10*f1e, marker='.', color='k', ecolor='k', linestyle='None')
            plt.show()
            plt.errorbar(sft, 1e10*p5y, yerr=1e10*f1e, marker='.', color='k', ecolor='k', linestyle='None')
            plt.show()

    def print_info(self, index=None):
        '''Print basic info of pulsar'''
        print("")
        print("Parameters in par file")
        print("Pulsar name:", self.psrn)
        print("Period epoch:", self.pepoch)
        print("Cadence", self.cadence)
        print("TOA length:", self.toanum)
        print("F0:", self.F0)
        print("F1:", self.F1)
        print("F2:", self.F2)
        print("")
        if isinstance(index, int) and 0 < index <= self.max_glitch:
            print("The {} glitch".format(index))
            print("Glitch epoch:", self.pglep[index-1])
            print("GLF0:", self.pglf0[index-1])
            print("GLF1:", self.pglf1[index-1])
            print("GLF2:", self.pglf2[index-1])
            print("GLF0D_1:", self.pglf0d[index-1], " - GLTD_1", self.pgltd[index-1])
            print("GLF0D2_1:", self.pglf0d2[index-1], " - GLTD2_1", self.pgltd2[index-1])
            print("Initial jump (GLFO(instant)_1):", self.pglf0ins[index-1])
            print("Decay jump (GLFO(tau_g)_1):", self.pglf0tg[index-1])
            print("tau_g:", self.taug[index-1])
            print("")
        else:
            for i in range(self.max_glitch):
                print("The {} glitch".format(i+1))
                print("Glitch epoch:", self.pglep[i])
                print("GLF0:", self.pglf0[i])
                print("GLF1:", self.pglf1[i])
                print("GLF2:", self.pglf2[i])
                print("GLF0D_1:", self.pglf0d[i], " - GLTD_1", self.pgltd[i])
                print("GLF0D2_1:", self.pglf0d2[i], " - GLTD2_1", self.pgltd2[i])
                print("Initial jump (GLFO(instant)_1):", self.pglf0ins[i])
                print("Decay jump (GLFO(tau_g)_1):", self.pglf0tg[i])
                print("tau_g:", self.taug[i])
                print("")

    def glitch_terms(self, t, gn=None):
        '''Calculate the glitch terms for MJD arrays t and the No.gn glitch in pulsar: x = time since glitch epoch in seconds'''
        glf0 = np.zeros_like(t)
        glf1 = np.zeros_like(t)
        glf2 = np.zeros_like(t)
        exp1 = np.zeros_like(t)
        exp2 = np.zeros_like(t)
        dglf1 = np.zeros_like(t)
        dglf2 = np.zeros_like(t)
        dexp = np.zeros_like(t)
        if isinstance(gn, int) and 0 <= gn < self.max_glitch:
            glep = self.pglep[gn]
            x = mjd2sec(t, glep)
            glf0[x > 0] += self.pglf0[gn]
            glf1[x > 0] += self.pglf1[gn] * x[x > 0]
            glf2[x > 0] += 0.5 * self.pglf2[gn] * x[x > 0]**2
            exp1[x > 0] += self.pglf0d[gn] * np.exp(-x[x > 0] / (self.pgltd[gn]*86400.0))
            exp2[x > 0] += self.pglf0d2[gn] * np.exp(-x[x > 0] / (self.pgltd2[gn]*86400.0))
            dglf1[x > 0] += self.pglf1[gn]
            dglf2[x > 0] += self.pglf2[gn] * x[x > 0]
            dexp[x > 0] += exp1[x > 0] / (self.pgltd[gn]*86400) + exp2[x > 0] / (self.pgltd2[gn]*86400)
        else:
            for i in range(self.max_glitch):
                glep = self.pglep[i]
                x = mjd2sec(t, glep)
                glf0[x > 0] += self.pglf0[i]
                glf1[x > 0] += self.pglf1[i] * x[x > 0]
                glf2[x > 0] += 0.5 * self.pglf2[i] * x[x > 0]**2
                exp1[x > 0] += self.pglf0d[i] * np.exp(-x[x > 0] / (self.pgltd[i]*86400.0))
                exp2[x > 0] += self.pglf0d2[i] * np.exp(-x[x > 0] / (self.pgltd2[i]*86400.0))
                dglf1[x > 0] += self.pglf1[i]
                dglf2[x > 0] += self.pglf2[i] * x[x > 0]
                dexp[x > 0] += exp1[x > 0] / (self.pgltd[i]*86400) + exp2[x > 0] / (self.pgltd2[i]*86400)
        return glf0, glf1, glf2, exp1, exp2, dglf1, dglf2, dexp

    def psr_taylor_terms(self, x):
        '''Calculate the pulsar taylor series terms for array x in second:
             x = time since period epoch in seconds'''
        tf1 = self.F1 * x
        tf2 = 0.5 * self.F2 * x * x
        tdf2 = self.F2 * x
        return tf1, tf2, tdf2

    def mask_glep(self, t, array):
        '''Mask data at GLEPs for MJD arrays t'''
        mask_index = []
        # for i, value in enumerate(t):
        for i in range(len(t)):   # using enumerate
            # for gi in range(self.max_glitch):
                # if t[i] <= self.pglep[gi] < t[i+1]:
            if any(t[i] <= glep < t[i+1] for glep in self.pglep):
                mask_index.append(i)   # or i+1
        mc = ma.array(array)
        mc[mask_index] = ma.masked
        return mc

    def pp_create_files(self):
        '''Call tempo2 to generate files for pulsar plots'''
        subprocess.call(["tempo2", "-output", "exportres", "-f", self.testpar, self.tim, "-nofit"])
        os.rename("out.res", "out2_{}.res".format(self.psrn))
        subprocess.call(["tempo2", "-output", "exportres", "-f", self.par, self.tim, "-writeres"])
        os.rename("param.labels", "param_{}.labels".format(self.psrn))
        os.rename("param.vals", "param_{}.vals".format(self.psrn))
        os.rename("cov.matrix", "cov_{}.matrix".format(self.psrn))
        os.rename("tnred.meta", "tnred_{}.meta".format(self.psrn))
        os.rename("out.res", "out_{}.res".format(self.psrn))
        os.rename("prefit.res", "prefit_{}.res".format(self.psrn))
        os.rename("postfit.res", "postfit_{}.res".format(self.psrn))
        os.rename("awhite.res", "awhite_{}.res".format(self.psrn))
        os.rename("design.matrix", "design_{}.matrix".format(self.psrn))
        os.rename("constraints.matrix", "constraints_{}.matrix".format(self.psrn))
        os.rename("adesign.matrix", "adesign_{}.matrix".format(self.psrn))

    def pp_calculate_data(self, start=None, finish=None):
        lab = np.loadtxt("param_{}.labels".format(self.psrn), dtype=np.str).T
        beta = np.loadtxt("param_{}.vals".format(self.psrn))
        meta = np.loadtxt("tnred_{}.meta".format(self.psrn), usecols=(1))
        omega, epoch = meta[0], meta[1]
        rx, ry, re = np.loadtxt("out_{}.res".format(self.psrn), usecols=(0, 5, 6), unpack=True)
        if isinstance(start, int) and isinstance(finish, int):
            t = np.linspace(start-0.5, finish+0.5, 1000)
        else:
            t = np.linspace(self.start-0.5, self.finish+0.5, 1000)
        y = np.zeros_like(t)
        cosidx=lab[1]=='param_red_cos'
        sinidx=lab[1]=='param_red_sin'
        maxwav=400
        nc=ns=0
        for i, (vcos, vsin) in enumerate(zip(cosidx, sinidx)):
            if vcos:
                nc+=1
                if nc > maxwav:
                    cosidx[i:] = False
            if vsin:
                ns+=1
                if ns > maxwav:
                    sinidx[i:] = False
            if nc>maxwav and ns>maxwav:
                break
        nwav = np.sum(sinidx)
        beta_mod = beta[np.logical_or(sinidx,cosidx)]
        M = np.zeros((2*nwav,len(t)))
        M2 = np.zeros((2*nwav,len(rx)))
        dM = np.zeros_like(M)
        ddM = np.zeros_like(M)
        with open("white_{}.asc".format(self.psrn),"w") as f:
            f.write("WAVE_OM {}\n".format(omega))
            f.write("WAVEEPOCH {}\n".format(epoch))
            for i in range(min(256,nwav)):
                f.write("WAVE{}  {}  {}\n".format(i+1, -beta_mod[i], -beta_mod[i+nwav]))
        #print("set up matricies")
        freqs=[]
        pwrs=np.power(beta_mod[:nwav],2) + np.power(beta_mod[nwav:],2)
        for i in range(nwav):
            omegai = omega*(i+1.0)
            M[i]        = np.sin(omegai * (t-epoch))
            M[i+nwav]   = np.cos(omegai * (t-epoch))
            freqs.append(365.25*omegai/2.0/np.pi)
            dM[i]      = -self.F0*omegai*M[i+nwav]/86400.0
            dM[i+nwav] = self.F0*omegai*M[i]/86400.0
            ddM[i]      = 1e15*self.F0*omegai*omegai*M[i]/(86400.0**2)
            ddM[i+nwav] = 1e15*self.F0*omegai*omegai*M[i+nwav]/(86400.0**2)                 
            M2[i]       = np.sin(omegai * (rx-epoch))
            M2[i+nwav]  = np.cos(omegai * (rx-epoch)) 
        #print("Do linear algebra")
        freqs=np.array(freqs)
        maxP=2*np.pi/omegai
        tt = mjd2sec(t, self.pepoch)
        M = M.T
        dM = dM.T
        ddM = ddM.T
        M2 = M2.T
        y = M.dot(beta_mod)
        yf = dM.dot(beta_mod)
        yd = ddM.dot(beta_mod)
        y_dat = M2.dot(beta_mod)
        yf2 = yf + (0.5*self.F2*tt*tt) #yf2 = yf + (0.5*F2*tt*tt + tt*F1 + F0)
        yd_model = np.zeros_like(yd)
        yd_model += 1e15*(self.F2*tt + self.F1)
        for i, ge in enumerate(self.pglep):
            gt = mjd2sec(t, ge)
            yf2[t>ge] +=  self.pglf0[i] + self.pglf1[i] * gt[t>ge] + 0.5 * self.pglf2[i] * gt[t>ge]**2 #bug
            yd_model[t>ge] += 1e15 * (self.pglf1[i] + self.pglf2[i] * gt[t>ge]) #bug
            if self.pglf0d[i] > 0:
                yf2[t>ge] += self.pglf0d[i] * np.exp(-(t[t>ge]-self.pglep[i])/self.pgltd[i])
                yf2[t>ge] += self.pglf0d2[i] * np.exp(-(t[t>ge]-self.pglep[i])/self.pgltd2[i])
                yd_model[t>ge] -= 1e15*self.pglf0d[i] * np.exp(-(t[t>ge]-self.pglep[i])/self.pgltd[i]) / (self.pgltd[i]*86400.0)
                yd_model[t>ge] -= 1e15*self.pglf0d2[i] * np.exp(-(t[t>ge]-self.pglep[i])/self.pgltd2[i]) / (self.pgltd2[i]*86400.0)
        yd2 = yd + yd_model
        ry += np.mean(y_dat-ry)
        with open("ifunc_{}.asc".format(self.psrn), "w") as f:
            for i, value in enumerate(t):
                f.write("{}  {}  {}\n".format(i+1, value, y[i])) #f.write("IFUNC{}  {}  {} {}\n".format(i+1, t[i], -y[i], 0))
        with open("deltanu_{}.asc".format(self.psrn), "w") as f:
            for i, value in enumerate(yf2):
                f.write("{} {}\n".format(t[i], value))
        with open("nudot_{}.asc".format(self.psrn), "w") as f:
            for i, value in enumerate(yd):
                f.write("{} {} {} {}\n".format(t[i], value, yd_model[i], yd2[i]))
        return rx, ry, re, y_dat, freqs, pwrs

    # def measure_prior
    # def set_prior


class Glitch(Pulsar):
    '''A Glitch subclass of Pulsar class consists glitch parameters and glitch model info'''
    def __init__(self, Pulsar, index):
        '''Initialize pulsar class with corresponding parameter file and TOA file'''
        self.parentpsr = Pulsar
        super().__init__(Pulsar.par, Pulsar.tim)
        if not (isinstance(index, int) and 0 < index <= Pulsar.max_glitch):
            index = Pulsar.max_glitch+1
            Pulsar.max_glitch += 1
        self.index = index-1
        self.inherit_pulsar()
        # self.create_new(self)

    def inherit_pulsar(self):
        '''Inherit glitch parameters from Pulsar info'''
        self.pglep = self.pglep[self.index]
        self.pglf0 = self.pglf0[self.index]
        self.pglf1 = self.pglf1[self.index]
        self.pglf2 = self.pglf2[self.index]
        self.pglf0d = self.pglf0d[self.index]
        self.pgltd = self.pgltd[self.index]
        self.pglf0d2 = self.pglf0d2[self.index]
        self.pgltd2 = self.pgltd2[self.index]
        self.taug = self.taug[self.index]
        self.pglf0ins = self.pglf0ins[self.index]
        self.pglf0tg = self.pglf0tg[self.index]
