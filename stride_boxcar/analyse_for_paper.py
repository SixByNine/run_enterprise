from astropy import coordinates as coord
from astropy import units as u
import argparse
import numpy as np
import scipy
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit # mod ly
import sys
from matplotlib.ticker import FormatStrFormatter
import numpy.ma as ma
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid.inset_locator import inset_axes

parser = argparse.ArgumentParser(description='Plot routine for fitting over a glitch.')
parser.add_argument('-p', '--parfile', help='Path to ephemeris', required=True)
parser.add_argument('-s', '--stride', help='Stride data text file', required=True)
args = parser.parse_args()
par = args.parfile
strdat = args.stride

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

#plt.rcParams["figure.figsize"] = [7.5,10.6]

pglep=np.zeros(100)
pglf0=np.zeros(100)
pglf1=np.zeros(100)

pglf0d=np.zeros(100)
pgltd=np.ones(100)
pglf0d2=np.zeros(100)
pgltd2=np.ones(100)
max_glitch=0

#Open par file and extract parameters
with open(par) as f:
    for line in f:
        line = line.strip()
        e=line.split()
        if e[0] == "PSRJ":
            psrn=e[1]
        if e[0].startswith("GLEP_"):
            i=int(e[0][5:])
            pglep[i-1] = float(e[1])
            max_glitch = max(i,max_glitch)
        if e[0].startswith("GLF0_"):
            i=int(e[0][5:])
            pglf0[i-1] = float(e[1])
        if e[0].startswith("GLF1_"):
            i=int(e[0][5:])
            pglf1[i-1] = float(e[1])
        if e[0].startswith("GLTD_"):
            i=int(e[0][5:])
            pgltd[i-1] = float(e[1])
        if e[0].startswith("GLF0D_"):
            i=int(e[0][6:])
            pglf0d[i-1] = float(e[1])
# mod ly
        if e[0].startswith("GLTD2_"):
            i=int(e[0][6:])
            pgltd2[i-1] = float(e[1])
        if e[0].startswith("GLF0D2_"):
            i=int(e[0][7:])
            pglf0d2[i-1] = float(e[1])
# mod ly
        if e[0] == "F0":
            F0=float(e[1])
        if e[0] == "PB":
            PB=float(e[1])
        if e[0] == "F1":
            F1=float(e[1])
        if e[0] == "F2":
            F2=float(e[1])
        if e[0] == "START":
            start=float(e[1])
        if e[0] == "FINISH":
            finish=float(e[1])
        if e[0] == "PEPOCH":
            pepoch=float(e[1])

print("")
print("Parameters in par file")
print("F0:", F0)
print("F1:", F1)
print("F2:", F2)
print("")

# for i in range(max_glitch):

print("Glitch epoch:", pglep[0])
print("GLF0:", pglf0[0])
print("GLF1:", pglf1[0])
print("")
print("GLF0D_1:", pglf0d[0], " - GLTD_1", pgltd[0], "GLF0D2_1:", pglf0d2[0], " - GLTD2_1", pgltd2[0])
print("Initial jump:", pglf0d[0]+pglf0d2[0]+pglf0[0])

#print("GLF0D_2:", pglf0d[1], " - T1", pgltd[1])
#print("GLF0D_3:", pglf0d[2], " - T1", pgltd[2])

#print("Initial jump:", np.sum(pglf0d)+pglf0[0])

dat="deltanu_{}.asc".format(psrn)
#par="test_final.par"

t, nu = np.loadtxt(dat,unpack=True)
t, no_nudot, nudot, nudot_model = np.loadtxt("nudot_{}.asc".format(psrn), unpack=True)
#this_mjd, this_dnu, this_dnu_err = np.loadtxt("mjd_dnu_dnuerr_nof2.txt", unpack=True)
panel1_mjd, panel1_nu, panel1_err = np.loadtxt("panel1_{}.txt".format(psrn), unpack=True)
panel2_mjd, panel2_nu, panel2_err = np.loadtxt("panel2_{}.txt".format(psrn), unpack=True)
panel3_mjd, panel3_nu, panel3_err = np.loadtxt("panel3_{}.txt".format(psrn), unpack=True)
str_mjd, str_nudot, str_err , str_nu2dot, str_2err = np.loadtxt(strdat, unpack=True, usecols=[0,3,4,5,6])
#print(np.max(str_err)) ; sys.exit(9)

glep = pglep[0]  # mod 
preepochs = []
for i in range(0, len(t)):
    if t[i] < glep:
        preepochs.append(t[i])
mask_len = len(preepochs) # mask length is the length of pre-glitch data

def glexp(xx,td,f0d):
    '''
    xx = time since glitch epoch
    td = decay timescale
    f0d = decay amplitude
    '''

    ee = np.zeros_like(xx)
    tau1=td*86400.0
    ee[xx>0] = f0d * np.exp(-xx[xx>0]/tau1)
    return ee

#Time since period epoch
x = (t-pepoch)*86400.0

#Time since glitch epoch
xx = (t-glep)*86400.0  # mod

#First derivative term of taylor series
tf1 = F1 * x

#second derivative term of taylor series
#f = f0 + f1 t + 0.5 f2 t^2
tf2 = 0.5 * x * x * F2

#Permanent change term (constant)
glf0=np.zeros_like(xx)
glf0[xx>0] = pglf0[0]
#glf0[t>glep] =  pglf0[0]

#GLF1 term
glf1 = np.zeros_like(xx)
glf1[xx>0] = xx[xx>0] * pglf1[0]

#transient terms
exp1=glexp(xx,pgltd[0],pglf0d[0])
if pglf0d2[0] != 0:
    exp2=glexp(xx,pgltd2[0],pglf0d2[0])
else:
    exp2=np.zeros_like(xx)
#if max_glitch>2:
#    exp3=glexp(xx,pgltd[2],pglf0d[2])
#else:
#    exp3=np.zeros_like(xx)

# mod ly
mdglf1 = np.zeros_like(xx)
mdglf1[xx>0] = pglf1[0]

strx = (str_mjd-glep)*86400
#str2nu = np.zeros_like(strx)
#str2nu[strx>0] = strx[strx>0] * str_nu2dot

strglf1 = np.zeros_like(strx)
strglf1[strx>0] = pglf1[0]

strexp1=glexp(strx,pgltd[0],pglf0d[0])
# mod ly

#frequency evolution with second derivative and spin-down change subtracted
numod = nu - tf2 # why only f2?
mc = ma.array(numod)
mc[mask_len] = ma.masked  # mc: nu subtract f2 with pre-glitch data masked
fig, ax = plt.subplots(nrows=5, ncols=1, figsize=(7.5, 12.8)) #10.6
#fig = plt.figure(figsize=(7.5, 10.6))
#gs = gridspec.GridSpec(4, 1)
#gs.update(wspace=0.025, hspace=0.0001)

plt.subplot(511)
plt.plot(t-glep,1e6*mc, 'k-', zorder=2)
plt.errorbar(panel1_mjd - glep, panel1_nu, yerr=panel1_err, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1)
#plt.errorbar(this_mjd - glep, this_dnu, yerr=this_dnu_err, marker='.', color='k', ecolor='k', linestyle='None', alpha=0.2, markersize=4)
plt.ylabel(r'$\delta \nu$ ($\mu$Hz)', fontsize=15)
plt.axvline(0, color='k', linestyle='dotted', alpha=0.3, linewidth=2)
#plt.axvline(0, color='k', linestyle='dashed', alpha=0.5)
#plt.xlim([-30, 50])
#plt.xlim([-10, 10])

frame = plt.gca()
frame.axes.xaxis.set_ticklabels([])

plt.subplot(512)
plt.plot(t-glep,1e6*(mc-glf0-glf1), 'k-', zorder=2)
plt.errorbar(panel2_mjd - glep, panel2_nu, yerr=panel2_err, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1, markersize=4)
#plt.xlim([-30, 30])
#plt.xlim([-30, 50])
#plt.xlim([-10, 10])
#plt.axvline(0, color='k', linestyle='dashed', alpha=0.5)
plt.axvline(0, color='k', linestyle='dotted', alpha=0.3, linewidth=2)
plt.ylabel(r'$\delta \nu$ ($\mu$Hz)', fontsize=15, labelpad=15)

frame = plt.gca()
frame.axes.xaxis.set_ticklabels([])

plt.subplot(513)
plt.ylabel(r'$\delta \nu$ ($\mu$Hz)', fontsize=15)
#plt.ylim([-0.3, 0.08])
#plt.xlim([-30, 50])
#plt.xlim([-10, 10])
plt.plot(t-glep,1e6*(mc-glf0-glf1-exp1), 'k-', zorder=2)
plt.errorbar(panel3_mjd - glep, panel3_nu, yerr=panel3_err, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1, markersize=4)
plt.axvline(0, color='k', linestyle='dotted', alpha=0.3, linewidth=2)
#plt.axvline(0, color='k', linestyle='dashed', alpha=0.5)
frame = plt.gca()
frame.axes.xaxis.set_ticklabels([])
#frame.axes.yaxis.set_ticklabels([])


#INSET
#inset_axes1 = inset_axes(ax[2], 
#                    width="40%", # width = 30% of parent_bbox
#                    height=0.8, # height : 1 inch
#                    loc=4, 
#                    borderpad=1
#                    )
#plt.plot(t-glep,1e6*(mc-glf0-glf1-exp1), 'k-', zorder=2)
#plt.errorbar(panel3_mjd - glep, panel3_nu, yerr=panel3_err, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1, markersize=4)
#plt.xlim([0, 50]) 
#frame = plt.gca()
#   #frame.axes.xaxis.set_ticklabels([])
#frame.axes.yaxis.set_ticklabels([])

#END INSET


# mod ly
#sf2 = x * F2
#nudot_mod = nudot_model - sf2
# mod ly
md = ma.array(nudot_model) 
md[mask_len] = ma.masked # md: nudot model with pre-glitch data masked


plt.subplot(514)
plt.plot(t-glep, md/1e5, 'k-', zorder=2)
plt.errorbar(str_mjd - glep, str_nudot/1e-10, yerr=str_err/1e-10, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1, markersize=4)
#plt.ylim([-3.725, -3.62])
#plt.xlim([-30, 50])
#plt.xlim([-10, 10])
plt.axvline(0, color='k', linestyle='dotted', alpha=0.3, linewidth=2)
plt.ylabel(r'$\dot{\nu}$ ($10^{-10}$ Hz s$^{-1}$)', fontsize=15)
plt.xlabel("Days since glitch epoch", fontsize=15)
plt.subplots_adjust(wspace=0, hspace=0.002)

frame = plt.gca()
frame.axes.xaxis.set_ticklabels([])


#INSET
#inset_axes2 = inset_axes(ax[3], 
#                    width="40%", # width = 30% of parent_bbox
#                    height=0.6, # height : 1 inch
#                    loc=1,
#                    borderpad=1)
#plt.plot(t-glep, md/1e5, 'k-', zorder=2)
#plt.errorbar(str_mjd - glep, str_nudot/1e-10, yerr=str_err/1e-10, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1, markersize=4)
#plt.xlim([0,50])
#frame = plt.gca()
##frame.axes.xaxis.set_ticklabels([])
#frame.axes.yaxis.set_ticklabels([])

#plt.plot(t-glep, md/1e5, 'k-', zorder=2)
#plt.errorbar(str_mjd - glep, str_nudot/1e-10, yerr=str_err/1e-10, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1, markersize=4)
#plt.xlim([-0.001, 5]) 
#plt.xticks([])
#plt.yticks([])
#END INSET


plt.subplot(515)
plt.plot(t-glep, md/1e5-(mdglf1-exp1/(pgltd[0]*86400))/1e-10, 'k-', zorder=2)
plt.errorbar(str_mjd - glep, (str_nudot-strglf1+strexp1/(pgltd[0]*86400))/1e-10, yerr=str_err/1e-10, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1, markersize=4)
#plt.ylim([-3.725, -3.62])
#plt.xlim([-30, 50])
#plt.xlim([-10, 10])
plt.axvline(0, color='k', linestyle='dotted', alpha=0.3, linewidth=2)
plt.ylabel(r'$\dot{\nu}$ ($10^{-10}$ Hz s$^{-1}$)', fontsize=15)
plt.xlabel("Days since glitch epoch", fontsize=15)
plt.subplots_adjust(wspace=0, hspace=0.002)

#frame = plt.gca()
#frame.axes.xaxis.set_ticklabels([])


plt.tight_layout()
plt.savefig("nu_nudot_gp_{}.pdf".format(psrn), format='pdf', dpi=400)
plt.show()
