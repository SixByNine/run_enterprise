#!/usr/bin/env python

from __future__ import print_function
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit # mod ly

parser = argparse.ArgumentParser(description='Generating data files fot plot routine of fitting over a glitch.')
parser.add_argument('-p', '--parfile', help='Path to ephemeris', required=True)
parser.add_argument('-s', '--stride', help='Stride data text file', required=True)
args = parser.parse_args()
par = args.parfile
strdat = args.stride

#from matplotlib.ticker import FormatStrFormatter
#plt.rcParams['pdf.fonttype'] = 42
#plt.rcParams['ps.fonttype'] = 42

#par="test_final.par"

t, f0, f0e, f1, f1e, f2, f2e, mjds, mjdf = np.loadtxt(strdat, unpack=True)

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


glep = pglep[0] # add

#Time since period epoch
x = (t-pepoch)*86400.0 # add

#Time since glitch epoch
xx = (t-glep)*86400.0 # add

#First derivative term of taylor series
tf1 = F1 * x

#second derivative term of taylor series
#f = f0 + f1 t + 0.5 f2 t^2
tf2 = 0.5 * x * x * F2

#Permanent change term (constant)
glf0 = np.zeros_like(xx)
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

def lin(x, a, b):
    return a*x+b

def qua(x, a, b, c):
    return a*x**2+b*x+c

#idxlen = len(t[t<glep])
#opt_pre, cov_pre = curve_fit(lin, t[:idxlen], f0[:idxlen], sigma=f0e[:idxlen])
#print(np.mean(f0), F0, opt_pre[0])
#slope = opt_pre[0]*(t-pepoch)

plt.errorbar(t, 1e6*(f0 - F0 - tf1 - tf2), yerr=1e6*f0e, marker='.', color='k', ecolor='k', linestyle='None')
plt.show()
plt.errorbar(t, 1e6*(f0 - F0 - tf1 - tf2 - glf0 - glf1), yerr=1e6*f0e, marker='.', color='k', ecolor='k', linestyle='None')
plt.show()
plt.errorbar(t, 1e6*(f0 - F0 - tf1 - tf2 - glf0 - glf1 - exp1), yerr=1e6*f0e, marker='.', color='k', ecolor='k', linestyle='None')
plt.show()

#for i in range(0, len(t)):
    #print(t[i], 1e6*(f0[i] - F0 - tf1[i] - tf2[i]), 1e6*f0e[i])
    #print(t[i], 1e6*(f0[i] - F0 - tf1[i] - tf2[i] - glf0[i] - glf1[i]), 1e6*f0e[i])
    #print(t[i], 1e6*(f0[i] - F0 - tf1[i] - tf2[i] - glf0[i] - glf1[i] - exp1[i]), 1e6*f0e[i])
# output to panel.txt

with open("panel1_{}.txt".format(psrn),"w") as file1:
    for i in range(0, len(t)):
        file1.write('%f   %e   %e   \n'%(t[i], 1e6*(f0[i] - F0 - tf1[i] - tf2[i]), 1e6*f0e[i]))
    file1.close()

with open("panel2_{}.txt".format(psrn),"w") as file2:
    for i in range(0, len(t)):
        file2.write('%f   %e   %e   \n'%(t[i], 1e6*(f0[i] - F0 - tf1[i] - tf2[i] - glf0[i] - glf1[i]), 1e6*f0e[i]))
    file2.close()

with open("panel3_{}.txt".format(psrn),"w") as file3:
    for i in range(0, len(t)):
        file3.write('%f   %e   %e   \n'%(t[i], 1e6*(f0[i] - F0 - tf1[i] - tf2[i] - glf0[i] - glf1[i] - exp1[i] -exp2[i]), 1e6*f0e[i]))
    file3.close()

