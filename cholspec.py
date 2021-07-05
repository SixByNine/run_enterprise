#!/usr/bin/env python


import argparse



parser = argparse.ArgumentParser(description="Run 'cholspec' spectral analalysis code")
parser.add_argument('par')
parser.add_argument('tim')
parser.add_argument("--fmax", default=None, type=float, help="Max freq in per day (default median spacing)")
parser.add_argument("--fmax-limit", default=None, type=float, help="use median spacing up to this limit")
parser.add_argument("--t2model", default=None, type=str, help="Use this T2model file")
parser.add_argument("--refmodel", default=None, type=str, help="Overplot this reference model")
parser.add_argument("--reftn", default=None, nargs=2, type=float, help="Overplot reference TN model (A,gamma)")
parser.add_argument("--plot-t2", default=None, type=str, help="cholSpectra output")
parser.add_argument("--simmodel", default=None, type=str, help="Overplot this simulated model")
parser.add_argument("--savefig", default=None, type=str, help="Save Figure")

args=parser.parse_args()

import matplotlib
if args.savefig:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from enterprise.pulsar import Pulsar

from spectrum import cholspec

from t2model import *

psr=Pulsar(args.par,args.tim,drop_t2pulsar=False)
is_binary=False
logA=gamma=F0=F1=PB=0
with open(args.par) as f:
    for line in f:
        split_line = line.split()
        if len(split_line)==0:
            continue
        if split_line[0] == 'TNRedAmp':
            logA = float(split_line[1])
        elif split_line[0] == 'TNRedGam':
            gamma = float(split_line[1])
        elif split_line[0] == 'F0':
            F0 = float(split_line[1])
        elif split_line[0] == 'F1':
            F1 = float(split_line[1])
        elif split_line[0] == 'PB':
            is_binary=True
            PB=float(split_line[1])

psr.sort_data()
toa = psr.toas/86400.0 # convert back to days for sanity
ntoa=len(toa)

toa_spacing = np.diff(toa)
if args.fmax==None:
    args.fmax = 0.5/np.median(toa_spacing)
    if not args.fmax_limit is None:
        args.fmax=min(args.fmax_limit,args.fmax)
tspan = np.amax(toa)-np.amin(toa)
fc_yr = 365.25/tspan / 4.0
# compute CVM
print("Tspan={:.1f} days, using fc = {:.4g} per year".format(tspan,fc_yr))
print("Compute CVM")

orig_model_psd = None
if args.t2model:
    print("Using t2model parameters")
    red_cvm,t2funcs = read_t2_model(args.t2model,psr)
elif logA != 0 and gamma != 0:
    print("Using TN parameters")
    red_cvm = cholspec.getC(toa, logA, gamma, fc_yr)
else:
    red_cvm = np.zeros((ntoa,ntoa))

print("Compute spectrum fmax={:.4g} per day".format(args.fmax))
f, psd_f_yr3, complex_spec, loglike = cholspec.cholspec(psr.residuals,toa,psr.toaerrs,red_cvm,f_max=args.fmax)

f *= 365.25


fig=plt.figure()
plt.loglog(f,psd_f_yr3,color='blue')

if not args.plot_t2 is None:
    x,y = np.loadtxt(args.plot_t2,usecols=(0,1),unpack=True)
    plt.loglog(x,y,color='black',label=args.plot_t2)

f2 = np.logspace(np.log10(np.amin(f)),np.log10(np.amax(f)),1024)
if logA!=0 and gamma!=0:
    tnModel_psd = cholspec.PSD(f2/365.25, fc_yr, logA, gamma)
    plt.loglog(f2,tnModel_psd,color='red',label='TN')

if args.t2model:
    for fn in t2funcs:
        plt.loglog(f2,fn(f2),color='green',label='T2Model')

if args.refmodel:
    _, reffuncs = read_t2_model(args.refmodel, psr)
    for fn in reffuncs:
        plt.loglog(f2, fn(f2), color='lightgreen', label='Ref T2Model')

if args.simmodel:
    _, simfuncs = read_t2_model(args.simmodel, psr)
    for fn in simfuncs:
        plt.loglog(f2, fn(f2), color='black',ls='--', label='Simulation Input')

if not args.reftn is None:
    tnModel_psd = cholspec.PSD(f2/365.25, fc_yr, *args.reftn)
    plt.loglog(f2, tnModel_psd, color='pink', label='Ref TN')



plt.ylim(np.amin(psd_f_yr3)*0.25,np.amax(psd_f_yr3)*4)

plt.legend()

if args.savefig:
    plt.savefig(args.savefig)
else:
    plt.show()
