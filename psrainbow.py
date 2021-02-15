#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

import argparse
import enterprise.pulsar

from spectrum import cholspec
from spectrum import fitspec


def rainbow(args):
    psr=enterprise.pulsar.Pulsar(args.par,args.tim,drop_t2pulsar=False)
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

    if logA==0 and gamma==0:
        print("No red noise model found")
        sys.exit(1)

    psr.sort_data()
    toa = psr.toas/86400.0 # convert back to days for sanity

    toa_spacing = np.diff(toa)
    if args.fmax==None:
        args.fmax = 1.0/np.median(toa_spacing)
    tspan = np.amax(toa)-np.amin(toa)
    fc_yr = 365.25/tspan / 4.0
    # compute CVM
    print("Tspan={:.1f} days, using fc = {:.4g} per year".format(tspan,fc_yr))
    print("Compute CVM")
    red_cvm = cholspec.getC(toa, logA, gamma, fc_yr)
    print("Compute spectrum fmax={:.4g} per day".format(args.fmax))
    f, psd_f_yr3, complex_spec, loglike = cholspec.cholspec(psr.residuals,toa,psr.toaerrs,red_cvm,f_max=args.fmax)
    orig_model_psd = cholspec.PSD(f,fc_yr,logA,gamma)

    ln_logz = fitspec.fit_lognorm(f, psd_f_yr3, complex_spec, orig_model_psd)
    red_qpnudot_logz = fitspec.fit_red_plus_qpnudot(f, psd_f_yr3, complex_spec, orig_model_psd)
    redpink_logz = fitspec.fit_redpink(f, psd_f_yr3, complex_spec, orig_model_psd)
    red_logz = fitspec.fit_red(f, psd_f_yr3, complex_spec, orig_model_psd)

    #print("delta logZ=",redpink_logz-red_logz)
    print("delta logZ=",red_qpnudot_logz-red_logz)






if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Run 'psrainbow' spectral analalysis code")
    parser.add_argument('par')
    parser.add_argument('tim')
    parser.add_argument("--fmax",default=None, help="Max freq in per day (default median spacing)")
    rainbow(parser.parse_args())

