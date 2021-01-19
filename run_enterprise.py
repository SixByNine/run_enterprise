#!/usr/bin/env python

from __future__ import division, print_function
import matplotlib

import derived_parameters

matplotlib.use('Agg')
import numpy as np

import model_components
import samplers

import matplotlib
matplotlib.use('Agg')

import argparse
import os,sys

from enterprise.pulsar import Pulsar
from enterprise.signals import gp_signals
from enterprise.signals import signal_base

import output_plots



def run_enterprise(args):

    print(vars(args))

    par = args.par
    tim=args.tim

    print("Read pulsar data")
    psr=Pulsar(par,tim,drop_t2pulsar=False)

    orig_toas=psr.t2pulsar.toas()
    issorted=np.all(orig_toas[:-1] <= orig_toas[1:])

    with open(par) as f:
        parfile = f.readlines()

    if args.outdir == None:
        outdir = 'chains/' + psr.name + "/"
    else:
        outdir = args.outdir


    ## @todo: Add supermodel stuff


    ### make model
    model = make_model(args,psr,parfile)
    pta = signal_base.PTA(model(psr))

    print("PSR Name: ", psr.name)
    print("")
    print("'least-squares' fit parameters:")
    for i,p in enumerate(psr.fitpars):
        print("  {:2d}   {}".format(i+1, p))

    print("MCMC fit parameters:")
    for i,p in enumerate(pta.param_names):
        print("  {:2d}   {:40s} ".format(i+1, p))

    ### run sampler
    sampler_info=None

    for sampler in samplers.all:
        if sampler.activated(args):
            sampler_info = sampler.run_solve(args,pta,outdir)
            break

    if sampler_info is None:
        print("No sampler seemed to be enabled")
        sys.exit(1)

    ### Get the posterior 'chain' (or samples)
    samples, log_prob, weights, xtra = sampler.get_posteriors(args,sampler_info)
    scaled_samples = np.copy(samples)

    parnames = []
    iparam = 0
    ### convert parameters
    for param in pta.params:
        if not param.to_par is None:
            name, scaled_samples[:,iparam] = param.to_par(param.name, samples[:,iparam])
        else:
            name = param.name
        parnames.append(name)
        iparam+=1
    parnames=np.array(parnames)
    ### Compute any derived parameters
    derived, derived_names = derived_parameters.compute_derived_parameters(args,psr,par,parfile,scaled_samples,log_prob,weights,parnames)

    ### Write results
    output_plots.write_results(args,psr,par,parfile,scaled_samples,log_prob,weights,parnames,derived,derived_names)

    ### Make plots
    output_plots.make_plots(args,psr,scaled_samples,log_prob,weights, parnames,raw_results=sampler_info)


def make_model(args,psr,parfile):
    model = gp_signals.TimingModel()

    for model_comp in model_components.all:
        m = model_comp.setup_model(args,psr,parfile)
        if m:
            model += m

    return model



if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="Run 'enterprise' on a single pulsar")
    parser.add_argument('par')
    parser.add_argument('tim')
    parser.add_argument('--outdir', '-o', type=str, help="Output directory for chains etc")

    for model_comp in model_components.all:
        grp = parser.add_argument_group(model_comp.name, model_comp.argdec)
        model_comp.setup_argparse(grp)


    grp = parser.add_argument_group("Sampling options")
    samplers.std_args(grp)
    for sampler in samplers.all:
        grp = parser.add_argument_group(sampler.name, sampler.argdec)
        sampler.setup_argparse(grp)

    grp = parser.add_argument_group("Output options")
    output_plots.setup_argparse(grp)

    args=parser.parse_args()

    run_enterprise(args)
