#!/usr/bin/env python

from __future__ import division, print_function
import matplotlib


matplotlib.use('Agg')
import numpy as np

import run_enterprise.model_components as model_components
import run_enterprise.samplers as samplers
import run_enterprise.qr_timing_model as qr_timing_model
import run_enterprise.output_plots as output_plots
import run_enterprise.time_domain_realisations as time_domain_realisations
from run_enterprise.modelfit import run_modelfit
from run_enterprise.psr_model import PsrModel

import argparse

from enterprise.pulsar import Pulsar
from enterprise.signals import gp_signals




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



    ## @todo: Add supermodel stuff


    ### make model
    print("Building model...")
    model = make_model(args,psr,parfile)


    psr_model = PsrModel(model, tim, par, psr)

    print("Run Model Fitting!")
    results, raw_results = run_modelfit(args, psr_model, parfile)
    if args.time_realisations > 0:
        time_domain_realisations.create_td_realisations(args,psr_model,raw_results)


def make_model(args,psr,parfile):
    if args.qr:
        model = qr_timing_model.QRTimingModel()
    else:
        model = gp_signals.TimingModel(use_svd=args.svd)

    for model_comp in model_components.all:
        m = model_comp.setup_model(args,psr,parfile)
        if m:
            model += m

    return model



