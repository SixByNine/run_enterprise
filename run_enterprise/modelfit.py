import sys
import numpy as np

import run_enterprise.samplers as samplers
import run_enterprise.output_plots as output_plots
import run_enterprise.derived_parameters as derived_parameters



def run_modelfit(args, psrmodel, parfile):
    psr = psrmodel.psr
    pta = psrmodel.get_pta()
    if args.outdir == None:
        outdir = 'chains/' + psr.name + "/"
    else:
        outdir = args.outdir

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
            sampler_info = sampler.run_solve(args,psrmodel,outdir)
            break

    if sampler_info is None:
        print("No sampler seemed to be enabled")
        sys.exit(1)

    ### Get the posterior 'chain' (or samples)
    samples, log_prob, weights, xtra = sampler.get_posteriors(args,sampler_info)
    raw_results=(samples, log_prob, weights)
    scaled_samples = np.copy(samples)

    parnames = []
    iparam = 0
    ### convert parameters
    for param in pta.params:
        # Test if the param has the to_par method
        if hasattr(param, 'to_par') and callable(param.to_par):
            name, scaled_samples[:,iparam] = param.to_par(param.name, samples[:,iparam])
        else:
            name = param.name
        parnames.append(name)
        iparam+=1
    parnames=np.array(parnames)

    par = psrmodel.par
    ### Compute any derived parameters
    derived, derived_names = derived_parameters.compute_derived_parameters(args,psr,par,parfile,scaled_samples,log_prob,weights,parnames)

    ### Write results
    results = output_plots.write_results(args,psr,par,parfile,scaled_samples,log_prob,weights,parnames,derived,derived_names)

    ### Make plots
    output_plots.make_plots(args,psr,scaled_samples,log_prob,weights, parnames,derived,derived_names,raw_results=sampler_info)

    return results, raw_results
