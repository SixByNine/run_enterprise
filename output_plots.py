import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import corner
import dynesty
from dynesty import plotting
from astropy import coordinates as coord
from astropy import units as units

def setup_argparse(parser):
    parser.add_argument('--white-corner', action='store_true', help='Make the efac/equad corner plots')
    parser.add_argument('--all-corner', '--corner-all', action='store_true', help='Make corner plots with all params')
    parser.add_argument('--plot-chain', action='store_true', help='Make a plot of the chains/posterior samples')
    parser.add_argument('--plot-derived', action='store_true', help='Include derived parameters in corner plots etc.')
    parser.add_argument('--burn', default=0.25, help='Fraction of chain to burn-in (MC only; default=0.25)')
    parser.add_argument('--truth-file', default=None, help='Truths values of parameters; default=None') # mod ly

def cornerplot(samples, labels, weights=None, truths=None, truth_color='cyan'): # mod ly

    if weights is None:
        weights = np.ones_like(samples[:,0])
    res={'samples': samples, 'weights':weights}
    fig,ax= plotting.cornerplot(res,truths=truths,labels=labels,truth_color=truth_color) # mod ly
    for aa in ax:
        for a in aa:
            if len(a.xaxis.get_label().get_text()) > 0:
                a.ticklabel_format(axis='x', style='sci', scilimits=(-3, 3), useMathText=True)
            if len(a.yaxis.get_label().get_text()) > 0:
                a.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useMathText=True)
    return fig


def make_plots(args, psr, samples, loglike, weights, pars, derived,derived_names,raw_results=None):
    print("Making plots....")
    if args.dynesty and args.dynesty_plots:
        try:
            fig,ax = plotting.runplot(raw_results)
            fig.savefig("%s.runplot.pdf" % psr.name)
        except Exception as e:
            print(e)

        try:
            fig,ax = plotting.traceplot(raw_results, labels=pars)
            fig.savefig("%s.traceplot.pdf" % psr.name)
        except Exception as e:
            print(e)

    if args.plot_derived and len(derived_names)>0:
        derived_names_x = []
        for n in derived_names:
            derived_names_x.append("*{}".format(n))
        pars = np.concatenate((pars,derived_names_x))
        samples = np.concatenate((samples,derived),axis=1)


    burn = int(len(samples) * args.burn)
    i_ef = np.array([i for i, v in enumerate(pars) if 'TNEF' in v], dtype=np.int)
    i_eq = np.array([i for i, v in enumerate(pars) if 'TNEQ' in v], dtype=np.int)

    m_ef = np.zeros(len(pars), dtype=np.bool)
    m_ef[i_ef] = True
    m_eq = np.zeros(len(pars), dtype=np.bool)
    m_eq[i_eq] = True

    n_w = np.logical_not(np.logical_or(m_ef, m_eq))
    imax = np.argmax(loglike)
    pmax = samples[imax, :]

    plt.rcParams.update({'font.size': 8})

# mod ly
    #print(pars[n_w])
    if args.truth_file is None:
        fig = cornerplot(samples[burn:, n_w], labels=pars[n_w], weights=weights, truths=pmax[n_w])
    else:
        d = {}
        with open(args.truth_file) as f:
            for line in f:
                key, val = line.rsplit(maxsplit=1)
                d[key] = float(val)
        fiducial= [None]*len(pars[n_w])
        for i in range(len(pars[n_w])):
            if pars[i].strip("*") in d: # strip * from derived parameters
                fiducial[i] = d[pars[i].strip("*")]
        fig = cornerplot(samples[burn:, n_w], labels=pars[n_w], weights=weights, truths=fiducial, truth_color='r')
# mod ly

    fig.savefig("%s.corner.pdf" % psr.name)

    # Make optional corner plots
    if args.all_corner:
        fig = cornerplot(samples[burn:, :], labels=pars, weights=weights,truths=pmax);
        fig.savefig("%s.corner_all.pdf" % psr.name)

    if args.white_corner:
        fig = cornerplot(samples[burn:, m_ef], labels=pars[m_ef], weights=weights,truths=pmax[m_ef]);
        fig.savefig("%s.corner_efac.pdf" % psr.name)
        fig = cornerplot(samples[burn:, m_eq], labels=pars[m_eq], weights=weights,truths=pmax[m_eq]);
        fig.savefig("%s.corner_equad.pdf" % psr.name)

    # Plot posterior if requested.
    if args.plot_chain:
        x = np.arange(samples.shape[0])
        with PdfPages("%s.chain.pdf" % psr.name) as pdf:
            for i in range(len(pars)):
                fig = plt.figure(figsize=(16, 8))
                plt.plot(x[:burn], samples[:burn, i], '.', color='gray')
                plt.plot(x[burn:], samples[burn:, i], '.', color='k')
                plt.title("%s" % pars[i])
                plt.ylabel("%s" % pars[i])
                pdf.savefig()
                plt.close()

    print("Plots saved")


def write_results(args, psr, par, parfile, samples, loglike, weights, pars, derived,derived_names):
    burn = int(len(samples) * args.burn)
    outpar = par + ".post"
    resfname = par + ".results"
    with open(outpar, "w") as outf:
        with open(resfname, "w") as resf:
            s = ("{:20s} {:10s} {:10s} {:10s} {:10s} {:10s} {:10s} {:10s} {:10s}".format("#param", " max-like", " mean",
                                                                                         " std", "2.5%", " 15.9%", " 50%",
                                                                                         " 84.1%", " 97.5%"))
            resf.write(s + "\n")
            print(s)



            imax = np.argmax(loglike)
            pmax = samples[imax, :len(pars)]

            ul68,ll68,ul,ll,means,stds,medians=get_stats(samples[burn:, :len(pars)], weights)



            for p, v, u, l, mean, median, sig, u68, l68 in zip(pars, pmax, ul, ll, means, medians, stds, ul68, ll68):
                s = "{:20s} {:< 10g} {:< 10g} {:< 10g} {:< 10g} {:< 10g}  {:< 10g} {:< 10g} {:< 10g} ".format(p, v, mean,
                                                                                                              sig,
                                                                                                              l, l68,
                                                                                                              median,
                                                                                                              u68, u)
                resf.write(s + "\n")
                print(s)


            for line in parfile:
                e = line.split()
                if e[0] in pars \
                        or (e[0]=="RAJ" and "dRAJ" in pars) \
                        or (e[0]=="DECJ" and "dDECJ" in pars):
                    outf.write("#"+line)
                else:
                    outf.write(line)
            for p,v in zip(pars,pmax):
                ## Still need a couple of special cases...
                if p=="dRAJ":
                    p="RAJ"
                    new_ra = psr.coord.ra + v*units.arcsecond
                    v=new_ra.to_string(unit=units.hourangle,sep=":",precision=8)
                if p == "dDECJ":
                    p = "DECJ"
                    new_dec = psr.coord.dec + v*units.arcsecond
                    v = new_dec.to_string(unit=units.degree, sep=":",precision=8)

                outf.write("%s %s\n" % (p, v))


            if len(derived) > 0:
                print("")
                print("Derived Params")
                print("")

                dimax = np.argmax(loglike)
                dpmax = derived[dimax, :]
                dul68,dll68,dul,dll,dmeans,dstds,dmedians=get_stats(derived[burn:, :], weights)

                for p, v, u, l, mean, median, sig, u68, l68 in zip(derived_names, dpmax, dul, dll, dmeans, dmedians, dstds, dul68,
                                                                   dll68):
                    s = "{:20s} {:< 10g} {:< 10g} {:< 10g} {:< 10g} {:< 10g}  {:< 10g} {:< 10g} {:< 10g} ".format(p, v,
                                                                                                                  mean,
                                                                                                                  sig,
                                                                                                                  l,
                                                                                                                  l68,
                                                                                                                  median,
                                                                                                                  u68,
                                                                                                                  u)
                    resf.write(s + "\n")
                    print(s)
                for p,v in zip(derived_names,pmax):
                    outf.write("%s %s\n" % (p, v))

    return pars, pmax, ul, ll, means, medians, stds, ul68, ll68

def get_stats(samps,weights):
    if weights is None:
        ul68 = np.percentile(samps, 84.1, axis=0)
        ll68 = np.percentile(samps, 15.9, axis=0)
        ul = np.percentile(samps, 97.5, axis=0)
        ll = np.percentile(samps, 2.5, axis=0)
        means = np.mean(samps, axis=0)
        stds = np.std(samps, axis=0)
        medians = np.median(samps, axis=0)
    else:
        ## use dynesty to resample
        ul68 = np.apply_along_axis(dynesty.utils.quantile, 0, samps, 0.841, weights)[0]
        ll68 = np.apply_along_axis(dynesty.utils.quantile, 0, samps, 0.159, weights)[0]
        ul = np.apply_along_axis(dynesty.utils.quantile, 0, samps, 0.975, weights)[0]
        ll = np.apply_along_axis(dynesty.utils.quantile, 0, samps, 0.025, weights)[0]
        medians = np.apply_along_axis(dynesty.utils.quantile, 0, samps, 0.500, weights)[0]
        means, covs = dynesty.utils.mean_and_cov(samps, weights=weights)
        stds = np.sqrt(np.diag(covs))

    return ul68,ll68,ul,ll,means,stds,medians