import dynesty
import numpy as np
from multiprocessing import Pool
import pickle, os
from dynesty import plotting
import matplotlib.pyplot as plt

name = "DynestySolver"
argdec = "Configure for nested sampling with DyNesty"


def setup_argparse(parser):
    parser.add_argument('--dynesty', action='store_true', help='Use dynesty sampler')
    parser.add_argument('--dynesty-dynamic', action='store_true', help='Use dynamic nested sampling')
    parser.add_argument('--dynesty-plots', action='store_true', help="make dynesty run plots")
    parser.add_argument('--dynesty-bound-eff', type=float, default=10.0, help="Efficiency to start bounding")
    parser.add_argument('--dynesty-bound',default="multi", help="Bounding method")
    parser.add_argument('--dynesty-sampler',default="auto", help="Sampling method")
    parser.add_argument('--dynesty-bootstrap',default=0,type=int,help="Bootstrap amount")


def activated(args):
    return args.dynesty

def run_solve(args,psrmodel,outdir):
    os.makedirs(outdir, exist_ok=True)
    filename = os.path.join(outdir, "dynesty_results.pkl")
    pta = psrmodel.get_pta()

    nlive = args.nlive
    print("Nlive: {} ndim: {}".format(nlive, len(pta.params)))
    if nlive < 2 * len(pta.params):
        nlive = 2 * len(pta.params)
        print("Warning... nlive too small, setting to {}".format(nlive))

    args.burn = 0  ## force no burn in for nested sampling
    if args.sample:
        psrmodel.pack_for_thread()
        with Pool(args.nthread, initializer=init, initargs=[psrmodel]) as tpool:
            print("Run Dynesty")
            if args.dynesty_dynamic:
                sampler = dynesty.DynamicNestedSampler(get_lnlikelihood, prior_transform, ndim=len(pta.params), pool=tpool,
                                                       queue_size=args.nthread, sample=args.dynesty_sampler, bound=args.dynesty_bound,
                                                       first_update={'min_ncall': 1000, 'min_eff': args.dynesty_bound_eff},
                                                       bootstrap=args.dynesty_bootstrap)
                sampler.run_nested(nlive_init=args.nlive,nlive_batch=args.nlive//2,dlogz_init=0.1)
                if args.nsample is not None:
                    while (sampler.results.samples.shape[0] < args.nsample):
                        print(f"Not enough samples... continue ({sampler.results.samples.shape[0]}/{args.nsample})")
                        sampler.add_batch(nlive=args.nlive)

            else:
                sampler = dynesty.NestedSampler(get_lnlikelihood, prior_transform, ndim=len(pta.params), nlive=nlive,
                                                pool=tpool, queue_size=args.nthread, sample=args.dynesty_sampler, bound=args.dynesty_bound,
                                                first_update={'min_ncall': 1000, 'min_eff': args.dynesty_bound_eff},bootstrap=args.dynesty_bootstrap)
                sampler.run_nested()

            print(sampler.results.summary())
            with open(filename, "wb") as f:
                pickle.dump(sampler.results, f)
            # print(sampler.results)
            return sampler.results
    else:
        with open(filename, 'rb') as f:
            return pickle.load(f)


def get_posteriors(args, results):
    try:
        weights = np.exp(results['logwt'] - results['logz'][-1])
    except:
        weights = results['weights']

    return results.samples, results.logz, weights, results

def init(psrmodel):
    global invTs
    global pta
    print("Init Thread...")
    invTs = []
    pta = psrmodel.get_pta()

    for p in pta.params:
        invTs.append(p.invT)
def prior_transform(u):
    global invTs
    """
    Surely there is a faster way... but I'm not sure how (mkeith jan 2021).
    """
    res = np.zeros_like(u)
    for i in range(len(u)):
        res[i] = invTs[i](u[i])
    return res


def get_lnlikelihood(p):
    global pta
    return pta.get_lnlikelihood(p)
