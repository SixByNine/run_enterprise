from . import xparameter as parameter
from enterprise.signals import gp_signals, utils
from enterprise.signals import deterministic_signals

from enterprise_extensions.chromatic import solar_wind as SW

name="SWGP version of Time Variable Solar wind"
argdec="Parameters for fitting dm variations due to solar wind"

def setup_argparse(parser):
    """
    :param parser:
    :return:
    """

    parser.add_argument('--swgp', action='store_true', help='Turn on SWGP')
    parser.add_argument('--swgp-const', action='store_true', help='Turn on constant SW from SWGP')
    parser.add_argument('--swgp-log10-A-dm-sw-min', type=float, default=-15,
                        help='Minimum log10_A_dm_sw for SWGP power-law prior')
    parser.add_argument('--swgp-log10-A-dm-sw-max', type=float, default=-8,
                        help='Maximum log10_A_dm_sw for SWGP power-law prior')
    parser.add_argument('--swgp-gamma-dm-sw-min', type=float, default=0,
                        help='Minimum gamma_dm_sw for SWGP power-law prior')
    parser.add_argument('--swgp-gamma-dm-sw-max', type=float, default=2,
                        help='Maximum gamma_dm_sw for SWGP power-law prior')
    parser.add_argument('--swgp-nmodes', type=int, default=15,
                        help='Number of Fourier modes for SWGP')
    



def setup_model(args, psr, parfile):
    components=[]
    if args.swgp_const:
        n_earth = parameter.Uniform(0,20,to_par=to_par)('n_earth')
        tmin = psr.toas.min()
        tmax = psr.toas.max()
        sw = SW.solar_wind(n_earth=n_earth,t_init=tmin,t_final=tmax)
        mean_sw = deterministic_signals.Deterministic(sw, name='mean_sw')
        components.append(mean_sw)
    if args.swgp:
        log10_A_dm_sw = parameter.Uniform(args.swgp_log10_A_dm_sw_min,
                                          args.swgp_log10_A_dm_sw_max,
                                          to_par=to_par)('log10_A_dm_sw')
        gamma_dm_sw = parameter.Uniform(args.swgp_gamma_dm_sw_min,
                                        args.swgp_gamma_dm_sw_max,
                                        to_par=to_par)('gamma_dm_sw')
        Tspan = psr.toas.max() - psr.toas.min()
        dm_sw_basis = SW.createfourierdesignmatrix_solar_dm(nmodes=args.swgp_nmodes,Tspan=Tspan)
        dm_sw_prior = utils.powerlaw(log10_A=log10_A_dm_sw, gamma=gamma_dm_sw)
        gp_sw = gp_signals.BasisGP(priorFunction=dm_sw_prior, basisFunction=dm_sw_basis, name='gp_sw')
        components.append(gp_sw)
    if len(components) > 0:
        model = components[0]
        for m in components[1:]:
            model += m
        return model
    else:
        return None




def to_par(self, p, chain):
   if "n_earth" in p:
       return "NE_SW", chain
   if "log10_A_dm_sw" in p:
       return "log10_A_dm_sw", chain
   if "gamma_dm_sw" in p:
       return "gamma_dm_sw", chain
   else:
       return None