from . import xparameter as parameter
import numpy as np
from enterprise.signals import gp_signals, signal_base, utils
import enterprise.constants as const
from enterprise.signals.selections import Selection
from enterprise.signals import selections

name="DM Variations"
argdec="Parameters for fitting dm variations"

def setup_argparse(parser):
    parser.add_argument('-D', '--dm', action='store_true', help='Enable DM variation search')
    parser.add_argument('--Adm-max', type=float, default=-12, help='Max log10A_DM')
    parser.add_argument('--Adm-min', type=float, default=-18, help='Min log10A_DM')
    parser.add_argument('--dm-ncoeff', type=int, default=None, help='Number of DM bins to use')
    parser.add_argument('--dm-prior-log', action='store_true',
                        help='Use uniform prior in log space for dm noise amplitude')




def setup_model(args, psr, parfile):
    dm=None
    if args.dm:
        if args.dm_ncoeff is None:
            nC = args.red_ncoeff
        else:
            nC = args.dm_ncoeff

        Tspan = psr.toas.max() - psr.toas.min()

        parfile.append("TNDMC {}\n".format(nC))

        # Note that enterprise models DM as time delay at 1400 MHz for some reason...
        # So we need to scale our values for the prior, and later unscale them
        dm_scale_factor = np.log10(1400.0**-2 / 2.41e-4)
        A_min = args.Adm_min + dm_scale_factor
        A_max = args.Adm_max + dm_scale_factor

        if args.dm_prior_log:
            log10_Adm = parameter.Uniform(A_min, A_max, to_par=to_par)('DM_A')
        else:
            log10_Adm = parameter.LinearExp(A_min, A_max, to_par=to_par)('DM_A')
        gamma_dm = parameter.Uniform(0,5,to_par=to_par)('DM_gamma')
        pldm = powerlaw_nogw(log10_A=log10_Adm, gamma=gamma_dm)
        dm = FourierBasisGP_DM(spectrum=pldm, components=nC, Tspan=Tspan)
    return dm

@signal_base.function
def powerlaw_nogw(f, log10_A=-16, gamma=5):
    df = np.diff(np.concatenate((np.array([0]), f[::2])))
    return ((10 ** log10_A) ** 2 *
            const.fyr ** (gamma - 3) * f ** (-gamma) * np.repeat(df, 2))

def FourierBasisGP_DM(spectrum, components=20,
                      selection=Selection(selections.no_selection),
                      Tspan=None, name=''):
    """Convenience function to return a BasisGP class with a
    fourier basis."""

    basis = utils.createfourierdesignmatrix_dm(nmodes=components, Tspan=Tspan,fref=1400.0)
    BaseClass = gp_signals.BasisGP(spectrum, basis, selection=selection, name=name)

    class FourierBasisGP_DM(BaseClass):
        signal_type = 'basis'
        signal_name = 'dm noise'
        signal_id = 'dm_noise_' + name if name else 'dm_noise'

    return FourierBasisGP_DM

def to_par(self,p,chain):
   dm_scale_factor = np.log10(1400.0 ** -2 / 2.41e-4)
   if "DM_A" in p:
      return "TNDMAmp", chain - dm_scale_factor
   if "DM_gamma" in p:
      return "TNDMGam", chain
   else:
      return None
