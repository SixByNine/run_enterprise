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





def setup_model(args, psr, parfile):
    dm=None
    if args.dm:
        if args.dm_ncoeff is None:
            nC = args.red_ncoeff
        else:
            nC = args.dm_ncoeff

        Tspan = psr.toas.max() - psr.toas.min()

        parfile.append("TNDMC {}\n".format(nC))

        log10_Adm = parameter.LinearExp(args.Adm_min,args.Adm_max,to_par=to_par)('DM_A')
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

    basis = utils.createfourierdesignmatrix_dm(nmodes=components, Tspan=Tspan)
    BaseClass = gp_signals.BasisGP(spectrum, basis, selection=selection, name=name)

    class FourierBasisGP_DM(BaseClass):
        signal_type = 'basis'
        signal_name = 'dm noise'
        signal_id = 'dm_noise_' + name if name else 'dm_noise'

    return FourierBasisGP_DM

def to_par(self,p,chain):
   if "DM_A" in p:
      return "TNDMAmp", chain
   if "DM_gamma" in p:
      return "TNDMGam", chain
   else:
      return None
