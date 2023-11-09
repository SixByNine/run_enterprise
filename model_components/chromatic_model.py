from . import xparameter as parameter
import numpy as np
from enterprise.signals import gp_signals, gp_bases, utils, signal_base
import enterprise.constants as const
from enterprise.signals.selections import Selection
from enterprise.signals import selections
from enterprise.signals import deterministic_signals

name="Chromatic (non-dm) Variations"
argdec="Parameters for fitting non-dm chromatic variations"

def setup_argparse(parser):
    parser.add_argument('--chrom', action='store_true', help='Enable Chrom variation search')
    parser.add_argument('--Achrom-max', type=float, default=-12, help='Max log10A_Chrom')
    parser.add_argument('--Achrom-min', type=float, default=-18, help='Min log10A_Chrom')
    parser.add_argument('--chrom-index', type=float, default=4, help='Chromatic index (default=4)')
    parser.add_argument('--chrom-gamma-max', type=float, default=5, help='Max gamma red')
    parser.add_argument('--chrom-gamma-min', type=float, default=0, help='Min gamma red')
    parser.add_argument('--chrom-ncoeff', type=int, default=None, help='Number of Chrom bins to use')
    parser.add_argument('--chrom-prior-log', action='store_true',
                        help='Use uniform prior in log space for chromatic noise amplitude')
    parser.add_argument('--chrom-tspan-mult', type=float, default=1, help='Multiplier for tspan for chromatic noise')
    parser.add_argument('--chrom-yr',action='store_true',help='Add annual chromatic sinusoid')
    parser.add_argument('--chrom-yr-amp-min',type=float,default=-20,help='min amp for annual chromatic sinusoid')
    parser.add_argument('--chrom-yr-amp-max',type=float,default=-5,help='max amp for annual chromatic sinusoid')
    parser.add_argument('--chrom-yr-idx-min',type=float,default=0,help='min index for annual chromatic sinusoid')
    parser.add_argument('--chrom-yr-idx-max',type=float,default=14,help='max index for annual chromatic sinuosid')
    parser.add_argument('--chrom-bump-start',default=[], nargs='+', type=float, help='Start epochs for chromatic bumps')
    parser.add_argument('--chrom-bump-end',default=[], nargs='+', type=float, help='End epochs for chromatic bumps')
    parser.add_argument('--chrom-bump-index-min',default=[], nargs='+', type=float, help='Indexes for chromatic bumps')
    parser.add_argument('--chrom-bump-index-max',default=[], nargs='+', type=float, help='Indexes for chromatic bumps')
    parser.add_argument('--chrom-bump-sigma-min',default=[], nargs='+', type=float, help='Sigma-min for chromatic bumps')
    parser.add_argument('--chrom-bump-sigma-max',default=[], nargs='+', type=float, help='Sigma-max for chromatic bumps')
    parser.add_argument('--chrom-bump-amp-min',default=[], nargs='+', type=float, help='Amp-min for chromatic bumps')
    parser.add_argument('--chrom-bump-amp-max',default=[], nargs='+', type=float, help='Amp-max for chromatic bumps')




def setup_model(args, psr, parfile):
    components = []
    if args.chrom:
        if args.chrom_ncoeff is None:
            nC = args.red_ncoeff
        else:
            nC = args.chrom_ncoeff

        Tspan = psr.toas.max() - psr.toas.min()
        Tspan *= args.chrom_tspan_mult
        nC = int(nC*args.chrom_tspan_mult)

        parfile.append("TNChromC {}\n".format(nC))
        parfile.append("TNChromIdx {}\n".format(args.chrom_index))

        A_min = args.Achrom_min
        A_max = args.Achrom_max

        if args.chrom_prior_log:
            log10_Achrom = parameter.Uniform(A_min, A_max, to_par=to_par)('Chrom_A')
        else:
            log10_Achrom = parameter.LinearExp(A_min, A_max, to_par=to_par)('Chrom_A')

        gamma_chrom = parameter.Uniform(args.chrom_gamma_min,args.chrom_gamma_max, to_par=to_par)('Chrom_gamma')
        plchrom = utils.powerlaw(log10_A=log10_Achrom, gamma=gamma_chrom)
        chrom = FourierBasisGP_Chrom(spectrum=plchrom, components=nC, Tspan=Tspan, index=args.chrom_index)
        components.append(chrom)

    if args.chrom_yr:
        # Yearly chromatic.
        log10_Amp_chrom1yr = parameter.Uniform(args.chrom_yr_amp_min, args.chrom_yr_amp_max,to_par=to_par)("CHROM_YR_AMP")
        phase_chrom1yr = parameter.Uniform(0, 2 * np.pi,to_par=to_par)("CHROM_YR_PHS")
        idx_chrom1yr = parameter.Uniform(args.chrom_yr_idx_min, args.chrom_yr_idx_max,to_par=to_par)("CHROM_YR_IDX")
        wf = chrom_yearly_sinusoid(log10_Amp=log10_Amp_chrom1yr, phase=phase_chrom1yr, idx=idx_chrom1yr)
        chrom1yr = deterministic_signals.Deterministic(wf, name="chrom1yr")
        components.append(chrom1yr)

    for ibump,(start,end) in enumerate(zip(args.chrom_bump_start,args.chrom_bump_end)):
        sign_param=parameter.Uniform(-1,1,to_par=to_par)(f"CHROM_BUMP_SIGN_{ibump+1}")
        t0_bump = parameter.Uniform(start,end,to_par=to_par)(f"CHROM_BUMP_EPOCH_{ibump+1}")
        while len(args.chrom_bump_sigma_min) <= ibump:
            args.chrom_bump_sigma_min.append(7)
        while len(args.chrom_bump_sigma_max) <= ibump:
            Tspan = psr.toas.max() - psr.toas.min()
            args.chrom_bump_sigma_max.append(Tspan/86400.0)
        while len(args.chrom_bump_amp_min) <= ibump:
            args.chrom_bump_amp_min.append(-10)
        while len(args.chrom_bump_amp_max) <= ibump:
            args.chrom_bump_amp_max.append(-1)
        while len(args.chrom_bump_index_min) <= ibump:
            args.chrom_bump_index_min.append(0)
        while len(args.chrom_bump_index_max) <= ibump:
            args.chrom_bump_index_max.append(14)

        sig_bump = parameter.Uniform(args.chrom_bump_sigma_min[ibump],
                                    args.chrom_bump_sigma_max[ibump],to_par=to_par)(f"CHROM_BUMP_SIG_{ibump+1}")
        amp_bump = parameter.Uniform(args.chrom_bump_amp_min[ibump],
                                     args.chrom_bump_amp_max[ibump],to_par=to_par)(f"CHROM_BUMP_AMP_{ibump+1}")
        idx_bump = parameter.Uniform(args.chrom_bump_index_min[ibump],
                                     args.chrom_bump_index_max[ibump],to_par=to_par)(f"CHROM_BUMP_IDX_{ibump+1}")
        wf = chrom_gaussian_bump(log10_Amp=amp_bump, t0=t0_bump,sigma=sig_bump,idx=idx_bump,sign_param=sign_param)
        components.append(deterministic_signals.Deterministic(wf))


    if len(components) > 0:
        model = components[0]
        for m in components[1:]:
            model += m
        return model
    else:
        return None

def FourierBasisGP_Chrom(spectrum, components=20,
                      selection=Selection(selections.no_selection),
                      Tspan=None, name='', index=4):
    """Convenience function to return a BasisGP class with a
    fourier basis."""

    basis = gp_bases.createfourierdesignmatrix_chromatic(nmodes=components, Tspan=Tspan,idx=index)
    BaseClass = gp_signals.BasisGP(spectrum, basis, selection=selection, name=name)

    class FourierBasisGP_DM(BaseClass):
        signal_type = 'basis'
        signal_name = 'chromatic noise'
        signal_id = 'chrom_noise_' + name if name else 'chrom_noise'

    return FourierBasisGP_DM




@signal_base.function
def chrom_gaussian_bump(toas, freqs, log10_Amp=-2.5, sign_param=1.0,
                    t0=53890, sigma=81, idx=2):
    """
    Chromatic time-domain Gaussian delay term in TOAs.
    Example: J1603-7202 in Lentati et al, MNRAS 458, 2016.
    """
    t0_sec=t0*86400.0
    sigma_sec=sigma*86400.0
    wf = 10**log10_Amp * np.exp(-(toas - t0_sec)**2/2/sigma_sec**2)
    return np.sign(sign_param) * wf * (1400 / freqs) ** idx



@signal_base.function
def chrom_yearly_sinusoid(toas, freqs, log10_Amp, phase, idx):
    """
    Chromatic annual sinusoid.
    :param log10_Amp: amplitude of sinusoid
    :param phase: initial phase of sinusoid
    :param idx: index of chromatic dependence
    :return wf: delay time-series [s]
    """

    wf = 10**log10_Amp * np.sin(2 * np.pi * const.fyr * toas + phase)
    return wf * (1400 / freqs) ** idx



def to_par(self,p,chain):
    if "Chrom_A" in p:
        return "TNChromAmp", chain
    elif "Chrom_gamma" in p:
        return "TNChromGam", chain
    elif "CHROM_YR" in p:
        return p,chain
    elif "CHROM_BUMP" in p:
        e=p.split("_")
        i=int(e[-1])
        if "EPOCH_" in p:
            p=f"GAUSEP_{i}"
        if "IDX_" in p:
            p=f"GAUSINDEX_{i}"
        if "SIG_" in p:
            p=f"GAUSSIG_{i}"
        return p, chain
    else:
        return None

