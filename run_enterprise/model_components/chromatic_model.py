from . import xparameter as parameter
import numpy as np
from enterprise.signals import gp_signals, gp_bases, utils, signal_base
import enterprise.constants as const
from enterprise.signals.selections import Selection
from enterprise.signals import selections
from enterprise.signals import deterministic_signals
from enterprise.signals.parameter import function

name="Chromatic (non-dm) Variations"
argdec="Parameters for fitting non-dm chromatic variations"

def setup_argparse(parser):
    parser.add_argument('--chrom', action='store_true', help='Enable Chrom variation search')
    parser.add_argument('--Achrom-max', type=float, default=-12, help='Max log10A_Chrom')
    parser.add_argument('--Achrom-min', type=float, default=-18, help='Min log10A_Chrom')
    parser.add_argument('--chrom-index', type=float, default=4, help='Chromatic index (default=4)')
    parser.add_argument('--chrom-index-max', type=float, default=None, help='Chromatic index max (fit for it!)')
    parser.add_argument('--chrom-index-min', type=float, default=None, help='Chromatic index min')
    parser.add_argument('--chrom-gamma-max', type=float, default=5, help='Max gamma red')
    parser.add_argument('--chrom-gamma-min', type=float, default=0, help='Min gamma red')
    parser.add_argument('--chrom-ncoeff', type=int, default=None, help='Number of Chrom bins to use')
    parser.add_argument('--chrom-prior-log', action='store_true',
                        help='Use uniform prior in log space for chromatic noise amplitude')
    parser.add_argument('--chrom-tspan-mult', type=float, default=1, help='Multiplier for tspan for chromatic noise')
    parser.add_argument('--chrom-nlog-below-tspan', type=int, default=0, help='Number of log-spaced chrom components below tspan')
    parser.add_argument('--chrom-nlog-factor', type=float, default=2.0, help='Factor for log-spaced chrom components below tspan') 
    parser.add_argument('--chrom-fit-cm0','--cm0', type=float, default=None, help='Fit for mean chromatic noise (CM)') 
    parser.add_argument('--chrom-fit-cm1','--cm1', type=float, default=None, help='Fit for linear chromatic noise (CM1)') 
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

    parser.add_argument('--chrom-exp-epoch', default=[], nargs='+', type=float, help='Epoch for chromatic exponential events')
    parser.add_argument('--chrom-exp-epoch-range', default=[], nargs='+', type=float, help='Epoch range for chromatic exponential events')
    parser.add_argument('--chrom-exp-amp-min', default=[], nargs='+', type=float, help='Amp-min for chromatic exponential events')
    parser.add_argument('--chrom-exp-amp-max', default=[], nargs='+', type=float, help='Amp-max for chromatic exponential events')
    parser.add_argument('--chrom-exp-tau-min', default=[], nargs='+', type=float, help='Tau-min for chromatic exponential events')
    parser.add_argument('--chrom-exp-tau-max', default=[], nargs='+', type=float, help='Tau-max for chromatic exponential events')
    parser.add_argument('--chrom-exp-index-min', default=[], nargs='+', type=float, help='Index-min for chromatic exponential events')
    parser.add_argument('--chrom-exp-index-max', default=[], nargs='+', type=float, help='Index-max for chromatic exponential event')
    



def setup_model(args, psr, parfile):
    components = []
    if args.chrom or args.chrom_fit_cm0 or args.chrom_fit_cm1:
        if args.chrom_index_max is not None:
            if args.chrom_index_min is None:
                args.chrom_index_min = 0
            chrom_index = parameter.Uniform(args.chrom_index_min,args.chrom_index_max,to_par=to_par)("Chrom_index")
        else:
            parfile.append("TNChromIdx {}\n".format(args.chrom_index))
            chrom_index=parameter.Constant(args.chrom_index)
        if args.chrom_fit_cm0 or args.chrom_fit_cm1:
            DMEPOCH=55000
            orig_cm1=0
            orig_cm0=0
            for line in parfile:
                e = line.strip().split()
                if len(e) > 1:
                    if e[0] == "CM":
                        orig_cm0 = float(e[1])
                    if e[0] == "CM1":
                        orig_cm1 = float(e[1])
                    elif e[0] == "DMEPOCH":
                        DMEPOCH = float(e[1])
            cm1 = parameter.Constant(0)
            cm0 = parameter.Constant(0)
            if args.chrom_fit_cm0:
                cm0 = parameter.Uniform(-args.chrom_fit_cm0+orig_cm0,args.chrom_fit_cm0+orig_cm0,to_par=to_par)("CM")
            if args.chrom_fit_cm1:
                cm1 = parameter.Uniform(-args.chrom_fit_cm1+orig_cm1,args.chrom_fit_cm1+orig_cm1,to_par=to_par)("CM1")
            
            components.append(deterministic_signals.Deterministic(fit_chrom_poly(cm=cm0,cm1=cm1,idx=chrom_index,dmepoch=parameter.Constant(DMEPOCH)),name="chrom_poly"))

        
        if args.chrom:
            if args.chrom_ncoeff is None:
                nC = args.red_ncoeff
            else:
                nC = args.chrom_ncoeff

            Tspan = psr.toas.max() - psr.toas.min()
            Tspan *= args.chrom_tspan_mult
            nC = int(nC*args.chrom_tspan_mult)

            parfile.append("TNChromC {}\n".format(nC))
            

            A_min = args.Achrom_min
            A_max = args.Achrom_max

            if args.chrom_prior_log:
                log10_Achrom = parameter.Uniform(A_min, A_max, to_par=to_par)('Chrom_A')
            else:
                log10_Achrom = parameter.LinearExp(A_min, A_max, to_par=to_par)('Chrom_A')

            gamma_chrom = parameter.Uniform(args.chrom_gamma_min,args.chrom_gamma_max, to_par=to_par)('Chrom_gamma')
            plchrom = utils.powerlaw(log10_A=log10_Achrom, gamma=gamma_chrom)

            if args.chrom_nlog_below_tspan > 0:
                freq0 = 1.0 / Tspan
                low_modes = sorted(freq0 * args.chrom_nlog_factor**-np.arange(1, args.chrom_nlog_below_tspan + 1)) # log spaced below Tspan
                hi_modes = freq0 * np.arange(1,nC+1)
                modes = np.concatenate((low_modes,hi_modes))
                print(f"Modes: {len(low_modes)} {len(hi_modes)}")
                print(modes)
                chrom = FourierBasisGP_Chrom(spectrum=plchrom, components=None, modes=modes, Tspan=Tspan, index=chrom_index)

                parfile.append(f"TNChromFLog {args.chrom_nlog_below_tspan}\n")
                parfile.append(f"TNChromFLog_factor {args.chrom_nlog_factor}\n")
            else:
                chrom = FourierBasisGP_Chrom(spectrum=plchrom, components=nC, Tspan=Tspan, index=chrom_index)

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
        components.append(deterministic_signals.Deterministic(wf, name=f"chrom_bump_{ibump+1}"))

    if args.chrom_exp_epoch:
        for iexp,epoch in enumerate(args.chrom_exp_epoch):
            # if we go past the end of array, just use the last value for the rest of the events
            while len(args.chrom_exp_amp_min) <= iexp:
                args.chrom_exp_amp_min.append(args.chrom_exp_amp_min[-1])
            while len(args.chrom_exp_amp_max) <= iexp:
                args.chrom_exp_amp_max.append(args.chrom_exp_amp_max[-1])
            while len(args.chrom_exp_tau_min) <= iexp:
                args.chrom_exp_tau_min.append(args.chrom_exp_tau_min[-1])
            while len(args.chrom_exp_tau_max) <= iexp:
                args.chrom_exp_tau_max.append(args.chrom_exp_tau_max[-1]) 
            while len(args.chrom_exp_index_min) <= iexp:
                args.chrom_exp_index_min.append(args.chrom_exp_index_min[-1])
            while len(args.chrom_exp_index_max) <= iexp:
                args.chrom_exp_index_max.append(args.chrom_exp_index_max[-1])
            if args.chrom_exp_tau_min[iexp] <= 0:
                raise ValueError(f"chrom_exp_tau_min must be > 0 for event {iexp+1}")
            if args.chrom_exp_tau_max[iexp] <= 0:
                raise ValueError(f"chrom_exp_tau_max must be > 0 for event {iexp+1}")
            if args.chrom_exp_amp_min[iexp] >= args.chrom_exp_amp_max[iexp]:
                raise ValueError(f"chrom_exp_amp_min must be < chrom_exp_amp_max for event {iexp+1}")
            if args.chrom_exp_index_min[iexp] >= args.chrom_exp_index_max[iexp]:
                raise ValueError(f"chrom_exp_index_min must be < chrom_exp_index_max for event {iexp+1}")
            if args.chrom_exp_tau_min[iexp] >= args.chrom_exp_tau_max[iexp]:
                raise ValueError(f"chrom_exp_tau_min must be < chrom_exp_tau_max for event {iexp+1}")
            if args.chrom_exp_amp_min[iexp] >= args.chrom_exp_amp_max[iexp]:
                raise ValueError(f"chrom_exp_amp_min must be < chrom_exp_amp_max for event {iexp+1}")
            if args.chrom_exp_index_min[iexp] >= args.chrom_exp_index_max[iexp]:
                raise ValueError(f"chrom_exp_index_min must be < chrom_exp_index_max for event {iexp+1}")
            if args.chrom_exp_tau_min[iexp] <= 0:
                raise ValueError(f"chrom_exp_tau_min must be > 0 for event {iexp+1}")
            if args.chrom_exp_tau_max[iexp] <= 0:
                raise ValueError(f"chrom_exp_tau_max must be > 0 for event {iexp+1}")
            if args.chrom_exp_amp_min[iexp] >= args.chrom_exp_amp_max[iexp]:
                raise ValueError(f"chrom_exp_amp_min must be < chrom_exp_amp_max for event {iexp+1}")
            if args.chrom_exp_index_min[iexp] >= args.chrom_exp_index_max[iexp]:
                raise ValueError(f"chrom_exp_index_min must be < chrom_exp_index_max for event {iexp+1}")
            if args.chrom_exp_tau_min[iexp] <= 0:
                raise ValueError(f"chrom_exp_tau_min must be > 0 for event {iexp+1}")
            if args.chrom_exp_epoch_range:
                epoch = parameter.Uniform(epoch-args.chrom_exp_epoch_range[iexp],epoch+args.chrom_exp_epoch_range[iexp],to_par=to_par)(f"EXPEP_{iexp+1}")
            else:
                epoch = parameter.Constant(epoch)
            amp_exp = parameter.Uniform(args.chrom_exp_amp_min[iexp],
                                        args.chrom_exp_amp_max[iexp],to_par=to_par)(f"EXPPH_{iexp+1}")
            tau_exp = parameter.Uniform(args.chrom_exp_tau_min[iexp],
                                         args.chrom_exp_tau_max[iexp],to_par=to_par)(f"EXPTAU_{iexp+1}")
            idx_exp = parameter.Uniform(args.chrom_exp_index_min[iexp],
                                         args.chrom_exp_index_max[iexp],to_par=to_par)(f"EXPIDX_{iexp+1}")
            wf = chrom_exponential(epoch=epoch, tau=tau_exp, amp=amp_exp, index=idx_exp)
            components.append(deterministic_signals.Deterministic(wf, name=f"chrom_exponential_{iexp+1}"))

    if len(components) > 0:
        model = components[0]
        for m in components[1:]:
            model += m
        return model
    else:
        return None

def FourierBasisGP_Chrom(spectrum, components=20,
                      selection=Selection(selections.no_selection),
                      Tspan=None, name='', modes=None, index=4):
    """Convenience function to return a BasisGP class with a
    fourier basis."""

    basis = createfourierdesignmatrix_chromatic(nmodes=components, Tspan=Tspan,idx=index,modes=modes)
    BaseClass = gp_signals.BasisGP(spectrum, basis, selection=selection, name=name)

    class FourierBasisGP_DM(BaseClass):
        signal_type = 'basis'
        signal_name = 'chromatic noise'
        signal_id = 'chrom_noise_' + name if name else 'chrom_noise'

    return FourierBasisGP_DM


# Note that this is copied from enterprise.signals.gp_bases because that lacks the ability to set modes. When enterprise is uipdated 
# this should be removed and the enterprise version used.
@function
def createfourierdesignmatrix_chromatic(toas, freqs, nmodes=30, Tspan=None, logf=False, fmin=None, fmax=None, idx=4,modes=None):

    """
    Construct Scattering-variation fourier design matrix.

    :param toas: vector of time series in seconds
    :param freqs: radio frequencies of observations [MHz]
    :param nmodes: number of fourier coefficients to use
    :param freq: option to output frequencies
    :param Tspan: option to some other Tspan
    :param logf: use log frequency spacing
    :param fmin: lower sampling frequency
    :param fmax: upper sampling frequency
    :param idx: Index of chromatic effects

    :return: F: Chromatic-variation fourier design matrix
    :return: f: Sampling frequencies
    """

    # get base fourier design matrix and frequencies
    F, Ffreqs = gp_bases.createfourierdesignmatrix_red(toas, nmodes=nmodes, Tspan=Tspan, logf=logf, fmin=fmin, fmax=fmax,modes=modes)

    # compute the DM-variation vectors
    Dm = (1400 / freqs) ** idx

    return F * Dm[:, None], Ffreqs


@signal_base.function
def chrom_exponential(toas, freqs, epoch, index, tau, amp):
    t = toas - epoch * 86400.0
    m = t >= 0
    res = np.zeros_like(toas)
    res[m] = amp * np.exp(-t[m] / (tau*86400.0)) * (1400 / freqs[m]) ** -index
    return res



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

@signal_base.function
def fit_chrom_poly(toas, freqs, dmepoch, cm, cm1,idx):
    x_yr = (toas/86400.0 - dmepoch)/365.25
    delta_cm = cm + x_yr*cm1
    t = delta_cm * (1400 / freqs) ** idx
    return t

def to_par(self,p,chain):
    if "Chrom_A" in p:
        return "TNChromAmp", chain
    elif "Chrom_gamma" in p:
        return "TNChromGam", chain
    elif "Chrom_index" in p:
        return "TNChromIdx", chain
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
    elif p=="CM":
        return "CM", chain
    elif p=="CM1":
        return "CM1", chain
    if "EXPIDX" in p:
        return p, -chain
    elif "EXP" in p:
        return p, chain
    else:
        return None

