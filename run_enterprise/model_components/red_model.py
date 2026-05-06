from html import parser

from . import xparameter as parameter
import numpy as np
from enterprise.signals import utils
from enterprise.signals import gp_signals
from enterprise.signals.parameter import function
import enterprise.constants as const
from enterprise.signals.selections import Selection
from enterprise.signals import selections


from enterprise.signals import deterministic_signals

import sys

name = "RedNoise"
argdec = "RedNoise Parameters. Note that the red noise model is default ENABLED"
Tspan=0

def setup_argparse(parser):
    parser.add_argument('--no-red-noise', dest='red', default=True, action='store_false',
                        help='Disable Power Law Red Noise search')
    parser.add_argument('--Ared-max', '-A', type=float, default=-12, help='Max log10A_Red')
    parser.add_argument('--Ared-min', type=float, default=-18, help='Min log10A_Red')
    parser.add_argument('--red-at', type=float, default=None, help='log10A measured at this timescale (in years)')
    parser.add_argument('--red-gamma-max', type=float, default=10, help='Max gamma red')
    parser.add_argument('--red-gamma-min', type=float, default=0, help='Min gamma red')
    parser.add_argument('--red-prior-log', action='store_true',
                        help='Use uniform prior in log space for red noise amplitude')
    parser.add_argument('--red-ncoeff', type=int, default=60, help='Number of red noise coefficients (nC)')
    parser.add_argument('--red-minf', type=float, default=None, help='Low frequency cut-off')
    parser.add_argument('--tspan-mult', type=float, default=1, help='Multiplier for tspan')
    parser.add_argument('--qp', action='store_true', help='Use QP nudot model')
    # parser.add_argument('--qp-f0-min', type=float, default=None)
    # parser.add_argument('--qp-f0-max', type=float, default=None)
    parser.add_argument('--qp-ratio-max', type=float, default=3.5)
    parser.add_argument('--qp-sigma-max', type=float, default=0.2)

    parser.add_argument('--qp-p-min-np', type=float, default=4.0)
    parser.add_argument('--qp-p-min', type=float, default=None)
    parser.add_argument('--qp-p-max', type=float, default=None)

    parser.add_argument('--plot-Ared-at-T', type=float, default=None, help="Plot derived Ared at freq 1/(T yr)")

    parser.add_argument('--double-powerlaw',action='store_true',help='Use double powerlaw model')
    parser.add_argument('--double-powerlaw-delta-gamma','--delta-gamma',type=float,default=None)
    parser.add_argument('--red-gamma2-max', type=float, default=3, help='Max gamma2 red for double pl')
    parser.add_argument('--red-gamma2-min', type=float, default=0, help='Min gamma2 red for double pl')

    parser.add_argument('--A2red-max', type=float, default=-12, help='Max log10A_Red for double pl')
    parser.add_argument('--A2red-min', type=float, default=-18, help='Min log10A_Red for double pl')

    parser.add_argument('--red-nlog-below-tspan', type=int, default=0, help='Number of log-spaced red noise components below tspan')
    parser.add_argument('--red-nlog-factor', type=float, default=2.0, help='Factor for log-spaced red noise components below tspan')
    parser.add_argument('--red-basis-subtract-quad', action='store_true', help='Subtract a quadratic from the red noise basis (equv to --red-basis-subtract-poly-order 2)')
    parser.add_argument('--red-basis-subtract-cubic', action='store_true', help='Subtract a cubic from the red noise basis (equv to --red-basis-subtract-poly-order 3)')
    parser.add_argument('--red-basis-subtract-poly-order', '--red-subpoly', type=int, default=None, help='Subtract an order N polynomial from red noise basis')
    parser.add_argument('--min-red-sin-scale', type=float, default=None, help='Limit red noise to reach a minimum to avoid numerical issues')

    parser.add_argument('--red-fit-cutoff', action='store_true', default=False,
                        help='Fit a cutoff frequency for the red noise model')
def setup_model(args, psr, parfile):
    global Tspan
    if args.red:
        if args.red_basis_subtract_cubic:
            args.red_basis_subtract_poly_order = 3
        elif args.red_basis_subtract_quad:
            args.red_basis_subtract_poly_order = 2

        # We scale Tspan if we are using tspan-mult
        nC = args.red_ncoeff
        Tspan = psr.toas.max() - psr.toas.min()

        Tspan *= args.tspan_mult
        nC *= int(args.tspan_mult)

        # Add the nC and flow parameters to the par file
        parfile.append("TNRedC {}\n".format(nC))
        parfile.append("TNRedFLow {}\n".format(np.log10(1.0 / args.tspan_mult)))

        if args.red_fit_cutoff:
            freq_cutoff = parameter.Uniform(10/Tspan, nC/Tspan, to_par=to_par)("freq_cutoff")
        else:
            freq_cutoff = parameter.Constant(10*nC/Tspan) # 10 times the max frequency
            parfile.append("TNRedC {}\n".format(nC))


        if args.min_red_sin_scale is not None:
            Aflat = np.sqrt((args.min_red_sin_scale**2 * 12.0 * np.pi ** 2  * Tspan )* const.fyr ** 3)
            log10Aflat = parameter.Constant(np.log10(Aflat))
            print(f"Min red noise scale: {args.min_red_sin_scale} -> {np.log10(Aflat)}")
        else:
            log10Aflat = parameter.Constant(-100)

        if args.red_at is None:
            red_log_sqrt_at = parameter.Constant(0)
        else:
            red_log_sqrt_at = parameter.Constant(0.5*np.log10(args.red_at))
            parfile.append("# red_at={}\n".format(args.red_at))

        # red noise (powerlaw with nC frequencies)
        if args.red_prior_log:
            log10_A = parameter.Uniform(args.Ared_min, args.Ared_max, to_par=to_par)
        else:
            log10_A = parameter.LinearExp(args.Ared_min, args.Ared_max, to_par=to_par)
        gamma = parameter.Uniform(args.red_gamma_min, args.red_gamma_max, to_par=to_par)

        log10_A2=None
        delta_gamma=None
        if args.double_powerlaw:
            if args.double_powerlaw_delta_gamma is not None:
                delta_gamma = parameter.Constant(args.double_powerlaw_delta_gamma, to_par=to_par)
                gamma2=None
            else:
                gamma2= parameter.Uniform(args.red_gamma2_min, args.red_gamma2_max, to_par=to_par)
                delta_gamma=None
            if args.red_prior_log:
                log10_A2 = parameter.Uniform(args.A2red_min, args.A2red_max, to_par=to_par)
            else:
                log10_A2 = parameter.LinearExp(args.A2red_min, args.A2red_max, to_par=to_par)
        if args.qp:
            tspan_days = (np.amax(psr.toas) - np.amin(psr.toas)) / 86400.0
            avg_cadence = tspan_days / len(psr.toas)

            # qp_f0_min = np.log10(2.0 / tspan_days) if args.qp_f0_min is None else args.qp_f0_min
            # qp_f0_max = np.log10(min(0.5 / avg_cadence, 1 / 20.0)) if args.qp_f0_max is None else args.qp_f0_max

            fmin = nC / tspan_days / args.tspan_mult
            qp_p_max = tspan_days / args.qp_p_min_np if args.qp_p_max is None else args.qp_p_max
            qp_p_min = max(2 * avg_cadence, 1.2 / fmin) if args.qp_p_min is None else args.qp_p_min
            qp_f0_max = np.log10(1.0 / qp_p_min)  ## for below calculation
            if 10 ** qp_f0_max > fmin:
                fstep = 1.0 / tspan_days / args.tspan_mult
                needC = (10 ** qp_f0_max) / fstep
                print("ERROR: Not enough components for QP model QPF0 max is {}, fmin is {}".format(10 ** qp_f0_max,
                                                                                                    fmin))
                print("       Would need nC={}".format(needC))
                sys.exit(1)

            log_qp_ratio = parameter.Uniform(-2, args.qp_ratio_max, to_par=to_par)
            # log_f0 = parameter.Uniform(qp_f0_min, qp_f0_max, to_par=to_par)
            qp_period = parameter.Uniform(qp_p_min, qp_p_max, to_par=to_par)
            sig = parameter.Uniform(1e-3, args.qp_sigma_max, to_par=to_par)
            lam = parameter.Uniform(0.01, 10, to_par=to_par)
            df = 86400.0 / (np.amax(psr.toas) - np.amin(psr.toas))
            df = parameter.Constant(df)
            if args.double_powerlaw:
                pl = double_powerlaw_qp(log10_A1=log10_A, gamma=gamma,freq_cutoff=freq_cutoff,
                                        log10_A2=log10_A2, delta_gamma=delta_gamma,
                                        log_qp_ratio=log_qp_ratio, p=qp_period, sig=sig, lam=lam,
                                 fmin=df,log_at=red_log_sqrt_at, log10_Aflat=log10Aflat,gamma2=gamma2)
            else:
                pl = powerlaw_qp(log10_A=log10_A, gamma=gamma,freq_cutoff=freq_cutoff, log_qp_ratio=log_qp_ratio, p=qp_period, sig=sig, lam=lam,
                                 fmin=df,log_at=red_log_sqrt_at, log10_Aflat=log10Aflat)
        else:
            if args.double_powerlaw:
                pl = double_powerlaw(log10_A1=log10_A, gamma=gamma, freq_cutoff=freq_cutoff,log10_A2 = log10_A2, delta_gamma=delta_gamma,
                                     log_at=red_log_sqrt_at, log10_Aflat=log10Aflat,gamma2=gamma2)
            else:
                pl = powerlaw(log10_A=log10_A, gamma=gamma,freq_cutoff=freq_cutoff,log_at=red_log_sqrt_at, log10_Aflat=log10Aflat)

        modes = None
        fmodes = None
        if not (args.red_minf is None):
            modes = 1.0 * np.arange(1, nC + 1) / Tspan
            sectoyr = 86400.0 * 365.25
            fmodes = modes[modes <= args.red_minf / sectoyr]
            modes = modes[modes > args.red_minf / sectoyr]
            print("Now have {} fourier components above {} yr^-1".format(len(modes), args.red_minf))

            if len(modes) < 1:
                print("Not enough modes! ", (1.0 * np.arange(1, nC + 1) / Tspan)[0])
                sys.exit(1)
            nC = len(modes)


        if args.red_basis_subtract_poly_order is not None:
             parfile.append("TNSubtractPoly 2\n")
             parfile.append(f"TNRedSubPolyOrd {args.red_basis_subtract_poly_order}\n")
             print(f"Subtracting order {args.red_basis_subtract_poly_order} polynomial from red noise basis")

        if args.red_nlog_below_tspan > 0:
            freq0 = 1.0 / Tspan
            low_modes = sorted(freq0 * args.red_nlog_factor**-np.arange(1, args.red_nlog_below_tspan + 1)) # log spaced below Tspan
            hi_modes = freq0 * np.arange(1,nC+1)
            modes = np.concatenate((low_modes,hi_modes))
            print(f"Modes: {len(low_modes)} {len(hi_modes)}")
            print(modes)
            rn = FourierBasisGP_SubPoly(spectrum=pl, components=None, Tspan=Tspan, modes=modes, poly_order=args.red_basis_subtract_poly_order)
            parfile.append(f"TNRedFLog {args.red_nlog_below_tspan}\n")
            parfile.append(f"TNRedFLog_factor {args.red_nlog_factor}\n")
        else:
            rn = FourierBasisGP_SubPoly(spectrum=pl, components=nC, Tspan=Tspan, modes=modes, poly_order=args.red_basis_subtract_poly_order)

        if not fmodes is None:
            for i, f in enumerate(fmodes):
                A = parameter.Uniform(-50, 50, to_par=to_par)("fsinA{}".format(i))
                B = parameter.Uniform(-50, 50, to_par=to_par)("fcosB{}".format(i))
                rn += deterministic_signals.Deterministic(sinusoid(A=A, B=B, om=2 * np.pi * f))



        return rn
    else:
        return None



@function 
def createfourierdesignmatrix_subpoly(toas, nmodes=30, Tspan=None, fref=1.0, modes=None, poly_order=3):
        
    T = Tspan if Tspan is not None else toas.max() - toas.min()

    if modes is None:
        modes = np.arange(1, nmodes + 1) / T
    else:
        nmodes = len(modes)
    f = modes

    N = len(toas)
    F = np.zeros((N, 2 * nmodes))

    # The sine/cosine modes
    F[:, ::2] = np.sin(2 * np.pi * toas[:, None] * f[None, :])
    F[:, 1::2] = np.cos(2 * np.pi * toas[:, None] * f[None, :])
    # fit and subtract the cubic term from each mode
    for imode in range(2*nmodes):
        coef = np.polyfit(toas, F[:, imode], poly_order)
        F[:, imode] -= np.polyval(coef, toas)

    Ffreqs = np.repeat(f, 2)

    return F, Ffreqs


def FourierBasisGP_SubPoly(spectrum, components=20,
                      selection=Selection(selections.no_selection),
                      Tspan=None, name="red_noise", modes=None, poly_order=3):
    """Convenience function to return a BasisGP class with a
    fourier basis."""

    if poly_order == None:
        return gp_signals.FourierBasisGP(spectrum, components=components, selection=selection, Tspan=Tspan, name=name, modes=modes)
    else:
        basis = createfourierdesignmatrix_subpoly(nmodes=components, Tspan=Tspan, modes=modes, poly_order=poly_order)
        # subtract a polynomial from each row of the deign matrix to remove the cubic term
        # this is equivalent to fitting for a polynomial and then subtracting it off the red noise

        BaseClass = gp_signals.BasisGP(spectrum, basis, selection=selection, name=name)

        class FourierBasisGP_SubPoly(BaseClass):
            signal_type = 'basis'
            signal_name = 'red noise'
            signal_id = name

        return FourierBasisGP_SubPoly

def to_par(self, p, chain):
    if "red_noise_log10_A1" in p:
        return "TNRedAmp", chain
    if "red_noise_log10_A2" in p:
        return "TNRedAmp2", chain
    if "red_noise_gamma2" in p:
        return "TNRedGam2", chain
    if "red_noise_delta_gamma" in p:
        return "TNRedDeltaGam", chain
    if "red_noise_log10_A" in p:
        return "TNRedAmp", chain
    if "red_noise_gamma" in p:
        return "TNRedGam", chain
    elif "_log_qp_ratio" in p:
        return "TN_QpRatio", chain
    elif "_log_f0" in p:
        return "TN_QpF0", chain
    elif "_p" in p:
        return "TN_QpPeriod", chain
    elif "_sig" in p:
        return "TN_QpSig", chain
    elif "_lam" in p:
        return "TN_QpLam", chain
    elif "fsin" in p:
        return p, chain
    elif "fcos" in p:
        return p, chain
    elif "freq_cutoff" in p:
        df = 1/Tspan
        return "TNRedC", np.floor(chain/df).astype(int)
    else:
        return None


from .chol_red_model import qp_term_cutoff


@function
def powerlaw(f, log10_A, gamma,freq_cutoff, components=2, log_at=0, log10_Aflat=-100):
    log10_A_1yr = log10_A - gamma * log_at
    df = np.diff(np.concatenate((np.array([0]), f[::components])))
    df2 = np.repeat(df, components)

    flat = ((10**log10_Aflat) ** 2) / 12.0 / np.pi ** 2 * const.fyr ** (-3) * df2
    red= ((10 ** log10_A_1yr) ** 2) / 12.0 / np.pi ** 2 * const.fyr ** (gamma - 3) * f ** (-gamma) * df2
    ret = flat+red
    return apply_cutoff(f,ret, freq_cutoff)

@function
def powerlaw_qp(f, log10_A, gamma,freq_cutoff, log_qp_ratio, p, sig, lam, fmin, components=2,log_at=0,log10_Aflat=-100):
    df = np.diff(np.concatenate((np.array([0]), f[::components])))
    df2 = np.repeat(df, components)

    log10_A_1yr = log10_A - gamma * log_at

    red = (10 ** log10_A_1yr) ** 2 / 12.0 / np.pi ** 2 * const.fyr ** (gamma - 3) * f ** (-gamma) * df2
    flat = ((10**log10_Aflat) ** 2) / 12.0 / np.pi ** 2 * const.fyr ** (-3) * df2

    f0 = 1.0 / p
    logPqp = log_qp_ratio + np.log10(
        (10 ** log10_A_1yr) ** 2 / 12.0 / np.pi ** 2 * const.fyr ** (gamma - 3) * (f0 * 365.25 * const.fyr) ** (
            -gamma) * df2)
    qp = qp_term_cutoff(f / const.fyr, logPqp, f0, sig, lam, fmin)
    ret = red + qp + flat
    return apply_cutoff(f,ret, freq_cutoff)


@function
def double_powerlaw(f, log10_A1, gamma,freq_cutoff, log10_A2, delta_gamma, components=2, log_at=0,log10_Aflat=-100,gamma2=None):
    gamma1=gamma
    if gamma2 is None:
        gamma2 = gamma - delta_gamma

    log10_A1_1yr = log10_A1  - gamma1 * log_at
    log10_A2_1yr = log10_A2  - gamma2 * log_at

    df = np.diff(np.concatenate((np.array([0]), f[::components])))
    df2 = np.repeat(df, components)
    flat = ((10**log10_Aflat) ** 2) / 12.0 / np.pi ** 2 * const.fyr ** (-3) * df2

    red1 = (10 ** log10_A1_1yr) ** 2 / 12.0 / np.pi ** 2 * const.fyr ** (gamma1 - 3) * f ** (-gamma1) * df2
    red2 = (10 ** log10_A2_1yr) ** 2 / 12.0 / np.pi ** 2 * const.fyr ** (gamma2 - 3) * f ** (-gamma2) * df2
    ret = red1 + red2 + flat
    return apply_cutoff(f,ret, freq_cutoff)


@function
def double_powerlaw_qp(f, log10_A1, gamma,freq_cutoff, log10_A2, delta_gamma, log_qp_ratio, p, sig, lam, fmin, components=2,log_at=0,log10_Aflat=-100,gamma2=None):
    gamma1=gamma
    if gamma2 is None:
        gamma2 = gamma - delta_gamma

    log10_A1_1yr = log10_A1 - gamma1 * log_at
    log10_A2_1yr = log10_A2 - gamma2 * log_at

    df = np.diff(np.concatenate((np.array([0]), f[::components])))
    df2 = np.repeat(df, components)
    red1 = (10 ** log10_A1_1yr) ** 2 / 12.0 / np.pi ** 2 * const.fyr ** (gamma1 - 3) * f ** (-gamma1) * df2
    red2 = (10 ** log10_A2_1yr) ** 2 / 12.0 / np.pi ** 2 * const.fyr ** (gamma2 - 3) * f ** (-gamma2) * df2
    red = red1 + red2
    flat = ((10**log10_Aflat) ** 2) / 12.0 / np.pi ** 2 * const.fyr ** (-3) * df2

    f0 = 1.0 / p
    logPqp = log_qp_ratio + np.log10(
        (10 ** log10_A1_1yr) ** 2 / 12.0 / np.pi ** 2 * const.fyr ** (gamma1 - 3) * (f0 * 365.25 * const.fyr) ** (
            -gamma1) * df2) + \
             np.log10(
                 (10 ** log10_A2_1yr) ** 2 / 12.0 / np.pi ** 2 * const.fyr ** (gamma2 - 3) * (f0 * 365.25 * const.fyr) ** (
                     -gamma2) * df2)
    qp = qp_term_cutoff(f / const.fyr, logPqp, f0, sig, lam, fmin)
    ret= red + qp + flat
    
    return apply_cutoff(f,ret, freq_cutoff)

@function
def sinusoid(toas, A, B, om):
    return B * np.sin(om * toas) + A * np.cos(om * toas)

def apply_cutoff(f,ret, freq_cutoff):
    ret[f>freq_cutoff] = 1e-40
    return ret