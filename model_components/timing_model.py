import numpy as np
import sys

from enterprise.signals import deterministic_signals
from enterprise.signals import signal_base

from astropy import coordinates as coord
from astropy import units as units

from . import xparameter as parameter

name = "BasicTimingModel"
argdec = "Basic pulsar spin and astrometric parameters."


def setup_argparse(parser):
    parser.add_argument('--f2', type=float, default=0, help='range of f2 to search')
    parser.add_argument('--pm', action='store_true', help='Fit for PMRA+PMDEC')
    parser.add_argument('--px', action='store_true', help='Fit for parallax')
    parser.add_argument('--px-range', type=float, default=10, help='Max parallax to search')
    parser.add_argument('--pm-angle', action='store_true', help='Fit for PM + angle rather than by PMRA/PMDEC')
    parser.add_argument('--pm-range', type=float, default=10, help='Search range for proper motion (deg/yr)')
    parser.add_argument('--pm-ecliptic', action='store_true', help='Generate ecliptic coords for proper motion')


def setup_model(args, psr, parfile):
    global orig_f2
    global orig_pmra
    global orig_pmdec
    for line in parfile:
        e = line.strip().split()
        if len(e) > 1:
            if e[0] == "F0":
                f0 = float(e[1])
            if e[0] == "F2":
                orig_f2 = float(e[1])
            elif e[0] == "PEPOCH":
                pepoch = float(e[1])
            elif e[0] == "POSEPOCH":
                posepoch = float(e[1])
            elif e[0] == "PMRA":
                orig_pmra = float(e[1])
            elif e[0] == "PMDEC":
                orig_pmdec = float(e[1])
            elif e[0] == "RAJ":
                psr_ra = e[1]
            elif e[0] == "DECJ":
                psr_dec = e[1]

    psr_coord = coord.SkyCoord(psr_ra, psr_dec, unit=(units.hourangle, units.degree))
    psr.coord = psr_coord
    psr.earthssb = psr.t2pulsar.earth_ssb
    psr.observatory_earth = psr.t2pulsar.observatory_earth

    if posepoch < 0:
        print("Assume POSEPOCH=PEPOCH")
        posepoch = pepoch

    if args.px and "PX" in psr.fitpars:
        print("ERROR: Can't fit for PX in both MCMC and least-squares")
        sys.exit(1)
    if args.pm and "PMRA" in psr.fitpars:
        print("ERROR: Can't fit for PMRA in both MCMC and least-squares")
        sys.exit(1)
    if args.pm and "PMDEC" in psr.fitpars:
        print("ERROR: Can't fit for PMDEC in both MCMC and least-squares")
        sys.exit(1)

    orig_toas = psr.t2pulsar.toas()
    issorted = np.all(orig_toas[:-1] <= orig_toas[1:])

    if (args.pm or args.px) and not issorted:
        print("ERROR: Toas must be sorted or wierd things happen with paralax and proper motion")
        sys.exit(1)

    components = []

    f2range = args.f2
    if f2range > 0:
        f2 = fit_f2(pepoch=parameter.Constant(pepoch), f0=parameter.Constant(f0),
                    valf2=parameter.Uniform(-f2range / 1e-27, f2range / 1e-27)('F2'))
        f2f = deterministic_signals.Deterministic(f2, name='f2')
        components.append(f2f)
    if args.pm:
        pmrange = args.pm_range
        pmra = parameter.Uniform(-pmrange + orig_pmra, pmrange + orig_pmra, to_par=to_par)
        pmdec = parameter.Uniform(-pmrange + orig_pmdec, pmrange + orig_pmdec, to_par=to_par)
        pm = deterministic_signals.Deterministic(
            fit_pm(posepoch=posepoch, pmra=pmra, pmdec=pmdec, orig_pmra=orig_pmra, orig_pmdec=orig_pmdec), name='pm')
        components.append(pm)

    if args.px:
        pxp = parameter.Uniform(0, args.px_range, to_par=to_par)
        px = deterministic_signals.Deterministic(fit_px(px=pxp), name='px')
        components.append(px)

    if args.pm_angle:
        pmrange = args.pm_range
        pmamp = parameter.Uniform(0, pmrange, to_par=to_par)
        pmangle = parameter.Uniform(0, 2 * np.pi, to_par=to_par)
        pm = deterministic_signals.Deterministic(fit_pm_angle(posepoch=posepoch, amp=pmamp, angle=pmangle), name='pm')
        components.append(pm)

    if len(components) > 0:
        model = components[0]
        for m in components[1:]:
            model += m
        return model
    else:
        return None


def to_par(self, p, chain):
    global orig_f2
    global orig_pmra
    global orig_pmdec
    if p == "f2":
        return "F2", chain*1e-27 + orig_f2
    if "px_px" in p:
        return "PX", chain
    if "pm_angle" in p:
        return "PM-angle", chain
    if "pm_amp" in p:
        return "PM", chain
    if "pm_pmra" in p:
        return 'PMRA', chain + orig_pmra
    if "pm_pmdec" in p:
        return 'PMDEC', chain + orig_pmdec
    else:
        return None


@signal_base.function
def fit_f2(toas, pepoch, f0, valf2):
    x = (toas - pepoch * 86400.0)
    return -1e-27 * valf2 * x * x * x / 6.0 / f0


@signal_base.function
def fit_pm_angle(toas, pos, earthssb, posepoch, amp, angle, orig_pmra=0, orig_pmdec=0):
    pmra = amp * np.cos(angle)
    pmdec = amp * np.sin(angle)
    return fit_pm(toas, pos, earthssb, posepoch, pmra, pmdec, orig_pmra, orig_pmdec)


@signal_base.function
def fit_px(toas, pos, earthssb, observatory_earth, px):
    # This is just the linearised parallax function from tempo2...
    pxconv = 1.74532925199432958E-2 / 3600.0e3
    AULTSC = 499.00478364

    rca = earthssb[:, 0:3] + observatory_earth[:, 0:3]
    rr = np.sum(rca * rca, axis=1)

    rcos1 = np.sum(pos * rca, axis=1)

    return px * 0.5 * pxconv * (rr - rcos1 * rcos1) / AULTSC


@signal_base.function
def fit_pm(toas, pos, earthssb, posepoch, pmra, pmdec, orig_pmra=0, orig_pmdec=0):
    # This is just the linearised proper motion function from tempo2...
    AULTSC = 499.00478364
    rce = earthssb[:, 0:3]
    re = np.sqrt(np.sum(rce * rce, axis=1))

    psrra = np.arctan2(pos[1], pos[0]) + 2 * np.pi
    psrdec = np.arcsin(pos[2] / np.sqrt(np.sum(pos * pos)))

    axy = rce[:, 2] / AULTSC
    s = axy / (re / AULTSC)
    deltae = np.arctan2(s, np.sqrt(1.0 - s * s))
    alphae = np.arctan2(rce[:, 1], rce[:, 0])

    # Convert radians/second to degrees per year
    # 60.0*60.0*1000.0*86400.0*365.25/24.0/3600.0
    rs2dy = (180.0 / np.pi) * 1314900000

    # This will be pmra in rad/s
    v_pmra = (pmra - orig_pmra) / rs2dy / np.cos(psrdec)
    v_pmdec = (pmdec - orig_pmdec) / rs2dy

    t0 = (toas / 86400.0 - posepoch)
    d_pmra = re * np.cos(deltae) * np.cos(psrdec) * np.sin(psrra - alphae) * t0
    d_pmdec = re * (np.cos(deltae) * np.sin(psrdec) * np.cos(psrra - alphae) - np.sin(deltae) * np.cos(psrdec)) * t0

    return v_pmra * d_pmra + v_pmdec * d_pmdec
