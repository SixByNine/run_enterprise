import numpy as np
from run_enterprise.spectrum import cholspec
from enterprise.pulsar import Pulsar
from scipy import linalg
import matplotlib.pyplot as plt

def qrfit(dm, b):
    """Least squares for data with known, whitened covariance."""
    q, r = np.linalg.qr(dm, mode="reduced")

    if np.linalg.matrix_rank(r) < dm.shape[1]:
        raise ValueError(
            "Fourier design is rank deficient after nuisance projection"
        )

    param = linalg.solve_triangular(r, q.T @ b)

    # Conditional coefficient covariance: (dm.T @ dm)^-1
    ri = linalg.solve_triangular(r, np.eye(r.shape[0]))
    param_cov = ri @ ri.T

    postfit = b - dm @ param
    chisq = postfit @ postfit

    return param, param_cov, chisq


def design_matrix(toas, freqssb, omega, exponents, reffreq=1400.0):
    """
    Design matrix for the GLS fit.
    
    toas: array of TOAs in days
    freqssb: ToA frequencies
    omega: fluctuation frequency of this mode.
    exponents: list of chromatic exponents to fit, e.g. [0,2]
    reffreq: reference frequency
    """
    exp=np.array(exponents)
    M = np.zeros((len(toas), len(exponents)*2))
    for i, exp in enumerate(exponents):
        M[:, i*2] = np.sin(omega * toas) * (freqssb / reffreq) ** -exp
        M[:, i*2 + 1] = np.cos(omega * toas) * (freqssb / reffreq) ** -exp
    return M

def margin_designmatrix(toas, freqssb, exponents, reffreq=1400.0):
    """
    Design matrix for marginalising over a constant term.
    Currently just marginalise over a cosntant term,
    but maybe should have polynomials as well.
    """
    exp=np.array(exponents)
    M = np.zeros((len(toas), len(exponents)))
    for i, exp in enumerate(exponents):
        M[:, i] = (freqssb / reffreq) ** -exp
    return M


def ptaspec(par, tim, reffreq=1400.0, fmax=None, fmax_limit=None):
    """
    Compute the PSD for a pulsar at reffreq with a par/tim file and a model.

    we use a GLS method fitting for sinusoids of varying chromaticity

    """

    psr=Pulsar(par,tim,drop_t2pulsar=False)
    is_binary=False
    redAmp=None
    redGam=None
    F0=F1=PB=0
    dmAmp=None
    dmGam=None
    chromAmp=None
    chromGam=None
    chromatic_index=0.0
    qp_period=None
    qp_ratio=None
    qp_logf0=None
    qp_sig=None
    qp_lam=None
    with open(par) as f:
        for line in f:
            split_line = line.split()
            if len(split_line)==0:
                continue
            if split_line[0] == 'TNRedAmp':
                redAmp = float(split_line[1])
            elif split_line[0] == 'TNRedGam':
                redGam = float(split_line[1])
            elif split_line[0] == 'TNChromAmp':
                chromAmp = float(split_line[1])
            elif split_line[0] == 'TNChromGam':
                chromGam = float(split_line[1])
            elif split_line[0] == 'TNChromIdx':
                chromatic_index = float(split_line[1])
            elif split_line[0] == 'TNDMAmp':
                # NOTE: This converts to amplitude at 1400 MHz, which is what we want for the PSD.
                dmAmp = float(split_line[1])+ np.log10(1400.0 ** -2 / 2.41e-4)
            elif split_line[0] == 'TNDMGam':
                dmGam = float(split_line[1])
            elif split_line[0] == 'F0':
                F0 = float(split_line[1])
            elif split_line[0] == 'F1':
                F1 = float(split_line[1])
            elif split_line[0] == 'PB':
                is_binary=True
                PB=float(split_line[1])
            elif split_line[0] == 'TN_QpPeriod':
                qp_period = float(split_line[1])
                qp_logf0 = -np.log10(qp_period)
            elif split_line[0] == 'TN_QpRatio':
                qp_ratio = float(split_line[1])
            elif split_line[0] == 'TN_QpLam':
                qp_lam = float(split_line[1])
            elif split_line[0] == 'TN_QpSig':
                qp_sig = float(split_line[1])


    psr.sort_data()
    toas = psr.toas/86400.0 # convert back to days for sanity
    ntoa=len(toas)

    toa_spacing = np.diff(toas)
    if fmax==None:
        fmax = 0.5/np.median(toa_spacing)
        if not fmax_limit is None:
            fmax=min(fmax_limit,fmax)
    tspan = np.amax(toas)-np.amin(toas)
    fc_yr = 365.25/tspan / 4.0


    Cs=[]
    exponents=[]
    print("Make covariance matrix")
    if redAmp is not None:
        print(f"Red A={redAmp}, gamma={redGam}, fc_yr={fc_yr}")
        Cred = cholspec.getC(toas, redAmp, redGam,fc_yr=fc_yr, log_qp_ratio=qp_ratio, log_f0=qp_logf0, sig=qp_sig, lam=qp_lam)
        exponents.append(0)
        Cs.append(Cred)
    if dmAmp is not None:
        print(f"DM A={dmAmp}, gamma={dmGam}, fc_yr={fc_yr}")
        Cdm = cholspec.getC(toas, dmAmp, dmGam, fc_yr=fc_yr, freqssb=psr.freqs, chromatic_index=2.0)
        exponents.append(2)
        Cs.append(Cdm)
    if chromAmp is not None:
        print(f"Chromatic A={chromAmp}, gamma={chromGam}, fc_yr={fc_yr}, chromatic_index={chromatic_index}")
        Cchrom = cholspec.getC(toas, chromAmp, chromGam, fc_yr=fc_yr, freqssb=psr.freqs, chromatic_index=chromatic_index)
        exponents.append(chromatic_index)
        Cs.append(Cchrom)

    if len(Cs)==0:
        raise ValueError("No red noise, DM noise, or chromatic noise parameters found in par file.")

    print("Done making covariance matrix.")
    C = np.sum(Cs, axis=0)
    white = np.diag(psr.toaerrs**2)
    C += white

    freqssb=psr.freqs 
    print(freqssb)
    M = margin_designmatrix(toas, freqssb, exponents, reffreq=reffreq)
    print("Design matrix shape:", M.shape)
    print("Whiten and project out nuisance parameters")
    # Now do linear algebra.
    L = linalg.cholesky(C, lower=True)
    def whiten(a):
        return linalg.solve_triangular(L, a, lower=True)
    # Whiten the nuisance design matrix.
    Mw = whiten(M)

    # SVD handles redundant nuisance columns safely.
    U, s, _ = linalg.svd(Mw, full_matrices=False)

    tol = (
        np.finfo(float).eps * max(Mw.shape) * s[0]
        if s.size else 0.0
    )
    n_nuisance = np.count_nonzero(s > tol)
    Qm = U[:, :n_nuisance]
    print("Number of nuisance parameters:", n_nuisance)
    def project(a):
        """Remove the whitened nuisance-column space."""
        return a - Qm @ (Qm.T @ a)

    # Prepare residuals once.
    y = project(whiten(np.asarray(psr.residuals, dtype=float)))

    f_Nyq = fmax
    delta_f = 1.0/(tspan) # frequency step, per day
    N_f = int(f_Nyq/delta_f) # number of spectral channels
    print("Computing {} frequency channels from {:.3g} to {:.3g} per day".format(N_f, delta_f, f_Nyq))
    freqs = np.linspace(delta_f,f_Nyq, N_f, endpoint=False) # per day
    omegas = 2 * np.pi * freqs
    psds=[]
    for i, omega in enumerate(omegas):
        # Design matrix for this mode.
        M = design_matrix(toas, freqssb, [omega], exponents, reffreq=reffreq)
        Mw = whiten(M)
        Z = project(Mw)

        # Least-squares fit for this mode.
        p, _, _ = qrfit(Z, y)
        # Compute the PSD for this mode.
        psd = np.sum(p ** 2) / 2.0
        psds.append(psd)



    plt.loglog(freqs, psds)
    plt.show()

    return freqs, np.array(psds)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("par", help="Par file")
    parser.add_argument("tim", help="Tim file")
    parser.add_argument("--reffreq", type=float, default=1400.0, help="Reference frequency in MHz")
    parser.add_argument("--fmax", type=float, default=1.0/180,help="Maximum frequency in per day")
    parser.add_argument("--fmax_limit", type=float, default=None,help="Limit on maximum frequency in per day")
    args = parser.parse_args()

    freqs, psds = ptaspec(args.par, args.tim,  reffreq=args.reffreq, fmax=args.fmax, fmax_limit=args.fmax_limit)