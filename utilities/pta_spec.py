import numpy as np
from run_enterprise.spectrum import cholspec
from enterprise.pulsar import Pulsar
from scipy import linalg
import matplotlib.pyplot as plt
from tqdm import tqdm

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

def margin_full_designmatrix(psr):
    return psr.t2pulsar.designmatrix(updatebats=False)

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


def component_correlations(pcov, component_a, component_b):
    """Return two canonical correlations, largest first."""
    ia = slice(2 * component_a, 2 * component_a + 2)
    ib = slice(2 * component_b, 2 * component_b + 2)
    
    Vaa = pcov[ia, ia]
    Vbb = pcov[ib, ib]
    Vab = pcov[ia, ib]
    
    La = linalg.cholesky(Vaa, lower=True)
    Lb = linalg.cholesky(Vbb, lower=True)
    
    # R = inv(La) @ Vab @ inv(Lb.T)
    tmp = linalg.solve_triangular(La, Vab, lower=True)
    R = linalg.solve_triangular(Lb, tmp.T, lower=True).T
    
    rho = linalg.svdvals(R)
    return np.clip(rho, 0.0, 1.0)

def ptaspec(par, tim, reffreq=1400.0, fmax=None, fmax_limit=None):
    """
    Compute the PSD for a pulsar at reffreq with a par/tim file and a model.

    we use a GLS method fitting for sinusoids of varying chromaticity

    """

    psr=Pulsar(par,tim,drop_t2pulsar=False,sort=False)
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
                # getC's PSD includes 1/(12*pi**2), unlike the DM model.
                dmAmp = (float(split_line[1])
                         + np.log10(1400.0 ** -2 / 2.41e-4)
                         + 0.5 * np.log10(12.0 * np.pi**2))
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


    order = np.argsort(psr.toas, kind="stable")
    toas = psr.toas[order]/86400.0 # convert back to days for sanity
    freqssb = psr.freqs[order]
    residuals = np.asarray(psr.residuals, dtype=float)[order]
    toaerrs = effective_toaerrs(psr.t2pulsar)[order]
    nuisance_matrix = margin_full_designmatrix(psr)[order, :]
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
        Cdm = cholspec.getC(toas, dmAmp, dmGam, fc_yr=fc_yr, freqssb=freqssb, chromatic_index=2.0)
        exponents.append(2)
        Cs.append(Cdm)
    if chromAmp is not None:
        print(f"Chromatic A={chromAmp}, gamma={chromGam}, fc_yr={fc_yr}, chromatic_index={chromatic_index}")
        Cchrom = cholspec.getC(toas, chromAmp, chromGam, fc_yr=fc_yr, freqssb=freqssb, chromatic_index=chromatic_index)
        exponents.append(chromatic_index)
        Cs.append(Cchrom)

    if len(Cs)==0:
        raise ValueError("No red noise, DM noise, or chromatic noise parameters found in par file.")

    print("Done making covariance matrix.")
    C = np.sum(Cs, axis=0)
    C_white = np.diag(toaerrs**2)
    C += C_white

    M = nuisance_matrix
    print("Design matrix shape:", M.shape)
    print("Whiten and project out nuisance parameters")
    # Now do linear algebra.
    L = linalg.cholesky(C, lower=True)
    def whiten(a):
        return linalg.solve_triangular(L, a, lower=True)
    # Whiten the nuisance design matrix.
    Mw = whiten(M)
    norms = np.linalg.norm(Mw, axis=0)
    keep = norms > 0
    Mw = Mw[:, keep] / norms[keep]
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
    y = project(whiten(residuals))

    f_Nyq = fmax
    delta_f = 1.0/(tspan) # frequency step, per day
    N_f = int(fmax / delta_f)
    freqs = np.arange(1, N_f + 1) * delta_f
    print("Computing {} frequency channels from {:.3g} to {:.3g} per day".format(N_f, freqs[0], freqs[-1]))
    omegas = 2 * np.pi * freqs
    total_psds=[]
    component_psds=[]
    coefficient_covariances=[]
    rho_max=[]
    rho_min=[]
    sec_per_year = 365.25 * 86400.0
    rng = np.random.default_rng(1234)
    white_quantiles = []
    conversion = 86400.0 / sec_per_year**3
    white_psd=[]
    expected_psd=[]

    for i, omega in enumerate(tqdm(omegas, desc="Processing modes")):
        # Design matrix for this mode.
        M = design_matrix(toas, freqssb, [omega], exponents, reffreq=reffreq)
        Mw = whiten(M)
        Z = project(Mw)

        # Least-squares fit for this mode.
        p, pcov, chisq = qrfit(Z, y)
        white_components, white_total, B = white_psd_estimator(
        Z, L, project, C_white, tspan)

        white_components *= conversion
        white_total *= conversion
        coeff = p.reshape(len(exponents), 2)
        component_psd = 0.5 * tspan * np.sum(coeff**2, axis=1)* conversion
        total_coeff = coeff.sum(axis=0)
        total_psd = 0.5 * tspan * np.sum(total_coeff**2)* conversion
        component_psds.append(component_psd)
        coefficient_covariances.append(pcov)
        total_psds.append(total_psd)
        rho = component_correlations(pcov, 0, 1)
        rho_max.append(rho[0])
        rho_min.append(rho[1])
        white_psd.append(white_total)
        # white_quantiles.append(
        # white_total_psd_quantiles(
        #     B, tspan, conversion=conversion, rng=rng
        # ))
        # Lsum = np.tile(np.eye(2), (1, len(exponents)))
        # expected_full_psd = (
        #     0.5 * tspan
        #     * np.trace(Lsum @ pcov @ Lsum.T)
        #     * conversion
        #     )
        
        # expected_psd.append(expected_full_psd)


    component_psds = np.asarray(component_psds)
    total_psds = np.asarray(total_psds)
    coefficient_covariances = np.asarray(coefficient_covariances)
    rho_max = np.asarray(rho_max)
    rho_min = np.asarray(rho_min)
    expected_psd = np.asarray(expected_psd)

    



    fig, ax = plt.subplots()
    ax.set_xscale('log')
    ax2=ax.twinx()
    ax.set_yscale('log')
    ax2.set_ylim(0, 1.0)
    ax.plot(freqs, total_psds, label="Total PSD", color='k', lw=2)
    ax2.plot(freqs, rho_max, color='k', alpha=0.5, label='Max Correlation',ls=':')

    colors={0:'blue', 2:'orange', chromatic_index:'green'}

    for exp in exponents:
        ax.plot(freqs, component_psds[:,exponents.index(exp)], label=r"$\alpha={}$".format(exp),ls='--',color=colors[exp])
    
    ax.plot(freqs, white_psd, label="White PSD", color='gray', ls='-.')

    # white_quantiles = np.asarray(white_quantiles)

    # ax.fill_between(
    #     freqs,
    #     white_quantiles[:, 0],
    #     white_quantiles[:, 2],
    #     color="gray",
    #     alpha=0.2,
    #     label="White-only 5–95% interval",
    # )

    # ax.plot(
    #     freqs,
    #     white_quantiles[:, 1],
    #     color="gray",
    #     ls=":",
    #     label="White-only median",
    # )

    # ax.plot(
    #     freqs,
    #     expected_psd,
    #     color="red",
    #     ls="-",
    #     label="Expected PSD",
    # )


    # save the yscale limits for the PSD plot
    ylims = ax.get_ylim()

    # Note that the model PSDs are two-sided, while the periodogram is one-sided.
    # The factor of 2 is included in the PSD functions, so we don't need to include it here.
    if redAmp is not None:
        red_model_psd = 2*cholspec.PSD(freqs, fc_yr, redAmp, redGam)
        ax.plot(freqs, red_model_psd, color="blue", ls="-", label="Red model PSD")
    if dmAmp is not None:
        dm_model_psd = 2*cholspec.PSD(freqs, fc_yr, dmAmp, dmGam)
        ax.plot(freqs, dm_model_psd, color="orange", ls="-", label="DM model PSD")
    if chromAmp is not None:
        chrom_model_psd = 2*cholspec.PSD(freqs, fc_yr, chromAmp, chromGam)
        ax.plot(freqs, chrom_model_psd, color="green", ls="-", label="Chrom model PSD")
    # restore the yscale limits for the PSD plot
    ax.set_ylim(ylims)
    
    ax.set_xlabel("Frequency (per day)")
    ax.set_ylabel("PSD")
    ax2.set_ylabel("Correlation")
    ax.set_title("PSD for {}".format(par))
    ax.legend()
    #ax2.legend()
    plt.show()

    return freqs, total_psds, component_psds, coefficient_covariances, rho_max, rho_min






import numpy as np
from scipy import linalg


def white_psd_estimator(Z, L, project, white_cov, span):
    """
    Expected white-noise contribution to a joint GLS periodogram.

    Parameters
    ----------
    Z : (n_toa, 2*n_component) array
        Projected, whitened Fourier design:
            Z = project(whiten(X))

        Columns must alternate sine/cosine for each component.

    L : (n_toa, n_toa) array
        Lower Cholesky factor of the FULL working covariance:
            C = L @ L.T

    project : callable
        Nuisance projection in whitened coordinates.
        Must accept a matrix with n_toa rows.

    white_cov : array
        Either:
        - shape (n_toa,): white-noise VARIANCES, not standard deviations;
        - shape (n_toa, n_toa): full white-noise covariance.

        This must be the white-noise contribution only, not full C.

    span : float
        Observing span in the same time units used for Fourier frequencies.
        For times in days, outputs have units residual_units**2 * days.

    Returns
    -------
    component_psd : (n_component,) array
        Expected white contribution to each component periodogram.

    total_psd : float
        Expected white contribution to the total reference-frequency
        periodogram, including cross-component coefficient covariance.

    B : (2*n_component, 2*n_component) array
        White-noise contribution to the fitted coefficient covariance.
    """
    Z = np.asarray(Z, dtype=float)
    white_cov = np.asarray(white_cov, dtype=float)

    n_toa, n_param = Z.shape

    if n_param == 0 or n_param % 2:
        raise ValueError("Z must contain sine/cosine column pairs")
    if not np.isfinite(span) or span <= 0:
        raise ValueError("span must be positive and finite")
    if n_toa < n_param:
        raise ValueError("Too few observations for the Fourier fit")

    Q, R = linalg.qr(Z, mode="economic")

    if np.linalg.matrix_rank(R) < n_param:
        raise ValueError("Projected Fourier design is rank deficient")

    # The coefficient estimator is:
    # A = inv(R) @ Q.T @ P @ inv(L)
    #
    # Construct its transpose-side factors using triangular solves,
    # without forming inv(L) or the dense projection matrix P.
    H = linalg.solve_triangular(
        L.T, project(Q), lower=False
    )
    A = linalg.solve_triangular(R, H.T, lower=False)

    if white_cov.ndim == 1:
        if white_cov.shape != (n_toa,):
            raise ValueError("White variances must have length n_toa")
        if np.any(~np.isfinite(white_cov)) or np.any(white_cov < 0):
            raise ValueError("White variances must be finite and nonnegative")

        B = (A * white_cov[None, :]) @ A.T

    elif white_cov.ndim == 2:
        if white_cov.shape != (n_toa, n_toa):
            raise ValueError("White covariance must have shape (n_toa, n_toa)")

        B = A @ white_cov @ A.T

    else:
        raise ValueError("white_cov must be a vector or matrix")

    B = 0.5 * (B + B.T)

    # Individual component periodograms: span/2 * (sin**2 + cos**2).
    component_psd = (
        0.5 * span * np.diag(B).reshape(-1, 2).sum(axis=1)
    )

    # At the reference frequency, sum all sine coefficients and
    # all cosine coefficients before calculating power.
    n_component = n_param // 2
    Lsum = np.tile(np.eye(2), (1, n_component))
    total_psd = 0.5 * span * np.trace(Lsum @ B @ Lsum.T)

    return component_psd, float(total_psd), B



def effective_toaerrs(psr):
    """Return TOA errors in seconds with TNEF/TNEQ applied."""
    err = np.asarray(psr.toaerrs).copy()

    if psr.noisemodel is None:
        return err*1e-6

    # Apply EFAC/TNEF first
    for name, noise in psr.noisemodel.items():
        if name.startswith("efac"):
            mask = psr.flagvals(noise.flag) == noise.flagval
            err[mask] *= noise.val

    # Then add EQUAD/TNEQ in quadrature
    for name, noise in psr.noisemodel.items():
        if name.startswith("equad"):
            mask = psr.flagvals(noise.flag) == noise.flagval
            err[mask] = np.sqrt(err[mask]**2 + noise.val**2)

    return err*1e-6  # convert to seconds

def check_white_psd(
    Z, whiten, project, white_variances, span,
    n_sim=1000, seed=1234,
):
    """
    Monte Carlo white-noise contribution using the actual GLS pipeline.

    white_variances: diagonal white-noise variances in residual units^2.
    """
    rng = np.random.default_rng(seed)

    white_variances = np.asarray(white_variances, dtype=float)
    n_toa = Z.shape[0]
    n_component = Z.shape[1] // 2

    # Independent white-noise realisations, one per column.
    noise = (
        np.sqrt(white_variances)[:, None]
        * rng.standard_normal((n_toa, n_sim))
    )

    # Apply exactly the same transformations as for the real residuals.
    simulated_y = project(whiten(noise))

    Q, R = linalg.qr(Z, mode="economic")
    p = linalg.solve_triangular(R, Q.T @ simulated_y)

    # Shape: component, quadrature, simulation.
    coeff = p.reshape(n_component, 2, n_sim)

    component_psds = 0.5 * span * np.sum(coeff**2, axis=1)
    total_psds = (
        0.5 * span * np.sum(coeff.sum(axis=0)**2, axis=0)
    )

    return component_psds.mean(axis=1), total_psds.mean()


def white_total_psd_quantiles(
    B, span, conversion=1.0, n_draw=20000, rng=None
):
    """
    White-only marginal PSD quantiles at one Fourier frequency.

    B: white contribution to the full fitted coefficient covariance.
    Returns the 5th, 50th and 95th percentiles.
    """
    if rng is None:
        rng = np.random.default_rng()

    n_component = B.shape[0] // 2
    Lsum = np.tile(np.eye(2), (1, n_component))

    Btotal = Lsum @ B @ Lsum.T
    Btotal = 0.5 * (Btotal + Btotal.T)

    eigenvalues = np.linalg.eigvalsh(Btotal)

    # Allow only tiny negative eigenvalues from roundoff.
    scale = np.max(np.abs(eigenvalues))
    if scale > 0 and eigenvalues.min() < -1e-10 * scale:
        raise ValueError("Total white coefficient covariance is not PSD")

    eigenvalues = np.maximum(eigenvalues, 0.0)

    # Rotation does not affect the sum of squared coefficients.
    z = rng.standard_normal((n_draw, 2))
    powers = (
        0.5 * span * conversion
        * np.sum(eigenvalues * z**2, axis=1)
    )

    return np.quantile(powers, [0.05, 0.50, 0.95])


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("par", help="Par file")
    parser.add_argument("tim", help="Tim file")
    parser.add_argument("--reffreq", type=float, default=1400.0, help="Reference frequency in MHz")
    parser.add_argument("--fmax", type=float, default=1.0/180,help="Maximum frequency in per day")
    parser.add_argument("--fmax_limit", type=float, default=None,help="Limit on maximum frequency in per day")
    args = parser.parse_args()

    ptaspec(args.par, args.tim,  reffreq=args.reffreq, fmax=args.fmax, fmax_limit=args.fmax_limit)