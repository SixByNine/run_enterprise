# run_enterprise
Single pulsar timing analysis using the [enterprise](https://doi.org/10.5281/zenodo.4059815) framework.
Sorry that there is no documentation! Maybe it will come one day.


You can cite this code and find a perminant URI for versions of it here:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5914351.svg)](https://doi.org/10.5281/zenodo.5914351)

Please cite also [enterprise](https://doi.org/10.5281/zenodo.4059815), [tempo2](https://bitbucket.org/psrsoft/tempo2), and your preferred sampler as appropriate.

## Usage

```
usage: run_enterprise [-h] [--outdir OUTDIR] [--no-red-noise] [--Ared-max ARED_MAX] [--Ared-min ARED_MIN] [--red-at RED_AT] [--red-gamma-max RED_GAMMA_MAX]
                      [--red-gamma-min RED_GAMMA_MIN] [--red-prior-log] [--red-ncoeff RED_NCOEFF] [--red-minf RED_MINF] [--tspan-mult TSPAN_MULT] [--qp]
                      [--qp-ratio-max QP_RATIO_MAX] [--qp-sigma-max QP_SIGMA_MAX] [--qp-p-min-np QP_P_MIN_NP] [--qp-p-min QP_P_MIN] [--qp-p-max QP_P_MAX]
                      [--plot-Ared-at-T PLOT_ARED_AT_T] [--double-powerlaw] [--double-powerlaw-delta-gamma DOUBLE_POWERLAW_DELTA_GAMMA] [--A2red-max A2RED_MAX]
                      [--A2red-min A2RED_MIN] [--red-nlog-below-tspan RED_NLOG_BELOW_TSPAN] [--red-nlog-factor RED_NLOG_FACTOR] [--min-red-sin-scale MIN_RED_SIN_SCALE]
                      [--no-white] [--jbo] [--be-flag BE_FLAG] [--white-prior-log] [--efac-max EFAC_MAX] [--efac-min EFAC_MIN] [--equad-max EQUAD_MAX]
                      [--equad-min EQUAD_MIN] [--ngecorr] [--ecorr] [--secorr] [--ecorr-max ECORR_MAX] [--ecorr-min ECORR_MIN] [-D] [--Adm-max ADM_MAX]
                      [--Adm-min ADM_MIN] [--dm-gamma-max DM_GAMMA_MAX] [--dm-gamma-min DM_GAMMA_MIN] [--dm-ncoeff DM_NCOEFF] [--dm-prior-log]
                      [--dm-tspan-mult DM_TSPAN_MULT] [--dm1 DM1] [--dm2 DM2] [--f2 F2] [--pm] [--px] [--px-range PX_RANGE] [--pm-angle] [--pm-range PM_RANGE]
                      [--pm-ecliptic] [--pos] [--pos-range POS_RANGE] [--legendre] [--leg-df0 LEG_DF0] [--leg-df1 LEG_DF1] [--wrap WRAP [WRAP ...]]
                      [--wrap-range WRAP_RANGE] [--tm-fit-file TM_FIT_FILE] [--fake-tm] [--qr] [--svd] [--glitch-all] [--glitches GLITCHES [GLITCHES ...]]
                      [--glitch-recovery GLITCH_RECOVERY [GLITCH_RECOVERY ...]] [--glitch-double-recovery GLITCH_DOUBLE_RECOVERY [GLITCH_DOUBLE_RECOVERY ...]]
                      [--glitch-triple-recovery GLITCH_TRIPLE_RECOVERY [GLITCH_TRIPLE_RECOVERY ...]] [--glitch-epoch-range GLITCH_EPOCH_RANGE]
                      [--glitch-td-min GLITCH_TD_MIN] [--glitch-td-max GLITCH_TD_MAX] [--glitch-f0-range GLITCH_F0_RANGE] [--glitch-f1-range GLITCH_F1_RANGE]
                      [--glitch-f2-range GLITCH_F2_RANGE] [--glitch-f0d-range GLITCH_F0D_RANGE] [--glf0d-prior-log] [--glitch-td-range GLITCH_TD_RANGE]
                      [--gltd-prior-linear] [--glitch-f0d-positive] [--glitch-td-split GLITCH_TD_SPLIT [GLITCH_TD_SPLIT ...]] [--glitch-alt-f0]
                      [--glitch-alt-f0t GLITCH_ALT_F0T [GLITCH_ALT_F0T ...]] [--alt-f0t-gltd] [--measured-prior] [--measured-without]
                      [--measured-sigma MEASURED_SIGMA [MEASURED_SIGMA ...]] [--auto-add] [--fit-planets] [--planets PLANETS] [--mass-max MASS_MAX [MASS_MAX ...]]
                      [--mass-min MASS_MIN [MASS_MIN ...]] [--mass-log-prior] [--period-max PERIOD_MAX [PERIOD_MAX ...]] [--period-min PERIOD_MIN [PERIOD_MIN ...]]
                      [--ecc-log-prior] [-sw] [--sw-sigma-max SW_SIGMA_MAX] [--sw-sigma-min SW_SIGMA_MIN] [--chrom] [--Achrom-max ACHROM_MAX] [--Achrom-min ACHROM_MIN]
                      [--chrom-index CHROM_INDEX] [--chrom-index-max CHROM_INDEX_MAX] [--chrom-index-min CHROM_INDEX_MIN] [--chrom-gamma-max CHROM_GAMMA_MAX]
                      [--chrom-gamma-min CHROM_GAMMA_MIN] [--chrom-ncoeff CHROM_NCOEFF] [--chrom-prior-log] [--chrom-tspan-mult CHROM_TSPAN_MULT] [--chrom-yr]
                      [--chrom-yr-amp-min CHROM_YR_AMP_MIN] [--chrom-yr-amp-max CHROM_YR_AMP_MAX] [--chrom-yr-idx-min CHROM_YR_IDX_MIN]
                      [--chrom-yr-idx-max CHROM_YR_IDX_MAX] [--chrom-bump-start CHROM_BUMP_START [CHROM_BUMP_START ...]]
                      [--chrom-bump-end CHROM_BUMP_END [CHROM_BUMP_END ...]] [--chrom-bump-index-min CHROM_BUMP_INDEX_MIN [CHROM_BUMP_INDEX_MIN ...]]
                      [--chrom-bump-index-max CHROM_BUMP_INDEX_MAX [CHROM_BUMP_INDEX_MAX ...]] [--chrom-bump-sigma-min CHROM_BUMP_SIGMA_MIN [CHROM_BUMP_SIGMA_MIN ...]]
                      [--chrom-bump-sigma-max CHROM_BUMP_SIGMA_MAX [CHROM_BUMP_SIGMA_MAX ...]] [--chrom-bump-amp-min CHROM_BUMP_AMP_MIN [CHROM_BUMP_AMP_MIN ...]]
                      [--chrom-bump-amp-max CHROM_BUMP_AMP_MAX [CHROM_BUMP_AMP_MAX ...]] [--nudot-epochs NUDOT_EPOCHS [NUDOT_EPOCHS ...]]
                      [--nudot-epoch-range NUDOT_EPOCH_RANGE [NUDOT_EPOCH_RANGE ...]] [--nudot-amplitude-range NUDOT_AMPLITUDE_RANGE [NUDOT_AMPLITUDE_RANGE ...]] [--cont]
                      [--nthread NTHREAD] [-N NSAMPLE] [-n] [--nlive NLIVE] [--nwalkers NWALKERS] [--emcee] [--dynesty] [--dynesty-dynamic] [--dynesty-plots]
                      [--dynesty-bound-eff DYNESTY_BOUND_EFF] [--dynesty-bound DYNESTY_BOUND] [--dynesty-sampler DYNESTY_SAMPLER] [--dynesty-bootstrap DYNESTY_BOOTSTRAP]
                      [--zeus] [--ptmcmc] [--multinest] [--multinest-prefix MULTINEST_PREFIX] [--multinest-ce-mode] [--multinest-disable-is]
                      [--multinest-eff MULTINEST_EFF] [--white-corner] [--all-corner] [--plot-chain] [--plot-derived] [--burn BURN] [--truth-file TRUTH_FILE]
                      [--dump-samples DUMP_SAMPLES] [--dump-samples-npy DUMP_SAMPLES_NPY] [--dump-chain] [--plotname PLOTNAME] [--smooth SMOOTH] [--skip-plots]
                      [--plotlibrary PLOTLIBRARY]
                      par tim

Run 'enterprise' on a single pulsar

positional arguments:
  par
  tim

optional arguments:
  -h, --help            show this help message and exit
  --outdir OUTDIR, -o OUTDIR
                        Output directory for chains etc

RedNoise:
  RedNoise Parameters. Note that the red noise model is default ENABLED

  --no-red-noise        Disable Power Law Red Noise search
  --Ared-max ARED_MAX, -A ARED_MAX
                        Max log10A_Red
  --Ared-min ARED_MIN   Min log10A_Red
  --red-at RED_AT       log10A measured at this timescale (in years)
  --red-gamma-max RED_GAMMA_MAX
                        Max gamma red
  --red-gamma-min RED_GAMMA_MIN
                        Min gamma red
  --red-prior-log       Use uniform prior in log space for red noise amplitude
  --red-ncoeff RED_NCOEFF
                        Number of red noise coefficients (nC)
  --red-minf RED_MINF   Low frequency cut-off
  --tspan-mult TSPAN_MULT
                        Multiplier for tspan
  --qp                  Use QP nudot model
  --qp-ratio-max QP_RATIO_MAX
  --qp-sigma-max QP_SIGMA_MAX
  --qp-p-min-np QP_P_MIN_NP
  --qp-p-min QP_P_MIN
  --qp-p-max QP_P_MAX
  --plot-Ared-at-T PLOT_ARED_AT_T
                        Plot derived Ared at freq 1/(T yr)
  --double-powerlaw     Use double powerlaw model
  --double-powerlaw-delta-gamma DOUBLE_POWERLAW_DELTA_GAMMA, --delta-gamma DOUBLE_POWERLAW_DELTA_GAMMA
  --A2red-max A2RED_MAX
                        Max log10A_Red for double pl
  --A2red-min A2RED_MIN
                        Min log10A_Red for double pl
  --red-nlog-below-tspan RED_NLOG_BELOW_TSPAN
                        Number of log-spaced red noise components below tspan
  --red-nlog-factor RED_NLOG_FACTOR
                        Factor for log-spaced red noise components below tspan
  --min-red-sin-scale MIN_RED_SIN_SCALE
                        Limit red noise to reach a minimum to avoid numerical issues

WhiteNoise:
  WhiteNoise Parameters. White noise defaults to Enabled.

  --no-white            Disable efac and equad
  --jbo, -j             Use -be flag for splitting backends
  --be-flag BE_FLAG, -f BE_FLAG
                        Use specified flag for splitting backends
  --white-prior-log     Use uniform prior in log space for Equad
  --efac-max EFAC_MAX   Max for efac prior
  --efac-min EFAC_MIN   Min for efac prior
  --equad-max EQUAD_MAX
                        Max for equad prior (default based on median error)
  --equad-min EQUAD_MIN
                        Min for equad prior (default based on median error)
  --ngecorr             Add ECORR for the nanograv backends
  --ecorr               Add ECORR for all
  --secorr              Add SECORR for all
  --ecorr-max ECORR_MAX
                        Max for ecorr prior
  --ecorr-min ECORR_MIN
                        Min for ecorr prior

DM Variations:
  Parameters for fitting dm variations

  -D, --dm              Enable DM variation search
  --Adm-max ADM_MAX, --dm-amp-max ADM_MAX
                        Max log10A_DM
  --Adm-min ADM_MIN, --dm-amp-min ADM_MIN
                        Min log10A_DM
  --dm-gamma-max DM_GAMMA_MAX
                        Max gamma red
  --dm-gamma-min DM_GAMMA_MIN
                        Min gamma red
  --dm-ncoeff DM_NCOEFF
                        Number of DM bins to use
  --dm-prior-log        Use uniform prior in log space for dm noise amplitude
  --dm-tspan-mult DM_TSPAN_MULT
                        Multiplier for tspan for dm
  --dm1 DM1             fit for DM1
  --dm2 DM2             fit for DM2

BasicTimingModel:
  Basic pulsar spin and astrometric parameters.

  --f2 F2               range of f2 to search
  --pm                  Fit for PMRA+PMDEC
  --px                  Fit for parallax
  --px-range PX_RANGE   Max parallax to search
  --pm-angle            Fit for PM + angle rather than by PMRA/PMDEC
  --pm-range PM_RANGE   Search range for proper motion (deg/yr)
  --pm-ecliptic         Generate ecliptic coords for proper motion
  --pos                 Fit for position (linear fit only)
  --pos-range POS_RANGE
                        Search range for position (arcsec)
  --legendre            Fit polynomial using Legendre series
  --leg-df0 LEG_DF0     Max offset in f0 parameters
  --leg-df1 LEG_DF1     Max offset in f1 parameters
  --wrap WRAP [WRAP ...]
                        Fit for missing phase wraps at this epoch
  --wrap-range WRAP_RANGE
                        Max number of wraps missing
  --tm-fit-file TM_FIT_FILE
                        Use setup file for timing model parameters ONLY WORKS WITH Multinest/MPI
  --fake-tm             Use a design matrix to fake the tm fit
  --qr                  use the QR decomposition for marginalisation
  --svd                 use the SVD decomposition for marginalisation

GlitchModel:
  Glitch and recovery parameters.

  --glitch-all, --gl-all
                        fit for all glitches
  --glitches GLITCHES [GLITCHES ...]
                        Select glitches to fit
  --glitch-recovery GLITCH_RECOVERY [GLITCH_RECOVERY ...]
                        fit for glitch recoveries on these glitches
  --glitch-double-recovery GLITCH_DOUBLE_RECOVERY [GLITCH_DOUBLE_RECOVERY ...]
                        fit for a second glitch recovery on these glitches
  --glitch-triple-recovery GLITCH_TRIPLE_RECOVERY [GLITCH_TRIPLE_RECOVERY ...]
                        fit for a third glitch recovery on these glitches
  --glitch-epoch-range GLITCH_EPOCH_RANGE, --glep-range GLITCH_EPOCH_RANGE
                        Window for glitch epoch fitting
  --glitch-td-min GLITCH_TD_MIN, --gltd-min GLITCH_TD_MIN
                        Min log10(td)
  --glitch-td-max GLITCH_TD_MAX, --gltd-max GLITCH_TD_MAX
                        Max log10(td)
  --glitch-f0-range GLITCH_F0_RANGE, --glf0-range GLITCH_F0_RANGE
                        Fractional change in GLF0
  --glitch-f1-range GLITCH_F1_RANGE, --glf1-range GLITCH_F1_RANGE
                        Fractional change in GLF1
  --glitch-f2-range GLITCH_F2_RANGE, --glf2-range GLITCH_F2_RANGE
                        Fractional change in GLF2
  --glitch-f0d-range GLITCH_F0D_RANGE, --glf0d-range GLITCH_F0D_RANGE
                        Fractional range of GLF0D compared to GLF0
  --glf0d-prior-log, --f0d-log
                        Use log prior for glf0d
  --glitch-td-range GLITCH_TD_RANGE, --gltd-range GLITCH_TD_RANGE
                        Fractional change in GLTD
  --gltd-prior-linear, --td-linear
                        Use linear prior on gltd
  --glitch-f0d-positive, --glf0d-positive
                        Allow only positive GLF0D
  --glitch-td-split GLITCH_TD_SPLIT [GLITCH_TD_SPLIT ...]
                        Where to split the td prior for multi exponentials
  --glitch-alt-f0       Use alternative parameterisation of glitches fitting for instantanious change in F0 rather than GLF0
  --glitch-alt-f0t GLITCH_ALT_F0T [GLITCH_ALT_F0T ...]
                        Replace GLF1 with change of spin frequency 'T' days after the glitches respectively
  --alt-f0t-gltd        Replace GLF1 with change of spin frequency 'GLTD' days after the glitch for all glitches
  --measured-prior      Use measured prior range for GLF0(instant) and GLF0(T=taug)
  --measured-without    Measure prior range without glitch for GLF0(instant) and GLF0(T=taug)
  --measured-sigma MEASURED_SIGMA [MEASURED_SIGMA ...], --sigma-range MEASURED_SIGMA [MEASURED_SIGMA ...]
                        Minus/Plus sigma range of GLF0(instant), and Minus/Plus sigma range of GLF0(T=taug) respectively
  --auto-add            Automatic add all existing glitches recoveries in the par file

Planet:
  Planet Orbital parameters.

  --fit-planets, -P     Fit for 1st planet orbit parameters
  --planets PLANETS     Number of planets to fit
  --mass-max MASS_MAX [MASS_MAX ...]
                        Max mass (Earth masses) prior for planets. Each planet fitted needs a value for this.
  --mass-min MASS_MIN [MASS_MIN ...]
                        Min mass (Earth masses) prior for planets. Each planet fitted needs a value for this.
  --mass-log-prior      Use log mass prior for all planets.
  --period-max PERIOD_MAX [PERIOD_MAX ...]
                        Max period (days) prior for planets. Each planet fitted needs a value for this.
  --period-min PERIOD_MIN [PERIOD_MIN ...]
                        Min period (days) prior for planets. Each planet fitted needs a value for this.
  --ecc-log-prior       Use log ecc prior for all planets.

Time Variable Solar wind:
  Parameters for fitting dm variations due to solar wind

  -sw, --solar-wind     Turn on solar wind
  --sw-sigma-max SW_SIGMA_MAX
                        Max solar wind
  --sw-sigma-min SW_SIGMA_MIN
                        Min solar wind

Chromatic (non-dm) Variations:
  Parameters for fitting non-dm chromatic variations

  --chrom               Enable Chrom variation search
  --Achrom-max ACHROM_MAX
                        Max log10A_Chrom
  --Achrom-min ACHROM_MIN
                        Min log10A_Chrom
  --chrom-index CHROM_INDEX
                        Chromatic index (default=4)
  --chrom-index-max CHROM_INDEX_MAX
                        Chromatic index max (fit for it!)
  --chrom-index-min CHROM_INDEX_MIN
                        Chromatic index min
  --chrom-gamma-max CHROM_GAMMA_MAX
                        Max gamma red
  --chrom-gamma-min CHROM_GAMMA_MIN
                        Min gamma red
  --chrom-ncoeff CHROM_NCOEFF
                        Number of Chrom bins to use
  --chrom-prior-log     Use uniform prior in log space for chromatic noise amplitude
  --chrom-tspan-mult CHROM_TSPAN_MULT
                        Multiplier for tspan for chromatic noise
  --chrom-yr            Add annual chromatic sinusoid
  --chrom-yr-amp-min CHROM_YR_AMP_MIN
                        min amp for annual chromatic sinusoid
  --chrom-yr-amp-max CHROM_YR_AMP_MAX
                        max amp for annual chromatic sinusoid
  --chrom-yr-idx-min CHROM_YR_IDX_MIN
                        min index for annual chromatic sinusoid
  --chrom-yr-idx-max CHROM_YR_IDX_MAX
                        max index for annual chromatic sinuosid
  --chrom-bump-start CHROM_BUMP_START [CHROM_BUMP_START ...]
                        Start epochs for chromatic bumps
  --chrom-bump-end CHROM_BUMP_END [CHROM_BUMP_END ...]
                        End epochs for chromatic bumps
  --chrom-bump-index-min CHROM_BUMP_INDEX_MIN [CHROM_BUMP_INDEX_MIN ...]
                        Indexes for chromatic bumps
  --chrom-bump-index-max CHROM_BUMP_INDEX_MAX [CHROM_BUMP_INDEX_MAX ...]
                        Indexes for chromatic bumps
  --chrom-bump-sigma-min CHROM_BUMP_SIGMA_MIN [CHROM_BUMP_SIGMA_MIN ...]
                        Sigma-min for chromatic bumps
  --chrom-bump-sigma-max CHROM_BUMP_SIGMA_MAX [CHROM_BUMP_SIGMA_MAX ...]
                        Sigma-max for chromatic bumps
  --chrom-bump-amp-min CHROM_BUMP_AMP_MIN [CHROM_BUMP_AMP_MIN ...]
                        Amp-min for chromatic bumps
  --chrom-bump-amp-max CHROM_BUMP_AMP_MAX [CHROM_BUMP_AMP_MAX ...]
                        Amp-max for chromatic bumps

NuNudotSwitch:
  Switching of frequency or frequency derivative.

  --nudot-epochs NUDOT_EPOCHS [NUDOT_EPOCHS ...]
                        Epochs for nudot switches
  --nudot-epoch-range NUDOT_EPOCH_RANGE [NUDOT_EPOCH_RANGE ...]
                        Epochs for nudot switches
  --nudot-amplitude-range NUDOT_AMPLITUDE_RANGE [NUDOT_AMPLITUDE_RANGE ...]
                        Prior range for nudot amplitudes

Sampling options:
  --cont                Continue existing run
  --nthread NTHREAD, -t NTHREAD
                        number of threads
  -N NSAMPLE, --nsample NSAMPLE
                        (max) number of samples
  -n, --no-sample       Disable the actual sampling...
  --nlive NLIVE         Number of live points (nested)
  --nwalkers NWALKERS   number of walkers (mcmc)

EmceeSolver:
  Configure for mcmc with EMCEE

  --emcee               Use emcee sampler

DynestySolver:
  Configure for nested sampling with DyNesty

  --dynesty             Use dynesty sampler
  --dynesty-dynamic     Use dynamic nested sampling
  --dynesty-plots       make dynesty run plots
  --dynesty-bound-eff DYNESTY_BOUND_EFF
                        Efficiency to start bounding
  --dynesty-bound DYNESTY_BOUND
                        Bounding method
  --dynesty-sampler DYNESTY_SAMPLER
                        Sampling method
  --dynesty-bootstrap DYNESTY_BOOTSTRAP
                        Bootstrap amount

ZeusSolver:
  Configure for mcmc with Zeus

  --zeus                Use zeus sampler

PTMCMCSolver:
  Configure for MCMC sampling with ptmcmcsampler

  --ptmcmc              Use ptmcmc sampler

MultiNestSolver:
  Configure for nested sampling with pymultinest

  --multinest           Use pymultinest sampler
  --multinest-prefix MULTINEST_PREFIX
                        Prefix for pymultinest runs
  --multinest-ce-mode   enable constant efficiency mode
  --multinest-disable-is
                        disable importance sampling
  --multinest-eff MULTINEST_EFF
                        Multinest sampling efficiency

Output options:
  --white-corner        Make the efac/equad corner plots
  --all-corner, --corner-all
                        Make corner plots with all params
  --plot-chain          Make a plot of the chains/posterior samples
  --plot-derived        Include derived parameters in corner plots etc.
  --burn BURN           Fraction of chain to burn-in (MC only; default=0.25)
  --truth-file TRUTH_FILE
                        Truths values of parameters; maxlike=maximum likelihood, default=None
  --dump-samples DUMP_SAMPLES
                        Dump N samples of final scaled parameters
  --dump-samples-npy DUMP_SAMPLES_NPY
                        Dump N samples of final scaled parameters as a numpy file
  --dump-chain          Dump chain after scaling to phyical parameters
  --plotname PLOTNAME   Set plot output file name stem
  --smooth SMOOTH       Factor for smoothing of corner plots (default=don't smooth)
  --skip-plots          Skip all plotting (for debugging I guess)
  --plotlibrary PLOTLIBRARY, --plot-library PLOTLIBRARY
                        Which plot library to use for cornerplots (chainconsumer, dynesty, corner)
```