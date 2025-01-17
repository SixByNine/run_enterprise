import numpy as np
import scipy.linalg as sl

import matplotlib.pyplot as plt


def setup_argparse(parser):
    parser.add_argument('--time-realisations', default=0, help='Make plots of time domain realisations', type=int)



def get_b(d, TNT, phiinv):
    # Taken from la_forge: https://github.com/nanograv/la_forge/tree/main/la_forge
    
    Sigma = TNT + (np.diag(phiinv) if phiinv.ndim == 1 else phiinv)
    try:
        u, s, _ = sl.svd(Sigma)
        mn = np.dot(u, np.dot(u.T, d)/s)
        Li = u * np.sqrt(1/s)
    except np.linalg.LinAlgError:
        Q, R = sl.qr(Sigma)
        Sigi = sl.solve(R, Q.T)
        mn = np.dot(Sigi, d)
        u, s, _ = sl.svd(Sigi)
        Li = u * np.sqrt(1/s)

    return mn + np.dot(Li, np.random.randn(Li.shape[0]))



def get_tdelay_from_chains(pta, psr, ch, pars, 
                           ipsr=0, 
                           nsamples=100, ch_idxs=None, 
                           plttypes=['achrom'],
                           signames=['all'], 
                           separe_signals=True,
                           verbose=True,
                           weights=None):    
    # Code replicated from ETPA gitlab

    if ch_idxs is None:
        if weights is None:
            ch_idxs = np.random.choice(len(ch),size=nsamples) 
        else:
            prob = weights/np.sum(weights)
            ch_idxs = np.random.choice(len(ch),size=nsamples, p=prob)

    if list(pars) != pta.param_names:
        print("WARNING: parameter names don't correspond with the ones in the PTA object created from strmodel.")
        print("\npars:")
        print("\n".join(pars))
        print("\nPTA param names:")
        print("\n".join(pta.param_names))
        
    

    sig_idxs = {}
    for s, idx in pta._signalcollections[ipsr]._idx.items():
        sig_idxs.update({s.signal_id:idx})
    
    if separe_signals:
        delays = np.zeros((len(signames), len(ch_idxs), len(psr.toas)))
    else:
        delays = np.zeros((len(ch_idxs), len(psr.toas)))

    pta_signals = pta._signalcollections[ipsr]._signals

    for i, ch_idx in enumerate(ch_idxs):
        if verbose:
            print("%02d / %02d"%(i+1, len(ch_idxs)), end="\r")
            
        #post_sample = {parn:val for parn, val in zip(pta.param_names, ch[ch_idx,:])}
        post_sample = pta.map_params(ch[ch_idx,:])

        # Get matrices, basis functions and indexes for the GP signals
        TNrs = pta.get_TNr(post_sample) # [len(basis functions)] * len(psrs)
        TNTs = pta.get_TNT(post_sample) # [len(basis function) x len(basis function)] * len(psrs)
        phiinvs = pta.get_phiinv(post_sample, logdet=False)  # [len(basis functions)] * len(psrs)
        Ts = pta.get_basis(post_sample) # [len(ToAs) x len(basis function)] * len(psrs)
        w = get_b(TNrs[ipsr], TNTs[ipsr], phiinvs[ipsr]) # basis functions
        ci=4
        if "Chrom_index" in pars:
            ci = post_sample["Chrom_index"]
        for j, signame in enumerate(signames):
            delay = 0

            if signame=='all':
                # Find all signals in pta object
                pta_sig = [s for s in pta_signals]
            else:
                # Find related signal in pta object
                #pta_sig = [s for s in pta_signals if signame in s.signal_id]

                # Strict choice
                pta_sig = [s for s in pta_signals if signame == s.signal_id]
            
            for s in pta_sig:
                # Add deterministic signals
                if s.signal_type == 'deterministic':
                    delay += s.get_delay(post_sample)
                    
                # Add GP signals
                if s.signal_type == 'basis':
                    idx = sig_idxs[s.signal_id]
                    delay += np.dot(Ts[ipsr][:, idx], w[idx])

            if plttypes[j]=="dm":
                delay *= psr.freqs**2*2.41e-4
            elif plttypes[j]=="sv":
                delay *= (psr.freqs/1400)**4
            elif plttypes[j]=="cm":
                delay *= (psr.freqs/1400)**ci
            
            if separe_signals:
                delays[j, i, :] = delay
            else:
                delays[i, :] += delay
    
    return delays


def plot_time_domain(ax, delays, psr, label=None,
                    plot_meds=True, color_meds='k',
                    plot_surf=True, color_surf=None,
                    plot_errs=False, color_errs=None, 
                    plot_real=False, color_real=None, ch_idx_real=None, 
                    plttype="achrom", rem_zeros=False,remove_quadratic=False):
    
    t = psr.toas / 86400

    if remove_quadratic:
        for i in range(len(delays)):
            delays[i] -= np.polyval(np.polyfit(t, delays[i], 2), t)

    maxs = np.quantile(delays, 0.16, axis=0)
    mins = np.quantile(delays, .84, axis=0)
    meds = np.median(delays, axis=0)    
    errl = meds - mins
    errh = maxs - meds

    i_label = False
    if plot_meds:
        label = label if not i_label else None
        if rem_zeros: 
            mask = np.where(np.abs(meds)>1e-15)
        else:
            mask = np.arange(len(meds))
        ax.plot(t[mask], meds[mask], '.', color=color_meds, ms=2,  rasterized=True, zorder=3, label=label)
        if label: i_label = True

    if plot_surf:
        label = label if not i_label else None
        if rem_zeros: 
            mask = np.where( (np.abs(mins)<1e-15) and (np.abs(maxs)<1e-15) )
        else:
            mask = np.arange(len(meds))
        ax.fill_between(t[mask], mins[mask], maxs[mask], where=maxs[mask]<=mins[mask], interpolate=True, alpha=.3, label=label,
                        color=color_surf, rasterized=True, zorder=0)
        if label: i_label = True
    
    if plot_errs:
        label = label if not i_label else None
        if rem_zeros: 
            mask = np.where(np.abs(meds)>1e-15)
        else:
            mask = np.arange(len(meds))
        ax.errorbar(t[mask], meds[mask], yerr=[errl[mask], errh[mask]], fmt='none', alpha=.4, capsize=3, label=label,
                    color=color_errs,  rasterized=True, zorder=1)
        if label: i_label = True
        
    if plot_real:
        if ch_idx_real is not None:
            delays = delays[ch_idx_real]
        for d in delays:
            label = label if not i_label else None
            mask = np.where(np.abs(d)>1e-15)
            ax.plot(t[mask], d[mask], alpha=.2, lw=.6, color=color_real, rasterized=True, zorder=2, label=label)
            if label: i_label = True


    ax.set_xlabel('MJD', fontsize=15)
    if plttype == "dm":
        ax.set_ylabel(r'$\Delta$ DM [$pc/cm^3$]', fontsize=15)
    elif plttype == "sv":
        ax.set_ylabel(r'Time delay (x $\mathrm{freqs}^4)$', fontsize=10)
    elif plttype == "cm":
        ax.set_ylabel(r'Time delay at 1400 MHz (s)', fontsize=10)
    else:
        ax.set_ylabel('Time delay [s]', fontsize=15)
        
    if i_label:
        ax.legend()
    ax.ticklabel_format(style='sci',scilimits=(0,0),axis='y')
    ax.axhline(0, ls=":")
    ax.grid(alpha=.2)


def create_td_realisations(args,psr_model, raw_results):
    samples, _, weights = raw_results
    pta = psr_model.get_pta()
    psr = psr_model.psr

    if args.plotname is None:
        plotname=psr.name
    else:
        plotname=psr.name+'_'+args.plotname
    print(f"Making {args.time_realisations} time-domain realisations...")

    pta_signals = set([s.signal_id for s in pta._signalcollections[0]._signals])
    signal = list(set(('red_noise','dm_noise','chrom_noise')).intersection(pta_signals))
    print("Making realisations for signals: ", signal, "out of ", pta_signals)
    types=[]
    for s in signal:
        if s in ['red_noise']:
            types.append('achrom')
        elif s in ['dm_noise']:
            types.append('dm')
        elif s in ['chrom_noise']:
            types.append('cm')
    
    delays = get_tdelay_from_chains(pta, psr, samples, pta.param_names, signames=signal, plttypes=types, weights=weights, separe_signals=True, nsamples=args.time_realisations)

    print(f"Made {delays.shape[0]}x{delays.shape[1]} timeseries of length {delays.shape[2]}")
    for i, (s,t) in enumerate(zip(signal,types)):
        # get new axis object
        print(f"Make {plotname}_realisations_{s}.pdf")
        fig, ax = plt.subplots(1, 1, figsize=(12, 6))
        ax.set_title(f"{psr.name} - {s}")
        plot_time_domain(ax,delays[i],psr,plot_real=True,color_real='gray',ch_idx_real=np.arange(min(128,args.time_realisations)),remove_quadratic=True,plttype=t)
        fig.savefig(f"{plotname}_realisations_{s}.pdf")
        plt.close(fig)
    
    print(f"Save {plotname}_realisations.npz")

    np.savez(f"{plotname}_realisations.npz", delays=delays,t=psr.toas/86400,signal=signal,types=types)