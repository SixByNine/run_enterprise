import numpy as np
from enterprise.signals import selections, white_signals, signal_base,utils
from . import xparameter as parameter


name = "WhiteNoise"
argdec = "WhiteNoise Parameters. White noise defaults to Enabled."


def setup_argparse(parser):
    parser.add_argument('--no-white', dest='white', default=True, action='store_false', help='Disable efac and equad')
    parser.add_argument('--jbo', '-j', action='store_true', help='Use -be flag for splitting backends')
    parser.add_argument('--be-flag', '-f', help='Use specified flag for splitting backends')
    parser.add_argument('--white-prior-log', action='store_true', help='Use uniform prior in log space for Equad')
    parser.add_argument('--efac-max', type=float, default=5, help='Max for efac prior')
    parser.add_argument('--efac-min', type=float, default=0.2, help='Min for efac prior')

    parser.add_argument('--equad-max', type=float, default=None, help='Max for equad prior (default based on median error)')
    parser.add_argument('--equad-min', type=float, default=None, help='Min for equad prior (default based on median error)')
    parser.add_argument('--ngecorr', action='store_true', help='Add ECORR for the nanograv backends')
    parser.add_argument('--ecorr', action='store_true', help='Add ECORR for all')
    parser.add_argument('--secorr', action='store_true', help='Add SECORR for all')

    parser.add_argument('--ecorr-max', type=float, default=-5, help='Max for ecorr prior')
    parser.add_argument('--ecorr-min', type=float, default=-9, help='Min for ecorr prior')


def setup_model(args, psr, parfile):
    selection = selections.Selection(selections.by_backend)
    selflag = '-f'
    if args.jbo:
        selection = selections.Selection(jbbackends)
        selflag = '-be'
    if args.be_flag:
        selflag = "-" + args.be_flag
        selection = selections.Selection(lambda flags: flagbackends(flags, args.be_flag))

    def to_par(self, p, chain):
        if "efac" in p:
            flag = p.split('_', 1)[1]
            flag = flag[:-5]
            return "TNEF %s %s" % (selflag, flag), chain
        if "equad" in p:
            flag = p.split('_', 1)[1]
            flag = flag[:-14]
            return "TNEQ %s %s" % (selflag, flag), chain
        if "secorr" in p:
            flag=p.split('_',1)[1]
            flag=flag[:-13]
            return "TNSECORR -be %s"%(flag), np.power(10,chain)*1e6 # Convert log-seconds to microseconds
        if "ecorr" in p:
            flag = p.split('_', 1)[1]
            flag = flag[:-12]
            return "TNECORR -be %s" % (flag), np.power(10, chain) * 1e6  # Convert log-seconds to microseconds

        else:
            return None
    s = selection(psr)

    if args.white:
        medianlogerr = np.log10(np.median(psr.toaerrs))
        efac = parameter.Uniform(args.efac_min, args.efac_max, to_par=to_par)
        if args.equad_min is None:
            equad_min = medianlogerr - 2
        else:
            equad_min = args.equad_min
        if args.equad_max is None:
            equad_max = medianlogerr + 2
        else:
            equad_max = args.equad_max

        if (args.white_prior_log):
            equad = parameter.Uniform(equad_min,equad_max, to_par=to_par)
        else:
            equad = parameter.LinearExp(equad_min,equad_max, to_par=to_par)
        ef = white_signals.MeasurementNoise(efac=efac, selection=selection)
        eq = white_signals.TNEquadNoise(log10_tnequad=equad, selection=selection)
        model = ef + eq
    else:
        ef = white_signals.MeasurementNoise(efac=parameter.Constant(1.0), selection=selection)
        model = ef


    if args.ecorr:
        if args.white_prior_log:
            ecorr = parameter.Uniform(args.ecorr_min, args.ecorr_max,to_par=to_par)
        else:
            ecorr = parameter.LinearExp(args.ecorr_min,args.ecorr_max,to_par=to_par)
        ec = white_signals.EcorrKernelNoise(log10_ecorr=ecorr, selection=selection)
        model += ec

    if args.secorr:
        if args.white_prior_log:
            ecorr = parameter.Uniform(args.ecorr_min, args.ecorr_max, to_par=to_par)
        else:
            ecorr = parameter.LinearExp(args.ecorr_min, args.ecorr_max, to_par=to_par)
        ec = SEcorrKernelNoise(log10_secorr=ecorr, selection=selection)
        model += ec

    if args.ngecorr:
        ngselection = selections.Selection(selections.nanograv_backends)
        if args.white_prior_log:
            ecorr = parameter.Uniform(-9,-5,to_par=to_par)
        else:
            ecorr = parameter.LinearExp(-9,-5,to_par=to_par)
        ec = white_signals.EcorrKernelNoise(log10_ecorr=ecorr, selection=ngselection)
        model += ec

    return model



def jbbackends(flags):
    backend_flags = flags['be']
    flagvals = np.unique(backend_flags)
    return {flagval: backend_flags == flagval for flagval in flagvals}


def flagbackends(flags, beflagval):
    backend_flags = flags[beflagval]
    flagvals = np.unique(backend_flags)
    return {flagval: backend_flags == flagval for flagval in flagvals}


def SEcorrKernelNoise(
    log10_secorr=parameter.Uniform(-10, -5),
    selection=selections.Selection(selections.no_selection),
    method="sherman-morrison",
    name="",
):
    import scipy.sparse
    if method not in ["sherman-morrison", "block", "sparse"]:
        msg = "SEcorrKernelNoise does not support method: {}".format(method)
        raise TypeError(msg)

    class SEcorrKernelNoise(signal_base.Signal):
        signal_type = "white noise"
        signal_name = "secorr_" + method
        signal_id = "_".join(["secorr", name, method]) if name else "_".join(["secorr", method])

        def __init__(self, psr):
            super(SEcorrKernelNoise, self).__init__(psr)
            self.name = self.psrname + "_" + self.signal_id


            sel = selection(psr)
            self._params, self._masks = sel("log10_secorr", log10_secorr)
            keys = sorted(self._masks.keys())
            masks = [self._masks[key] for key in keys]

            Umats = []
            for key, mask in zip(keys, masks):
                Umats.append(utils.create_quantization_matrix(psr.toas[mask], nmin=2)[0])

            nepoch = sum(U.shape[1] for U in Umats)
            U = np.zeros((len(psr.toas), nepoch))
            self._slices = {}
            netot = 0
            for ct, (key, mask) in enumerate(zip(keys, masks)):
                nn = Umats[ct].shape[1]
                U[mask, netot: nn + netot] = Umats[ct]
                self._slices.update({key: utils.quant2ind(U[:, netot: nn + netot])})
                netot += nn

            # Get TOBS
            self.secorr_scales=dict()
            for ct, (key,mask)in enumerate(zip(keys, masks)):
                tobs = psr.flags['length'][mask]
                observation_times = np.ones(len(self._slices[key]))
                for islice, slice in enumerate(self._slices[key]):
                    observation_times[islice] = tobs[slice][0]
                self.secorr_scales[key] = observation_times/3600.0 ## secorr scales by sqrt(tobs/1h)... we are looking at cvm, so no squaring


            # initialize sparse matrix
            self._setup(psr)
        @property
        def ndiag_params(self):
            """Get any varying ndiag parameters."""
            return [pp.name for pp in self.params]

        @signal_base.cache_call("ndiag_params")
        def get_ndiag(self, params):
            if method == "sherman-morrison":
                return self._get_ndiag_sherman_morrison(params)
            elif method == "sparse":
                return self._get_ndiag_sparse(params)
            elif method == "block":
                return self._get_ndiag_block(params)

        def _setup(self, psr):
            if method == "sparse":
                self._setup_sparse(psr)

        def _setup_sparse(self, psr):
            Ns = scipy.sparse.csc_matrix((len(psr.toas), len(psr.toas)))
            for key, slices in self._slices.items():
                for slc in slices:
                    if slc.stop - slc.start > 1:
                        Ns[slc, slc] = 1.0
            self._Ns = signal_base.csc_matrix_alt(Ns)

        def _get_ndiag_sparse(self, params):
            for p in self._params:
                for slc in self._slices[p]:
                    if slc.stop - slc.start > 1:
                        self._Ns[slc, slc] = 10 ** (2 * self.get(p, params))
            return self._Ns

        def _get_ndiag_sherman_morrison(self, params):
            slices, jvec = self._get_jvecs(params)
            return signal_base.ShermanMorrison(jvec, slices)

        def _get_ndiag_block(self, params):
            slices, jvec = self._get_jvecs(params)
            blocks = []
            for jv, slc in zip(jvec, slices):
                nb = slc.stop - slc.start
                blocks.append(np.ones((nb, nb)) * jv)
            return signal_base.BlockMatrix(blocks, slices)

        def _get_jvecs(self, params):
            slices = sum([self._slices[key] for key in sorted(self._slices.keys())], [])
            jvec = np.concatenate(
                [
                    self.secorr_scales[key] * 10 ** (2 * self.get(key, params))
                    for key in sorted(self._slices.keys())
                ]
            )
            return (slices, jvec)

    return SEcorrKernelNoise

