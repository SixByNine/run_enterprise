import numpy as np
import sys

from enterprise.signals import deterministic_signals
from enterprise.signals import signal_base

from astropy import coordinates as coord
from astropy import units as units

from . import xparameter as parameter

name = "GlitchModel"
argdec = "Glitch and recovery parameters."


def setup_argparse(parser):
    parser.add_argument('--glitch-all','--gl-all', action='store_true', help='fit for all glitches')
    parser.add_argument('--glitches', type=int, nargs='+', help='Select glitches to fit')
    parser.add_argument('--glitch-recovery', type=int, nargs='+', help='fit for glitch recoveries on these glitches')
    parser.add_argument('--glitch-double-recovery', type=int, nargs='+', help='fit for a second glitch recovery on these glitches')
    parser.add_argument('--glitch-epoch-range', '--glep-range', type=float, default=0, help='Window for glitch epoch fitting')
    parser.add_argument('--glitch-td-min','--gltd-min', type=float, default=0, help='Min log(td)')
    parser.add_argument('--glitch-td-max', '--gltd-max', type=float, default=4, help='Max log10(td)')
    parser.add_argument('--glitch-f0-range', '--glf0-range', type=float, default=0.3, help="Fractional change in glF0")
    parser.add_argument('--glitch-f1-range', '--glf1-range', type=float, default=0.3, help="Fractional change in glF1")
    parser.add_argument('--glitch-f0d-range', '--glf0d-range', type=float, default=0.3, help="Fractional range of f0d compared to glF0")
    parser.add_argument('--glitch-f0d-positive', '--glf0d-positive', action='store_true' , help="Allow only positive GLF0D")


def setup_model(args, psr, parfile):
    global glitches
    glitches = {}
    for line in parfile:
        e = line.strip().split()
        if len(e) > 1:
            if e[0] == "F0":
                f0 = float(e[1])
            elif e[0] == "PEPOCH":
                pepoch = float(e[1])
            elif e[0].startswith("GLEP"):
                glitches[int(e[0][5:])] = {'EP': float(e[1])}
            elif e[0].startswith("GLF0D"):
                glitches[int(e[0][6:])]['F0D'] = float(e[1])
            elif e[0].startswith("GLTD"):
                glitches[int(e[0][5:])]['TD'] = float(e[1])
            elif e[0].startswith("GLF0_"):
                glitches[int(e[0][5:])]['F0'] = float(e[1])
            elif e[0].startswith("GLF1_"):
                glitches[int(e[0][5:])]['F1'] = float(e[1])
            elif e[0].startswith("GLF2_"):
                glitches[int(e[0][5:])]['F2'] = float(e[1])

    if args.glitch_recovery:
        for gl in glitches:
            if "GLF0D_{}".format(gl) in psr.fitpars or "GLTD_{}".format(gl) in psr.fitpars:
                print("ERROR: Can't fit for GLF0D/GLTD in both MCMC and least-squares")
                sys.exit(1)
            if "GLF0_{}".format(gl) in psr.fitpars or "GLF1_{}".format(gl) in psr.fitpars:
                print("ERROR: Can't fit for GLF0/GLF1 in both MCMC and least-squares")
                sys.exit(1)

    components = []
    if args.glitches is None:
        args.glitches=[]
    if args.glitch_recovery is None:
        args.glitch_recovery = []

    for gl in glitches:
        if args.glitch_all or gl in args.glitches or gl in args.glitch_recovery:
            ## pre-compute the model for the existing glitch
            glf0=glf1=glf2=glf0d=gltd=0
            glep = glitches[gl]['EP']
            if 'F0' in glitches[gl]:
                glf0 = glitches[gl]['F0']
            if 'F1' in glitches[gl]:
                glf1 = glitches[gl]['F1']
            if 'F2' in glitches[gl]:
                glf2 = glitches[gl]['F2']
            if 'F0D' in glitches[gl]:
                glf0d = glitches[gl]['F0D']
            if 'TD' in glitches[gl]:
                gltd = glitches[gl]['TD']
            # This is the 'old model'
            glitches[gl]['old_model'] = model_glitch(psr.toas, f0, glep, glf0,glf1,glf2,glf0d,gltd)


            if args.glitch_epoch_range > 0:
                glep = parameter.Uniform(glitches[gl]['EP']-args.glitch_epoch_range, glitches[gl]['EP']+args.glitch_epoch_range, to_par=to_par)("GLEP_{}".format(gl))
            else:
                glep=parameter.Constant(glitches[gl]['EP'])

            glf0_range = sorted([glf0 - args.glitch_f0_range*glf0, glf0 + args.glitch_f0_range*glf0])
            glf1_range = sorted([glf1 - args.glitch_f1_range*glf0, glf1 + args.glitch_f1_range*glf1])
            glf0d_range = [0 if args.glitch_f0d_positive else -args.glitch_f0d_range*abs(glf0), args.glitch_f0d_range*abs(glf0)]

            glf0 = parameter.Uniform(glf0_range[0]*1e6, glf0_range[1]*1e6,to_par=to_par)("GLF0_{}".format(gl))  # units of 10^-6
            glf1 = parameter.Uniform(glf1_range[0]*1e12,glf1_range[1]*1e12,to_par=to_par)("GLF1_{}".format(gl))  # units of 10^-12
            if gl in args.glitch_recovery or gl in args.glitch_double_recovery:
                f0d = parameter.Uniform(glf0d_range[0],glf0d_range[1],to_par=to_par)("GLF0D_{}".format(gl))  # units of 10^-6
                td = parameter.Uniform(args.glitch_td_min, args.glitch_td_max,to_par=to_par)("GLTD_{}".format(gl))
                if gl in args.glitch_double_recovery:
                    f0d_2 = parameter.Uniform(glf0d_range[0], glf0d_range[1], to_par=to_par)(
                        "GLF0D2_{}".format(gl))  # units of 10^-6
                    td_2 = parameter.Uniform(args.glitch_td_min, args.glitch_td_max, to_par=to_par)("GLTD2_{}".format(gl))
            else:
                f0d = parameter.Constant(0)
                td = parameter.Constant(0)

            if gl in args.glitch_double_recovery:
                ## At some point we should update this to more than 2 exponentials.
                glexp = model_glitch_residual_multiexp(f0=parameter.Constant(f0), glep=glep,
                                              glf0=glf0, glf1=glf1, glf0d=f0d, log10_td=td,
                                              glitch_id=parameter.Constant(gl),
                                                       glf0d2 = f0d_2, log10_td2=td_2)
            else:
                glexp = model_glitch_residual(f0=parameter.Constant(f0), glep=glep,
                                              glf0=glf0, glf1=glf1, glf0d=f0d, log10_td=td,
                                              glitch_id=parameter.Constant(gl))
            glexp_f = deterministic_signals.Deterministic(glexp, name="gl{}:".format(gl))
            components.append(glexp_f)

    if len(components) > 0:
        model = components[0]
        for m in components[1:]:
            model += m
        return model
    else:
        return None


def to_par(self, p, chain):
    if p.startswith("GL"):
        pname=p
        if "GLF0" in pname:
            return pname, chain * 1e-6
        elif "GLF1" in pname:
            return pname, chain * 1e-12
        elif "GLF2" in pname:
            return pname, chain * 1e-18
        elif "GLTD" in pname:
            return pname, 10**chain
        else:
            return pname,chain
    else:
        return p,chain


def model_glitch(toas, f0, glep, glf0=0, glf1=0, glf2=0, glf0d=0, gltd=0):
    t = toas - glep * 86400.0
    m = t >= 0
    glif = np.zeros_like(toas)
    glif[m] = glf0 * t[m] + 0.5 * glf1 * (t[m]) ** 2 + glf2 * t[m] ** 3 / 6.0
    phs = -glif

    expf = np.zeros_like(toas) + 1.0  # Cancel recovery when t<0
    expf[m] = np.exp(-t[m] / (gltd * 86400.0))  # No +1.0!
    phs += - glf0d * gltd * 86400 * (1.0 - expf)

    return phs / f0


@signal_base.function
def model_glitch_residual(toas, f0, glep, glf0, glf1=0, glf2=0, glf0d=0, log10_td=0, glitch_id=1):
    global glitches
    # remove old glitch
    return model_glitch(toas, f0=f0, glep=glep, glf0=glf0 * 1e-6, glf1=glf1 * 1e-12, glf2=glf2 * 1e-18,  #
                        glf0d=glf0d * 1e-6, gltd=10**log10_td) - glitches[int(glitch_id)]['old_model']

@signal_base.function
def model_glitch_residual_multiexp(toas, f0, glep, glf0, glf1=0, glf2=0, glf0d=0, log10_td=0, glitch_id=1, glf0d2=0, log10_td2=0):
    global glitches
    # remove old glitch
    return model_glitch(toas, f0=f0, glep=glep, glf0=glf0 * 1e-6, glf1=glf1 * 1e-12, glf2=glf2 * 1e-18,  #
                        glf0d=glf0d * 1e-6, gltd=10**log10_td) \
           - model_glitch(toas, f0=f0, glep=glep,glf0d=glf0d2, gltd=log10_td2) \
           - glitches[int(glitch_id)]['old_model']
