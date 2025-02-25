from . import xparameter as parameter
import numpy as np
from enterprise.signals import utils
from enterprise.signals import gp_signals
from enterprise.signals.parameter import function
import enterprise.constants as const

from enterprise.signals import deterministic_signals

import sys

name = "NuNudotSwitch"
argdec = "Switching of frequency or frequency derivative."


def setup_argparse(parser):
    parser.add_argument('--nudot-epochs', default=[], nargs='+', type=float, help='Epochs for nudot switches')
    parser.add_argument('--nudot-epoch-range', default=[], nargs='+', type=float, help='Epochs for nudot switches')
    parser.add_argument('--nudot-amplitude-range', default=[], nargs='+', type=float, help='Prior range for nudot amplitudes')
    parser.add_argument('--nudot-amplitude-centre', default=[], nargs='+', type=float, help='Prior centre for nudot amplitudes')
    parser.add_argument('--nudot-transition-time-range', default=[], nargs='+', type=float, help='Transition time for nudot switches')


def setup_model(args, psr, parfile):
    pass

    for line in parfile:
        e = line.strip().split()
        if len(e) > 1:
            if e[0] == "F0":
                f0 = float(e[1])
    components = []
    for i, epoch in enumerate(args.nudot_epochs):
        num=i+1
        depoch = args.nudot_epoch_range[i]
        amp_range = args.nudot_amplitude_range[i]
        amp_centre=0
        if args.nudot_amplitude_centre:
            amp_centre=args.nudot_amplitude_centre[i]
        epoch_par = parameter.Uniform(epoch - depoch, epoch + depoch, to_par=to_par)(f"NUDOT_EPOCH_{num}")
        amp_par = parameter.Uniform(amp_centre-amp_range, amp_centre + amp_range, to_par=to_par)(f"NUDOT_AMP_{num}")
        if args.nudot_transition_time_range and args.nudot_transition_time_range[i] > 0:
            transition_time=parameter.Uniform(0,args.nudot_transition_time_range[i],to_par=to_par)(f"NUDOT_TIME_{num}")
            components.append(deterministic_signals.Deterministic(nudot_transition(f0=parameter.Constant(f0),
                                                                                   delta_nudot=amp_par,
                                                                                   epoch=epoch_par,transition_time=transition_time),name=f"nudotswitch{num}"))
        else:
            components.append(deterministic_signals.Deterministic(nudot_step(f0=parameter.Constant(f0),
                                                                               delta_nudot=amp_par,
                                                                               epoch=epoch_par),name=f"nudotswitch{num}"))
    if len(components) > 0:
        model = components[0]
        for m in components[1:]:
            model += m
        return model
    else:
        return None


@function
def nudot_step(toas, f0, delta_nudot, epoch):
    t = toas - epoch * 86400.0
    m = t >= 0
    out = np.zeros_like(toas)
    out[m] = 0.5 * delta_nudot * t[m] ** 2
    return -out / f0

@function
def nudot_transition(toas, f0, delta_nudot, epoch, transition_time):
    t = toas - epoch * 86400.0
    t2 = toas - (epoch + transition_time) * 86400.0
    w=transition_time*86400.0

    m = np.logical_and(t >= 0,t<w)
    m2 = t2 >= 0

    out = np.zeros_like(toas)
    
    #\frac{A\left(x-w\right)^{2}}{6w}+\frac{Aw\left(x-w\right)}{2}+\frac{Aw^{2}}{6}\left\{x>w\right\}

    out[m] = (delta_nudot * t[m]**3)/6/w 

    c =  delta_nudot * w**2 / 6

    out[m2] = delta_nudot*(t2[m2]**2)/2 + delta_nudot*w*t2[m2]/2 + c
    return -out / f0

def to_par(self, p, chain):
    return p, chain
