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


def setup_model(args, psr, parfile):
    pass

    for line in parfile:
        e = line.strip().split()
        if len(e) > 1:
            if e[0] == "F0":
                f0 = float(e[1])
    components = []
    for i, epoch in enumerate(args.nudot_epochs):
        depoch = args.nudot_epoch_range[i]
        amp_range = args.nudot_amplitude_range[i]
        epoch_par = parameter.Uniform(epoch - depoch, epoch + depoch, to_par=to_par)(f"nudot_epoch{i}")
        amp_par = parameter.Uniform(-amp_range, amp_range, to_par=to_par)(f"nudot_amp{i}")
        components.append(deterministic_signals.Deterministic(nudot_transition(f0=parameter.Constant(f0),
                                                                               delta_nudot=amp_par,
                                                                               epoch=epoch_par)))
    if len(components) > 0:
        model = components[0]
        for m in components[1:]:
            model += m
        return model
    else:
        return None


@function
def nudot_transition(toas, f0, delta_nudot, epoch):
    t = toas - epoch * 86400.0
    m = t >= 0
    out = np.zeros_like(toas)
    out[m] = 0.5 * delta_nudot * t[m] ** 2
    return -out / f0


def to_par(self, p, chain):
    return p, chain
