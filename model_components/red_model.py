from . import xparameter as parameter
import numpy as np
from enterprise.signals import utils
from enterprise.signals import gp_signals




name="RedNoise"
argdec="RedNoise Parameters. Note that the red noise model is default ENABLED"

def setup_argparse(parser):
   parser.add_argument('--no-red-noise',dest='red',default=True,action='store_false', help='Disable Power Law Red Noise search')
   parser.add_argument('--Ared-max','-A',type=float,default=-12,help='Max log10A_Red')
   parser.add_argument('--Ared-min',type=float,default=-18,help='Min log10A_Red')
   parser.add_argument('--red-gamma-max',type=float,default=10,help='Max gamma red')
   parser.add_argument('--red-gamma-min',type=float,default=0,help='Min gamma red')
   parser.add_argument('--red-prior-log',action='store_true',help='Use uniform prior in log space for red noise amplitude')
   parser.add_argument('--red-ncoeff',type=int,default=60,help='Number of red noise coefficients (nC)')
   parser.add_argument('--tspan-mult', type=float, default=2, help='Multiplier for tspan')


def setup_model(args, psr, parfile):
   if args.red:
      # We scale Tspan if we are using tspan-mult
      nC = args.red_ncoeff
      Tspan = psr.toas.max() - psr.toas.min()


      Tspan *= args.tspan_mult
      nC *= int(args.tspan_mult)

      # Add the nC and flow parameters to the par file
      parfile.append("TNRedC {}\n".format(nC))
      parfile.append("TNRedFLow {}\n".format(np.log10(1.0/args.tspan_mult)))

      # red noise (powerlaw with nC frequencies)
      if args.red_prior_log:
         log10_A = parameter.Uniform(args.Ared_min, args.Ared_max,to_par=to_par)
      else:
         log10_A = parameter.LinearExp(args.Ared_min, args.Ared_max,to_par=to_par)
      gamma = parameter.Uniform(args.red_gamma_min, args.red_gamma_max,to_par=to_par)

      pl = utils.powerlaw(log10_A=log10_A, gamma=gamma)
      rn = gp_signals.FourierBasisGP(spectrum=pl, components=nC, Tspan=Tspan)
      return rn
   else:
      return None

def to_par(self,p,chain):
   if "red_noise_log10_A" in p:
      return "TNRedAmp", chain
   if "red_noise_gamma" in p:
      return "TNRedGam", chain
   else:
      return None
