
all=[]

try:
    from  . import emcee_solver
    all.append(emcee_solver)
except ImportError:
    print("Cannot import emcee")

try:
    from . import dynesty_solver
    all.append(dynesty_solver)
except ImportError:
    print("Cannot import dynesty")

try:
    from . import zeus_solver
    all.append(zeus_solver)
except ImportError:
    print("Cannot import dynesty")



def std_args(parser):
    parser.add_argument('--cont', action='store_true', help='Continue existing run')
    parser.add_argument('--nthread', '-t', type=int, default=1, help="number of threads")
    parser.add_argument('-N','--nsample', type=float, default=1e6, help='(max) number of samples')
    parser.add_argument('-n','--no-sample',dest='sample',default=True,action='store_false', help='Disable the actual sampling...')
