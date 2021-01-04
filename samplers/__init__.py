from  . import emcee_solver, dynesty_solver

all = [emcee_solver,dynesty_solver]


def std_args(parser):
    parser.add_argument('--cont', action='store_true', help='Continue existing run')
    parser.add_argument('--nthread', '-t', type=int, default=1, help="number of threads")
    parser.add_argument('-N','--nsample', type=float, default=1e6, help='(max) number of samples')
    parser.add_argument('-n','--no-sample',dest='sample',default=True,action='store_false', help='Disable the actual sampling...')
