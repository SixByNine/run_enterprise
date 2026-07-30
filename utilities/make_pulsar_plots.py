#!/usr/bin/env python
# Import the libraries we need.
import numpy as np
import matplotlib
import argparse

matplotlib.use('Agg')

import matplotlib.pyplot as plt

from scipy import linalg
import sys
import subprocess


def subtract_legendre_series_fit_on_ref(
    basis_ref,
    x_ref,
    basis_target,
    x_target,
    order,
    domain=None,
    d_target=None,
    dd_target=None,
    d_scale=1.0,
    dd_scale=1.0,
):
    if order is None:
        return

    if order < 0:
        return

    if domain is None:
        x_min = np.min(x_ref)
        x_max = np.max(x_ref)
    else:
        x_min, x_max = domain

    if x_max == x_min:
        return

    x_ref_scaled = 2.0 * (x_ref - x_min) / (x_max - x_min) - 1.0
    x_target_scaled = 2.0 * (x_target - x_min) / (x_max - x_min) - 1.0

    # Build Legendre basis rows P_0..P_order on the reference and target abscissae.
    v_ref = np.zeros((order + 1, len(x_ref)))
    v_target = np.zeros((order + 1, len(x_target)))
    dv_target = None
    ddv_target = None

    need_d = d_target is not None
    need_dd = dd_target is not None
    dxi_dt = 2.0 / (x_max - x_min)
    if need_d:
        dv_target = np.zeros_like(v_target)
    if need_dd:
        ddv_target = np.zeros_like(v_target)

    for n in range(order + 1):
        poly = np.polynomial.legendre.Legendre.basis(n)
        v_ref[n] = poly(x_ref_scaled)
        v_target[n] = poly(x_target_scaled)

        if need_d:
            dv_target[n] = poly.deriv(1)(x_target_scaled) * dxi_dt
        if need_dd:
            ddv_target[n] = poly.deriv(2)(x_target_scaled) * (dxi_dt * dxi_dt)

    # Least-squares fit of polynomial coefficients on basis_ref for every row.
    # Solve V_ref^T * C^T ~= basis_ref^T.
    coeff_t, _, _, _ = np.linalg.lstsq(v_ref.T, basis_ref.T, rcond=None)
    coeff = coeff_t.T

    # Subtract the same fitted polynomial series from both reference and target bases.
    basis_ref -= coeff.dot(v_ref)
    basis_target -= coeff.dot(v_target)

    # dM and ddM represent scaled first/second derivatives of the basis rows.
    if need_d:
        d_target += d_scale * coeff.dot(dv_target)
    if need_dd:
        dd_target += dd_scale * coeff.dot(ddv_target)

def polymodel(x,f0,meta):
    res=np.zeros_like(x)
    print(meta)

    for i in range(len(meta)):
        k=meta[i,0]
        if k==0:
            res += meta[i,1]*x*86400.0/f0
        if k==1:
            res += meta[i,1]*x*x/2
        if k==2:
            res += (meta[i,1]/6.0)*x*x*x/1.0e9
    return res


def polymodel_nu(x,f0,meta):
    res=np.zeros_like(x)
    print(meta)
    d2s=86400.0
    s2d=1.0/86400.0
    for i in range(len(meta)):
        k=meta[i,0]
        if k==0:
            res += meta[i,1]*86400.0/f0
        if k==1:
            res += meta[i,1]*x
        if k==2:
            res += (meta[i,1]/2.0)*x*x/1.0e9
    return -f0*res*s2d

def polymodel_nudot(x,f0,meta):
    res=np.zeros_like(x)
    print(meta)

    d2s=86400.0
    s2d=1.0/86400.0
    for i in range(len(meta)):
        k=meta[i,0]
        if k==0:
            pass
        if k==1:
            res += meta[i,1]
        if k==2:
            res += (meta[i,1])*x/1.0e9
    return -f0*res*s2d*s2d

def nudot_transition(toas, f0, delta_nudot, epoch, transition_time):
    t = toas - epoch * 86400.0
    t2 = toas - (epoch + transition_time) * 86400.0
    w=transition_time*86400.0
    m = np.logical_and(t >= 0,t<w)
    m2 = t2 >= 0
    out = np.zeros_like(toas)
    out[m] = (delta_nudot * t[m]**3)/6/w 
    c =  delta_nudot * w**2 / 6
    out[m2] = delta_nudot*(t2[m2]**2)/2 + delta_nudot*w*t2[m2]/2 + c
    return -out / f0

parser=argparse.ArgumentParser(description="Make plots from tempo2 outputs")
parser.add_argument('par',type=str)
parser.add_argument('tim',type=str)
parser.add_argument('par2',type=str,nargs='?',default=None)
parser.add_argument('--noplots',dest="plots",action='store_false')
parser.add_argument('--noerr',action='store_true')
parser.add_argument('--perturb',action='store_true')
parser.add_argument('--only-low-modes',action='store_true')
parser.add_argument('--deproject-cvm',action='store_true')
parser.add_argument('-s','--start',type=float,default=None)
parser.add_argument('-f','--finish',type=float,default=None)
parser.add_argument('-N',type=int,default=1000)
args=parser.parse_args()

par=args.par
tim=args.tim
par2=args.par2
if par2 is None:
    par2=par

do_cov_method=False
use_P_in_cvm=False


F0=0
start =46500
finish=58200
psrn="??"

glep=np.zeros(100)
glf0=np.zeros(100)
glf1=np.zeros(100)

glf0d=np.zeros(100)
gltd=np.zeros(100)

glf0d2=np.zeros(100)
gltd2=np.zeros(100)

nudot_tran_epoch=np.zeros(100)
nudot_tran_amp=np.zeros(100)
nudot_tran_time=np.zeros(100)

max_glitch=0
max_nudot_tran=0
PB=0

inpar=[]

TNRedAmp=None
TNRedGam=None
TNRedLog=0
TNRedLog_factor=2.0
TNRedC=0
TNRedGam2= None
TNRedAmp2=None
TNRedSubPolyOrd=None

F2=0

parameters_to_zero=[]

with open(par) as f:
    for line in f:
        if line.startswith("TNRed"):
            e=line.split()
            if e[0]=="TNRedAmp":
                TNRedAmp=float(e[1])
            if e[0] == "TNRedGam":
                TNRedGam = float(e[1])
            if e[0] == "TNRedGam2":
                TNRedGam2 = float(e[1])
            if e[0] == "TNRedAmp2":
                TNRedAmp2 = float(e[1])
            if e[0] == "TNRedFLog":
                TNRedLog = int(e[1])
            if e[0] == "TNRedC":
                TNRedC = int(e[1])
            if e[0] == "TNRedFLog_factor":
                TNRedLog_factor = float(e[1])
            if e[0] == "TNRedSubPolyOrd":
                TNRedSubPolyOrd = int(e[1])
        else:
            inpar.append(line)
        line = line.strip()
        e=line.split()
        if e[0] == "PSRJ":
            psrn=e[1]
        if e[0].startswith("GLEP_"):
            i=int(e[0][5:])
            glep[i-1] = float(e[1])
            max_glitch = max(i,max_glitch)
        if e[0].startswith("GLF0_"):
            i=int(e[0][5:])
            glf0[i-1] = float(e[1])
            parameters_to_zero.append(e[0])
        if e[0].startswith("GLTD_"):
            i=int(e[0][5:])
            gltd[i-1] = float(e[1])
        if e[0].startswith("GLF0D_"):
            i=int(e[0][6:])
            glf0d[i-1] = float(e[1])
            parameters_to_zero.append(e[0])
        if e[0].startswith("GLTD2_"):
            i=int(e[0][6:])
            gltd2[i-1] = float(e[1])
        if e[0].startswith("GLF0D2_"):
            i=int(e[0][7:])
            glf0d2[i-1] = float(e[1])
            parameters_to_zero.append(e[0])
        if e[0].startswith("GLF1_"):
            i=int(e[0][5:])
            glf1[i-1] = float(e[1])
            parameters_to_zero.append(e[0])
        if e[0].startswith("NUDOT_EPOCH_"):
            i=int(e[0][12:])
            nudot_tran_epoch[i-1] = float(e[1])
            max_nudot_tran = max(i,max_nudot_tran)
        if e[0].startswith("NUDOT_AMP_"):
            i=int(e[0][10:])
            nudot_tran_amp[i-1] = float(e[1])
        if e[0].startswith("NUDOT_TIME_"):
            i=int(e[0][11:])
            nudot_tran_time[i-1] = float(e[1])
        if e[0] == "F0":
            F0=float(e[1])
        if e[0] == "PB":
            PB=float(e[1])
        if e[0] == "F1":
            F1=float(e[1])
        if e[0] == "F2":
            F2=float(e[1])
            parameters_to_zero.append(e[0])
        if e[0] == "START":
            start=float(e[1])
        if e[0] == "FINISH":
            finish=float(e[1])
        if e[0] == "PEPOCH":
            pepoch=float(e[1])
            
glep=glep[:max_glitch]
nudot_tran_epoch=nudot_tran_epoch[:max_nudot_tran]


cmd=["tempo2","-output","exportres","-f",par2,tim,"-writeres","-nofit",'-npsr','1','-nobs','50000']
if args.par2 is None:
    for p in parameters_to_zero:
        cmd.append("-set")
        cmd.append(p)
        cmd.append("0.0")
print(" ".join(cmd))
subprocess.call(cmd)

rx2,ry2,re2 = np.loadtxt("out.res",usecols=(0,5,6),unpack=True)

subprocess.call(["tempo2","-output","exportres","-f",par,tim,"-writeres",'-npsr','1','-nobs','50000'])
#subprocess.run(["ls","-l"])

lab = np.loadtxt("param.labels",dtype=str).T
beta = np.loadtxt("param.vals")
cvm  = np.loadtxt("cov.matrix")



meta = np.loadtxt("tnred.meta",usecols=(1,))
try:
    subtract_meta = np.loadtxt("tnred_subtract.meta",usecols=(1,2,3))
except FileNotFoundError:
    subtract_meta=np.zeros((1,3))


omega=meta[0]
epoch=meta[1]

rx,ry,re = np.loadtxt("out.res",usecols=(0,5,6),unpack=True)

try:
    t2red,t2dm,t2chrom = np.loadtxt("out.res",usecols=(9,10,11),unpack=True)
except ValueError:
    print("No tempo2 red noise, dm or chromaticity parameters found in out.res")
    t2red = np.zeros_like(rx)
    t2dm = np.zeros_like(rx)
    t2chrom = np.zeros_like(rx)

ry -= t2dm
ry -= t2chrom

ry2-=t2dm
ry2-=t2chrom


dat_t=rx
dat_e=re
NSAMPLES=1000

if args.start is not None:
    start=args.start
if args.finish is not None:
    finish=args.finish
if args.N is not None:
    NSAMPLES = args.N

t=np.linspace(start-30,finish+30,NSAMPLES)
y=np.zeros_like(t)

cosidx=lab[1]=='param_red_cos'
sinidx=lab[1]=='param_red_sin'

dm_cos_idx = lab[1]=='param_red_dm_cos'
dm_sin_idx = lab[1]=='param_red_dm_sin'

all_cos_idx = np.logical_or(cosidx,dm_cos_idx)
all_sin_idx = np.logical_or(sinidx,dm_sin_idx)
        
nwav = np.sum(sinidx)
ndm_wav = np.sum(dm_sin_idx)




beta_all_mod = beta[np.logical_or(all_sin_idx,all_cos_idx)]
cvm_all_mod = cvm[np.logical_or(all_sin_idx,all_cos_idx)][:,np.logical_or(all_sin_idx,all_cos_idx)]



if args.perturb:
    print("Perturbing the model")
    #pert=np.random.normal(np.zeros_like(beta_mod), np.sqrt(np.diag(cvm_mod)))
    try:
        pert=np.random.multivariate_normal(np.zeros_like(beta_all_mod), cvm_all_mod, size=1,check_valid='raise',tol=1e-6)[0]
    except ValueError as e:
        print(f"\n\n\n\n******\n{par} Error in multivariate normal sampling: {e}",file=sys.stderr)
        print("******\n\n\n\n",file=sys.stderr)
        pert=np.random.multivariate_normal(np.zeros_like(beta_all_mod), cvm_all_mod, size=1)[0]
    beta_all_mod += pert
    beta[np.logical_or(all_sin_idx,all_cos_idx)] = beta_all_mod


beta_mod = beta[np.logical_or(sinidx,cosidx)]
cvm_mod = cvm[np.logical_or(sinidx,cosidx)][:,np.logical_or(sinidx,cosidx)]

beta_dm_mod = beta[np.logical_or(dm_sin_idx,dm_cos_idx)]
cvm_dm_mod = cvm[np.logical_or(dm_sin_idx,dm_cos_idx)][:,np.logical_or(dm_sin_idx,dm_cos_idx)]


cvm_mod_orig = np.copy(cvm_mod)

W=np.power(dat_e,2)  # np.zeros_like(dat_e)+np.power(0.5e-3,2)#

frqidx=lab[1]=='param_f'

if args.deproject_cvm:
    CCof = cvm[np.logical_or(sinidx,cosidx),:][:,frqidx]
    CCf  = cvm[frqidx,:][:,frqidx]
    LLf = linalg.cholesky(CCf)
    LLfinv = linalg.inv(LLf)

    A=CCof.dot(LLfinv)
    cvm_mod_unf = cvm_mod_orig - np.dot(A,A.T)

    cvm_mod=np.diag(np.diag(cvm_mod_unf))
    #cvm_mod = cvm_mod_unf



s2d2 = 1.0/(86400.0*86400.0)
s2d = 1.0/86400.0

M = np.zeros((2*nwav,len(y)))
M2 = np.zeros((2*nwav,len(dat_t)))

dM = np.zeros_like(M)
ddM = np.zeros_like(M)

with open("white.par","w") as f:
    f.writelines(inpar)
    f.write("WAVE_OM {}\n".format(omega))
    f.write("WAVEEPOCH {}\n".format(epoch))
    for i in range(min(256,nwav-TNRedLog)): # Currently fitwaves can't handle the log-spaced waves
        f.write("WAVE{}  {}  {}\n".format(i+1,-beta_mod[i],-beta_mod[i+nwav]))




print("set up matricies")

freqs=[]
pwrs=np.power(beta_mod[:nwav],2) + np.power(beta_mod[nwav:],2)
pwrs2=np.zeros_like(pwrs)
if TNRedAmp2 is not None:
    pwrs3=np.zeros_like(pwrs)
else:
    pwrs3=None

hi_modes = omega * np.arange(1,nwav+1-TNRedLog) # linear spaced above Tspan

if TNRedLog > 0:
    low_modes = omega * TNRedLog_factor**-np.arange(1, TNRedLog + 1) # log spaced below Tspan
    modes = np.concatenate((hi_modes,low_modes))
    low_dfs = -np.diff(np.concatenate([low_modes,[0]])) / 2.0 / np.pi / 86400.0  # df in per second
    hi_dfs = np.diff(np.concatenate([low_modes[-1:], hi_modes])) / 2.0 / np.pi / 86400.0  # df in per second
    dfs = np.concatenate([hi_dfs, low_dfs])

else:
    modes = hi_modes
    dfs = np.diff(np.concatenate([[0], modes])) / 2.0 / np.pi / 86400.0  # df in per second





for i, omegai in enumerate(modes):
    freqs.append(365.25*omegai/2.0/np.pi)
    M[i]        = np.sin(omegai * (t-epoch))
    M[i+nwav]   = np.cos(omegai * (t-epoch))

    M2[i]       = np.sin(omegai * (dat_t-epoch))
    M2[i+nwav]  = np.cos(omegai * (dat_t-epoch))

    
    dM[i]      = -F0*omegai*s2d*M[i+nwav]
    dM[i+nwav] = F0*omegai*s2d*M[i]
    
    ddM[i]      = F0*omegai*omegai*s2d2*M[i]
    ddM[i+nwav] = F0*omegai*omegai*s2d2*M[i+nwav]

                         
if TNRedSubPolyOrd is not None:
    print(f"Subtracting Legendre series up to order {TNRedSubPolyOrd} from red-noise basis")
    # Fit on the data sampling (M2/dat_t) and apply the same subtraction to M.
    leg_domain = (np.min(dat_t), np.max(dat_t))
    subtract_legendre_series_fit_on_ref(
        M2,
        dat_t,
        M,
        t,
        TNRedSubPolyOrd,
        domain=leg_domain,
        d_target=dM,
        dd_target=ddM,
        d_scale=F0 * s2d,
        dd_scale=F0 * s2d2,
    )




if do_cov_method:
    dt_D_D = np.zeros((len(dat_t), len(dat_t)))
    dt_D_S = np.zeros((len(dat_t), len(t)))
    dt_S_S = np.zeros((len(t), len(t)))
    C_D_D = np.zeros_like(dt_D_D)
    C_D_S = np.zeros_like(dt_D_S)
    C_S_S = np.zeros_like(dt_S_S)

    d2C_D_S = np.zeros_like(dt_D_S)
    d4C_S_S = np.zeros_like(dt_S_S)

    for i, it in enumerate(dat_t):
        dt_D_D[i] = np.abs(dat_t - it)*86400.0
        dt_D_S[i] = np.abs(t - it)*86400.0

    for i, it in enumerate(t):
        dt_S_S[i] = np.abs(t-it)*86400.0

    print(TNRedAmp, TNRedGam,P0)
    print("Make Cov Matricies")
    for i, om in enumerate(modes):
        omegai = om / 86400.0  # Want omega in per second not per day
        df = dfs[i]  # (omega / 2.0 / np.pi) / 86400.0  # df in per second as well.
        P = P0 * (freqs[i])**-TNRedGam # freqs[] is in per year, so no need to divide

        C_D_D += P*np.cos(omegai*dt_D_D)*df
        C_D_S += P*np.cos(omegai*dt_D_S)*df
        C_S_S += P*np.cos(omegai*dt_S_S)*df

        d2C_D_S += F0 * omegai * omegai * P*np.cos(omegai*dt_D_S)*df
        d4C_S_S += ((F0*omegai * omegai)**2) *P*np.cos(omegai*dt_S_S)*df


        #cvm_mod[i,i] = P*df
        #cvm_mod[i+nwav,i+nwav] = P*df

    Cn_D_D = C_D_D + np.diag(W)
    print("Invert C_D_D")
    Cn_D_D_inv = linalg.inv(Cn_D_D)
    C_S_D = C_D_S.T
    d2C_S_D = d2C_D_S.T

print("Do linear algebra")

freqs=np.array(freqs)

fyr = 1/(365.25*86400.0)
P0 = (10**(2*TNRedAmp) / (12*np.pi**2)) * fyr**-3
P0_2 = (10**(2*TNRedAmp2) / (12*np.pi**2)) * fyr**-3 if TNRedAmp2 is not None else None

if use_P_in_cvm:
    print("Use PSD directly in place of tempo2 cvm")
    cvm_mod = np.zeros_like(cvm_mod)

keep_mask = np.ones_like(beta_mod, dtype=bool)
for i, om in enumerate(modes):
    omegai = om / 86400.0  # Want omega in per second not per day
    df = dfs[i] #(omega / 2.0 / np.pi) / 86400.0  # df in per second as well.
    P = P0 * (freqs[i]) ** -TNRedGam  # freqs[] is in per year, so no need to divide
    pwrs2[i] = P*df
    if TNRedAmp2 is not None:
        P2 = P0_2 * (freqs[i]) ** -TNRedGam2
        pwrs3[i] = P2*df
        smallmode=pwrs2[i] < 0.1*pwrs3[i]
        if (om > omega) and (args.only_low_modes and smallmode):
            print(f"Mode deleted >> {i} {om} {omega} {smallmode}")
            keep_mask[i] = False
            keep_mask[i+nwav] = False
           

for i, om in enumerate(modes):
    omegai = om / 86400.0  # Want omega in per second not per day
    df = dfs[i] #(omega / 2.0 / np.pi) / 86400.0  # df in per second as well.
    P = P0 * (freqs[i]) ** -TNRedGam  # freqs[] is in per year, so no need to divide
    if TNRedAmp2 is not None:
        P2 = P0_2 * (freqs[i]) ** -TNRedGam2
    if use_P_in_cvm:
        cvm_mod[i,i] = P*df
        cvm_mod[i+nwav,i+nwav] = P*df
        if TNRedAmp2 is not None:
            cvm_mod[i,i] += P2*df
            cvm_mod[i+nwav,i+nwav] += P2*df

beta_mod[~keep_mask] = 0.0

beta_mod = beta_mod[keep_mask]
cvm_mod = cvm_mod[np.ix_(keep_mask, keep_mask)]
M = M[keep_mask,:]
M2 = M2[keep_mask,:]
dM = dM[keep_mask,:]
ddM = ddM[keep_mask,:]
orig_freqs= freqs.copy()
freqs = freqs[keep_mask[:nwav]]


orig_nwav = nwav
nwav = len(beta_mod) // 2




pwrs=np.power(beta_mod[:nwav],2) + np.power(beta_mod[nwav:],2)

                         
#print(cvm.shape,cvm_mod.shape)
M = M.T
ddM = ddM.T
M2 = M2.T
dM = dM.T

y=M.dot(beta_mod) - polymodel(t-pepoch,F0,subtract_meta)

y_dat = M2.dot(beta_mod) - polymodel(dat_t-pepoch,F0,subtract_meta)

yd = (ddM.dot(beta_mod) - polymodel_nudot(t-pepoch,F0,subtract_meta))*1e15

if do_cov_method:
    ALT_y = C_S_D.dot(Cn_D_D_inv.dot(ry))
    ALT_y_C = C_S_S - C_S_D.dot(Cn_D_D_inv).dot(C_D_S)
    print(np.diag(dt_S_S)[:10])
    print(np.diag(C_S_S)[:10])


    ALT_y_err = np.sqrt(np.diag(ALT_y_C))

    ALT_dd = 1e15*d2C_S_D.dot(Cn_D_D_inv.dot(ry))
    ALT_dd_C = d4C_S_S - d2C_S_D.dot(Cn_D_D_inv).dot(d2C_D_S)
    ALT_dd_err = 1e15*np.sqrt(np.diag(ALT_dd_C))
    print(ALT_y_err[0],ALT_dd_err[0])

"""
    if (t < mjd[0]){
        // we are before the first jump
        // so our gradient is just the zeroth offset.
        return yoffs[0];
    } else if(t > mjd[N-1]){
        return yoffs[N-1];
    } else{
        // find the pair we are between...
        for (int ioff =0;ioff<N;ioff++){
            if(t >= mjd[ioff] && t < mjd[ioff+1]){
                double x1 = mjd[ioff];
                double x2 = mjd[ioff+1];
                double x = (t-x1)/(x2-x1);
                double y1=yoffs[ioff];
                double y2=yoffs[ioff+1];
                return (y2-y1)*x + y1;
            }
        }
    }
"""

def ifunc(t,mjd,yoffs):
    if t < mjd[0]:
        return yoffs[0]
    elif t > mjd[-1]:
        return yoffs[-1]
    else:
        for i in range(len(mjd)):
            if t > mjd[i] and t < mjd[i+1]:
                x1=mjd[i]
                x2=mjd[i+1]
                x=(t-x1)/(x2-x1)
                y1=yoffs[i]
                y2=yoffs[i+1]
                return (y2-y1)*x + y1


Coo=M.dot(cvm_mod).dot(M.T)



Cf = M2.dot(cvm_mod).dot(M2.T) 
Cof = M.dot(cvm_mod).dot(M2.T)

#Cf = correlation_tools.cov_nearest(Cf,n_fact=1000,method="clipped") 

#print(linalg.eigvals(Cf))
#origCf = np.copy(Cf)
#p_m = np.trace(Cf)/float(len(dat_e))
#d=10000000.0
#b=0.0
#a=d
#while True:
#    try:
#        #print(numpy.linalg.cond(Cf,'fro'))
#        print("T",a)
#        Lf = linalg.cholesky(Cf)
#        
#        break
#    except:
#        a-=2
#        b = np.sqrt(d*d-a*a)
#        Cf = np.power(a/d,2)*origCf + np.power(b/d,2) * p_m*np.eye(len(dat_e))
#        #print(np.power(a/d,2),p_m,Cf[0][0],origCf[0][0])
#        #print(origCf[0][0],np.power(a/d,2)*origCf[0][0],np.power(b/d,2) * p_m)
#print(a)


if args.noerr:
    ey=np.zeros_like(y)
else:
    Cf+= np.diag(W)


    Lf = linalg.cholesky(Cf)
    Lf_inv = linalg.inv(Lf)
    Cf_inv = linalg.inv(Cf)
    #Co = Coo - np.dot(np.power(a/d,2)*Cof,np.dot(Cf_inv,np.power(a/d,2)*Cof.T))
    #Co = Coo - np.dot(Cof,np.dot(Cf_inv,Cof.T))

    A=Cof.dot(Lf_inv)
    Co = Coo - np.dot(A,A.T)





    print("Negative Co: ",np.sum(np.diag(Co)<0))
    #print("a^2/d^2 = ",np.power(a/d,2))


    #Co *= np.power(b/d,2)
    ey = np.sqrt(np.diag(Co))



if args.noerr:
    ed=np.zeros_like(y)
else:
    ddCoo    =ddM.dot(cvm_mod).dot(ddM.T)                         
    ddCof    = ddM.dot(cvm_mod).dot(M2.T)


    A=ddCof.dot(Lf_inv)
    ddCo = ddCoo - np.dot(A,A.T)
    print("Negative ddCo: ",np.sum(np.diag(ddCo)<0))


    ed=np.sqrt(np.diag(ddCo))*1e15

tt=(t-pepoch)*86400.0
yd_model=np.zeros_like(yd)
yd_model += 1e15*(F2*tt + F1)

i=0
for ge in glep:
    print("Glitch {} = {} >> {}".format(i,ge,1e15*glf1[i]))
    yd_model[t>ge] += 1e15*glf1[i]
    if gltd[i] > 0:
        print("Exponential")
        yd_model[t>ge] -= 1e15 * glf0d[i] * np.exp(-(t[t>ge]-glep[i])/gltd[i]) / (gltd[i]*86400.0)
    if gltd2[i] > 0:
        print("2xExponential")
        yd_model[t>ge] -= 1e15 * glf0d2[i] * np.exp(-(t[t>ge]-glep[i])/gltd2[i]) / (gltd2[i]*86400.0)
    i+=1

for i,nep in enumerate(nudot_tran_epoch):
    print(f"Nudot transition {nep} {nudot_tran_amp[i]} {nudot_tran_time[i]}")
    nep2 = nep + nudot_tran_time[i]
    dt=(t-nep)
    m=(t>nep) & (t < nep2)
    yd_model[m] += 1e15 * (dt[m]/nudot_tran_time[i]) * nudot_tran_amp[i]
    yd_model[t >= nep2] += 1e15 * nudot_tran_amp[i]



yd2 = yd + yd_model

if args.noerr:
    ef=np.zeros_like(y)
else:
    dCoo    = dM.dot(cvm_mod).dot(dM.T)                         
    dCof    = dM.dot(cvm_mod).dot(M2.T)        
    A=dCof.dot(Lf_inv)
    dCo = dCoo - np.dot(A,A.T)
    print("Negative dCo: ",np.sum(np.diag(dCo)<0))

    if np.sum(np.diag(ddCo)<0)>0 or np.sum(np.diag(dCo)<0)>0 or np.sum(np.diag(Co)<0)>0:
        print("ERROR: there are some negative variances!")
    ef= np.sqrt(np.diag(dCo))

yf = dM.dot(beta_mod) - polymodel_nu(t-pepoch,F0,subtract_meta)


#yf2 = yf + (0.5*F2*tt*tt + tt*F1 + F0)
yf2 = yf + (0.5*F2*tt*tt)

i=0
for ge in glep:
    gt=(t-ge)*86400.0
    yf2[t>ge] += glf1[i] * gt[t>ge] + glf0[i]
    if gltd[i] > 0:
        yf2[t>ge] += glf0d[i] * np.exp(-(t[t>ge]-glep[i])/gltd[i])

    if gltd2[i] > 0:
        yf2[t>ge] += glf0d2[i] * np.exp(-(t[t>ge]-glep[i])/gltd2[i])
    i+=1
for i,nep in enumerate(nudot_tran_epoch):
    nep2 = nep + nudot_tran_time[i]
    x=(t-nep)*86400.0
    x2=(t-nep2)*86400.0
    w=nudot_tran_time[i]*86400.0
    m=(t>nep) & (t < nep2)
    yf2[m] += nudot_tran_amp[i] * x[m]**2 / 2/w
    yf2[t >= nep2] += nudot_tran_amp[i] * x2[t >= nep2] + nudot_tran_amp[i]*w/2


ry += np.mean(y_dat-ry)


with open("white_ifunc.par","w") as f:
    f.writelines(inpar)
    f.write("SIFUNC 2 0\n")
    for i in range(len(t)):
        f.write("IFUNC{}  {}  {} {}\n".format(i+1,t[i],-y[i],ey[i]))



with open("resid.asc","w") as f:
    for i in range(len(dat_t)):
        f.write("{} {} {} {} {}\n".format(dat_t[i],ry[i],ry2[i],y_dat[i],dat_e[i]))


with open("nudot.asc","w") as f:
    for i in range(len(yd)):
        f.write("{} {} {} {} {}\n".format(t[i],yd[i],yd_model[i],yd2[i],ed[i]))

with open("deltanu.asc","w") as f:
    for i in range(len(yf2)):
        f.write("{} {} {}\n".format(t[i],yf2[i],y[i]))


if args.plots:


    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(311)

    plt.errorbar(dat_t,ry2,yerr=re2,color='k',marker='.',ls='None',ms=3.0,alpha=0.7)

    plt.title("PSR "+psrn)
    plt.xlabel("MJD")
    plt.ylabel("residual (s)")

    for ge in glep:
        plt.axvline(ge,linestyle="--",color='purple',alpha=0.7)

    ax2 = ax.twinx()
    ax2.set_ybound(np.array(ax.get_ybound())*F0)
    ax2.set_ylabel("Residual (turns)")


    ax = fig.add_subplot(312)
    plt.plot(t,y,color='green')
    plt.errorbar(dat_t,ry,yerr=dat_e,color='k',marker='.',ls='None',ms=3.0,alpha=0.7)

    if do_cov_method:
        plt.plot(t,ALT_y,color='red')
        plt.fill_between(t,ALT_y-ALT_y_err,ALT_y+ALT_y_err,color='red',alpha=0.5)


    plt.fill_between(t,y-ey,y+ey,color='green',alpha=0.5)
    #plt.title("PSR "+psrn)
    plt.xlabel("MJD")
    plt.ylabel("residual (s)")

    for ge in glep:
        plt.axvline(ge,linestyle="--",color='purple',alpha=0.7)

    ax2 = ax.twinx()
    ax2.set_ybound(np.array(ax.get_ybound())*F0)
    ax2.set_ylabel("Residual (turns)")

    ax = fig.add_subplot(313)
    plt.plot(t,y-y,color='green')
    plt.errorbar(dat_t,ry-y_dat,yerr=dat_e,color='k',marker='.',ls='None',ms=3.0,alpha=0.7)

    plt.fill_between(t,-ey,+ey,color='green',alpha=0.5)
    if do_cov_method:
        plt.fill_between(t, -ALT_y_err, ALT_y_err, color='red', alpha=0.5)

    #plt.title("PSR "+psrn)
    plt.xlabel("MJD")
    plt.ylabel("Residual - Model (s)")

    for ge in glep:
        plt.axvline(ge,linestyle="--",color='purple',alpha=0.7)
    ax2 = ax.twinx()
    ax2.set_ybound(np.array(ax.get_ybound())*F0)
    ax2.set_ylabel("Residual - Model (turns)")


    plt.savefig("residuals_{}.pdf".format(psrn))

    plt.figure(figsize=(16,9))
    plt.plot(t,yd,color='blue')
    plt.fill_between(t,yd-ed,yd+ed,color='blue',alpha=0.5)
    if do_cov_method:
        plt.plot(t,ALT_dd,color='brown')
        plt.fill_between(t,ALT_dd-ALT_dd_err,ALT_dd+ALT_dd_err,color='brown',alpha=0.5)
    for ge in glep:
        plt.axvline(ge,linestyle="--",color='purple',alpha=0.7)

    plt.title("PSR "+psrn)
    plt.xlabel("MJD")
    plt.ylabel("$\\dot{\\nu}$ ($10^{-15}$ Hz$^2$)")

    plt.savefig("nudot_{}.pdf".format(psrn))




    fig = plt.figure(figsize=(16,9))
    fig.suptitle("PSR "+psrn)

    ax = fig.add_subplot(321)

    plt.errorbar(dat_t,ry2,yerr=re2,color='k',marker='.',ls='None',ms=3.0,alpha=0.7)

    #plt.xlabel("MJD")
    plt.ylabel("residual (s)")

    for ge in glep:
        plt.axvline(ge,linestyle="--",color='purple',alpha=0.7)

        
    ax3 = ax.twiny()
    ax3.set_xbound((np.array(ax.get_xbound())-53005)/365.25+2004)
    ax3.xaxis.tick_top()
    ax3.set_xlabel("Year")
    ax3.xaxis.set_tick_params(direction='inout',labeltop=True)
    ax.xaxis.set_tick_params(labelbottom=False,direction='in')



    ax2 = ax.twinx()
    ax2.set_ybound(np.array(ax.get_ybound())*F0)
    ax2.set_ylabel("Residual (turns)")



    ax = fig.add_subplot(323)
    plt.plot(t,y,color='green')
    plt.errorbar(dat_t,ry,yerr=dat_e,color='k',marker='.',ls='None',ms=3.0,alpha=0.7)

    plt.fill_between(t,y-ey,y+ey,color='green',alpha=0.5)
    #plt.title("PSR "+psrn)
    #plt.xlabel("MJD")
    plt.ylabel("residual (s)")

    for ge in glep:
        plt.axvline(ge,linestyle="--",color='purple',alpha=0.7)

        
    #ax3 = ax.twiny()
    #ax3.set_xbound((np.array(ax.get_xbound())-53005)/365.25+2004)
    #ax3.xaxis.tick_top()
    #ax3.xaxis.set_tick_params(direction='in',labeltop=False)
    ax.xaxis.set_tick_params(labelbottom=False,direction='in')
        
    ax2 = ax.twinx()
    ax2.set_ybound(np.array(ax.get_ybound())*F0)
    ax2.set_ylabel("Residual (turns)")


    ax = fig.add_subplot(325)
    plt.plot(t,y-y,color='green')
    plt.errorbar(dat_t,ry-y_dat,yerr=dat_e,color='k',marker='.',ls='None',ms=3.0,alpha=0.7)

    plt.fill_between(t,-ey,+ey,color='green',alpha=0.5)
    #plt.title("PSR "+psrn)
    plt.xlabel("MJD")
    plt.ylabel("Residual - Model (s)")

    for ge in glep:
        plt.axvline(ge,linestyle="--",color='purple',alpha=0.7)

    #ax3 = ax.twiny()
    #ax3.set_xbound((np.array(ax.get_xbound())-53005)/365.25+2004)
    #ax3.xaxis.tick_top()
    #ax3.xaxis.set_tick_params(direction='in',labeltop=False)
    ax.xaxis.set_tick_params(labelbottom=True,direction='inout')    
        
    ax2 = ax.twinx()
    ax2.set_ybound(np.array(ax.get_ybound())*F0)
    ax2.set_ylabel("Residual - Model (turns)")





    ax = fig.add_subplot(322)
    plt.plot(t,1e6*yf2,color='orange')
    plt.fill_between(t,1e6*(yf2-ef),1e6*(yf2+ef),color='orange',alpha=0.5)
    for ge in glep:
        plt.axvline(ge,linestyle="--",color='purple',alpha=0.7)

    plt.xlabel("MJD")
    plt.ylabel("$\\Delta{\\nu}$ ($\\mathrm{\\mu}$Hz)")
    ax.yaxis.set_label_position("right")
    ax.yaxis.set_tick_params(labelleft=False,labelright=True,right=True,left=False, direction='in')
    ax3 = ax.twiny()
    ax3.set_xbound((np.array(ax.get_xbound())-53005)/365.25+2004)
    ax3.xaxis.tick_top()
    ax3.set_xlabel("Year")


    ax = fig.add_subplot(324)
    plt.plot(t,yd_model,color='lightblue',ls='--')
    plt.plot(t,yd2,color='blue')
    plt.fill_between(t,yd2-ed,yd2+ed,color='blue',alpha=0.5)

    for ge in glep:
        plt.axvline(ge,linestyle="--",color='purple',alpha=0.7)

    plt.xlabel("MJD")
    plt.ylabel("$\\dot{\\nu}$ ($10^{-15}$ Hz$^2$)")
    ax.yaxis.set_label_position("right")
    ax.yaxis.set_tick_params(labelleft=False,labelright=True,right=True,left=False, direction='in')


    ax = fig.add_subplot(326)
    plt.plot(t,yd,color='blue')
    plt.fill_between(t,yd-ed,yd+ed,color='blue',alpha=0.5)
    for ge in glep:
        plt.axvline(ge,linestyle="--",color='purple',alpha=0.7)

    plt.xlabel("MJD")
    plt.ylabel("$\\Delta\\dot{\\nu}$ ($10^{-15}$ Hz$^2$)")
    ax.yaxis.set_label_position("right")
    ax.yaxis.set_tick_params(labelleft=False,labelright=True,right=True,left=False, direction='in')


    plt.subplots_adjust(hspace=0,wspace=0.15)
    plt.figtext(x=0.47,y=0.94,s="$P$:{:.0f} ms".format(1000.0/F0),horizontalalignment='center')
    plt.figtext(x=0.53,y=0.94,s="$\\dot{{P}}$:{:.0g}".format(-F1/F0/F0),horizontalalignment='center')
    if PB > 0:
        plt.figtext(x=0.59,y=0.94,s="$P_B$:{:.1g}".format(PB),horizontalalignment='center')

    plt.savefig("combined_{}.pdf".format(psrn))


    plt.figure()
    plt.figtext(x=0.47,y=0.94,s="$P$:{:.0f} ms".format(1000.0/F0),horizontalalignment='center')
    plt.figtext(x=0.53,y=0.94,s="$\\dot{{P}}$:{:.0g}".format(-F1/F0/F0),horizontalalignment='center')
    if PB > 0:
        plt.figtext(x=0.59,y=0.94,s="$P_B$:{:.1g}".format(PB),horizontalalignment='center')


    idx = np.argsort(freqs)
    plt.loglog(freqs[idx],pwrs[idx])
    idx = np.argsort(orig_freqs)
    plt.loglog(orig_freqs[idx],pwrs2[idx])
    if pwrs3 is not None:
        plt.loglog(orig_freqs[idx],pwrs3[idx],ls='--')
    plt.title("PSR "+psrn)
    plt.xlabel("Freq (yr^-1)")
    plt.ylabel("Power")
    plt.savefig("pwrspec_{}.pdf".format(psrn))

