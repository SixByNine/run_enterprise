import numpy as np
from scipy import linalg
from scipy import interpolate as interp
import sys

"""
Power spectral density computation codes. Based on code by Iuliana Nitu 2020.
Similar approach to cholSpec tempo2 plugin.
"""

def qrfit(dm, b):
    """Least Squares fitting using the QR decomposition"""
    q, r = np.linalg.qr(dm, mode='reduced')
    p = np.dot(q.T, b)

    ri = linalg.inv(r)
    param = np.dot(linalg.inv(r), p)
    newcvm = ri.dot(ri.T)
    model = param.dot(dm.T)
    postfit = b - model
    chisq = np.sum(np.power(postfit, 2))
    newcvm *= (chisq / float(len(b) - len(param)))
    return param, newcvm, chisq

def cholspec(residuals, times, errors, cov, f_max=2.0):
    C = cov + np.diag(errors**2)


    sec_per_year = 365.25 * 60. * 60. * 24.
    sec_per_day = 60. * 60. * 24.

    L = np.linalg.cholesky(C)
    Linv = linalg.inv(L) # Not efficient, but ok for now


    w = Linv.dot(residuals)  # This is the "whitened residuals"
    # First we will fit just a mean value to get the comparative log-likelihood
    # Might want to fit pulsar model here too?
    M = np.zeros((1, len(times)))  # Design matrix with 1 parameter

    M[0] = 1.0  # Set all elements to 1.0

    wM = Linv.dot(M.T)  # "Whitened" design matrix
    beta, cvm, chisq = qrfit(wM, w)  # Solve least-squares for whitened design matrix and whitened data.
    default_ll = -0.5 * chisq  # Assume log-likelihood is just chi-square

    total_time = (np.amax(times)-np.amin(times)) # days
    avg_delta_times = total_time/(len(times)-1)

    #f_Nyq = 1./(2.*avg_delta_times) # per day; max freq
    #print ('MAX F = ',f_Nyq)
    f_Nyq = f_max
    #delta_f = 1.0/(10.*total_time) # frequency step, per day
    delta_f = 1.0/(total_time) # frequency step, per day
    N_f = int(f_Nyq/delta_f) # number of spectral channels
    f = np.linspace(delta_f,f_Nyq, N_f, endpoint=False) # per day

    # Create the output arrays
    dpower = np.zeros(N_f)
    complex_spec = np.zeros(N_f,dtype=complex)
    loglike = np.zeros(N_f)


    angf = f*2.*np.pi # rad/day
    wt_matrix = np.outer (angf,times)
    sin_wt_matrix = np.sin(wt_matrix)
    cos_wt_matrix = np.cos(wt_matrix)
    # sin and cos have shapes = len(f),len(times)
    print ('chisq = '+ str(chisq))
    print ('redchisq = '+str(chisq/len(times)))
    print ('\n')

    for ifreq in range(len(f)):
        if ifreq%100==0:
            print("\r{:03d}/{:03d}          ".format(ifreq,N_f),end="")
        sys.stdout.flush()

        M = np.zeros((3,len(times))) # Design Matrix of three parameters: sin, cos, mean
        M[0] = sin_wt_matrix[ifreq] # Elements of DM for sin
        M[1] = cos_wt_matrix[ifreq] # Elements of DM for cos
        M[2] = 1.0            # Elements of DM for mean value.

        # Whiten and solve the least-squares problem
        wM = Linv.dot(M.T)
        beta, cvm, chisq = qrfit(wM,w)
        ll = -0.5*chisq
        delta_ll = ll - default_ll

        #Output results
        complex_spec[ifreq] = beta[0]*1.j + beta[1]
        dpower[ifreq] = beta[0]**2 + beta[1]**2 # sec^2
        loglike[ifreq] = delta_ll

    print("Done.")
    psd_f_yr3 = dpower*sec_per_day/delta_f/(sec_per_year**3) # yr^3
    return f, psd_f_yr3, complex_spec, loglike



def PSD(freqs_days,fc_yr,logA,gamma):
    days_per_yr=365.25
    # logA is log10(A)
    Amp2 = np.power(np.power(10.,logA),2)
    fc_days = fc_yr/days_per_yr

    psd = 1./(12.* np.power(np.pi,2))* Amp2 * np.power(np.power((freqs_days/fc_days),2)+1.,-gamma/2.) * (fc_yr)**(-gamma) # yr^3

    return psd


def getC(times,logA,gamma,fc_yr):
    days_per_yr = 365.25
    sec_per_year = 365.25 * 60. * 60. * 24.
    sec_per_day = 60. * 60. * 24.
    # note: 'times' can be irregularly sampled; also negative
    total_time = (np.amax(times)-np.amin(times)) # days
    ndays_time = int(total_time+1)
    #fc_days = 1./total_time
    #fc_days = 0.01/days_per_year
    fc_days = fc_yr/days_per_yr
    # create regular sequence of times..
    #.. to include 'times' after the iFFT
    npts_regtimes = 128
    delta_regtimes = 1 #day
    while npts_regtimes*delta_regtimes < (ndays_time+1)*2 or npts_regtimes*delta_regtimes < (2./fc_days):
        npts_regtimes *= 2

    # this creates an array of times to use in covFunc...
    # ... regtimes =  [0,1,...,npts_regtimes-1]*delta_regtimes
    # ... in days

    # which then gives the FT frequencies...
    # ... [0,1/npts_regtimes,...,(npts_regtimes/2)/npts_regtimes]]*(1/delta_regtimes)
    # ... in days^-1

    nu = np.fft.rfftfreq(npts_regtimes,delta_regtimes) #day^-1
    # len(nu) == npts_regtimes/2+1
    # delta_nu == 1./(npts_regtimes*delta_regtimes)
    #delta_nu = (nu[-1] - nu[0])/nu.shape[0] #day^-1

    # calculate power spectral density
    psd_nu = PSD(nu, fc_yr, logA, gamma) # yr^3

    # get the covariance function:
    covFunc = np.fft.irfft(psd_nu*sec_per_year**3/(delta_regtimes*sec_per_day)) #sec^2
    # len(covFunc) == (len(nu)-1)*2 == npts_regtimes

    # need covFunc to work for abs(times[j]-times[i]), any j,i
    # ndays_time-1 <max (abs(times[j]-times[i])) =< ndays_time

    # so interpolate it for the regular times up to ndays_time:
    covFunc = interp.interp1d(np.arange(ndays_time+1), covFunc[:ndays_time+1])

    C = np.empty([len(times),len(times)])
    #alltimes = np.empty([len(times),len(times)])
    for i in range(len(times)):
        try:
            C[i] = covFunc(np.abs(times-times[i]))
            #alltimes[i] = np.abs(times-times[i])
        except ValueError as e:
            print (e)
            print (np.abs(times-times[i]))
    return C # have covariance matrix!