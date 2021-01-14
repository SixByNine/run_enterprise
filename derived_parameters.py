import numpy as np
import numpy.polynomial as polynomial


def compute_derived_parameters(args,psr,par,parfile,scaled_samples,log_prob,weights,parnames):
    derived=[]
    derived_names=[]
    if args.legendre:
        l1=l2=l3=0
        if "L1" in parnames:
            l1=scaled_samples[:,parnames=="L1"][:,0]
        else:
            l1 = np.zeros(scaled_samples.shape[0])
        if "L2" in parnames:
            l2 = scaled_samples[:, parnames=="L2"][:,0]
        else:
            l2 = np.zeros(scaled_samples.shape[0])
        if "L3" in parnames:
            l3 = scaled_samples[:, parnames=="L3"][:,0]
        else:
            l3=np.zeros_like(l1)
        domain = [psr.toas[0],psr.toas[-1]]

        for line in parfile:
            e = line.strip().split()
            if len(e) > 1:
                if e[0] == "F0":
                    f0 = float(e[1])


        print("XXX ",scaled_samples.shape,l1.shape,l2.shape,l3.shape)

        xxx = polynomial.legendre.Legendre([0,0,0,6],domain=domain,window=[-1,1])
        poly = xxx.convert(kind=polynomial.polynomial.Polynomial, domain=domain,window=domain)
        print(poly.coef*f0*np.array([1,1,2,6]))
        nsamples=len(l1)
        df0=np.zeros_like(l1)
        df1=np.zeros_like(l1)
        df2=np.zeros_like(l1)
        for s in range(nsamples):
            leg = polynomial.legendre.Legendre([0,l1[s],l2[s],l3[s]],domain=domain,window=[-1,1])
            poly = leg.convert(kind=polynomial.polynomial.Polynomial, domain=domain,window=domain)
            df0[s] = poly.coef[1]*f0
            df1[s] = 2*poly.coef[2]*f0
            df2[s] = 6*poly.coef[3]*f0

        derived.append(df0)
        derived.append(df1)
        derived.append(df2)
        derived_names.append("dF0")
        derived_names.append("dF1")
        derived_names.append("F2")

    return np.array(derived).T,derived_names




