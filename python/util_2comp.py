import numpy as np
import math
import pyfits
import os

def par_struc_2comp():
    
# ----- name of results file
    fname = 'planck_2comp.fits'
    
# ----- Planck+DIRBE best-fit two-component model global parameters
    f1 = 0.0485
    q1_over_q2 = 8.219
    beta1 = 1.63
    beta2 = 2.82
    
    nu0 = 2997.92458 # GHz
    
    nu_ref = 545. # GHz
    
# ----- conversion factor from 545 GHz optical depth to E(B-V)
    tau2ebv = 2.46e3
    
    offs_tau_ebv = 0.0006 # mag E(B-V)
    
# ----- HEALPix Nside of results summary file
    nside = 2048
    
    par = {'fname'        : fname,
           'f1'           : f1,
           'q1_over_q2'   : q1_over_q2,
           'beta1'        : beta1,
           'beta2'        : beta2,
           'tau2ebv'      : tau2ebv,
           'offs_tau_ebv' : offs_tau_ebv,
           'nu0'          : nu0,
           'nu_ref'       : nu_ref,
           'nside'        : nside        }
    
    return par

def zeta(s, n=10):

# this function needs to be validated against IDL version!!!!

# ----- compute only first n terms in infinite sum
#        n = 10 seems about right for s ~ 5.7

  z = np.sum(1./((np.arange(n)+1)**s))

  return z


def Z(alpha=None, n=10):
# ----- is there a reasonable default value for alpha ?
    z = zeta(4+alpha, n=n)*math.gamma(4+alpha)
    return z

def get_t1(T2, q1_over_q2=par_struc_2comp()['q1_over_q2'],
           beta1=par_struc_2comp()['beta1'], beta2=par_struc_2comp()['beta2']):
    
    hnu0_over_kb = 143.977300455 # h*nu_0/k_B, MKS
    
    fac = ((1./q1_over_q2)*Z(beta2)/Z(beta1)*((hnu0_over_kb)**(beta1-beta2)))**(1/(4+beta1))
    pow = (4+beta2)/(4+beta1)
    
    T1 = fac*(T2**pow)
    return T1

def i_2comp(nu, T2, f1=par_struc_2comp()['f1'],
            beta1=par_struc_2comp()['beta1'], beta2=par_struc_2comp()['beta2'],
            q1_over_q2=par_struc_2comp()['q1_over_q2']):

# ----- inten1, inten2 returned as part of a tuple, and aren't optional outputs
#       as in the IDL version

    par = par_struc_2comp()

    nu0 = par['nu_ref']

    hk = 0.0479924335 

    T1 = get_t1(T2, q1_over_q2=q1_over_q2, beta1=beta1, beta2=beta2)

    inten1 = f1*q1_over_q2*((nu/nu0)**(3+beta1))*(1/(np.exp(hk*nu/T1)-1))
    inten2 = (1-f1)*((nu/nu0)**(3+beta2))*(1./(np.exp(hk*nu/T2)-1))

    return (inten1+inten2), inten1, inten2

def pred_spec(nu, T2, nu_ref, i_ref):

    m_ref, _, __ = i_2comp(nu_ref, T2)
    m_nu, _, __  = i_2comp(nu, T2)

    pred = i_ref*(m_nu/m_ref)

    return pred

def getval_2comp(nu=par_struc_2comp()['nu_ref'],
                 ind=None, ebv=False):

    par = par_struc_2comp()

    fname = os.path.join(os.environ['ETC_2COMP'], par['fname'])

# ----- is there a way to only read in the necessary field(s) as opposed
#       to reading in entire ~1.4G file
    hdus = pyfits.open(fname)
    tab = hdus[1].data # should cache this so it's only read once

# ----- still need to test case of single-element ind special case of ind = 0
    if ind is not None:
        tab = tab[ind] # checks on input ind ?

    if not ebv:
        iref = tab['m545'] # memory waste
        T2 = tab['T2']     # memory waste
# ----- if nu keyword doesn't specify frequency, then assume 545 GHz desired
        vals = (iref if (nu == par['nu_ref']) else pred_spec(nu, T2, par['nu_ref'], iref))
    else:
        vals = par['tau2ebv']*tab['tau545'] + par['offs_tau_ebv']

    return vals
