"""Planck-based two-component emission/extinction prediction utilities."""

import numpy as np
import math
import pyfits
import os

def par_struc_2comp():
    """Return dictionary containing various two-component parameters"""

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
    tau2ebv = 2.6244472e3
    
    offs_tau_ebv = -0.00260618 # mag E(B-V)
    
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
  """
  Calculate Riemann zeta function

  Inputs:
      s - scalar argument of Riemann zeta function
  Keyword Inputs:
      n - compute only first n terms in infinite sum, default 10

  Comments:
      intended for real arguments
  """

# ----- this function needs to be furthervalidated against IDL version

# ----- compute only first n terms in infinite sum
#        n = 10 seems about right for s ~ 5.7

  z = np.sum(1./((np.arange(n)+1)**s))

  return z


def Z(alpha=None, n=10):
    """
    Compute Z(alpha), as defined in FDS99 equation 16

    Inputs:
        alpha - see FDS99 equation 16

    Keyword Inputs:
        n - number of terms in Riemann zeta sum, default 10
    """
    
# ----- is there a reasonable default value for alpha ?
    z = zeta(4+alpha, n=n)*math.gamma(4+alpha)
    return z

def get_t1(T2, q1_over_q2=par_struc_2comp()['q1_over_q2'],
           beta1=par_struc_2comp()['beta1'], beta2=par_struc_2comp()['beta2']):
    """
    compute cold dust temperature based on hot dust temperature

    Inputs:
        T2 - hot dust temperature, K

    Keyword Inputs:
        q1_over_q2 - two-component parameter q1/q2, default 8.219
        beta1 - cold dust emissivity power law index, default 1.63
        beta2 - hot dust emissivity power law index, default 2.82
    """
    
    hnu0_over_kb = 143.977300455 # h*nu_0/k_B, MKS
    
    fac = ((1./q1_over_q2)*Z(beta2)/Z(beta1)*((hnu0_over_kb)**(beta1-beta2)))**(1/(4+beta1))
    pow = (4+beta2)/(4+beta1)
    
    T1 = fac*(T2**pow)
    return T1

def i_2comp(nu, T2, f1=par_struc_2comp()['f1'],
            beta1=par_struc_2comp()['beta1'], beta2=par_struc_2comp()['beta2'],
            q1_over_q2=par_struc_2comp()['q1_over_q2']):
    """
    Evaluate two-component dust spectrum
    
    Inputs:
        nu - frequencies at which to evaluate the model, GHz
        T2 - temperature of hot dust component, K

    Keyword Inputs:
        f1 - default to 0.0485
        beta1 - cold component emissivity power law index, default to 1.63
        beta2 - hot component emissivity power law index, default to 2.82
        q1_over_q2 - two-component parameter q1/q2, default to 8.219
    """
    par = par_struc_2comp()

    nu0 = par['nu0']

    hk = 0.0479924335 

    T1 = get_t1(T2, q1_over_q2=q1_over_q2, beta1=beta1, beta2=beta2)

    inten1 = f1*q1_over_q2*((nu/nu0)**(3+beta1))*(1/(np.exp(hk*nu/T1)-1))
    inten2 = (1-f1)*((nu/nu0)**(3+beta2))*(1./(np.exp(hk*nu/T2)-1))

    return (inten1+inten2), inten1, inten2

def pred_spec(nu, T2, nu_ref, i_ref, sig_ref):
    """
    Make two-component emission predictions

    Inputs:
         nu     - frequency or frequencies at which emission predictions  are
                  desired, GHz
         T2     - hot dust temperature, K
         nu_ref - reference frequency for input i_ref, GHz
         i_ref  - monochromatic two-component model intensity at reference
                  frequency
    """
    
    m_ref, _, __ = i_2comp(nu_ref, T2)
    m_nu, _, __  = i_2comp(nu, T2)

    fac = (m_nu/m_ref)
    pred = i_ref*fac

    sig_m = sig_ref*fac

    return pred, sig_m

def getval_2comp(nu=par_struc_2comp()['nu_ref'], ind=None, ebv=False, unc=False):
    """
    Predict emission or extinction with Planck-based two-component model
    
    Keyword Inputs:
        nu  - if retrieving emission predictions, gives the frequency or
              frequencies at which to evaluate two-component model,
              default 545 GHz
        ind - HEALPix indices for which predictions desired, if not set then
              full-sky predictions are returned
        ebv - if True, retrieve reddening predictions instead of emission
              predictions
    """

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
        sig_ref = tab['sig_m545']
# ----- if nu keyword doesn't specify frequency, then assume 545 GHz desired
        vals, sig_vals = pred_spec(nu, T2, par['nu_ref'], iref, sig_ref)
    else:
        vals = par['tau2ebv']*tab['tau545'] + par['offs_tau_ebv']
        if unc: sig_vals = par['tau2ebv']*tab['sig_tau545']

    if unc:
        return vals, sig_vals
    else:
        return vals
