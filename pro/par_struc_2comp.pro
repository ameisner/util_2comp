;+
; NAME:
;   par_struc_2comp
;
; PURPOSE:
;   serve as a repository for various two-component parameters
;
; CALLING SEQUENCE:
;   par = par_struc_2comp()
;
; INPUTS:
;   
; OPTIONAL INPUTS:
;   
; KEYWORDS:
;   
; OUTPUTS:
;   par - structure holding various Planck-based two-component parameters
;
; OPTIONAL OUTPUTS:
;   
; EXAMPLES:
;   
; COMMENTS:
;   
; REVISION HISTORY:
;   2014-Jul-31 - Aaron Meisner
;----------------------------------------------------------------------
function par_struc_2comp

; ----- what about the directory that this file is in?
  fname = 'planck_2comp.fits'

  f1 = 0.0485
  q1_over_q2 = 8.219
  beta1 = 1.63
  beta2 = 2.82

  nu0 = 2997.92458 ; GHz

  nu_ref = 545.

  tau2ebv = 2.46e3
  offs_tau_ebv = 0.0006 ; mag E(B-V)
  
; ----- HEALPix Nside of results summary file
  nside = 2048

  par = {fname        : fname,        $ 
         f1           : f1,           $ 
         q1_over_q2   : q1_over_q2,   $
         beta1        : beta1,        $
         beta2        : beta2,        $
         tau2ebv      : tau2ebv,      $
         offs_tau_ebv : offs_tau_ebv, $
         nu0          : nu0,          $
         nu_ref       : nu_ref,       $ 
         nside        : nside          }

  return, par

end
