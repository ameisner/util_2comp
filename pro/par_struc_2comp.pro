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

; ----- name of results file
  fname = 'planck_2comp.fits'

; ----- Planck+DIRBE best-fit two-component model global parameters
  f1 = 0.0485d
  q1_over_q2 = 8.219d
  beta1 = 1.63d
  beta2 = 2.82d

  nu0 = 2997.92458d ; GHz

  nu_ref = 545.d ; GHz

; ----- boundaries for recommended range of two-component model applicability
  nu_min = 100  ; GHz
  nu_max = 3000 ; GHz

; ----- conversion factor from 545 GHz optical depth to E(B-V)
  tau2ebv = 2.6244472e3

  offs_tau_ebv = -0.00260618 ; mag E(B-V)
  
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
         nu_min       : nu_min,       $
         nu_max       : nu_max,       $ 
         nside        : nside          }

  return, par

end
