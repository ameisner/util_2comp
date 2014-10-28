;+
; NAME:
;   getval_2comp
;
; PURPOSE:
;   predict emission or extinction with Planck-based two-component model
;
; CALLING SEQUENCE:
;   vals = getval_2comp(nu=, ind=, ebv=, unc=)
;
; INPUTS:
;   
; OPTIONAL INPUTS:
;   nu  - if retrieving emission predictions, gives the frequency or 
;         frequencies at which to evaluate two-component model,
;         default 545 GHz
;   ind - HEALPix indices (nested order) for which predictions are desired, if
;         not set then full-sky predictions are returned
;
; KEYWORDS:
;   ebv - if set, retrieve reddening predictions instead of emission 
;         predictions
;
; OUTPUTS:
;   vals - by default, predicted emission values in MJy/sr but if ebv keyword 
;          set, output is reddening in mag E(B-V)
;
; OPTIONAL OUTPUTS:
;   unc - retrieve 1 sigma uncertainties on output values
;
; EXAMPLES:
;   
; COMMENTS:
;
; REVISION HISTORY:
;   2014-Aug-3 - Aaron Meisner
;----------------------------------------------------------------------
function getval_2comp, nu=nu, ind=ind, ebv=ebv, unc=unc

  par = par_struc_2comp()
  check_freq, nu

  if (n_elements(ind) EQ 0) then ind = lindgen(12L*par.nside*par.nside)

  fname = concat_dir(getenv('ETC_2COMP'), par.fname)

; ----- could try caching to reduce read-in times
  str = mrdfits(fname, 1, rows=ind)

  if ~keyword_set(ebv) then begin
; ----- if nu keyword doesn't specify frequency, then assume 545 GHz desired
      if ~keyword_set(nu) then nu = par.nu_ref
      vals = ((n_elements(nu) EQ 1) && (nu EQ par.nu_ref)) ?  str.m545 : $ 
          pred_spec(nu, double(str.T2), par.nu_ref, double(str.m545), $
                    double(str.sig_m545), sig_m)
      if arg_present(unc) then begin
          unc = ((n_elements(nu) EQ 1) && (nu EQ par.nu_ref)) ? $ 
              str.sig_m545 : sig_m
      endif
  endif else begin
      vals = par.tau2ebv*str.tau545 + par.offs_tau_ebv
      if arg_present(unc) then unc = par.tau2ebv*str.sig_tau545
  endelse

  return, vals
end
