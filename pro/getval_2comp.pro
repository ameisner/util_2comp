;+
; NAME:
;   getval_2comp
;
; PURPOSE:
;   predict emission or extinction with Planck-based two-component model
;
; CALLING SEQUENCE:
;   vals = getval_2comp(nu=, ind=, ebv=)
;
; INPUTS:
;   
; OPTIONAL INPUTS:
;   nu  - if retrieving emission predictions, gives the frequency or 
;         frequencies at which to evaluate two-component model,
;         default 545 GHz
;   ind - HEALPix indices for which predictions desired, if not set then
;         full-sky predictions are returned
;
; KEYWORDS:
;   ebv - if set, retrieve reddening predictions instead of emission 
;         predictions
;
; OUTPUTS:
;   vals - by default, predicted emission values in MJy/sr but if ebv keyword 
;          set, output is reddening in mag E(B-v)
;
; OPTIONAL OUTPUTS:
;   
; EXAMPLES:
;   
; COMMENTS:
;
; REVISION HISTORY:
;   2014-Aug-3 - Aaron Meisner
;----------------------------------------------------------------------
function getval_2comp, nu=nu, ind=ind, ebv=ebv

  par = par_struc_2comp()

; ----- still need to test case of single-element ind special case of ind = 0
  if (n_elements(ind) EQ 0) then ind = lindgen(12L*par.nside*par.nside)

  fname = concat_dir(getenv('ETC_2COMP'), par.fname)

; --- is there a way to only read in the necessary field(s) as opposed
;     to reading in entire ~1.4G file
  str = mrdfits(fname, 1) ; should cache this so it's only read once

  if ~keyword_set(ebv) then begin
      iref = str[ind].m545 ; memory waste
      T2 = str[ind].T2     ; memory waste
; ----- if nu keyword doesn't specify frequency, then assume 545 GHz desired
      if ~keyword_set(nu) then nu = par.nu_ref
      vals = (nu EQ par.nu_ref) ?  iref : pred_spec(nu, T2, par.nu_ref, iref)
  endif else begin
      vals = par.tau2ebv*str[ind].tau545 + par.offs_tau_ebv
  endelse

  return, vals
end
