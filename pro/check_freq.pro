;+
; NAME:
;   check_freq
;
; PURPOSE:
;   check user input frequencies and issue warning when appropriate
;
; CALLING SEQUENCE:
;   check_freq, nu
;
; INPUTS:
;   nu - frequency or array of frequencies, assumed to be in GHz
;
; OPTIONAL INPUTS:
;   
; KEYWORDS:
;   
; OUTPUTS:
;   
; OPTIONAL OUTPUTS:
;   
; EXAMPLES:
;   
; COMMENTS:
;   
; REVISION HISTORY:
;   2014-Oct-19 - Aaron Meisner
;----------------------------------------------------------------------
pro check_freq, nu

  if ~keyword_set(nu) then return

  par = par_struc_2comp()
  if (max(nu) GT par.nu_max) OR (min(nu) LT par.nu_min) then begin
      print, $ 
          'WARNING : frequency out of recommended range of model applicability'
      print, $ 
          '          note that input frequencies are expected to be in GHz'
  endif

end
