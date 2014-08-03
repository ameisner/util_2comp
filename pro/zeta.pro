;+
; NAME:
;   zeta
;
; PURPOSE:
;   calculate Riemann zeta function
;
; CALLING SEQUENCE:
;   z = zeta(s, n=)
;
; INPUTS:
;   s - scalar, argument of zeta function
;
; OPTIONAL INPUTS:
;   n - compute only first n terms in infinite sum, default 10
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
;   intended for real arguments
;
; REVISION HISTORY:
;   2014-Feb-25 - Aaron Meisner
;----------------------------------------------------------------------
function zeta, s, n=n

; ----- compute only first n terms in infinite sum
  if ~keyword_set(n) then n = 10 ; seems about right for s ~ 5.7

; ----- should eventually rewrite to be vectorized at some point
  z = total(1./((findgen(n)+1)^s))

  return, z
end
