;+
; NAME:
;   z_fds
;
; PURPOSE:
;   compute Z(alpha), as defined in FDS99 equation 16
;
; CALLING SEQUENCE:
;   z_fds = z_fds(alpha, n=)
;
; INPUTS:
;   alpha - alpha, see FDS99 equation 16
;
; OPTIONAL INPUTS:
;   n - include first n terms of infinite sum when computing zeta(4+alpha)
;
; KEYWORDS:
;   
; OUTPUTS:
;   z - Z(alpha), see FDS99 equation 16
;
; OPTIONAL OUTPUTS:
;   
; EXAMPLES:
;   
; COMMENTS:
;   
; REVISION HISTORY:
;   2014-Feb-25 - Aaron Meisner
;----------------------------------------------------------------------
function z_fds, alpha, n=n

  z = zeta(4.d + alpha, n=n)*gamma(4.d + alpha)

  return, z
end
