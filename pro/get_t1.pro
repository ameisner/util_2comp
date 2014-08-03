;+
; NAME:
;   get_t1
;
; PURPOSE:
;   compute cold dust temperature based on hot dust temperature
;
; CALLING SEQUENCE:
;   T1 = get_t1(T2)
;
; INPUTS:
;   T2 - temperature of hot dust component, K
;
; OPTIONAL INPUTS:
;   q1_over_q2 - two-component parameter q1/q2
;   beta1 - cold dust emissivity power law index
;   beta2 - hot dust emissivity power law index
;
; KEYWORDS:
;   
; OUTPUTS:
;   T1 - temperature of cold dust component, K
;
; OPTIONAL OUTPUTS:
;   
; EXAMPLES:
;   
; COMMENTS:
;   FDS99 equation (19)
;
; REVISION HISTORY:
;   2014-Aug-3 - Aaron Meisner
;----------------------------------------------------------------------
function get_t1, T2, q1_over_q2=q1_over_q2, beta1=beta1, beta2=beta2

  if ~keyword_set(q1_over_q2) then q1_over_q2 = 13.
  if ~keyword_set(beta1) then beta1 = 1.67
  if ~keyword_set(beta2) then beta2 = 2.7

  hnu0_over_kb = 143.977300455 ; h*nu_0/k_B, MKS

  fac = ((1./q1_over_q2)*Z(beta2)/Z(beta1)*((hnu0_over_kb)^(beta1-beta2)))^(1/(4+beta1))
  pow = (4+beta2)/(4+beta1)

  T1 = fac*(T2^pow)
  return, T1

end
