;+
; NAME:
;   i_2comp
;
; PURPOSE:
;   evaluate two-component dust spectrum
;
; CALLING SEQUENCE:
;   inten = i_2comp(nu, T2, inten1=, inten2=, f1=, beta1=, beta2=,
;       q1_over_q2=)
;
; INPUTS:
;   nu - frequencies at which to evaluate the model, GHz
;   T2 - temperature of hot dust component, K
;
; OPTIONAL INPUTS:
;   f1 - default to 0.0485
;   beta1 - cold component emissivity power law index, default to 1.63
;   beta2 - hot component emissivity power law index, default to 2.82
;   q1_over_q2 - two-component parameter q1/q2, default to 8.219
;
; OUTPUTS:
;   inten - two-component intensity I_nu evaluated at nu
;
; OPTIONAL OUTPUTS:
;   inten1 - contribution of cold dust to inten
;   inten2 - contribution of hot dust to inten
;
; EXAMPLES:
;   
; COMMENTS:
;
;   intended for either scalar nu and vector T2, or scalar T2 and vector
;   nu
;
;   also works for nu and T2 both vectors of same length, but this
;   seems like a much less common use case ...
;
;
; REVISION HISTORY:
;   2014-Aug-3 - Aaron Meisner
;----------------------------------------------------------------------
function i_2comp, nu, T2, inten1=inten1, inten2=inten2, f1=f1, $ 
                      beta1=beta1, beta2=beta2, q1_over_q2=q1_over_q2

  par = par_struc_2comp()

  if ~keyword_set(f1) then f1 = par.f1
  if ~keyword_set(beta1) then beta1 = par.beta1
  if ~keyword_set(beta2) then beta2 = par.beta2
  if ~keyword_set(q1_over_q2) then q1_over_q2 = par.q1_over_q2

  nu0 = par.nu0

  hk = 0.0479924335d

  T1 = get_t1(T2, q1_over_q2=q1_over_q2, beta1=beta1, beta2=beta2)

  inten1 = f1*q1_over_q2*((nu/nu0)^(3+beta1))*(1./(exp(hk*nu/T1)-1))
  inten2 = (1-f1)*((nu/nu0)^(3+beta2))*(1./(exp(hk*nu/T2)-1))

  return, (inten1+inten2)

end
