;+
; NAME:
;   pred_spec
;
; PURPOSE:
;   make two-component emission predictions
;
; CALLING SEQUENCE:
;   pred = pred_spec(nu, T2, nu_ref, i_ref)
;
; INPUTS:
;   nu     - frequency or frequencies at which emission predictions  are 
;            desired, GHz
;   T2     - hot dust temperature, K
;   nu_ref - reference frequency for input i_ref, GHz
;   i_ref  - monochromatic two-component model intensity at reference frequency
;
; OPTIONAL INPUTS:
;   
; KEYWORDS:
;   
; OUTPUTS:
;   pred - emission predictions, same units as i_ref
;
; OPTIONAL OUTPUTS:
;   
; EXAMPLES:
;   
; COMMENTS:
;   intended for either vector nu and scalar T2 or vector T2 and
;   scalar nu
;
; REVISION HISTORY:
;   2014-Aug-3 - Aaron Meisner
;----------------------------------------------------------------------
function pred_spec, nu, T2, nu_ref, i_ref

  m_ref = i_2comp(nu_ref, T2)
  m_nu = i_2comp(nu, T2)

  pred = i_ref*(m_nu/m_ref)

  return, pred
end
