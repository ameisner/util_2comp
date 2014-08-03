;+
; NAME:
;   pred_spec
;
; PURPOSE:
;   make two-component emission predictions
;
; CALLING SEQUENCE:
;   
; INPUTS:
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
;   2014-Aug-3 - Aaron Meisner
;----------------------------------------------------------------------
function pred_spec, nu, T2, nu_ref, i_ref

  m_ref = i_2comp(nu_ref, T2)
  m_nu = i_2comp(nu, T2)

  pred = i_ref*(m_nu/m_ref)

  return, pred
end
