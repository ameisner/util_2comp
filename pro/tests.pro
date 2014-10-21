;+
; NAME:
;   tests
;
; PURPOSE:
;   unit tests for util_2comp IDL routines
;
; CALLING SEQUENCE:
;   tests
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
;   2014-Oct-19 - Aaron Meisner
;----------------------------------------------------------------------
pro test_em_ref

; test that reference frequency emission prediction is returned properly

  par = par_struc_2comp()
  pred = getval_2comp()
  fname = concat_dir(getenv('ETC_2COMP'), par.fname)
  str = mrdfits(fname, 1)
  assert, n_elements(pred) EQ n_elements(str.m545)
  assert, total(str.m545 NE pred) EQ 0

end

pro test_unc_ref

; test that reference frequency emission uncertainty is returned properly

  par = par_struc_2comp()
  pred = getval_2comp(unc=unc)
  fname = concat_dir(getenv('ETC_2COMP'), par.fname)
  str = mrdfits(fname, 1)

  assert, n_elements(pred) EQ n_elements(str.m545)
  assert, total(str.m545 NE pred) EQ 0

  assert, n_elements(unc) EQ n_elements(str.sig_m545)
  assert, total(str.sig_m545 NE unc) EQ 0

end

pro test_em_ref_partial

; test reference frequency emission prediction again, but for some
; subregion of sky

  pix = lindgen(2048L*2048)*12

  par = par_struc_2comp()
  pred = getval_2comp(ind=pix)

  fname = concat_dir(getenv('ETC_2COMP'), par.fname)
  str = mrdfits(fname, 1)

  assert, n_elements(pred) EQ n_elements(pix)
  assert, total(pred NE str[pix].m545) EQ 0

end

pro test_unc_ref_partial

; test reference frequency emission uncertainty for some subregion of sky

  pix = lindgen(2048L*2048)*12

  par = par_struc_2comp()
  pred = getval_2comp(ind=pix, unc=unc)

  fname = concat_dir(getenv('ETC_2COMP'), par.fname)
  str = mrdfits(fname, 1)

  assert, n_elements(pred) EQ n_elements(pix)
  assert, total(pred NE str[pix].m545) EQ 0

  assert, n_elements(unc) EQ n_elements(pix)
  assert, total(unc NE str[pix].sig_m545) EQ 0

end

pro test_em_nu

; test emission prediction at non-reference frequency

  nu = 353. ; arb.

  par = par_struc_2comp()
  pred = getval_2comp(nu=nu)

  fname = concat_dir(getenv('ETC_2COMP'), par.fname)
  str = mrdfits(fname, 1)

  assert, n_elements(pred) EQ n_elements(str.m545)
  assert, total(pred EQ str.m545) EQ 0

end

pro test_unc_nu

; test that uncertainty and emission are scaled appropriately to 
; non-reference frequency

end

pro test_em_nu_partial

; test emission prediction again at non-reference frequency, but for some
; subregion of sky

  pix = lindgen(2048L*2048)*12

  nu_test = 857. ; arb

  pred = getval_2comp(ind=pix, nu=nu_test)

  pred_full = getval_2comp(nu=nu_test)

  assert, n_elements(pred) EQ n_elements(pix)
  assert, total(pred NE pred_full[pix]) EQ 0

end

pro test_one_pix_nu_ref

; test the case of a single pixel emission prediction at reference frequency

  pix = 10 ; arb

  pred = getval_2comp(ind=pix)

  pred_full = getval_2comp()

  assert, n_elements(pred) EQ n_elements(pix)
  assert, pred EQ pred_full[pix]

end

pro test_one_pix_one_freq

; emission prediction for a single pixel at a single non-reference frequency

  pix = 10000000 ; arb

  nu_test = 353. ; arb

  pred = getval_2comp(ind=pix, nu=nu_test)

  pred_full = getval_2comp(nu=nu_test)

  assert, n_elements(pred) EQ n_elements(pix)
  assert, pred EQ pred_full[pix]

end

pro test_rat_em_unc

; test that ratio of predicted emission to its uncertainty scales properly

  pred_ref = getval_2comp(unc=unc_ref)

  nu_test = 1250. ; arb
  pred = getval_2comp(nu=nu_test, unc=unc_test)

  tol = 1e-6

  assert, max(abs((unc_test/pred)-(unc_ref/pred_ref))) LT tol

end

; repeat all unit tests for reddening rather than emission ?

; tests for case in which single pixel requested, but for multiple freqs
; (this only can apply to emission predictions)

; tests for case in which frequencies and pixels are arrays of same length?

pro tests

; run unit tests

  test_em_ref
  test_em_ref_partial
  test_em_nu_partial
  test_one_pix_nu_ref
  test_one_pix_one_freq
  test_unc_ref
  test_unc_ref_partial
  test_em_nu
  test_rat_em_unc

end
