#!/usr/bin/env python
"""
Unit tests for util_2comp Python routines

"""

import util_2comp
import pyfits
import numpy as np
import os

def test_em_ref():
    """test that reference frequency emission prediction is returned properly"""
    print test_em_ref.__doc__

    par = util_2comp.par_struc_2comp()
    pred = util_2comp.getval_2comp()
    fname = os.path.join(os.environ['ETC_2COMP'], par['fname'])
    hdus = pyfits.open(fname)
    str = hdus[1].data
    assert len(pred) == len(str['m545'])
    assert np.sum(str['m545'] != pred) == 0

def test_unc_ref():
    """test that reference frequency emission uncertainty is returned properly"""
    print test_unc_ref.__doc__

    par = util_2comp.par_struc_2comp()
    pred, unc = util_2comp.getval_2comp(unc=True)
    fname = os.path.join(os.environ['ETC_2COMP'], par['fname'])
    hdus = pyfits.open(fname)
    str = hdus[1].data

    assert len(pred) == len(str['m545'])
    assert np.sum(str['m545'] != pred) == 0

    assert len(unc) == len(str['sig_m545'])
    assert np.sum(str['sig_m545'] != unc) == 0

def test_em_ref_partial():
    """test reference frequency emission prediction again, but for some subregion of sky"""
    print test_em_ref_partial.__doc__

    pix = np.arange(0,12*2048*2048, 12)

    par = util_2comp.par_struc_2comp()
    pred = util_2comp.getval_2comp(ind=pix)

    fname = os.path.join(os.environ['ETC_2COMP'], par['fname'])
    hdus = pyfits.open(fname)
    str = hdus[1].data

    assert len(pred) == len(pix)

    m545 = str['m545']
    assert np.sum(pred != m545[pix]) == 0

def test_unc_ref_partial():
    """test reference frequency emission uncertainty for some subregion of sky"""
    print test_unc_ref_partial.__doc__

    pix = np.arange(0,12*2048*2048, 12)

    par = util_2comp.par_struc_2comp()
    pred, unc = util_2comp.getval_2comp(ind=pix, unc=True)

    fname = os.path.join(os.environ['ETC_2COMP'], par['fname'])
    hdus = pyfits.open(fname)
    str = hdus[1].data

    assert len(pred) == len(pix)
    m545 = str['m545']
    assert np.sum(pred != m545[pix]) == 0

    assert len(unc) == len(pix)
    sig_m545 = str['sig_m545']
    assert np.sum(unc != sig_m545[pix]) == 0

def test_em_nu():
    """test emission prediction at non-reference frequency"""
    print test_em_nu.__doc__

    nu = 353. # arb.

    par = util_2comp.par_struc_2comp()
    pred = util_2comp.getval_2comp(nu=nu)

    fname = os.path.join(os.environ['ETC_2COMP'], par['fname'])
    hdus = pyfits.open(fname)
    str = hdus[1].data

    assert len(pred) == len(str['m545'])
    assert np.sum(pred == str['m545']) == 0

def test_em_nu_partial():
    """test emission prediction again at non-reference frequency, but for some subregion of sky"""
    print test_em_nu_partial.__doc__

    pix = np.arange(0,12*2048*2048, 12)

    nu_test = 857. # arb

    pred = util_2comp.getval_2comp(ind=pix, nu=nu_test)

    pred_full = util_2comp.getval_2comp(nu=nu_test)

    assert len(pred) == len(pix)
    assert np.sum(pred != pred_full[pix]) == 0

def test_one_pix_nu_ref():
    """test the case of a single pixel emission prediction at reference frequency"""
    print test_one_pix_nu_ref.__doc__

    pix = 10 # arb

    pred = util_2comp.getval_2comp(ind=pix)

    pred_full = util_2comp.getval_2comp()

    assert pred.size == 1
    assert pred == pred_full[pix]

def test_one_pix_one_freq():
    """emission prediction for a single pixel at a single non-reference frequency"""
    print test_one_pix_one_freq.__doc__

    pix = 10000000 # arb

    nu_test = 353. # arb

    pred = util_2comp.getval_2comp(ind=pix, nu=nu_test)

    pred_full = util_2comp.getval_2comp(nu=nu_test)

    assert pred.size == 1
    assert pred.astype('float32') == pred_full[pix]


def test_rat_em_unc():
    """test that ratio of predicted emission to its uncertainty scales properly"""
    print test_rat_em_unc.__doc__

    pred_ref, unc_ref = util_2comp.getval_2comp(unc=True)

    nu_test = 1250. # arb
    pred, unc_test = util_2comp.getval_2comp(nu=nu_test, unc=True)

    tol = 1e-6

    assert np.max(np.abs((unc_test/pred)-(unc_ref/pred_ref))) < tol

def test_multifreq():
    """test case of single-pixel, multi-frequency emission predictions"""
    print test_multifreq.__doc__

    nu = 150. + np.arange(2800) # arb
    ind = 1000000 # arb

    pred, unc = util_2comp.getval_2comp(nu=nu, ind=ind, unc=True)

    assert len(pred) == len(nu)
    assert len(unc) == len(nu)

def test_index_0():
    """case of ind=0"""
    print test_index_0.__doc__

    pix = 0

    em = util_2comp.getval_2comp(ind=pix)
    em_all = util_2comp.getval_2comp()

    assert em.size == 1
    assert em_all[pix] == em

    ebv = util_2comp.getval_2comp(ind=pix, ebv=True)
    ebv_all = util_2comp.getval_2comp(ebv=True)

    assert ebv.size == 1
    assert ebv_all[pix] == ebv.astype('float32')

def test_ebv():
    """test full-sky reddening query"""
    print test_ebv.__doc__

    par = util_2comp.par_struc_2comp()
    ebv = util_2comp.getval_2comp(ebv=True)

    assert len(ebv) == (12*par['nside']*par['nside'])

def test_ebv_partial():
    """test reddening query for some subset of pixels"""
    print test_ebv_partial.__doc__

    pix = np.arange(0,12*2048*2048, 2048)
    ebv = util_2comp.getval_2comp(ebv=True, ind=pix)
    ebv_full = util_2comp.getval_2comp(ebv=True)

    assert len(pix) == len(ebv)
    assert np.sum(ebv != ebv_full[pix]) == 0

if __name__ == '__main__':
    """run the above tests"""

    test_em_ref()
    test_em_ref_partial()
    test_em_nu_partial()
    test_one_pix_nu_ref()
    test_one_pix_one_freq()
    test_unc_ref()
    test_unc_ref_partial()
    test_em_nu()
    test_rat_em_unc()
    test_multifreq()
    test_ebv()
    test_ebv_partial()
    test_index_0()
