"""Selection, observable information and statistics are separate claims."""
import numpy as np

from . import mechanism


def test_local_charge_and_statistics_identities():
    result=mechanism.report()
    assert result['disjoint_site_creators_commute']
    assert result['two_site_created_charges']==[-2,0,2]
    assert not result['internal_spin_half_implies_Lorentz_spin_half']


def test_protected_information_needs_a_sensitive_readout():
    # The admitted identity-repair branches preserve every state but give
    # the same complete record distribution, including the zero branch.
    identity=np.eye(2)
    branches=(identity/np.sqrt(3),np.sqrt(2/3)*identity,np.zeros((2,2)))
    plus=np.array([[.5,.5],[.5,.5]])
    minus=np.array([[.5,-.5],[-.5,.5]])
    for a in branches:
        for b in branches:
            history=b@a
            assert np.isclose(np.trace(history@plus@history.T),
                              np.trace(history@minus@history.T))
    pointer=plus
    assert np.isclose(np.trace(pointer@plus),1)
    assert np.isclose(np.trace(pointer@minus),0)


def test_charge_conservation_allows_redistribution_inside_a_sector():
    # Independent finite control: two distinct configurations of charge+1
    # mix although the charge remains exactly protected.
    charge=np.diag([0.,1.,1.])
    h=np.array([[0.,0.,0.],[0.,0.,1.],[0.,1.,0.]])
    assert np.array_equal(charge@h,h@charge)
    initial=np.array([0.,1.,0.]); final=h@initial
    assert np.dot(initial,final)==0
    assert np.dot(initial,charge@initial)==np.dot(final,charge@final)==1
