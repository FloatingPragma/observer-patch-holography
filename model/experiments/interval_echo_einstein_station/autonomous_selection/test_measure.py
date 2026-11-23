import json
from math import comb
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.executable_model import federation
from ..equilibrium_refinement.measure import apparatus, generator
from . import ancestral, measure, open_host


def test_unmonitored_spectrum_and_equal_energy_preparations():
    car = federation(0)
    _,contacts,_,words = apparatus(car,[0])
    h = measure.pair_hamiltonian(car,contacts,words)
    np.testing.assert_array_equal(h.toarray(),generator(car,words,contacts))
    indices,weights=measure.launch_pairs(h.diagonal())
    np.testing.assert_allclose(h.diagonal()[indices] @ weights,h.diagonal().mean())
    vals,vecs=np.linalg.eigh(h.toarray())
    direct = abs(expm(-.5j*h.toarray())[:,indices.ravel()])**2
    spectral = abs(vecs @ (np.exp(-.5j*vals)[:,None]*vecs[indices.ravel()].T))**2
    np.testing.assert_allclose(direct,spectral,atol=1e-12)
    # All within-block coherences of a degenerate identity Hamiltonian survive dephasing.
    modes,_=np.linalg.qr(np.array([[3,1],[1,-3]],float))
    avg,n=measure.diagonal_record_laws(np.array([4,4]),modes,np.array([0,1]),1e-12)
    assert n==1
    np.testing.assert_allclose(avg,np.eye(2),atol=1e-14)


def test_native_open_board_stationary_obstruction_and_dilution():
    fresh = open_host.run()
    directory=Path(__file__).parent
    saved=json.loads((directory/'open_host.json').read_text())
    for row,comparison in zip(fresh['rows'],saved['rows'],strict=True):
        assert row['full_generator_eigenvector_residual']<1e-10
        for reading,stored in zip(row['preparations'],comparison['preparations'],strict=True):
            assert reading['inclusive_counts']==stored['inclusive_counts']
            np.testing.assert_allclose(reading['signed_slope_nats'],
                                       np.log((9-reading['alpha'])/(9+reading['alpha'])), atol=2e-15)
        assert abs(reading['measured_slope_nats'])>5*reading['standard_error']
    report=ancestral.run()
    stored=json.loads((directory/'ancestral.json').read_text())
    for row,previous in zip(report['rows'],stored['rows'],strict=True):
        last=row['preparations'][-1]
        assert last['counts']==previous['preparations'][-1]['counts']
        assert last['transferred_probability']>.0005
        assert abs(last['signed_entropy_slope']-last['measured_slope'])<.001


def test_archived_endpoint_scores_and_lowest_cluster():
    rows=json.loads(Path(__file__).with_name('results.json').read_text())['rows']
    assert len(rows)==12
    for r in rows:
        flight=r['flight']
        assert flight['uniform_tv_lower_bound_for_any_lowest_cluster_state']>.92
        assert flight['maximum_cluster_eigenvector_residual']<1e-10
        reference=[]
        n_ports=sum(len(c)==1 for c in r['cut_words'])
        for word in r['interior_words']:
            # Multiplicity is a choice of the remaining cut ports.
            reference.append(comb(n_ports,2-len(word)))
        prob=np.array(reference)/r['pair_words']
        assert abs(prob.sum()-1)<1e-13
        for group in (flight['appointment_average'],flight['spectral_average']):
            for setting in group:
                reading=np.asarray(setting['coarse_interior_counts'])/measure.SHOTS
                np.testing.assert_allclose(measure.tv(reading,prob),setting['sampled_bulk_tv'])
                assert setting['sampled_bulk_tv']>.3


def test_half_filled_exact_cut_nullspace_and_accessible_sharp_record():
    from . import dense
    measured=dense.run()
    stored=json.loads(Path(__file__).with_name('dense.json').read_text())
    assert measured['n_words']==924
    assert measured['integer_polynomial']==stored['integer_polynomial']
    assert len(measured['integer_polynomial'])==5
    assert measured['exact_bound_L']==stored['exact_bound_L']
    for s,t in zip(measured['stationary_states'],stored['stationary_states'],strict=True):
        assert s['coarse_counts']==t['coarse_counts']
    assert measured['stationary_states'][-1]['normalized_signed_bulk_slope']<-.21
    assert abs(measured['stationary_states'][1]['normalized_signed_bulk_slope'])<1e-13
    assert measured['sharp_launch_energy_projector_difference_audit']<1e-9
    assert measured['guaranteed_time_averaged_terminal_score_difference']>.000198
