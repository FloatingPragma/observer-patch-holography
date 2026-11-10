"""Independent exact residue assembly, source controls and rejected transfers."""
import copy
import json
from pathlib import Path
import pytest
import sympy as s
from . import model, residue, bounds, transfer


def test_actual_native_operator_and_retained_reference():
    result = model.operator_certificate()
    assert result['complete_Fock_operator_norm'] == '60'
    assert result['complete_cubic_Fourier_word_checks'] == 576
    assert result['retained_reference_rotation_translation_checks'] == 72
    assert result['actual_full_occupied_source_translation_controls'] == 9
    assert not result['continuous_spin_symmetry_at_nonzero_hopping_assumed']
    for bad in (-1, 3, True, 0.0):
        with pytest.raises(ValueError):
            model.native_spin(bad)
    for bad in (-1, s.sqrt(2), 0.5):
        with pytest.raises(ValueError):
            model.hamiltonian_column(None, bad)


def test_independent_tensor_assembly_and_exact_residue_without_eigenvectors():
    # Tensor construction uses site zero as the most significant bit,
    # independently of the production bit-column and magnetic-block routines.
    sigma = (s.Matrix([[0,1],[1,0]]), s.Matrix([[0,-s.I],[s.I,0]]), s.diag(1,-1))
    spin = [[s.ImmutableSparseMatrix(s.kronecker_product(*(
        sigma[a]/2 if target==site else s.eye(2) for target in range(8))))
        for a in range(3)] for site in range(8)]
    zero = s.SparseMatrix(256,256,{})
    h = sum((spin[u][a]*spin[v][a] for u in range(8) for v in range(u+1,8)
        if (u^v) in (1,2,4) for a in range(3)), zero)
    n = sum(((-1)**u.bit_count()*spin[u][2] for u in range(8)), zero)
    indices = [i for i in range(256) if i.bit_count()==4]
    a = 4*h.extract(indices,indices); n = n.extract(indices,indices)
    cubic = a**3+20*a**2-16*a-576*s.eye(70)
    w = s.Matrix.hstack(*cubic.nullspace()); g = w.T*w
    small = g.inv()*w.T*a*w
    t, = (a+16*s.eye(70)).nullspace()
    v = g.inv()*w.T*n*t
    line = v*v.T*g/(t.T*t)[0]
    total = g.inv()*w.T*n**2*w
    assert [s.trace(small**k*line) for k in range(3)] == [s.Rational(v,3) for v in (23,-436,8368)]
    assert [s.trace(small**k*total) for k in range(3)] == [s.Rational(v,3) for v in (32,-448,8704)]
    exact = residue.exact_data()
    assert exact['line_moments'] == [s.trace(small**k*line) for k in range(3)]
    assert exact['total_moments'] == [s.trace(small**k*total) for k in range(3)]
    assert s.cancel(exact['line']+exact['outside']-exact['success']) == 0
    assert residue.coordinate_join()['actual_source_Fourier_phase_checks'] == 9


def test_full_source_probabilities_and_complete_original_preparation():
    result = transfer.certificate()
    assert s.Rational(result['source']['hopping_coefficient']) == s.Rational(1,2**34)
    line = result['actual_interacting_energy_line']
    assert s.Rational(line['excitation_energy_interval_in_source_units'][0]) > 0
    for name in ('first_line', 'outside_line', 'total_pump_success', 'pump_failure'):
        lo, hi = map(s.Rational,line[name+'_probability_conditioned_on_ground_preparation'])
        assert 0 < lo <= hi < 1
    records = result['complete_preparation_and_detector_records']
    assert len(records['records']) == records['positive_record_count'] == 8
    dimension = (5*2**30)**8*119**24
    assert s.Rational(records['ground_preparation_success_probability']) == s.Rational(1,dimension)
    assert s.Rational(records['records'][0]['probability']) == 1-s.Rational(1,dimension)
    for row in records['records'][2:]:
        assert sum(row['source_momentum_in_pi_units']) == 1
        assert row['joint_momentum_in_pi_units'] == [0,0,0]
    assert not records['conditioned_initial_filter_failures_omitted']


def test_invalid_cluster_and_source_bounds_are_rejected(monkeypatch):
    for args in ((0,360,s.Rational(1,2**34),1,1),
                 (1,360,1,1,1), (1,1,s.Rational(1,10),s.Rational(1,1000),1)):
        with pytest.raises(ValueError):
            bounds.residue_errors(*args)
    for args in (((-1,1),0), ((0,2),0), ((0,1),-1), ((0,1),0.1)):
        with pytest.raises(ValueError):
            bounds.probability_interval(*args)
    original = transfer.geometry_data
    def zero_exchange():
        geometry, spin = copy.deepcopy(original())
        geometry['actual_exchange_interval'] = ['0','1']
        return geometry,spin
    monkeypatch.setattr(transfer,'geometry_data',zero_exchange)
    with pytest.raises(ValueError,match='positive ordered'):
        transfer.data()


def test_geometry_rank_and_dependency_custody_reject_mutation(monkeypatch,tmp_path):
    original_root = transfer.ROOT
    receipt = json.loads((original_root/transfer.GEOMETRY_RECEIPT).read_text())
    receipt['full_interacting_transfer']['first_excited_cluster_rank'] = 2
    path = tmp_path/transfer.GEOMETRY_RECEIPT
    path.parent.mkdir(parents=True); path.write_text(json.dumps(receipt))
    monkeypatch.setattr(transfer,'ROOT',tmp_path)
    with pytest.raises(ValueError,match='rank-three'):
        transfer.geometry_data()
    from . import certificate as custody
    monkeypatch.setattr(custody,'DEPENDENCIES',(*custody.DEPENDENCIES,'missing-collective-spin.json'))
    with pytest.raises(FileNotFoundError,match='missing declared collective spin'):
        custody.input_hashes()
