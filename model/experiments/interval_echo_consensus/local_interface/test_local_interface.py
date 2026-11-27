"""Independent full-register and physically accessible record checks."""
import copy
import json

import numpy as np
import pytest

from . import measure as m
from model.candidate_laws.contrast_memory_spherical_v1 import law


def independent_branches(adaptive):
    u = law.echo()
    out = {}
    for r, (a,b) in enumerate(law.SRC):
        p = np.diag(np.tile(np.arange(9) == r,11))
        f = u if adaptive and a<b else np.eye(99)
        for s in range(9):
            q = np.diag(np.tile(np.arange(9) == s,11))
            out[r,s] = q@u@f@p@u
    return out


@pytest.mark.parametrize('adaptive',[False,True])
def test_all_native_branches_on_full_dirty_collar_alphabet(adaptive):
    actual=m.branch_maps(adaptive)
    expected=independent_branches(adaptive)
    gram=np.zeros((99,99),complex)
    charge=np.tile([a+b for a,b in law.SRC],11)
    for key,k in actual.items():
        np.testing.assert_allclose(k,expected[key],atol=1e-14,rtol=0)
        assert np.max(np.abs(k*(charge[:,None]-charge[None,:])))<1e-14
        gram+=k.conj().T@k
    np.testing.assert_allclose(gram,np.eye(99),atol=1e-14,rtol=0)
    # Rectangular blank-input maps also sum to the9-input identity.
    blank=sum((k[:,:9].conj().T@k[:,:9] for k in actual.values()),np.zeros((9,9),complex))
    np.testing.assert_allclose(blank,np.eye(9),atol=1e-14,rtol=0)


def test_fixed_charge_reference_every_remote_matrix_entry_and_spectator():
    c=m.balanced_reference()
    q=np.tile([a+b for a,b in law.SRC],11)
    assert np.linalg.matrix_rank(c)==99
    assert np.max(np.abs(c*(q[:,None]+q[None,:]-4)))==0
    # Complex spectator is retained, correlated with a full-rank charged reference.
    rng=np.random.default_rng(901244)
    c3=np.stack([c,c*(1j*np.exp(1j*rng.normal(size=(99,99))))],axis=2)/np.sqrt(2)
    joined=c3.reshape(99,198)
    remote=joined.T@joined.conj()
    after=np.zeros((198,198),complex)
    for k in independent_branches(True).values():
        out=k@joined
        after+=out.T@out.conj()
    np.testing.assert_allclose(after,remote,atol=1e-14,rtol=0)


def test_collar_cannot_be_dropped_from_a_reusable_local_interface():
    u=law.echo()
    blank=np.zeros(99);blank[law.SRC.index((1,1))]=1
    dirty=np.zeros(99);dirty[9*(1+law.LABELS.index((0,-2)))+law.SRC.index((1,1))]=1
    def material(v):return (np.abs(v.reshape(11,9))**2).sum(axis=0)
    np.testing.assert_array_equal(material(blank),material(dirty))
    before,after=material(u@blank),material(u@dirty)
    assert before[law.SRC.index((0,2))]==0
    assert abs(after[law.SRC.index((0,2))]-.5)<1e-14


def test_scheduler_consults_only_accessible_material_record_and_changes_it():
    controls=m.controls()
    assert controls['feedback_record_total_variation']>.1
    assert controls['distinct_collars_shared_port_order_record_tv']>.1
    assert abs(controls['open_outside_load1_probability']-.5)<1e-13
    # Native U is not an inverse of a measurement: first record is retained.
    v=np.zeros(99);v[law.SRC.index((0,1))]=1
    p=m.probabilities(v,True)
    assert abs(p.sum()-1)<1e-14
    assert abs(p.sum(axis=1)[law.SRC.index((0,1))]-.5)<1e-14
    assert abs(p.sum(axis=1)[law.SRC.index((1,0))]-.5)<1e-14


def test_charge_split_is_direct_sum_not_fixed_sector_tensor_product():
    geo=m.geometry()
    d=geo['local_charge_dimensions']
    assert sum(d)==99 and sum(d[q]*d[4-q] for q in range(5))==2299
    assert geo['four_port_two_collar_sector_dimensions']['4']==2299
    assert sum(geo['four_port_two_collar_sector_dimensions'].values())==99**2
    assert [1,5] in geo['boundary_seams']
    assert geo['reference_qutrits']==geo['reference_collars']==180


def test_exact_certificate_and_receipt_tamper():
    assert all(m.exact_certificate().values())
    actual=m.report()
    saved=json.loads(m.RECEIPT.read_text())
    m.compare(saved,actual)
    changed=copy.deepcopy(saved)
    changed['public_records']['all_pairs']=80
    with pytest.raises(ValueError):m.compare(changed,actual)
    changed=copy.deepcopy(saved)
    changed['reference']['normalized_remote_marginal_error']+=.01
    with pytest.raises(ValueError):m.compare(changed,actual)


def test_record_selects_actual_service_support_and_union_complement():
    evidence=m.selected_support_audit()
    assert all(evidence['exact_per_first_record_completeness'].values())
    assert evidence['complement_marginal_max_error']<1e-13
    assert evidence['selected_B_conditional_complement_trace_distance']>.01
    assert abs(evidence['wrong_support_record_total_variation']-.5)<1e-13
    records=np.array(evidence['source_selected_records'])
    assert abs(records[m.LOADS.index((0,1)),m.LOADS.index((1,1))]-.5)<1e-13
    # Independently reconstruct factored branches from the owning native echo.
    u=law.echo(); identity=np.eye(99)
    for r,t,support,left,right in m.selected_support_branches():
        p=np.diag(np.tile(np.arange(9)==r,11))
        q=np.diag(np.tile(np.arange(9)==t,11))
        if m.LOADS[r][0]<m.LOADS[r][1]:
            assert support==(2,3)
            np.testing.assert_allclose(m.numeric(left),p@u,atol=1e-14,rtol=0)
            np.testing.assert_allclose(m.numeric(right),q@u,atol=1e-14,rtol=0)
        else:
            assert support==(0,1)
            np.testing.assert_allclose(m.numeric(left),q@u@p@u,atol=1e-14,rtol=0)
            np.testing.assert_array_equal(m.numeric(right),identity)


def test_archive_checks_current_runtime_sources_and_each_physical_output():
    import copy
    import json
    from .provenance import PRODUCER, REVIEWED_PRODUCER, compare_replay

    saved = json.loads(m.RECEIPT.read_text())
    live = copy.deepcopy(saved)
    live['source_sha256'][PRODUCER] = REVIEWED_PRODUCER
    live['source_sha256']['model/candidate_laws/interval_echo/SPEC.md'] = 'changed prose'
    compare_replay(saved, live)
    bad_archive = copy.deepcopy(saved)
    bad_archive['source_sha256']['model/candidate_laws/interval_echo/SPEC.md'] = 'invented original'
    with pytest.raises(ValueError, match='differs from the deposit'):
        compare_replay(bad_archive, live)
    # The message experiment archived a different authentic revision of the same SPEC.
    # Existence in the shared archive does not establish enrollment in this experiment.
    from model.experiments.interval_echo_message.channel import SPEC as message_protocol
    from model.experiments.record_provenance import verify_deposit
    spec_path = 'model/candidate_laws/interval_echo/SPEC.md'
    other = message_protocol['sources'][spec_path]
    assert other != saved['source_sha256'][spec_path]
    verify_deposit({spec_path: other})
    bad_archive['source_sha256'][spec_path] = other
    with pytest.raises(ValueError, match='source manifest differs'):
        compare_replay(bad_archive, live)
    for path in ('model/candidate_laws/contrast_memory_spherical_v1/law.py', PRODUCER):
        bad_input = copy.deepcopy(live)
        bad_input['source_sha256'][path] = 'changed executable'
        with pytest.raises(ValueError):
            compare_replay(saved, bad_input)
    bad_result = copy.deepcopy(live)
    bad_result['geometry']['ports'] -= 1
    with pytest.raises(ValueError):
        compare_replay(saved, bad_result)
