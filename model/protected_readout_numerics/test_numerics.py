from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import json
import numpy as np
import pytest
from . import model,kernel,arithmetic,interval_replay
from model.protected_particle_readout import packet


@pytest.fixture(scope='module')
def reference():
    return model.retained()


def test_exact_kernel_and_literal_control_are_distinct():
    x=kernel.certificate()
    assert len(x['exact_template_bond_pairs'])==12
    assert x['all1820_eta_image_columns_annihilated']
    assert x['exact_target_template_all_terminal_native_probabilities']==0
    assert not x['literal_finite_control_word_terminal_probability_is_forced_zero']
    assert 0<F(x['literal_all_terminal_native_probability_upper'])<F(7318,10**24)


def test_original_raw_arrays_certified_with_full_complex_phase(reference):
    value,arrays=reference
    rows=model.validate(value,arrays)
    assert all(F(r['complete_state_error_upper'])<F(36,10**12)for r in rows)
    changed=dict(arrays)
    changed['analyzer_0_full_success_state']=arrays['analyzer_0_full_success_state'].copy()
    changed['analyzer_0_full_success_state']*=np.exp(1e-5j)
    with pytest.raises(ValueError,match='full state exceeds'):
        model.validate(value,changed)


def test_backend_bounds_are_outward_and_recomputed():
    bound,primitives,budget=model.backend_certificate()
    assert F(bound['global_prefix_state_error_upper'])<F(budget['complete_state_error_upper'])
    assert all(F(x['complete_M6_error_upper'])<=F(budget['local_matrix_error_cap'])for x in primitives)
    target=arithmetic.ident(2)
    assert arithmetic.matrix_error(np.eye(2),target)==0
    damaged=np.eye(2,dtype=complex);damaged[1,0]=1e-4j
    assert arithmetic.matrix_error(damaged,target)>F(9,10**5)


def test_arbitrary_small_probabilities_are_not_silently_accepted(reference):
    value,arrays=reference
    for kind in ('native','outside','impossible_export'):
        changed=deepcopy(value)
        if kind=='native':
            rows=changed['executed_native_readout']['all130_repair_records'];rows[0]['probability']=1e-8
            changed['executed_native_readout']['terminal_native_repair_probability']=sum(r['probability']for r in rows if r['record'][0]=='native')
        elif kind=='outside':changed['executed_native_readout']['all130_repair_records'][8]['probability']=.999
        else:changed['executed_native_readout']['pointer_settings'][0]['pointer_probabilities'][2]=1e-20
        with pytest.raises(ValueError):model.validate(changed,arrays)


def test_complete_failure_arrays_and_padding_remain_required(reference):
    value,arrays=reference
    for change in ('padding','large','missing'):
        altered=dict(arrays)
        if change=='missing':altered.pop('copy_failure_weights')
        else:
            altered['copy_failure_weights']=arrays['copy_failure_weights'].copy()
            altered['copy_failure_weights'][1,2 if change=='padding'else 0]=1e-25 if change=='padding'else 1e-8
        with pytest.raises(ValueError):model.validate(value,altered)


def test_science_bounds_scope_and_unlisted_metadata_remain_exact(reference):
    value,_=reference;_,_,pins=model.specification()
    current=deepcopy(value);current['source_bindings']=pins
    model.compare_science(value,current)
    for change in ('status','bound','source','scope','sample'):
        altered=deepcopy(current)
        if change=='status':altered['physical_particle_identification']=True
        elif change=='bound':altered['exact_native_error_certificate']['full_Fock_Strang_error']='1/2'
        elif change=='source':altered['source_bindings']['model/protected_particle_readout/packet.py']='0'*64
        elif change=='scope':altered['executed_native_readout']['resolved_all_eight_frequencies_from_these_two_settings']=True
        else:altered['executed_native_readout']['pointer_settings'][0]['simulated_counts']['pointer_zero']+=1
        with pytest.raises(ValueError):model.compare_science(value,altered)


def test_original_checker_keeps_its_probability_rejection(reference):
    value,_=reference;altered=deepcopy(value)
    altered['executed_native_readout']['all130_repair_records'][0]['probability']=7.367274713651182e-35
    with pytest.raises(ValueError,match='numerical evidence changed'):packet.compare(value,altered)
    spec=json.loads((model.HERE/'contract.json').read_text())
    assert model.digest(model.ROOT/spec['original_producer'])==spec['original_producer_sha256']


def test_legacy_projection_preserves_science_and_only_named_source_transition():
    _,old,pins=model.specification();current=deepcopy(old);current['source_bindings']=pins
    model.verify_legacy_projection(old,current,pins)
    for kind in ('number','status','unknown_source'):
        changed=deepcopy(current)
        if kind=='number':changed['executed_native_readout']['terminal_native_repair_probability']=0.
        elif kind=='status':changed['integrated_gate_accepted']=True
        else:changed['source_bindings']['unlisted.py']='0'*64
        with pytest.raises(ValueError):model.verify_legacy_projection(old,changed,pins)
    projected=model.live_projection();projected['source_bindings']=old['source_bindings']
    assert projected==old
