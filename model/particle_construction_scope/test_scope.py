"""Hostile controls for the conclusions used in the short app replacements."""
import copy
import json

import pytest

from .exact_controls import summarize
from .replay import (OUTPUTS, RETAINED, ROOT, adapt_source, differences,
                     nonfinite_leaves, project_path_outputs, retained_manifest)


def fixtures():
    fresh = {name+'/'+output: json.loads((RETAINED/name/(output+'.raw')).read_text())
             for name, outputs in OUTPUTS.items() for output in outputs}
    return fresh, {'native_V': {'exact_equal': True}}


def test_report_findings_keep_failed_and_conditional_domains():
    fresh, exact = fixtures()
    result = summarize(fresh, exact)
    assert result['CM2']['all_four_passes'] == 0
    assert result['HT']['growing_register_free_scalar_limit_excluded'] is False
    assert result['source_derived_SM_GeV_mass'] is False
    assert result['DR2']['clock_choice_gap_coupling_and_dressed_compression_declared'] is True


@pytest.mark.parametrize('changed', ['native_pass', 'fit_pass', 'hold_weight', 'analyzer_rank', 'second_tick'])
def test_unlisted_scientific_promotion_is_rejected(changed):
    fresh, exact = fixtures()
    fresh = copy.deepcopy(fresh)
    if changed == 'native_pass':
        fresh['CM/r2_compare.json']['summary']['native_all_four_pass'] = 1
    elif changed == 'fit_pass':
        fresh['CM/r2_compare.json']['fit']['A']['smallest_four'][0]['verdict']['all_four_pass'] = True
    elif changed == 'hold_weight':
        fresh['DR/register.json']['rows']['P1_link_kept_(S)']['weights']['H'] = []
    elif changed == 'analyzer_rank':
        fresh['SIM/derive/lnat3/lnat3_t4.json']['Phi_Schmidt_rank_C|E'] = 1
    elif changed == 'second_tick':
        fresh['SIM/derive/lnat3/lnat3_t4.json']['multi_tick_unread_K_n']['n=2']['K_n']['b (R)'] = '1/12'
    with pytest.raises(ValueError):
        summarize(fresh, exact)


def test_leaf_comparison_has_no_silent_float_or_boolean_equivalence():
    assert differences({'status': False}, {'status': 0})
    assert differences({'value': 1.0}, {'value': 1.0+1e-15})
    assert differences({'source': 'A'}, {'source': 'B'})
    assert differences({'a': [1,2]}, {'a': [1,2]}) == []


def test_undefined_cp_phase_is_explicit_not_a_zero_or_match():
    value = {'sin_delta': float('nan')}
    assert nonfinite_leaves(value) == [{'path': ['sin_delta'], 'value': 'nan'}]
    assert differences(value, {'sin_delta': 0.0})
    assert differences(value, {'sin_delta': float('nan')}) == []


def test_retained_original_bytes_and_typed_path_adapter():
    manifest = retained_manifest()
    assert len(manifest['files']) == 21
    name = 'CM/cm_r2_compare.py'
    raw = (RETAINED/(name+'.raw')).read_text()
    adapted, replacements = adapt_source(name, raw, ROOT)
    assert len(replacements) == 1
    assert str(ROOT) in adapted
    with pytest.raises(ValueError):
        adapt_source(name, raw.replace('/repair-law/reverse-engineering-reality', '/unregistered-source-root'), ROOT)


def test_output_path_adapter_rejects_other_source_and_preserves_science():
    old = {'RER_path_used':'historical-path','physical_mass':False}
    new = {'RER_path_used':str(ROOT),'physical_mass':True}
    a,b,record = project_path_outputs('SIM/derive/lnat3/lnat3_t1.json',old,new,ROOT)
    assert record and differences(a,b) == [{'path':['physical_mass'],'old':False,'fresh':True}]
    with pytest.raises(ValueError):
        project_path_outputs('SIM/derive/lnat3/lnat3_t1.json',old,
                             {'RER_path_used':'unreviewed-other-source','physical_mass':False},ROOT)
