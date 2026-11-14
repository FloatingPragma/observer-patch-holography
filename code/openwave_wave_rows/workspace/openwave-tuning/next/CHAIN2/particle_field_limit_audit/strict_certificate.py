"""Corrected preparation envelopes over every supported selector setting.

The frozen producer and its native words are unchanged. This facade replaces
only the undersized plus_i-only preparation error allowance and its totals.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE = ROOT/'openwave-tuning/next/CHAIN2/refinement_apparatus'
sys.dont_write_bytecode = True
for path, expected in json.loads((HERE/'INPUTS.json').read_text())['files_sha256'].items():
    if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != expected:
        raise ValueError('strict facade source pin changed: '+path)
sys.path.insert(0, str(SOURCE))
import apparatus as b

SETTINGS = ('plus', 'minus', 'plus_i', 'minus_i', 'zero', 'one')


def preparation_bounds(size, precision):
    values = {}
    for setting in SETTINGS:
        bound = F()
        for kind, targets, angle in b.selector_template(size, setting)+b.preparation_template(size):
            event = b.a.certified_pulse(kind, tuple(targets), angle, precision)
            b.a.verify_pulse(event)
            bound += F(event['full_response_operator_error_exact'])
        values[setting] = bound
    # External selector has no prescribed selector prefix. Its reset list
    # must exclude the selector; the conservative carrier count can remain.
    values['external_selector'] = sum((F(b.a.certified_pulse(kind, tuple(targets), angle, precision)
                                        ['full_response_operator_error_exact'])
                                       for kind, targets, angle in b.preparation_template(size)), F())
    return values


def require_preparation_envelope(size, precision, claimed):
    values = preparation_bounds(size, precision)
    if F(claimed) < max(values.values()):
        raise ValueError('preparation allowance undercounts a supported selector setting')
    return True


def _replace_preparation(row, size, precision):
    values = preparation_bounds(size, precision)
    original = F(row['delta_preparation_exact'])
    corrected = max(values.values())
    assert corrected >= original
    row['original_delta_preparation_exact'] = str(original)
    row['preparation_setting_errors_exact'] = {key: str(value) for key, value in values.items()}
    row['preparation_budget_increase_exact'] = str(corrected-original)
    row['delta_preparation_exact'] = str(corrected)
    row['preparation_envelope_scope'] = 'all six prescribed selector settings and unchanged external selector/reference'
    row['certificate_facade'] = 'oph.chain2.particle_field_limit_audit.strict_preparation.v1'
    row['producer_native_words_modified'] = False
    return corrected-original


def family_certificate(size, steps, tau0=F(1)):
    row = b.family_certificate(size, steps, tau0)
    gap = _replace_preparation(row, size, row['precision_bits'])
    total = F(row['uncapped_TV_upper_exact'])+gap
    row['uncapped_TV_upper_exact'] = str(total)
    row['complete_record_TV_upper_exact'] = str(min(F(1), total))
    if F(row['delta_preparation_exact']) > F(row['h_exact'])**2:
        raise ValueError('member precision does not cover all-setting h^2 preparation budget')
    require_preparation_envelope(size, row['precision_bits'], row['delta_preparation_exact'])
    return row


def development_budget(size, steps, precision=64):
    row = b.development_budget(size, steps, precision)
    _replace_preparation(row, size, precision)
    total = (F(row['delta_preparation_exact'])+F(row['n_delta_W_exact'])+
             F(row['clock_and_detector_exact'])+F(row['reset']['timeout_bound_exact']))
    row['native_vs_ideal_complete_law_TV_exact'] = str(min(F(1), total))
    require_preparation_envelope(size, precision, row['delta_preparation_exact'])
    return row


def native_macro(size=2, steps=2, setting='plus', precision=64):
    from generate import build
    row = build(size, steps, setting, precision)
    row['budget'] = development_budget(size, steps, precision)
    row['budget_facade'] = 'oph.chain2.particle_field_limit_audit.strict_preparation.v1'
    return row
