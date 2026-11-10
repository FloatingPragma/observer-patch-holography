"""Typed NumPy data evidence retains the existing fail-closed claim schema."""
import pytest

from tools import check_claim_registry as claims


def test_npy_is_numeric_data_without_changing_other_artifact_kinds():
    assert claims.evidence_medium('model/affine_gauss_full_spectrum/trial_state.npy') == 'numeric_data'
    assert claims.evidence_medium('model/witness.npy#array') == 'numeric_data'
    assert claims.evidence_medium('model/witness.npz') == 'emitted_record'
    assert claims.evidence_medium('model/witness.py') == 'executed_program'
    assert claims.evidence_medium('model/witness.json.raw') is None
    assert claims.evidence_medium('model/witness.bin') is None


def test_numeric_data_does_not_allow_a_paper_proof_declaration():
    row = {'claim_id': 'numeric-fixture', 'claim_class': 'conditional_implication',
           'evidence': ['model/witness.npy']}
    assert claims.check_evidence_medium(row) is False
    with pytest.raises(SystemExit, match='admissible only on a theorem-asserting row'):
        claims.check_evidence_medium(dict(row, proof_medium='paper_prose'))


@pytest.mark.parametrize('path', ['model/witness.raw', 'model/witness.npy.raw', 'model/witness'])
def test_unclassified_binary_data_stays_rejected(path):
    with pytest.raises(SystemExit, match='evidence artifact medium is unclassified'):
        claims.check_evidence_medium({'claim_id': 'numeric-fixture',
                                     'claim_class': 'emitted_artifact', 'evidence': [path]})
