"""Pin the actual native source property used by the correlated reset lemma."""
from pathlib import Path
import hashlib
import json

import composition
from model.native_clock import configuration

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def run():
    paths = [Path(__file__), HERE / 'RESET_DOMAIN.md', HERE / 'COMPOSITION.json',
             Path(configuration.__file__)]
    pins = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in paths}
    audit = configuration.source_audit()
    assert set(audit) == {'C', 'E'}
    assert all(row['columns_per_seam'] == 78 and row['exact_gram'] == 'I'
               and row['current_inverse_verified'] for row in audit.values())
    record = {
        'schema': 'oph.local12.correlated_reset_domain.v1',
        'checks': {'both_actual_fine_source_laws_audited': True,
                   'all_78_input_columns_per_seam': True,
                   'exact_gram_and_partial_injections': True,
                   'signed_current_predecessor_verified': True,
                   'source_bytes_unchanged': True},
        'actual_native_source_audit': audit,
        'proof': 'RESET_DOMAIN.md',
        'domain': 'full native total-two sector and arbitrary quantum references',
        'success_Kraus_form': 'a_w |target><x(w)|, since a conclusive tape resolves its original column',
        'success_effect': 'diag(s_x); lambda_min=min(s_x)',
        'bound': 'conditional first-failure union bound; no initial independence premise',
        'private_tape_boundary': 'fine tapes remain in the complete instrument; marginalizing the ideal resolved label/tape gives blank tensor original reference marginal',
        'source_sha256': pins,
        'official_rows_changed': []}
    assert pins == {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in paths}
    (HERE / 'RESET_DOMAIN.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({'checks': record['checks'], 'actual_native_source_audit': audit}))


if __name__ == '__main__':
    run()
