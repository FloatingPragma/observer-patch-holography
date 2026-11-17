"""Bind the partial OpenWave results on native records to their evidence."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
try:
    from tools.pin_scope import bound_pins
except ImportError:  # run or imported from tools/ as a top-level module
    from pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / 'code/openwave_partial_records'
OUTPUT = ROOT / 'claims/openwave_partial_records.json'
PROOF = 'code/openwave_partial_records/PROOF.md'

ROWS = [
    {'row': 3, 'criterion': 'de Broglie clock (Zitterbewegung)',
     'anchor': 'row-3-de-broglie-clock',
     'providers': ['compact_internal_clock', 'compact_internal_clock_root_review'],
     'result': ('A native internal reference on the compact charged state of the sealed 240-port source: '
                'the walk W = S(2AA*-I) built from the actual native event arcs rotates span{u,w} by '
                'theta = acos(209/210) = 0.0976287750 per reflection word, and the native triangle cycle v '
                '(Wv = v) makes an observable theta harmonic in ordinary coordinate records of (u+v)/sqrt2, '
                'with all 524,288 source and coordinate outcomes retained (303,660 native calls).'),
     'scope': ('The reflection walk and the coherent preparation are declared memory services on the '
               'sealed source, and theta is measured in reflection words.')},
    {'row': 4, 'criterion': 'Particle stability (Derrick escape)',
     'anchor': 'row-4-particle-stability',
     'providers': ['compact_internal_clock', 'compact_internal_clock_root_review'],
     'result': ('The exact invariant span{u, Su, v} of the same walk stays in one twelve-port cell '
                '(ports 0 to 11) at every integer time, by an exact integer identity.'),
     'scope': 'The statement covers the invariant span of the declared reflection walk on the sealed source.'},
    {'row': 6, 'criterion': 'Angular momentum J (spin hbar/2)',
     'anchor': 'row-6-spin',
     'providers': ['compact_field_spin_join__existing_audit', 'compact_field_spin_root_review'],
     'result': ('The exact CAR field embedding J|p,s> = c^dagger_(p,s)|vacuum> satisfies '
                'V_h J = J(K_h tensor I_2) for every seam, coin and current of the current-only record law '
                '(crossing sign +1, 604,800 columns per law, 302,400 native calls); it carries a conserved '
                'internal field spin together with the compact energy certificate E6 = 6.014331228348673 '
                'on the same state.'),
     'scope': ('The embedding acts on the positive one-particle sector of the C record law, and the spin '
               'is the internal CAR field spin of that sector.')},
    {'row': 8, 'criterion': 'Antimatter + annihilation',
     'anchor': 'row-8-annihilation',
     'providers': ['neutral_source_completion', 'neutral_boundary_root_review', 'neutral_radiation_join'],
     'result': ('A signed pair (j,-j) meets and is removed by native contact with current j under both '
                'coins on all 2,520 contact branches; the adjacent contact state has a complete 17-outcome '
                'energy pointer with mean contact energy 11.051928582384567 in transfer units (exact '
                'enclosure), and the contact current releases a record field E = 2j, B = 0 that propagates '
                'by the exact discrete wave equation with local flux continuity (360 of 360 contact flights '
                'launched and absorbed, 147,840 wave equalities).'),
     'scope': ('The contact energy and the ring field norm are each measured on their own register and '
               'in their own units.')},
    {'row': 20, 'criterion': 'Electric force (Coulomb 1/r, electrostatics)',
     'anchor': 'row-20-electric-force',
     'providers': ['row20_residence'],
     'result': ('Residence records of one unit, reconstructed from its signed native currents on the '
                'declared degree-six cubic routing, give the far-field 1/r form: a + b mean(r^-p) over '
                'shells 3 to 14 fits p = 1.035144 with approximate 95% profile interval '
                '[0.983376, 1.086373] (5,040 trajectories, 13,213,995 currents, no censoring).'),
     'scope': 'The form is measured on the declared routing, and the strength alpha is a separate question.'},
    {'row': 31, 'criterion': 'Orbital quantization (atomic structure)',
     'anchor': 'row-31-orbital-levels',
     'providers': ['row31_probe_binding', 'row31_coulomb_limit', 'unified_wave_submission'],
     'result': ('Native first-hit feedback (hit probabilities 5/22, 1/11 and 1/22 from the six-neighbour '
                'native walk) creates an attractive internal-band well whose lowest squared energies '
                '2.2720093211, 4.2188978370 and 4.6123309041 lie below the band thresholds 5-sqrt5, 6 and '
                '5+sqrt5 (305,157 native calls, byte-identical replay); its controlled limit gives '
                'hydrogen-like standing shells e_n = -(lambda c_*)^2/(4n^2), c_* = 1/(4 pi G(0)), with n^2 '
                'orbital states per shell.'),
     'scope': ('Routing, target and feedback are declared, and the shell spectrum is a controlled-limit '
               'theorem with a finite native calibration of the well.')},
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate():
    spec = importlib.util.spec_from_file_location('partial_records_verify', PACKAGE / 'verify.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    verification = module.verify()
    rows = []
    for row in ROWS:
        entry = {k: v for k, v in row.items() if k != 'anchor'}
        entry.update(status='partial', official_status='untested', proof=f"{PROOF}#{row['anchor']}",
                     sealed=verification['sealed'][f"row_{row['row']}"])
        rows.append(entry)
    return {
        'schema': 'oph.openwave.partial_records.v1',
        'scope': ('Partial OpenWave results carried by native records. Each row states the positive '
                  'finite result and one scope sentence; every row is partial and officially untested.'),
        'rows': rows,
        'evidence': {'proof': PROOF, 'readme': 'code/openwave_partial_records/README.md',
                     'bindings': 'code/openwave_partial_records/EVIDENCE.json',
                     'recompute': 'python3 code/openwave_partial_records/verify.py',
                     'archived_files_verified': verification['archived_files_verified'],
                     'providers': verification['providers'],
                     'recomputed': verification['recomputed']},
        'official_OpenWave_awards_added': 0,
        'source_sha256': bound_pins({str(p.relative_to(ROOT)): sha(p) for p in
                          [Path(__file__), PACKAGE / 'verify.py', PACKAGE / 'recompute.py',
                           PACKAGE / 'PROOF.md', PACKAGE / 'EVIDENCE.json',
                           ROOT / 'code/openwave_static_signed_source/evidence/openwave/shared/'
                                  'localized_rest_join/evidence/graph.json']}),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    text = json.dumps(generate(), indent=2, sort_keys=True) + '\n'
    if args.write:
        OUTPUT.write_text(text)
    elif not OUTPUT.is_file() or OUTPUT.read_text() != text:
        raise SystemExit('Stale OpenWave partial-records claim ledger')
    print('PASS_OPENWAVE_PARTIAL_RECORDS_CLAIM')


if __name__ == '__main__':
    main()
