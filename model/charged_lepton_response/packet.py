"""Recompute scoped charged-field certificates with explicit calibration ancestry."""
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from . import (car_koide, phase_readout, phase_selection, relativistic_mass,
              mass_response_bridge, calibrated_sanity, p_dependence, candidate_map)

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PARENT = 'model/standard_model_source_map/particle_map.json'
PARENTS = (PARENT, candidate_map.DEFECTS,
           'model/mcpr_native_memory_bridge/results.json',
           'model/pixel_mass_closure/results.json')
OUTPUT = HERE/'results.json'
MARKDOWN = HERE/'RESULTS.md'
CANDIDATES = HERE/'candidate_map.json'


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def source_bindings():
    parent = strict_json.load(ROOT/PARENT)
    if parent['physical_particle_map'] != []:
        raise ValueError('parent physical admission or provenance requires review')
    for receipt in PARENTS:
        value = strict_json.load(ROOT/receipt)
        pins = value.get('source_bindings', value.get('source_sha256'))
        if not pins:
            raise ValueError('missing parent source manifest: '+receipt)
        for name, expected in pins.items():
            path = (ROOT/name).resolve()
            if not path.is_relative_to(ROOT) or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                raise ValueError('stale charged-field parent: '+name)
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/charged_lepton_response/PROOF.md',
        'source_receipts': list(PARENTS),
        'independent_checks': ['model/charged_lepton_response/'+name for name in
            ('README.md', 'CAR_KOIDE.md', 'PHASE_AND_PREMISES.md', 'RELATIVISTIC_MASS.md',
             'MASS_RESPONSE_BRIDGE.md', 'CALIBRATED_SANITY.md', 'SOURCE_PHASE_SELECTION.md',
             'P_DEPENDENCE.md', 'CANDIDATE_MAP.md', 'calibration_reference.json')]+
            list(p_dependence.SOURCE_INPUTS)+
            ['Lean/Screen/SeamU1HolonomyClassification.lean',
             'Lean/Screen/HolonomyInterference.lean',
             'Lean/ObserverPatchHolography/KoideCirculant.lean',
             'model/mcpr_native_memory_bridge/PROOF.md',
             'model/pixel_mass_closure/PROOF.md',
             'code/particles/leptons/derive_charged_source_law_rigidity_conditional.py',
             'code/particles/leptons/derive_charged_mcpr_completion_conditional.py']})
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


@lru_cache(maxsize=1)
def _science_json(source_key):
    calibrated = calibrated_sanity.certificate()
    pixel = p_dependence.certificate()
    candidates = candidate_map.summarize(strict_json.load(ROOT/PARENT), calibrated,
        pixel, strict_json.load(ROOT/candidate_map.DEFECTS))
    return canonical({'full_CAR_realization': car_koide.certificate(),
        'phase_and_action_readout': phase_readout.certificate(),
        'source_phase_selection_audit': phase_selection.certificate(),
        'free_relativistic_field': relativistic_mass.certificate(),
        'conditional_mass_response': mass_response_bridge.certificate(),
        'calibrated_sanity_test': calibrated,
        'conditional_P_dependence': pixel,
        'candidate_particle_map': candidates})


def build():
    pins = source_bindings()
    science = json.loads(_science_json(tuple(pins.items())))
    value = {'schema': 'oph.charged_lepton_response.v1',
        'source_bindings': pins, 'certificates': science,
        'source_selected_phase': None, 'physical_particle_map': [],
        'physical_mass_predictions': [],
        'calibrated_action_comparison_is_independent_prediction': False,
        'parent_source_bindings_authenticated': True,
        'parent_scientific_experiments_cold_replayed_here': False}
    if source_bindings() != pins:
        raise RuntimeError('source changed during charged-field export')
    return value


def validate(value):
    if canonical(value) != canonical(build()):
        raise ValueError('charged-field evidence, scope or provenance changed')


def render(value):
    from fractions import Fraction as F
    calibrated = value['certificates']['calibrated_sanity_test']
    ratios = calibrated['comparison']['source_action_ratios']
    reference = calibrated['reference']
    def reference_text(name):
        row = reference[name]
        return row['value']+' ± '+row['standard_uncertainty']
    tau_inside = calibrated['comparison']['held_out_tau_ratio_inside_reported_one_sigma_window']
    def middle(pair): return float(sum(map(F, pair))/2)
    lines = ['# Charged-field readout results', '',
        'The configured coherent program has a full-CAR realization and a controlled free Dirac limit.',
        'All lattice valleys and quantum archives are retained. The interacting chiral limit and measured-mass attachment are unproved.', '',
        'The phase below is calibrated to the muon/electron ratio. Tau is excluded from that calibration.', '',
        '| Quantity | Action result | Experimental reference | Role |',
        '| --- | ---: | ---: | --- |',
        f"| muon/electron | {middle(ratios[1]):.12f} | {reference_text('muon_electron_ratio')} | Calibration |",
        f"| tau/electron | {middle(ratios[2]):.12f} | {reference_text('tau_electron_ratio')} | Conditional comparison |", '',
        f"Reference: [{reference['edition']}]({reference['source']}).",
        ('The held-out tau ratio lies inside this reported one-standard-uncertainty interval.'
         if tau_inside else
         'The held-out tau ratio does not lie wholly inside this reported one-standard-uncertainty interval.'),
        'No joint statistical significance or actual interaction-correction bound is claimed.', '',
        'The exact rational phase, coefficient intervals and input ancestry are in [results.json](results.json).',
        'The [calibration proof](CALIBRATED_SANITY.md) and [dependency map](PROOF.md) state the physical requirements.',
        'The [source-phase audit](SOURCE_PHASE_SELECTION.md) distinguishes C/E duties from orientation reversal and the conditional pixel law.',
        'Physically identified particles: **none**. Source-selected physical phase: **unproved**.', '']
    lines += ['## Conditional C-duty spectrum', '',
        'The complete memory construction realizes the supplied pixel-dependent law.',
        'Its four coordinate branches give the following ratios. Each lightest coefficient is normalized to one.', '',
        '| Declared branch | Middle/light | Heavy/light |',
        '| --- | ---: | ---: |']
    for row in value['certificates']['conditional_P_dependence']['branches']:
        ratios = row['source']['charged_action_ratios']
        lines.append(f"| {row['branch']} | {middle(ratios[1]):.12f} | {middle(ratios[2]):.12f} |")
    inside = [row['branch'] for row in value['certificates']['conditional_P_dependence']['branches']
              if all(row['comparison'][name]['source_interval_inside_reference_one_sigma']
                     for name in ('muon_electron_ratio', 'tau_electron_ratio'))]
    comparison_text = ('Only the measured-alpha endpoint lies inside both quoted individual ratio intervals.'
                       if inside == ['measured_endpoint'] else
                       'Branches wholly inside both quoted individual ratio intervals: '
                       +(', '.join(inside) if inside else 'none')+'.')
    lines += ['', comparison_text,
        'The measured-alpha endpoint inherits a measured fine-structure input. Other branches retain their specified coordinates or conditional census and scale law.',
        'The [capacity proof](P_DEPENDENCE.md) gives interval bounds and input ancestry.', '',
        'The [candidate map](CANDIDATE_MAP.md) covers eighteen target channels and preserves the separate eleven-entry protected-object inventory.',
        'Its [machine-readable map](candidate_map.json) assigns charged spectral modes as hypotheses; physical masses and defect-to-species links are unassigned.', '']
    return '\n'.join(lines)
