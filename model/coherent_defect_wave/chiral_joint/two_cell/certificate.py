"""Replay all six finite-source checks with their explicit physical boundary."""
import hashlib
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate import action
from .basis import complete_metadata
from .hamiltonian import DIMENSION,PARAMETERS,structural_checks
from .spectrum import sector_certificate
from .instrument import certificate_data as instrument_data
from .readout import complete_experiment,unconditional_experiment,reflected_experiments
from .symmetries import certificate as symmetry_data

ROOT=Path(__file__).resolve().parents[4]
LEAN_INPUTS=('Lean/QFT/EffectiveQuantumComparison.lean',
    'Lean/QFT/PeriodicMomentum.lean','Lean/QFT/FiniteUnitaryScatteringNoGo.lean',
    'Lean/Geometry/MassShellKinematics.lean','Lean/Geometry/LorentzOverlapCocycle.lean',
    'Lean/Geometry/ObserverRestSpace.lean','Lean/Geometry/ObserverFrameHyperboloid.lean',
    'Lean/Geometry/CelestialNullCone.lean','Lean/Geometry/CanonicalLorentzModule.lean',
    'Lean/lean-toolchain','Lean/lakefile.lean','Lean/lake-manifest.json')


def input_hashes():
    here=Path(__file__).parent
    seeds={str(p.relative_to(ROOT)) for p in here.iterdir()
           if p.suffix in ('.py','.md')}
    seeds.add(str((here/'eigenbasis_witnesses.npz').relative_to(ROOT)))
    seeds.update(action.SOURCE_FILES)
    seeds.update(('model/requirements.txt',
        'model/coherent_defect_wave/chiral_source/PROOF.md',
        'model/coherent_defect_wave/chiral_higgs/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/PROOF.md',
        'paper/tex_fragments/EFFECTIVE_QUANTUM_COMPARISON.tex'))
    seeds.update(LEAN_INPUTS)
    missing=sorted(name for name in seeds if not (ROOT/name).is_file())
    if missing:raise ValueError('missing two-cell source dependency: '+', '.join(missing))
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(set(import_closure(seeds))|seeds)}


def certificate():
    from .candidate_data import certificate_data as candidate_data
    from .refinement import certificate_data as refinement_data
    from .propagation import certificate_data as propagation_data
    return {'schema':'oph.coherent_defect_wave.chiral_joint.two_cell.v1',
        'source':{'active_dimension':DIMENSION,'parameters':PARAMETERS,
                  'state':'complete neutral Gauss ground, with retained preparation reference',
                  'clock':'the same uniform open two-cell Hamiltonian',
                  'comparison_source':'explicit separate two-link periodic completion'},
        'Gauss_projection':[complete_metadata(q) for q in (0,-1,1,(1,0))],
        'Hamiltonian':structural_checks(),
        'full_spectra':[sector_certificate(q) for q in (0,-1,1)],
        'instrument':instrument_data(),
        'conditional_readout':complete_experiment(),
        'original_attempt_readout':unconditional_experiment(),
        'reflected_experiment':reflected_experiments(),
        'symmetries':symmetry_data(),
        'propagation':propagation_data(),
        'candidates':candidate_data(),
        'refinement':refinement_data(),
        'scope':{'entire_reduced_source_projected':True,
                 'same_H_ground_across_every_conserved_sector':True,
                 'complete_compensated_energy_records':True,
                 'physical_rest_mass_or_GeV_identification':False,
                 'interacting_relativistic_limit_proved':False,
                 'observed_species_identified':False,
                 'finite_gap_is_Clay_mass_gap':False},
        'input_sha256':input_hashes()}
