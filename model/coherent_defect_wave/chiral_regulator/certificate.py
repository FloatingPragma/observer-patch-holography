"""Exact spectral source checks and reproducible evidence custody."""
import hashlib
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate.action import SOURCE_FILES
from .model import (PROGRAM_ID,spectral_certificate,position_certificate,
                    filled_reference_certificate,instrument_certificate,tail_bound)

ROOT=Path(__file__).resolve().parents[3]


def certificate():
    folder=Path(__file__).parent
    paths={str(p.relative_to(ROOT)) for p in folder.glob('*.py') if not p.name.startswith('test_')}
    paths.update(SOURCE_FILES)
    paths.update(('model/coherent_defect_wave/chiral_regulator/PROOF.md',
                  'model/coherent_defect_wave/chiral_source/PROOF.md',
                  'model/coherent_defect_wave/field_identity_audit/PROOF.md',
                  'model/coherent_defect_wave/mass_instrument/PROOF.md',
                  'model/standard_model_source_map/KINETIC_ATTACHMENT.md',
                  'model/common_source/chiral_candidate/REGULATOR_ROUTES.md',
                  'model/interacting_extension/chiral_lattice/PROOF.md',
                  'Lean/Screen/WeylYukawaConventions.lean','Lean/QFT/PeriodicMomentum.lean',
                  'model/requirements.txt'))
    missing=[p for p in paths if not (ROOT/p).is_file()]
    if missing:raise FileNotFoundError(missing)
    hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
            for p in sorted(import_closure(paths)) if (ROOT/p).is_file()}
    return {'schema':'oph.coherent_defect_wave.chiral_regulator.v1','program':PROGRAM_ID,
            'symbol_and_chirality':spectral_certificate(),'position_source':position_certificate(),
            'reference':filled_reference_certificate(),'complete_free_instrument':instrument_certificate(),
            'free_continuum':{'retained_window_symbol_and_covariance_error':'0',
                 'field_tail_bound':'(1+Lambda_N^2)^(-s/2)*norm_Hs(f)',
                 'Lambda_N':'pi*(N+1)/L','exact_tail_control_at_cutoff3_order2_norm10':str(tail_bound(3,2,10)),
                 'continuum_rest_mass':'0','all_field_word_convergence':True,
                 'finite_nonlocality_retained':True},
            'scope':{'full_finite_spectral_left_Weyl_source':True,'all_cutoff_modes_retained':True,
                 'positive_free_particle_hole_reference':True,'controlled_free_local_continuum':True,
                 'finite_strict_causal_cone':False,'local_interacting_gauge_continuum':False,
                 'physical_Standard_Model_identification':False,'GeV_mass_prediction':False},
            'input_sha256':hashes}
