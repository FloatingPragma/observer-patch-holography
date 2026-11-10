"""Reproducible custody for the finite cubic CAR readout."""
import hashlib
from fractions import Fraction as F
from pathlib import Path

from model.particle_common_contract.provenance import import_closure
from .lattice import State,packet,tick,translate,plus,ZERO
from .readout import experiment,higher_momentum_experiment

ROOT=Path(__file__).resolve().parents[3]


def input_hashes():
    here=Path(__file__).parent
    paths={p for p in here.glob('*.py') if not p.name.startswith('test_')}
    paths.update(here/name for name in ('PROOF.md','README.md'))
    paths.update(ROOT/name for name in (
        'model/requirements.txt','model/particle_common_contract/provenance.py',
        'model/coherent_defect_wave/species/PROOF.md',
        'model/coherent_defect_wave/species/spatial_fock.py',
        'model/coherent_defect_wave/spatial/PROOF.md',
        'evidence/native_geometric_source_20260925/spec.json',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py'))
    names=import_closure({str(p.relative_to(ROOT)) for p in paths})
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(names) if (ROOT/name).is_file()}


def geometric_certificate():
    initial=packet(); advanced=tick(initial,4)
    for axis in range(3):
        displacement=tuple(2 if i==axis else 0 for i in range(3))
        if tick(translate(initial,4,displacement),4)!=translate(advanced,4,displacement):
            raise ArithmeticError('full cubic tick does not commute with its supercell translation')
    # Whole-cell momentum must translate the link field as well as the occupation.
    full=translate(advanced,4,(2,0,0))
    bare=State({(((p[0]+2)%4,p[1],p[2]),i,field):v
                for (p,i,field),v in advanced.amplitudes},advanced.norm_denominator)
    difference=full.as_dict()
    for key,(a,b) in bare.amplitudes:
        difference[key]=plus(difference.get(key,ZERO),(-a,-b))
    witness=State(difference,full.norm_denominator).norm_squared()
    if witness!=F(1,225):
        raise ArithmeticError('the bare-wave translation countercontrol changed')
    small=tick(packet(origin=(2,2,2)),8)
    large=tick(packet(origin=(4,4,4)),12)
    if translate(large,12,(-2,-2,-2))!=small:
        raise ArithmeticError('local source execution depends on a remote boundary')
    return {'complete_cell_translations_checked_axes':[0,1,2],
            'translation_stride_cells':2,
            'bare_wave_vs_complete_translation_squared_difference':'1/225',
            'fixed_spacing_one_tick_volume_comparison':{
                'cube_sides':[8,12],'seed_origins':[[2,2,2],[4,4,4]],
                'translated_complete_sparse_vectors_equal':True,
                'output_amplitude_support':len(small.amplitudes)},
            'causal_bounds_per_tick':{'maximum_coordinate_displacement':2,
                                       'maximum_Manhattan_displacement':6},
            'continuum_limit_inferred_from_fixed_spacing_test':False}


def certificate():
    from .lattice import local_column,mode
    index=mode(7,3,0,0)
    diagonal=dict(local_column(1,index))[1,index]
    if diagonal!=F(-29,30):
        raise ArithmeticError('source diagonal used by the Bloch expectation formula changed')
    return {'schema':'oph.coherent_defect_wave.cubic_readout.v1',
            'program':'parity-controlled spatial CAR law on the cubic matching circuit',
            'clock':'one onsite U_F layer followed by six root-SWAP hopping layers',
            'reference':'empty-wave product of positive-Gauss face-0 link pairs',
            'reference_is_unique_physical_vacuum':False,
            'spatial_interpretation':{'wave_occupation_moves_between_cells':True,
                                      'defect_face_can_change_within_its_carrier':True,
                                      'regional_topological_charge_moves_between_cells':False,
                                      'carrier_face_and_spatial_cell_are_distinct_positions':True},
            'local_encoding':{'gauge_orbit_count':5,'wave_modes_per_cell':7200,
                              'physical_one_wave_local_dimension':36000,
                              'no_link_register_is_discarded':True,
                              'five_to_power_volume_stationary_empty_wave_basis_states':True},
            'seeded_native_U_diagonal':str(diagonal),
            'experiment':experiment(),'higher_momentum_experiment':higher_momentum_experiment(),
            'geometry':geometric_certificate(),
            'input_sha256':input_hashes()}
