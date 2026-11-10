"""Complete periodic-source spectrum and compensated charged readout receipt."""
import hashlib
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate import action
from .basis import DIMENSION,metadata
from .spectrum import sector_certificate
from .symmetry import certificate_data as symmetry_data
from .readout import charged_spectroscopy

ROOT=Path(__file__).resolve().parents[4]


def input_hashes():
    here=Path(__file__).parent
    seeds={str(p.relative_to(ROOT)) for p in here.iterdir() if p.suffix in ('.py','.md')}
    seeds.add(str((here/'eigenbasis_witnesses.npz').relative_to(ROOT)))
    seeds.update(action.SOURCE_FILES)
    seeds.update(('model/requirements.txt',
        'model/coherent_defect_wave/chiral_source/PROOF.md',
        'model/coherent_defect_wave/chiral_higgs/PROOF.md',
        'model/coherent_defect_wave/chiral_higgs/GAUSS_APPARATUS.md',
        'model/coherent_defect_wave/chiral_joint/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/two_cell/BASIS_PROOF.md',
        'model/coherent_defect_wave/chiral_joint/two_cell/SPECTRUM_PROOF.md',
        'model/coherent_defect_wave/chiral_joint/two_cell/INSTRUMENT_PROOF.md',
        'model/coherent_defect_wave/chiral_joint/two_cell/READOUT_PROOF.md',
        'model/coherent_defect_wave/chiral_joint/two_cell/SYMMETRIES.md',
        'model/coherent_defect_wave/chiral_joint/two_cell/REFINEMENT_PROOF.md'))
    missing=sorted(p for p in seeds if not (ROOT/p).is_file())
    if missing:raise ValueError('missing periodic-source input: '+', '.join(missing))
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(set(import_closure(seeds))|seeds)}


def certificate():
    return {'schema':'oph.coherent_defect_wave.chiral_joint.periodic_spectrum.v1',
        'source':{'active_dimension':DIMENSION,'geometry':'two cells and two oriented quantum links',
            'Hamiltonian':'two_cell.symmetries.periodic_column, both quantum Higgs memories',
            'Yukawa_e':'1','Yukawa_u':'0','Yukawa_d':'0',
            'scalar_and_link_energy_coefficients':'1','hopping_magnitude':'1/4',
            'operator_norm_upper_bound':'30','clock':'exp(-i t H)',
            'energy_reference':'actual complete neutral periodic ground E0',
            'supplied_architecture_parameters_are_not_empirical_mass_calibrations':True},
        'complete_Gauss_sectors':[metadata(q,execute_kernel=True) for q in (0,-1,1)],
        'complete_physical_spectra':[sector_certificate(q) for q in (0,-1,1)],
        'same_boundary_conjugation':symmetry_data(),
        'compensated_charged_spectroscopy':charged_spectroscopy(),
        'scope':{'complete_reduced_source_physical_ground_certified':True,
            'all_occupation_sectors_included':True,
            'full_R15_family_ground_calculated':False,
            'same_source_ground_charge_energy_and_joint_phase_records':True,
            'ground_preparation_is_declared_finite_controller':True,
            'all_preparation_and_measurement_failures_retained':True,
            'charged_phase_degeneracy_comes_from_retained_reference_position':True,
            'source_only_particle_momentum_or_relativistic_dispersion_proved':False,
            'GeV_masses_or_observed_species_identified':False,
            'interacting_continuum_limit_proved':False},
        'input_sha256':input_hashes()}
