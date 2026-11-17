"""Retained complete-field clock execution and independent controls."""
from __future__ import annotations

from fractions import Fraction as F
import ast
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from . import clock
from ..interacting_field import source
from ..observer.source_scalar import compiler
from ..local_repair_family import inference

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
COUNTS='model/local_repair_family/record_clock_results.json'


def retained_counts():
    """Only complete seam counts cross into the rate decoder."""
    raw=json.loads((ROOT/COUNTS).read_text())
    counts={tuple(b['seam']):(b['trials'],b['coin1_count']) for b in raw['acquisition']['blocks']}
    inference._counts(counts)
    return counts


def _probabilities(outputs):
    return [float(np.vdot(v,v).real) for v in outputs]


def _reference_density(state):
    # Columns are an untouched reference basis; this is Tr_field |psi><psi|.
    return state.T@state.conj()


def source_bindings():
    paths={p.relative_to(ROOT).as_posix() for p in HERE.glob('*.py')}
    paths.add('model/band_clock/README.md')
    paths.update(compiler.source_summary()['source_sha256'])
    paths.update(inference.source_bindings())
    paths.update((COUNTS,'model/particle_clock/clock.py',
        'model/bound_pair/fiber.py','model/bound_pair/protection.py','model/bound_pair/spectrum.py',
        'model/interacting_field/source.py','model/interacting_field/lattice.py',
        'model/interacting_field/bounds.py','model/interacting_field/experiment.py',
        'model/interacting_field/replay.py',
        'model/common_source/chiral_candidate/gaussian_source.py',
        'model/common_source/mass_response/model.py','model/observer/source_scalar/compiler.py'))
    # Python executes parent package initializers and imported local modules.
    # Follow both, including imports inside functions; archive/dynamic sources
    # remain supplied by the existing compiler manifest above.
    pending = list(paths)
    inspected = set()
    while pending:
        name = pending.pop()
        if name in inspected:
            continue
        inspected.add(name)
        path = ROOT/name
        discovered = set()
        for parent in path.relative_to(ROOT).parents:
            initializer = parent/'__init__.py'
            if (ROOT/initializer).is_file():
                discovered.add(initializer.as_posix())
        if path.suffix == '.py':
            module = list(path.relative_to(ROOT).with_suffix('').parts)
            package = module if path.name == '__init__.py' else module[:-1]
            for node in ast.walk(ast.parse(path.read_text())):
                modules = []
                if isinstance(node, ast.Import):
                    modules = [alias.name.split('.') for alias in node.names]
                elif isinstance(node, ast.ImportFrom):
                    base = (package[:len(package)-node.level+1] if node.level else [])
                    base += node.module.split('.') if node.module else []
                    modules = [base]+[base+alias.name.split('.') for alias in node.names
                                      if alias.name != '*']
                for imported in modules:
                    if imported and imported[0] in ('model', 'tools'):
                        candidates = (Path(*imported).with_suffix('.py'),
                                      Path(*imported)/'__init__.py')
                        discovered.update(p.as_posix() for p in candidates if (ROOT/p).is_file())
        for item in discovered-paths:
            paths.add(item)
            pending.append(item)
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def run():
    counts=retained_counts(); space=clock.ClockField()
    initial,band,energies,coefficients=clock.packet_with_reference()
    decision,outputs=clock.execute(counts,initial,steps=32,space=space)
    if decision['status']!='accept':
        raise ArithmeticError('retained complete source counts failed the fixed rate gate')
    duration=decision['Hamiltonian_model_duration']; beta=math.log(25/8)
    executed=source.evolve(space,initial,duration,32)
    exact_learned=clock.exact_state(space,initial,duration)
    exact_true=clock.exact_state(space,initial,beta)
    band_exact=band@(coefficients*np.exp(-1j*duration*energies[:,None]))
    pointer_reference=clock.pointers(space,exact_true)
    pointer_initial=clock.pointers(space,initial)
    output_values=tuple(outputs.values())
    # The discarded-contact control acts on the SAME full input and Hilbert space.
    free=clock.ClockField(interaction=0.)
    free_state=clock.exact_state(free,initial,duration)
    free_outputs=clock.pointers(free,free_state)
    before=np.sum(abs(initial[space.contact_indices])**2,axis=1)
    after=np.sum(abs(executed[space.contact_indices])**2,axis=1)
    separated0=float(np.sum(abs(initial[space.relative_distance>0])**2))
    separated1=float(np.sum(abs(executed[space.relative_distance>0])**2))
    motion_tv=(float(np.sum(abs(after-before)))+abs(separated1-separated0))/2
    resources=source.resources(space,duration,32,segments=1)
    _,pointer_info=source.occupation()
    for key in source.COUNT_KEYS:
        resources[key]+=2*pointer_info['resources'][key]
    resources['floating_local_residual_sum_with_routing']+=2*pointer_info['source_target_residual']
    resources['floating_residual_telescoping_estimate']=math.expm1(resources['floating_local_residual_sum_with_routing'])
    resources.update({'fresh_pointer_blanks':2,'committed_pointer_bits':2,
                      'field_mode_register_size':space.modes,'untouched_reference_dimension':2,
                      'bound_packet_preparation_cost_included':False,
                      'field_Gibbs_preparation_cost_included':False,
                      'record_acquisition_attempts_separately_retained':30*(1<<20),
                      'scope':'worst-branch full-field evolution plus two occupation pointers'})
    reduced=sum((_reference_density(v) for v in output_values),np.zeros((2,2),complex))
    cq=clock.cq_distance(output_values,pointer_reference)
    budget=clock.error_budget()
    if (cq>float(F(budget['complete_average_halfdiamond_upper']))
            or np.linalg.norm(exact_learned-band_exact)>2e-12
            or abs(sum(_probabilities(output_values))-1)>2e-8
            or motion_tv<1e-4):
        raise ArithmeticError('complete moving-band source history failed its fixed gates')
    return {'schema':'oph.band_clock.v1','source_bindings':source_bindings(),
            'scope':'SUPPLIED_FULL_FIELD_REFERENCE_AND_RECORD_RATE_TRANSFER_SERVICE',
            'action':{'same_periodic_Weyl_contact_action':True,'side':4,
                'full_N2_dimension':space.dimension,'bound_band_dimension':64,
                'hopping_matchings':len(space.groups),'hopping_bonds':sum(map(len,space.groups)),
                'coordinate_speed':'1/8','repulsion':'1',
                'scattering_states_and_all_eight_Weyl_tastes_retained':True,
                'band_Hamiltonian_fitted_instead_of_local_action':False},
            'rate_decision':decision,'uniform_error_budget':budget,
            'common_reference':clock.reference_report(),
            'source_execution':{
                'initial_field_reference_norm':float(np.vdot(initial,initial).real),
                'source_norm_residual':abs(float(np.vdot(executed,executed).real)-1),
                'source_vs_same_rate_exact_vector_norm':float(np.linalg.norm(executed-exact_learned)),
                'independent_position_vs_momentum_evolution_norm':float(np.linalg.norm(exact_learned-band_exact)),
                'complete_CQ_half_trace_vs_true_rate':cq,
                'complete_pointer_probabilities':_probabilities(output_values),
                'true_rate_pointer_probabilities':_probabilities(pointer_reference),
                'untouched_reference_marginal_residual':float(np.linalg.norm(reduced-_reference_density(initial))),
                'complete_contact_site_plus_separated_motion_TV':motion_tv,
                'contact_site_probabilities_before':before.tolist(),
                'contact_site_probabilities_after':after.tolist(),
                'separated_probability_before':separated0,'separated_probability_after':separated1,
                'bound_band_leakage_from_split':float(np.linalg.norm(executed-band@(band.conj().T@executed))**2),
                'zero_clock_complete_CQ_contrast':clock.cq_distance(output_values,pointer_initial),
                'deleted_contact_complete_CQ_contrast':clock.cq_distance(output_values,free_outputs),
                'all_four_pointer_outcomes_retained':True,
                'proper_code_export_failure_slots':'zero for ideal invariant-code words; numerical local residuals are diagnostic',
                'requested_resources':resources},
            'unfulfilled_joins':{
                'native_product_reference_intertwines_full_field':False,
                'new_full_Gibbs_reference_prepared_or_selected':False,
                'bound_packet_source_preparation_executed':False,
                'native_repair_is_Hamiltonian_flow':False,
                'physical_time_energy_mass_calibration':False,
                'continuum_chiral_QFT_or_physical_species_identified':False,
                'CP04_or_CP06_accepted':False}}
