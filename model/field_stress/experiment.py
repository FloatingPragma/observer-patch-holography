"""Reproduce the same-field energy/entropy component, with source custody."""
from __future__ import annotations
from fractions import Fraction as F
import hashlib
from pathlib import Path
import math
import numpy as np
from ..interacting_field import source, bounds
from ..interacting_field.lattice import WeylLattice
from ..common_source.chiral_candidate import gaussian_source as gaussian
from .energy import certify, commutator, SCALE, sum_exact
from . import entropy

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def dynamics(lattice, local, currents):
    p, q = lattice.origin+1,lattice.neighbor
    initial = lattice.basis_state(lattice.origin,lattice.origin+1)
    state = source.hop(lattice,initial,p,q,-1j,math.pi/4)
    preparation = gaussian.local_gate('hopping',-1j,math.pi/4)
    swap = gaussian.local_gate('fswap')
    routing = 2*(q-p-1)
    prep_counts = {key:preparation.resources[key]+routing*swap.resources[key] for key in source.COUNT_KEYS}
    prep_counts.update({'completed_hopping_blocks':1,'routing_FSWAP_calls':routing,
                        'floating_local_residual_sum_with_routing':source.local_error(preparation)+routing*source.local_error(swap),
                        'floating_residual_estimate_outward_certified':False})
    origin_site = lattice.site_index[(0,0,0)]
    current_sum = sum_exact([current if x == origin_site else
                            commutator(local[origin_site],local[x]).times_i()
                            for (x,y),current in currents.items() if x == origin_site or y == origin_site],lattice.dimension)
    initial_outflow = float(np.vdot(state,current_sum.matrix(SCALE**2)@state).real)
    initial_energy = float(np.vdot(state,lattice.hamiltonian@state).real)
    duration = F(1,8)
    exact = lattice.exact_evolve(state,float(duration))
    exact_energy = float(np.vdot(exact,lattice.hamiltonian@exact).real)
    rows = []
    b = bounds.norm_sum(F(1,8),F(1))
    for steps in (4,16):
        changed = source.evolve(lattice,state,float(duration),steps)
        error = F(bounds.history(((duration,steps),))['complete_history_half_trace_upper'])
        resources = source.resources(lattice,float(duration),steps,segments=1)
        # Keep per-call numerical residuals elsewhere; aggregate source errors
        # are platform-sensitive and are not upgraded to rigorous error bars.
        counts = {key:value for key,value in resources.items() if isinstance(value,int) and not isinstance(value,bool)}
        energy = float(np.vdot(changed,lattice.hamiltonian@changed).real)
        rows.append({'steps':steps,'ideal_unitary_error_upper':str(error),
                     'ideal_energy_drift_upper':str(2*b*error),
                     'source_energy':energy,'source_energy_drift':abs(energy-initial_energy),
                     'source_state_to_sparse_exact_norm':float(np.linalg.norm(changed-exact)),
                     'source_norm_residual':abs(float(np.vdot(changed,changed).real)-1),
                     'requested_evolution_resources':counts,
                     'floating_local_hopping_max_residual':resources['maximum_local_hopping_target_residual'],
                     'floating_contact_max_residual':resources['contact_target_residual'],
                     'outward_certified_floating_or_hardware_error':False})
    return {'model_time_exact':str(duration),'initial_state':'source-prepared real superposition of an origin contact pair and its adjacent spin-flipped separated pair',
            'preparation_resources':prep_counts,'initial_origin_outward_energy_flux':initial_outflow,
            'initial_total_energy':initial_energy,'sparse_exact_energy_drift':abs(exact_energy-initial_energy),
            'source_runs':rows,'unitary_dynamics_only':True,
            'detector_backaction_and_external_control_work_not_silently_included':True,
            'physical_time_or_energy_calibrated':False}


def source_bindings():
    paths = [p.relative_to(ROOT).as_posix() for p in HERE.glob('*.py')]+[
        'model/field_stress/README.md',
        'model/field_regions/regions.py',
        'model/interacting_field/lattice.py','model/interacting_field/source.py',
        'model/interacting_field/bounds.py','model/interacting_field/DERIVATION.md',
        'model/common_source/chiral_candidate/gaussian_source.py',
        'model/common_source/chiral_candidate/kinetic.py',
        'model/observer/source_scalar/compiler.py',
        'model/observer/scalar_response/experiment.py',
        'Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean',
        'Lean/ObserverPatchHolography/EinsteinBranch/Composition.lean',
        'Lean/Thermodynamics/EinsteinPremiseLink.lean']
    return {path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def run():
    lattice = WeylLattice(1)
    certificate, local, currents = certify(lattice)
    free_certificate, _, free_currents = certify(WeylLattice(1,interaction=0))
    changed = sum(not (value-free_currents[key]).is_zero() if key in free_currents else not value.is_zero()
                  for key,value in currents.items())
    certificate['interaction_deletion_changed_current_pairs'] = changed
    certificate['free_nonzero_unordered_currents'] = free_certificate['nonzero_unordered_currents']
    dynamic = dynamics(lattice,local,currents)
    thermal = entropy.run(lattice)
    okay = (certificate['operator_identities_exact'] and changed > 0
            and abs(dynamic['initial_origin_outward_energy_flux']) > 1e-3
            and dynamic['sparse_exact_energy_drift'] < 2e-10
            and thermal['first_law_derivative_residual'] < 2e-10
            and thermal['exact_identity_numerical_residual'] < 2e-10
            and thermal['variation_reference_commutator_Hilbert_Schmidt'] > 1e-6
            and thermal['regional_modular_bare_energy_commutator_Hilbert_Schmidt'] > 1e-5)
    for row in dynamic['source_runs']:
        okay &= row['source_energy_drift'] < float(F(row['ideal_energy_drift_upper']))+2e-10
    return {'schema':'oph.same_field_energy_entropy.v1','all_pass':bool(okay),
            'scope':'exact finite energy continuity and same-action quantum entropy; conditional gravity premises remain open',
            'action':'H=dGamma[-(1/8) sigma.sin(k)]+sum_x n_x0 n_x1, complete N=2 sector',
            'energy_continuity':certificate,'compiled_energy_response':dynamic,
            'same_action_entropy':thermal,
            'premise_attachment':{
                'PR-23':'full-sector finite Gibbs first law on the actual action; geometric cap kernel and native full-field reference identification absent',
                'PR-24':'local model-time energy continuity; no common Lorentz stress tensor or half-sided inclusion',
                'PR-25':'no null stress/geometry tomography or universal source map',
                'PR-26':'quantum entropy tangent and exact relative-entropy remainder; no fixed-volume area or generalized-entropy stationarity',
                'PR-27':'no small geodesic balls, Bianchi geometry or cofinal controlled geometric tail',
                'PR-28':'no physical gravitational tower or independent Newton calibration'},
            'EinsteinRegisterData_inhabited':False,'integrated_gate_accepted':False,
            'physical_comparison_count_added':0,'source_bindings':source_bindings()}
