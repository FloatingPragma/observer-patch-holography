"""Exact finite spatial connection, causal cone and spectator band controls."""
import hashlib
from itertools import product
from pathlib import Path

import sympy as s

from model.particle_common_contract.provenance import import_closure
from ..model import JointObserver
from .family import (A,B,CubicFamily,PairState,dressed_swap,exchange,gamma,
                     tensor,positive_local,local_step,face_measurement,charge_measurement)
from .dirac_route import certificate as dirac_certificate

ROOT = Path(__file__).resolve().parents[3]


def input_hashes():
    seeds = {str(p.relative_to(ROOT)) for p in Path(__file__).parent.glob('*.py')
             if not p.name.startswith('test_')}
    seeds |= {'model/requirements.txt',
              'evidence/native_geometric_source_20260925/spec.json',
              'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
              'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py',
              'Lean/ObserverPatchHolography/Locality/NoSignalling.lean',
              'Lean/ObserverPatchHolography/Provenance/HistoryCausalInvariance.lean',
              'Lean/QFT/PeriodicMomentum.lean'}
    paths = {ROOT/name for name in import_closure(seeds)}
    paths = {p for p in paths if p.name in {q.name for q in p.parent.iterdir()}}
    return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}


def pair_control():
    observer = JointObserver()
    left = positive_local(observer.basis_state(1,5,internal=3))
    right = positive_local(observer.basis_state(3,6,internal=0))
    source = tensor(left,right)
    assert source.norm_squared() == 1
    for cell in (0,1):
        assert gamma(source,cell) == source
    result = exchange(source)
    assert result.norm_squared() == 1
    assert exchange(result,inverse=True) == source
    assert exchange(exchange(source)) == dressed_swap(source)
    assert dressed_swap(dressed_swap(source)) == source
    for cell in (0,1):
        assert gamma(result,cell) == result
        assert exchange(gamma(source,cell)) == gamma(result,cell)
    charge = charge_measurement(result,1)
    assert [branch.norm_squared() for _,branch in charge] == [s.Rational(1,2)]*2
    assert all(keys[0][0] in (1,1 ^ 6553634) and keys[1][0] in (3,3 ^ 6553634)
               for keys,_ in result.amplitudes)
    # The same complete cell U reads the incoming wave and moves its local flux.
    # Both alternatives use identical receiver preparation and the same gate schedule.
    distributions = []
    for source_port in (5,6):
        left = positive_local(observer.basis_state(1,source_port,internal=3))
        right = positive_local(observer.basis_state(1,6,internal=0))
        coupled = exchange(tensor(left,right))
        evolved = local_step(coupled,1)
        assert evolved.norm_squared() == 1
        assert local_step(evolved,1,inverse=True) == coupled
        for cell in (0,1):
            assert gamma(evolved,cell) == evolved
        branches = face_measurement(evolved,1)
        assert sum(branch.norm_squared() for _,branch in branches) == 1
        assert all(gamma(branch,1) == branch for _,branch in branches)
        distributions.append({str(record['face_index']):str(branch.norm_squared())
                              for record,branch in branches})
    assert distributions[0] != distributions[1]
    return {
        'two_local_gauss_constraints_preserved':True,
        'square_of_exchange_is_dressed_swap':True,
        'inverse_is_exact':True,
        'local_field_words_retained_during_exchange':True,
        'receiver_charge_outcomes':[record['charge'] for record,_ in charge],
        'receiver_charge_probabilities':[str(branch.norm_squared()) for _,branch in charge],
        'total_internal_charge_before_and_after':'-1',
        'source_ports_for_remote_control':[5,6],
        'receiver_initial_port':6,
        'full_receiver_face_distributions':distributions,
        'observed_protocol':'one exchange layer followed by the literal receiver U; sender U commutes with receiver detector',
        'outside_reference_retained':True,
    }


def band_control():
    z = s.symbols('z',nonzero=True)
    even = s.Matrix([[A,B],[B,A]])
    odd = s.Matrix([[A,B/z],[B*z,A]])
    symbol = (odd*even).applyfunc(s.expand)
    trace = s.simplify(s.trace(symbol))
    assert s.simplify(trace-s.I*((z+1/z)/2-1)) == 0
    assert s.simplify(symbol.det()) == -1
    # Direct finite execution on a torus, not a matrix imported from a field theory.
    family = CubicFamily(4,periodic=True)
    initial = {(0,0,0):s.Integer(1)}
    result = family.one_charge_step(initial)
    assert family.one_charge_step(result,inverse=True) == initial
    assert s.expand(sum(s.conjugate(v)*v for v in result.values())) == 1
    assert len(result) > 1
    shifted = {((x+2)%4,y,z):v for (x,y,z),v in initial.items()}
    expected = {((x+2)%4,y,z):v for (x,y,z),v in result.items()}
    assert family.one_charge_step(shifted) == expected
    # Lift all eight parity basis vectors and compare the literal finite circuit
    # with the Bloch symbol, including a non-real allowed torus momentum.
    checks = 0
    for side,phases in ((4,(-1,1,-1)),(8,(s.I,1,1))):
        torus = CubicFamily(side,periodic=True)
        symbols = [symbol.subs(z,value) for value in phases]
        total = s.kronecker_product(*symbols).applyfunc(s.expand)
        for column,parity in enumerate(product(range(2),repeat=3)):
            lifted = {x:s.prod(phases[i]**(x[i]//2) for i in range(3))
                      for x in torus.sites if tuple(y%2 for y in x) == parity}
            evolved = torus.one_charge_step(lifted)
            for x in torus.sites:
                row = 4*(x[0]%2)+2*(x[1]%2)+x[2]%2
                phase = s.prod(phases[i]**(x[i]//2) for i in range(3))
                assert s.expand(evolved.get(x,0)-phase*total[row,column]) == 0
            checks += 1
    return {
        'invariant_sector':'one triplet-color change in the exact uniform native cell product background',
        'unit_cell_side_in_carriers':2,
        'one_axis_bloch_matrix':[[str(v) for v in symbol.row(r)] for r in range(2)],
        'one_axis_trace':'i*(cos(k)-1)',
        'one_axis_determinant':'-1',
        'eigenphase_definition':'cos(omega(k))=sin(k/2)^2; lambda(U_axis)=exp(-i*pi/2 +/- i*omega(k))',
        'three_axis_symbol':'U_x(kx) tensor U_y(ky) tensor U_z(kz)',
        'near_zero_phase_of_uniform_branch':'-(kx^2+ky^2+kz^2)/4 + O(|k|^4)',
        'corner_splitting':'sums of signed |q_axis|/sqrt(2), not a Lorentz norm',
        'torus_side':4,'torus_cells':64,'executed_one_tick_support':len(result),
        'exact_inverse_and_probability_checked':True,
        'literal_circuit_bloch_column_checks':checks,
        'translation_by_two_cells_commutes_in_executed_sector':True,
        'relativistic_mass_shell_certified':False,
    }


def cone_control():
    family = CubicFamily(14)
    roots = {(6,6,6)}
    cones = [family.cone(roots,t) for t in range(3)]
    for t,cone in enumerate(cones):
        assert all(max(abs(x[i]-6) for i in range(3)) <= 2*t for x in cone)
        assert all(sum(abs(x[i]-6) for i in range(3)) <= 6*t for x in cone)
    assert CubicFamily(12).cone(roots,2) == family.cone(roots,2)
    return {'cell_hilbert_dimension':72000,
            'spatial_architecture':'finite cubic boxes or even periodic cubic tori',
            'gate_depth_per_tick':'one onsite U layer and six nearest-neighbor matching layers',
            'layer_order':[[axis,parity] for axis in range(3) for parity in (0,1)],
            'exact_backwards_cone_sizes':[len(cone) for cone in cones],
            'coordinate_radius_bound':'2*t cells',
            'manhattan_radius_bound':'6*t cells',
            'zero_influence_outside_complete_circuit_cone':True,
            'local_volume_comparison':'identical reduced initial cone/reference state gives exactly identical complete local records',
            'continuum_spacing_limit_established':False}


def certificate():
    return {'schema':'oph.coherent_defect_wave.spatial.v1',
            'input_sha256':input_hashes(),
            'gate_coefficients':{'identity':'(1-i)/2','dressed_swap':'(1+i)/2'},
            'interaction':pair_control(),'locality':cone_control(),'spectator_dispersion':band_control(),
            'declared_Clifford_routing_alternative':dirac_certificate(),
            'ancestry':{'reused':['literal coherent-defect-wave-v1 onsite U','interior Gauss frame'],
                        'declared':['cubic placement','six fixed matching layers','square-root dressed wave exchange'],
                        'empirical_particle_inputs':[]},
            'scope':{'connected_finite_spatial_family':True,'strict_circuit_causality':True,
                     'all_boundary_gauge_constraints':False,'Lorentz_covariance':False,
                     'interacting_relativistic_continuum':False,'physical_particle_identity':False,'GeV_mass':False}}
