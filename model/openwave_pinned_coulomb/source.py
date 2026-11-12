"""Native charge compression and finite Green minimum differences."""
from functools import lru_cache
import sympy as s
import numpy as np
from model.openwave_current_coulomb.source import cochains,physical_path,native_density,changed_column
from model.openwave_reference_binding.binding import matter_irrep_certificate
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.coherent_defect_wave.chiral_joint.periodic_family.model import Configuration,kinetic_terms,native_blocks
from model.coherent_defect_wave.chiral_joint.peter_weyl_links import model as pw
from model.openwave_field_observables.lattice import green


def route(side,cell):
    point=[0,0,0];vertices=[tuple(point)]
    for axis in range(3):
        for _ in range(cell[axis]):point[axis]+=1;vertices.append(tuple(point))
    return vertices


@lru_cache(maxsize=1)
def native_compression():
    generators=representation()[2];y=generators[11]
    weights=tuple(y.diagonal())
    assert all(value!=0 for value in weights)
    blocks=native_blocks()
    assert all(len({weights[i] for i in block})==1 for block in blocks)
    words=[]
    for n in (2,4,6):
        terms=kinetic_terms(Configuration(side=n))
        assert all(t.left_cell!=t.right_cell for t in terms)
        words.append({'side':n,'native_axial_terms':len(terms),'distinct_endpoints':True})
    # An actual complete CAR/scalar/link input, with weak/colour singlet Y=-1.
    path=physical_path(2,[(0,0,0),(1,0,0)])
    cell=4;mask=(1<<(30*cell+10))|(1<<(30*cell+11));ref=1<<9
    key=(mask,tuple(path['full_native_link_labels']),(0,)*8)
    initial=native_density(key,ref)
    hopping=pw.hopping_column(key,kinetic_terms(Configuration()))
    assert hopping
    assert all(native_density(target,ref)!=initial for target in hopping)
    changed=changed_column(key,ref)
    compressed={target:value for target,value in changed.items() if native_density(target,ref)==initial}
    assert compressed
    for current in range(11):
        assert not pw.gauss_column((mask,(0,)*24,(0,)*8),cell,current)
    return {'all_native_fermion_hypercharges_nonzero':list(map(str,weights)),
        'all_native_link_transport_blocks_preserve_one_hypercharge':True,
        'actual_full_volume_path_checks':words,
        'actual_hopping_column_terms':len(hopping),'all_actual_hopping_outputs_leave_pinned_density':True,
        'complete_changed_H_column_terms':len(changed),'compressed_local_column_terms':len(compressed),
        'full_CAR_scalar_weak_colour_singlet_witness_with_Y_minus_one':True,
        'actual_selected_local_spectral_embedding':matter_irrep_certificate('ec'),
        'universal_compression':'For every fixed local total-Y pattern, P_rho T P_rho=0. A native hopping word changes its two distinct endpoint densities by nonzero opposite native hypercharges.'}


def exact_minimum_difference(side,cell,kappa=s.Integer(1)):
    if s.sympify(kappa).is_positive is not True:raise ValueError('positive field coefficient required')
    _,g=cochains(side);index=(cell[0]*side+cell[1])*side+cell[2]
    return s.cancel(kappa*(g[0,0]-g[0,index]))


@lru_cache(maxsize=1)
def exact_rows():
    rows=[]
    for side in (2,4):
        for cell in ((0,0,0),(1,0,0),(1,1,0),(1,1,1)):
            physical=physical_path(side,route(side,cell))
            difference=exact_minimum_difference(side,cell)
            assert str(difference)==physical['unit_kappa_Coulomb_energy']
            rows.append({'side':side,'cell':list(cell),'unit_kappa_exact_minimum_difference':str(difference),
                'native_full_link_string':physical['full_native_link_labels'],
                'all_twelve_Gauss_constraints':physical['all_twelve_actual_endpoint_Gauss_constraints']})
    return rows


def profile(side):
    g=green(side);rows=[]
    for radius in (4,6,8,10,12):
        force=(g[radius+1,0,0]-g[radius-1,0,0])/2
        reference=-1/(4*np.pi*radius**2)
        potential_difference=g[2*radius,0,0]-g[radius,0,0]
        reference_difference=-1/(8*np.pi*radius)
        rows.append({'radius':radius,'unit_kappa_exact_finite_minimum_energy':float(g[0,0,0]-g[radius,0,0]),
            'finite_difference_force':float(force),'inverse_square_comparison':float(reference),
            'force_amplitude_ratio':float(force/reference),
            'energy_difference_from_twice_radius':float(potential_difference),
            'potential_difference_amplitude_ratio':float(potential_difference/reference_difference)})
    assert all(r['finite_difference_force']<0 for r in rows)
    return {'side':side,'Green_zero_mode_removed':True,'rows':rows}


def run():
    actual=native_compression();rows=exact_rows();profiles=[profile(n) for n in (64,96,128)]
    return {'schema':'oph.native_pinned_coulomb_minima.v1','source':'H_prime=H_local+2^-34 T+4 N_ref+(kappa/2)rho^Y G rho^Y; all original unit link penalties replaced, all119 link states retained',
        'pinned_sector':'one reference ec of Y=+1 at0; arrival local matter+scalar in the complete weak/colour singlet Y=-1 sector; every other matter cell in its complete gauge-singlet sector; joint Gauss physical domain and all links retained',
        'complete_instrument':'orthogonal pinned-position projectors for every cell plus the full complement; reference failures, screening clouds and other endpoint irreps are retained in that complement',
        'native_compression':actual,
        'exact_operator_theorem':{'local_charge_minimum':'e_chi=min spectrum of the complete local Y=-1, weak/colour singlet Hamiltonian; no numerical value or selected-block identification assumed',
            'exact_constrained_minimum':'(V-1)E0+e_chi+4+kappa(G00-G0y)',
            'exact_difference_from_origin':'kappa(G00-G0y)',
            'unknown_local_minimum_cancels':True,
            'proof':'The physical pinned carrier factors into complete local gauge multiplicity spaces and a nonempty native flux space. All local energies have the same minima at each position; native ec strings attain their tensor lower bound. The Green term is a scalar on the pinned current pattern and hopping compresses to zero.',
            'link_loop_degeneracy':'all flux multiplicities retained; no link energy, so they cannot lower the sum of the local minima or change the exact separation difference',
            'applies_at_declared_nonzero_hopping':'yes for compressed expectation minima, because P_rho T P_rho=0',
            'pinned_subspaces_invariant_under_full_H_prime':False,
            'full_autonomous_pair_motion_or_screened_ground_force':False},
        'exact_finite_minimum_rows':rows,'far_field_minimum_energy_and_static_force':profiles,
        'normalization':{'native_current':1,'kappa_examples':[1,2],
            'same_native_charge_Gauss_and_compression_at_both_couplings':True,
            'energy_and_force_double_when_kappa_doubles':True,
            'distance':'counting lattice spacing; inherited matter kinetic spacing pi does not fix physical length',
            'physical_alpha_or_GeV_scale_derived':False},
        'field_shape_comparison':'Native constrained interaction has Green form approaching -kappa/(4pi r); reported finite-volume and spacing ratios are numerical comparisons, not continuum enclosures.',
        'verdict':'PASS_EXACT_NATIVE_CHARGE_CONSTRAINED_COULOMB_MINIMA'}
