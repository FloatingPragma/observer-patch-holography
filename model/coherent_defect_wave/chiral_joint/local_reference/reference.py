"""Exact symmetries and native occupations of the local ground reference."""
from functools import lru_cache
import sympy as s

from ..higgs_composites import model as car
from ..neutral_local_spectrum.symmetries import conjugate_polynomial
from ...chiral_higgs.model import higgs_memory
from .spectrum import certificate as spectral_certificate


@lru_cache(maxsize=1)
def certificate():
    spectrum=spectral_certificate()
    assert spectrum['global_physical_ground_rank']==2
    assert spectrum['global_physical_ground_twice_spin']==1
    parts,couplings=car.quantum_yukawa()
    h=car.add(*parts.values(),car.monomial(matrix=higgs_memory()[2]))
    h=car.clean({key:value.subs({y:1 for y in couplings.values()}) for key,value in h.items()})
    assert conjugate_polynomial(h)==h
    for index in range(12):
        gauge=car.gauge_current(index)
        assert not car.commutator(gauge,h)
        assert conjugate_polynomial(gauge)==car.scale(gauge,-1)
    for axis in range(3):
        spin=car.spin_current(axis)
        assert not car.commutator(spin,h)
        assert conjugate_polynomial(spin)==car.scale(spin,-1)
    groups={'uc':(0,1,4),'Q':(2,3,5,6,7,8),'ec':(9,),
            'L':(10,11),'dc':(12,13,14)}
    assert sorted(j for group in groups.values() for j in group)==list(range(15))
    records={}
    for name,group in groups.items():
        number=car.one_body(s.diag(*(int(i%15 in group) for i in range(30))))
        dimension=len(group)
        assert conjugate_polynomial(number)==car.add(
            car.scale(car.monomial(),2*dimension),car.scale(number,-1))
        assert all(not car.commutator(number,car.spin_current(axis)) for axis in range(3))
        records[name]={'native_irrep_dimension':dimension,'native_modes':2*dimension,
            'particle_hole_number_identity':f'N_{name} -> {2*dimension}-N_{name}',
            'exact_ground_mean_occupation':str(dimension)}
    return {'density':'P_ground/2','rank':2,'stationary':True,'rotation_invariant':True,
        'physical_Gauss_invariant':True,'full_native_H_charge_conjugation_identity':True,
        'charge_conjugation_reverses_all_spin_generators':True,
        'native_irrep_ground_occupations':records,'exact_total_mean_native_occupation':'15',
        'mean_occupations_also_hold_on_each_pure_ground_spin_state':True,
        'fermion_parity':'odd','raw_L':'1','raw_B3':'0','centered_raw_L_minus_one':'0',
        'raw_charge_identified_with_measured_particle_number':False,
        'bare_one_fermion_reference':False,
        'many_cell_interacting_stationary_reference_supplied':False}
