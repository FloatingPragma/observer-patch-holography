"""A reservoir, filled reference and complete energy/charge detector.

The matrix execution uses all sixteen states of one four-mode momentum
block and both states of its particle reservoir. The other source modes
are retained in their filled reference. No success-only normalization is
substituted for the complete experiment.
"""
from functools import lru_cache
from itertools import combinations
import sympy as s

from ..many_body_route.spectrum import symbol, energy, projectors, source_mass


def clean(matrix):
    return matrix.applyfunc(s.simplify)


@lru_cache(maxsize=8)
def annihilator(index, modes=4):
    if type(modes) is not int or modes < 1 or type(index) is not int or index not in range(modes):
        raise ValueError('a mode in a positive finite CAR domain is required')
    out=s.zeros(2**modes)
    for mask in range(2**modes):
        if mask & (1 << index):
            out[mask ^ (1 << index),mask]=(-1)**((mask & ((1 << index)-1)).bit_count())
    return s.ImmutableMatrix(out)


def source_block(spacing,momentum):
    h=symbol(spacing,momentum);e=energy(spacing,momentum)
    _,negative=projectors(spacing,momentum)
    modes=[annihilator(i) for i in range(4)]
    number=s.diag(*(mask.bit_count() for mask in range(16)))
    hamiltonian=clean(sum((h[i,j]*modes[i].H*modes[j]
                          for i in range(4) for j in range(4)),s.zeros(16))+2*e*s.eye(16))
    rho=s.zeros(16)
    pairs=list(combinations(range(4),2))
    for left in pairs:
        row=sum(1 << i for i in left)
        for right in pairs:
            column=sum(1 << i for i in right)
            rho[row,column]=negative.extract(left,right).det()
    rho=clean(rho)
    assert s.trace(rho)==1 and clean(rho*rho-rho)==s.zeros(16)
    assert clean(hamiltonian*rho)==s.zeros(16)
    assert clean(number*rho-2*rho)==s.zeros(16)
    return hamiltonian,e,number,rho


def energy_projectors(hamiltonian,e):
    # The full four-mode excitation spectrum is {0,E,2E,3E,4E}.
    reduced=clean(hamiltonian/e)
    result=[]
    for value in range(5):
        projector=s.eye(16)
        for other in range(5):
            if other!=value:
                projector=clean(projector*(reduced-other*s.eye(16))/(value-other))
        result.append(projector)
    assert clean(sum(result,s.zeros(16)))==s.eye(16)
    for j,p in enumerate(result):
        assert clean(p*p-p)==s.zeros(16) and p.H==p
        assert clean(hamiltonian*p-j*e*p)==s.zeros(16)
    return result


@lru_cache(maxsize=1)
def reservoir_swap():
    """Full fermionic mode permutation: source mode zero <-> reservoir four."""
    permutation=(4,1,2,3,0);u=s.zeros(32)
    for mask in range(32):
        occupied=[permutation[i] for i in range(5) if mask & (1 << i)]
        inversions=sum(occupied[i]>occupied[j] for i in range(len(occupied))
                       for j in range(i+1,len(occupied)))
        target=sum(1 << i for i in occupied)
        u[target,mask]=(-1)**inversions
    assert u.H*u==s.eye(32) and u*u==s.eye(32)
    for j in range(5):
        assert u.H*annihilator(j,5)*u==annihilator(permutation[j],5)
    return s.ImmutableMatrix(u)


def reservoir_branches(rho,initial):
    if initial not in (0,1):raise ValueError('occupied or empty reservoir required')
    total=s.zeros(32)
    for row in range(16):
        for col in range(16):total[row+16*initial,col+16*initial]=rho[row,col]
    out=reservoir_swap()*total*reservoir_swap().H
    return [clean(out.extract(range(16*flag,16*(flag+1)),range(16*flag,16*(flag+1))))
            for flag in (0,1)]


def complete_experiment(spacing,momentum):
    h,e,n,rho=source_block(spacing,momentum);ep=energy_projectors(h,e)
    c=annihilator(0);occupation=c.H*c
    w=s.simplify(s.trace(occupation*rho))
    joint=[]
    for energy_index,projector in enumerate(ep):
        for q in range(-2,3):
            charge=s.diag(*(int(mask.bit_count()-2==q) for mask in range(16)))
            assert clean(projector*charge-charge*projector)==s.zeros(16)
            effect=clean(projector*charge)
            # Remove impossible labels, never an effect with zero probability
            # merely in the selected reference.
            if effect!=s.zeros(16):joint.append((energy_index,q,effect))
    assert len(joint)==9 and clean(sum((p for _,_,p in joint),s.zeros(16)))==s.eye(16)
    records=[]
    for initial in (1,0):
        setting='add' if initial else 'remove'
        outputs=reservoir_branches(rho,initial)
        expected=(c.H*rho*c,occupation*rho*occupation) if initial else (
            (s.eye(16)-occupation)*rho*(s.eye(16)-occupation),c*rho*c.H)
        for observed,actual,target in zip((0,1),outputs,expected,strict=True):
            assert clean(actual-target)==s.zeros(16)
            for energy_index,q,effect in joint:
                prob=s.simplify(s.trace(effect*actual)/2)
                assert prob.is_nonnegative is True
                records.append({'setting':setting,'reservoir_out':observed,
                                'success':observed!=initial,'relative_pair_charge':q,
                                'energy':str(s.simplify(energy_index*e)),
                                'battery_final_level':4-energy_index,
                                'probability':str(prob)})
    assert s.simplify(sum((s.sympify(r['probability']) for r in records),s.Integer(0)))==1
    assert len(records)==36
    expected_nonzero=[('add',0,1,e,(1-w)/2),('add',1,0,0,w*w/2),
                      ('add',1,0,2*e,w*(1-w)/2),('remove',1,-1,e,w/2),
                      ('remove',0,0,0,(1-w)**2/2),('remove',0,0,2*e,w*(1-w)/2)]
    for setting,flag,q,which,p in expected_nonzero:
        row=next(r for r in records if r['setting']==setting and r['reservoir_out']==flag
                 and r['relative_pair_charge']==q and s.simplify(s.sympify(r['energy'])-which)==0)
        assert s.simplify(s.sympify(row['probability'])-p)==0
    return {'spacing':str(spacing),'physical_momentum':list(map(str,momentum)),
            'one_particle_energy':str(e),'negative_band_probe_weight':str(w),
            'full_source_block_dimension':16,'with_reservoir_dimension':32,
            'energy_conserving_apparatus':{'battery_levels':9,'initial_level':4,
                                          'total_dimension':288,
                                          'initial_battery_energy':str(4*e),
                                          'complete_shell_isometries_instantiated':True,
                                          'full_288_matrix_materialized':False},
            'complete_records':records,'nonzero_record_count':sum(r['probability']!='0' for r in records),
            'filled_negative_band_reference':True,'reservoir_number_retained':True,
            'failed_preparation_can_create_neutral_pair':True,
            'source_charge_is_not_identified_electric_charge':True}


def mass_relation(spacing,momentum):
    a=s.sympify(spacing);p=tuple(map(s.sympify,momentum));e=energy(a,p)
    kinetic=s.simplify(sum((s.sin(a*x)**2/a**2 for x in p),s.Integer(0)))
    regulator=s.simplify(sum((1-s.cos(a*x) for x in p),s.Integer(0))/a)
    # This is an inversion of the *declared* exact regulator, not a fitted correction.
    exact_mass=s.simplify(s.sqrt(s.simplify(e*e-kinetic))-regulator)
    assert s.simplify(exact_mass-source_mass())==0
    return {'spacing':str(a),'physical_momentum':list(map(str,p)),'energy':str(e),
            'lattice_kinetic_squared':str(kinetic),'Wilson_mass_term':str(regulator),
            'inverted_source_parameter':str(exact_mass),
            'inversion_uses_declared_regulator':True}


def certificate():
    experiments=[complete_experiment(s.Integer(1),p) for p in
                 ((0,0,0),(0,0,s.pi/2))]
    relations=[mass_relation(s.Integer(1),p) for p in
               ((0,0,0),(0,0,s.pi/2),(s.pi/2,s.pi/2,0),(s.pi,s.pi,s.pi))]
    return {'schema':'oph.coherent_defect_wave.mass_instrument.v1',
            'source':'coherent-graded-flux-wilson-v1',
            'experiments':experiments,'exact_regulator_relations':relations,
            'continuum_rest_energy_in_source_units':str(source_mass()),
            'neutral_pair_threshold_in_source_units':str(2*source_mass()),
            'unit_map':{'rest_mass_GeV':'2 * source_energy_unit_in_GeV',
                        'source_energy_unit_in_GeV':None,'measured_particle_input_used':False},
            'scope':{'complete_source_mass_readout':True,
                     'source_parameter_selected_by_declared_face_count_law':True,
                     'original_Floquet_clock_logarithm':False,
                     'physical_Standard_Model_identification':False,
                     'interacting_relativistic_mass':False,'GeV_mass_derived':False}}
