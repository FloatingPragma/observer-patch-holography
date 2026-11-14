"""Exact current, anomaly and mass conditions for the fifteen-state attachment.

Only representation and invariant tensors are consumed from the public
action package. Its trial coefficients, Higgs vacuum and mass kernel are
excluded from this construction.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations,combinations_with_replacement
import sympy as s

from model.common_source.chiral_candidate import action


def unit(n,i,j):
    result=s.zeros(n);result[i,j]=1
    return result


def clean(matrix):
    return s.ImmutableSparseMatrix(matrix.applyfunc(s.expand))


@lru_cache(maxsize=1)
def representation():
    subsets,retained_g5,_,operators,vertices=action._representation()
    basis=[]
    for i,j in combinations(range(3),2):
        basis.extend(((unit(5,i,j)+unit(5,j,i))/2,
                      -s.I*(unit(5,i,j)-unit(5,j,i))/2))
    basis.extend((s.diag(1,-1,0,0,0)/2,s.diag(0,1,-1,0,0)/2))
    basis.extend(((unit(5,3,4)+unit(5,4,3))/2,
                  -s.I*(unit(5,3,4)-unit(5,4,3))/2,
                  s.diag(0,0,0,1,-1)/2))
    basis.append(s.diag(-s.Rational(1,3),-s.Rational(1,3),-s.Rational(1,3),
                        s.Rational(1,2),s.Rational(1,2)))
    canonical=tuple(clean(t) for t in basis)
    currents=tuple(action._exterior_derivative(t,subsets) for t in canonical)
    assert len(subsets)==15 and all(len(row) in (2,4) for row in subsets)
    assert all(t==t.H and s.trace(t)==0 and t[:3,3:]==s.zeros(3,2) for t in canonical)
    assert all(s.trace(g)==0 and g[:3,3:]==s.zeros(3,2) and g[3:,:3]==s.zeros(2,3)
               for g in retained_g5)
    # The retained port basis and this rational current basis span the same algebra.
    native_columns=s.Matrix.hstack(*(s.Matrix(g).reshape(25,1) for g in retained_g5))
    assert native_columns.to_DM(extension=True).to_field().rank()==12
    assert currents[-1]==operators['q']/6
    return subsets,canonical,currents,vertices


def anomaly_certificate():
    _,_,currents,_=representation()
    assert all(s.trace(t)==0 for t in currents)
    checks=0
    for i,j,k in combinations_with_replacement(range(12),3):
        value=s.trace(currents[i]*(currents[j]*currents[k]+currents[k]*currents[j]))
        assert s.simplify(value)==0
        checks+=1
    return {'symmetric_cubic_current_trace_checks':checks,
            'all_Tr_Ta_anticommutator_Tb_Tc':'0','mixed_gravity_linear_traces':'0',
            'weak_doublets_counting_color':4,'weak_doublet_count_even':True,
            'nonperturbative_chiral_gauge_measure_constructed':False}


def yukawa_certificate():
    subsets,g5,currents,vertices=representation()
    y=currents[-1];weights=list(y.diagonal())
    pair_sums=[s.simplify(weights[i]+weights[j]) for i in range(15) for j in range(i,15)]
    assert all(z!=0 for z in pair_sums)
    admissible=Counter()
    for i,j in combinations_with_replacement(range(15),2):
        if abs(weights[i]+weights[j])==s.Rational(1,2):
            admissible[tuple(sorted((6*weights[i],6*weights[j])))]+=1
    assert set(admissible)=={(-4,1),(1,2),(-3,6)}
    checks=0;determinant_checks=0;rows=[]
    for name,vertex in sorted(vertices.items()):
        components=vertex['components'];conjugate=vertex['conjugate_Higgs']
        for t,g in zip(currents,g5,strict=True):
            matter=-s.I*t;higgs=-s.I*g[3:,3:]
            if conjugate:
                higgs=-higgs.T
            for k,c in enumerate(components):
                assert c==c.T
                residual=matter.T*c+c*matter+sum((higgs[j,k]*components[j] for j in range(2)),s.zeros(15))
                assert clean(residual)==s.zeros(15)
                checks+=1
        for j in range(5):
            degree=action._exterior_derivative(unit(5,j,j),subsets)
            higgs_degree=unit(5,j,j)[3:,3:]
            if conjugate:
                higgs_degree=-higgs_degree.T
            for k,c in enumerate(components):
                charged=degree.T*c+c*degree+sum(
                    (higgs_degree[l,k]*components[l] for l in range(2)),s.zeros(15))
                assert clean(charged-c)==s.zeros(15)
                determinant_checks+=1
        rows.append({'channel':name,'six_Y_pair':list(vertex['origin']),
                     'conjugate_Higgs':conjugate,'neutral_component_rank':components[1].rank(),
                     'independent_complex_Yukawa_coefficient':True})
    # Neutral Higgs component, one family, independent positive test coefficients.
    yu,yd,ye,v=s.symbols('y_u y_d y_e v',positive=True)
    kernel=clean(v*(yu*vertices['up']['components'][1]+yd*vertices['down']['components'][1]
                    +ye*vertices['charged_lepton']['components'][1]))
    square=clean(kernel.H*kernel)
    assert square==s.diag(*square.diagonal())
    spectrum=Counter(map(str,square.diagonal()))
    expected=Counter({'v**2*y_u**2':6,'v**2*y_d**2':6,'v**2*y_e**2':2,'0':1})
    assert spectrum==expected
    null=[i for i,value in enumerate(square.diagonal()) if value==0]
    assert len(null)==1
    q=action.electric_charge(families=False)
    assert q[null[0],null[0]]==0 and weights[null[0]]==-s.Rational(1,2)
    # All allowed coefficients are compatible with gauge covariance. It does not fix a ratio.
    first=kernel.subs({yu:1,yd:1,ye:1,v:1})
    second=kernel.subs({yu:1,yd:2,ye:3,v:1})
    assert s.trace(first.H*first)==14 and s.trace(second.H*second)==48
    return {'unbroken_symmetric_mass_entries_checked':120,
            'unbroken_hypercharge_invariant_Weyl_mass_dimension':0,
            'minimal_Higgs_hypercharge':'1/2','full_current_Yukawa_identities':checks,
            'full_U3_times_U2_phase_characters':[1,1,1,1,1],
            'determinant_character_component_checks':determinant_checks,
            'common_U5_phase_charge_of_each_Yukawa_vertex':5,
            'constant_Yukawa_invariant_under_full_U3_times_U2':False,
            'allowed_channels':rows,
            'neutral_K_dagger_K_spectrum':dict(spectrum),
            'massless_neutral_index':null[0],'massless_neutral_auxiliary_subset':list(subsets[null[0]]),
            'one_family_independent_complex_coefficients':3,
            'three_family_independent_complex_matrix_entries':27,
            'two_gauge_compatible_mass_square_trace_controls':['14','48'],
            'charged_mass_ratio_determined_by_current_symmetry':False,
            'native_defect_program_supplies_Higgs_vacuum_or_Yukawa_coefficients':False,
            'physical_pole_masses_identified':False}


def certificate():
    subsets,_,currents,_=representation()
    y=currents[-1];chirality_rows=2*15
    return {'auxiliary_internal_dimension':15,'auxiliary_basis_subsets':[list(x) for x in subsets],
            'one_family_hypercharge_multiplicities':dict(Counter(map(str,y.diagonal()))),
            'pure_internal_tracefree_current_algebra_dimension':12,
            'current_source':'exterior representation of the actual block-internal commutant',
            'one_declared_Weyl_family_Cauchy_dimension':chirality_rows,
            'Dirac_coin_tensor_same_R15_Cauchy_dimension':60,
            'Dirac_coin_tensor_same_R15_is_chiral_SM_family':False,
            'auxiliary_even_masks_are_physical_fermionic_creators':False,
            'representation_input_files':list(action.SOURCE_FILES),
            'trial_source_coefficients_or_mass_kernel_consumed':False,
            'anomalies':anomaly_certificate(),'mass_and_Yukawa':yukawa_certificate(),
            'physical_species_identified':False,'source_GeV_masses':{}}
