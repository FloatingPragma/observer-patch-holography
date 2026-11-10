"""Actual local-ground spin coefficients and the full CAR two-hop moment."""
from collections import defaultdict
from functools import lru_cache
import sympy as s
from ..local_reference import physical,algebra
from ..local_reference.spectrum import certificate as local_spectrum
from ..local_reference.reference import certificate as reference_certificate
from ..higgs_composites.model import normal_order
from ..family_cell.quotient import evaluate,divide,outward
from ..periodic_family.model import native_blocks
from ..peter_weyl_links.clebsch import link
from ...chiral_source.model import PAULI


@lru_cache(maxsize=1)
def actual_spin_coefficients():
    keys,basis,g,k,h=physical.physical_block(3,1)
    native=s.zeros(30);native[9,9]=s.Rational(1,2);native[24,24]=-s.Rational(1,2)
    constant,matrix=algebra.transformed_one_body(native);assert constant==0
    operator=physical.tensor_map(keys,keys,lepton=algebra.one_body_fock(matrix.extract(algebra.LEPTONS,algebra.LEPTONS)))
    compressed=s.ImmutableSparseMatrix(g.inv()*basis.H*physical.metric_on(keys)*operator*basis)
    assert g*compressed==compressed.H*g
    x=s.Symbol('e');polynomial=s.Poly(h.charpoly(x).as_expr(),x).sqf_part()
    factor=s.Poly((x*(x-1))**3-313*(x*(x-1))**2+31708*x*(x-1)-1041468,x)
    assert polynomial.rem(factor).is_zero
    numerator=s.Integer(0);power=s.eye(h.rows)
    for i in range(polynomial.degree()):
        coefficient=s.Poly(sum(polynomial.nth(j)*x**(j-i-1)
                              for j in range(i+1,polynomial.degree()+1)),x).rem(factor)
        numerator+=2*coefficient.as_expr()*s.trace(power*compressed)
        power=power*h
    numerator=s.Poly(numerator,x).rem(factor);denominator=polynomial.diff().rem(factor)
    expected_n=14*x**5-2786*x**4+2540*x**3+572294*x**2-416014*x-29053944
    expected_d=2*x**5+3370*x**4-9060*x**3-275038*x**2+485662*x-4370808
    assert s.Poly(numerator.as_expr()*expected_d-denominator.as_expr()*expected_n,x).rem(factor).is_zero
    roots=s.polys.polytools.intervals(factor,eps=s.Rational(1,10**80))
    assert len(roots)==6 and all(m==1 for interval,m in roots)
    root=tuple(map(s.Rational,roots[0][0]))
    e0lo,e0hi=map(s.Rational,local_spectrum()['global_physical_ground_energy_interval'])
    assert e0lo<=root[0]<root[1]<=e0hi
    beta=divide(evaluate(numerator,root),evaluate(denominator,root))
    chi_n=s.Poly(3*numerator.as_expr()**2-2*numerator.as_expr()*denominator.as_expr()+denominator.as_expr()**2,x).rem(factor)
    chi_d=s.Poly(2*denominator.as_expr()**2,x).rem(factor)
    chi=divide(evaluate(chi_n,root),evaluate(chi_d,root))
    assert beta[1]<0 and chi[0]>s.Rational(1,3)
    reference=reference_certificate();assert reference['exact_total_mean_native_occupation']=='15'
    full=algebra.certificate()
    assert full['quark_collective_color_and_spin_commutators']==11
    # Multiplicity one of the maximal Sp4 module makes its commuting SU(2)
    # quark-spin representation trivial. Every family quark-spin vector
    # consequently has zero compression to this module.
    from ..local_reference.characters import certificate as character_certificate
    maximal=character_certificate()['maximal_module']
    assert maximal['multiplicity']==1 and maximal['dimension']==140
    for block in native_blocks():
        if 9 in block or 10 in block:continue
        projector=s.diag(*(int(a in block) for a in range(15)))
        for i,j in ((0,1),(1,2),(2,0)):
            total=s.kronecker_product(PAULI[i]/2,s.diag(*(int(a not in (9,10,11)) for a in range(15))))
            family=s.kronecker_product(PAULI[j]/2,projector)
            expected=s.I*s.kronecker_product(PAULI[3-i-j]/2,projector)
            assert total*family-family*total==expected
    return {'ground_energy_polynomial':str(factor.as_expr()),'selected_root':'smallest real root',
            'ground_energy_isolating_interval':[str(v) for v in root],
            'actual_Sec_z_orthogonal_compression_used':True,
            'beta_ec_definition':'P_ground S_ec P_ground = beta_ec S_ground',
            'beta_ec_numerator':str(numerator.as_expr()),'beta_ec_denominator':str(denominator.as_expr()),
            'beta_ec_interval':outward(beta,48),'beta_L':'1-beta_ec',
            'all_three_quark_family_spin_coefficients':'0',
            'quark_spin_zero_from_actual_commutants_and_unique_maximal_module':True,
            'chi_definition':'beta_ec^2+(1-beta_ec)^2/2',
            'chi_interval':outward(chi,48),
            'chi_numerator':str(chi_n.as_expr()),'chi_denominator':str(chi_d.as_expr()),
            'strict_source_independent_chi_lower_bound':'1/3'}


def project_word(word,coefficients):
    """Full CAR word projected through the proved one-body ground tensor."""
    blocks=native_blocks();which={a:i for i,block in enumerate(blocks) for a in block}
    out=s.zeros(4)
    for normal,sign in normal_order(word):
        inversions=sum(int(a//30>b//30) for i,(a,_) in enumerate(normal) for b,_ in normal[i+1:])
        factors=[]
        for site in range(2):
            local=tuple((mode%30,create) for mode,create in normal if mode//30==site)
            if not local:factors.append(s.eye(2));continue
            assert len(local)==2 and local[0][1] and not local[1][1]
            spin,a=divmod(local[0][0],15);other,b=divmod(local[1][0],15)
            if a!=b:factors.append(s.zeros(2));continue
            family=which[a];dimension=len(blocks[family])
            tensor=s.eye(2)*s.Rational(int(spin==other),2)
            tensor+=coefficients[family]/dimension*sum((axis[other,spin]*axis/2 for axis in PAULI),s.zeros(2))
            factors.append(tensor)
        out+=sign*(-1)**inversions*s.kronecker_product(*factors)
    return out


@lru_cache(maxsize=1)
def second_moment_certificate():
    blocks=native_blocks();beta=s.symbols('beta0:5',real=True);chi=sum(beta[i]**2/len(block) for i,block in enumerate(blocks))
    vacuum=[]
    for block in blocks:
        for a in block:
            for b in block:
                u=link(a,b);vacuum.extend((u[:,0],u.H[:,0]))
    v=s.ImmutableSparseMatrix.hstack(*vacuum);gram=v.H*v
    expected=[]
    for block in blocks:expected.extend([s.Rational(1,len(block))]*(2*len(block)**2))
    assert gram==s.diag(*expected) and v.cols==118
    axes=[]
    for axis,sigma in enumerate(PAULI):
        result=s.zeros(4);terms=0
        for family,block in enumerate(blocks):
            d=len(block)
            for a in block:
                for b in block:
                    for (spin,t),left in sigma.todok().items():
                        for (u,w),right in sigma.todok().items():
                            words=(((30+15*t+b,True),(15*spin+a,False),(15*u+a,True),(30+15*w+b,False)),
                                   ((15*spin+a,True),(30+15*t+b,False),(30+15*w+b,True),(15*u+a,False)))
                            factors=(s.conjugate(left)*right,left*s.conjugate(right))
                            for word,factor in zip(words,factors,strict=True):
                                result+=factor*project_word(word,beta)/(16*d);terms+=1
        result=result.applyfunc(s.expand)
        exchange=sum(((1 if i==axis else -1)*s.kronecker_product(a/2,a/2) for i,a in enumerate(PAULI)),s.zeros(4))
        target=(15*s.eye(4)-4*chi*exchange)/16
        assert (result-target).applyfunc(s.expand)==s.zeros(4)
        axes.append({'axis':axis,'full_CAR_word_controls':terms,'exact_ground_spin_matrix_identity':True})
    actual=actual_spin_coefficients();lo,hi=map(s.Rational,actual['chi_interval'])
    return {'full_PW_vacuum_link_Gram_dimension':118,'vacuum_forward_and_adjoint_channels_orthogonal':True,
            'all_native_blocks_and_two_spin_components_retained':True,'axis_checks':axes,
            'ground_one_body_tensor':'P c_(s,a)^dagger c_(t,b) P = delta_ab [delta_st I/2 + (beta_family/d_family) sum_i sigma_i,ts S_i]',
            'actual_edge_coefficient_absolute_squared':'1/16',
            'exact_P0_Tedge_squared_P0':'[15 I - 4 chi sum_i R_axis,ii S_x,i S_y,i]/16',
            'R_axis':'diagonal +1 along the edge axis, -1 on the other two spin axes',
            'actual_local_ground_spin_coefficients':actual,
            'z_edge_antialigned_minus_aligned_second_moment_interval':[str(lo/8),str(hi/8)],
            'strict_non_scalar_contrast_lower_bound':'1/24',
            'witness_uses_actual_quantum_onsite_ground':True,
            'moment_identified_with_on_shell_energy_denominator_weighted_exchange':False}
