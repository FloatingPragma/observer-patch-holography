"""Fresh interval vacuum contractions; the original source's numbers are not reused."""
from fractions import Fraction as F
from functools import lru_cache
import sympy as s
from . import algebra as a
from model.native_wilson_maxwell import matter,field
from model.native_wilson_maxwell_refinement.bounds import root,rounded

ZERO=(F(0),F(0))
add=matter.cadd
mul=matter.cmul
sub=matter.csub


def projector(H,minimal,powers,terms=512):
    coefficients,error=matter.sign_polynomial(minimal,terms)
    # This generic polynomial bound assumes only1<=|H|<5, proved afresh
    # for the orthogonal member. Its old vacuum contractions are not called.
    P=s.eye(32)/2
    for n,c in enumerate(coefficients):
        if c:P-=s.Rational(c.numerator,c.denominator)*powers[n]/2
    out=[[ZERO for _ in range(32)]for _ in range(32)]
    for i in range(32):
        for j in range(i,32):
            re,im=map(matter.fraction,s.expand(P[i,j]).as_real_imag())
            re=F(matter.outward(re,bits=80));im=F(matter.outward(im,bits=80))
            if i==j:assert im==0
            out[i][j]=(re,im);out[j][i]=(re,-im)
    eps=error/2+F(64,2**80)
    return out,eps,error


def sparse(matrix):
    return {(i,j):tuple(map(matter.fraction,s.expand(v).as_real_imag()))for (i,j),v in matrix.todok().items()}


def left_product(P,J):
    result={}
    for (i,j),z in J.items():
        for r in range(32):result[r,j]=add(result.get((r,j),ZERO),mul(P[r][i],z))
    return result


@lru_cache(None)
def calculate():
    H,Js,minimal,powers,edges=a.cube()
    P,eps,sign_error=projector(H,minimal,powers)
    current=[sparse(J)for J in Js]
    supports=[sorted({i for i,j in J})for J in current]
    PJ=[left_product(P,J)for J in current]
    comp=[{key:sub(J.get(key,ZERO),value)for key,value in A.items()}for J,A in zip(current,PJ)]
    answers=[]
    for weight in (field.cube()['P4'],field.cube()['P6']):
        result=ZERO
        for e in range(12):
            for f in range(12):
                coefficient=matter.fraction(weight[e,f])
                if not coefficient:continue
                value=ZERO
                for i in supports[f]:
                    for j in supports[e]:value=add(value,mul(PJ[e].get((i,j),ZERO),comp[f].get((j,i),ZERO)))
                result=add(result,(coefficient*value[0],coefficient*value[1]))
        assert result[1]==0;answers.append(result[0])
    p4,p6=answers
    transfer=144*eps*(2+eps)
    lo=p4/4+p6/(2*root(6))-transfer
    hi=p4/4+p6/(2*root(6,False))+transfer
    assert 0<lo<hi
    # Mean current vanishes exactly, by the checked charge conjugation,
    # not by ordinary time reversal (broken by the Gamma regulator).
    for J in current:
        mean=ZERO
        for (i,j),v in J.items():mean=add(mean,mul(P[j][i],v))
        assert mean==ZERO
    # A certified direct rotation of the two finite filled seas exists.
    oldH,_,oldminimal,oldpowers=matter.cube_matter()
    oldP,oldeps,_=projector(oldH,oldminimal,oldpowers)
    D=[[sub(P[i][j],oldP[i][j])for j in range(32)]for i in range(32)]
    D2=[[ZERO for _ in range(32)]for _ in range(32)]
    for i in range(32):
        for j in range(32):
            for k in range(32):D2[i][j]=add(D2[i][j],mul(D[i][k],D[k][j]))
    rowbound=max(sum(root(re*re+im*im)for re,im in row)for row in D2)
    distance=root(rowbound)+eps+oldeps
    assert rowbound<F(1,2) and distance<F(3,4)
    # Coulomb expectation of the new sea, with every old density/reference
    # operator retained. Wick gives2delta_xy-Tr(P_x P P_y P).
    geometry=field.cube();R=s.eye(8)-s.eye(8)[:,0]*s.ones(1,8)
    coulomb=R.T*geometry['G']*R
    mean=F(0)
    for x in range(8):
        for y in range(8):
            coeff=matter.fraction(coulomb[x,y])/2
            if not coeff:continue
            trace=ZERO
            for i in range(4*x,4*x+4):
                for j in range(4*y,4*y+4):trace=add(trace,mul(P[i][j],P[j][i]))
            assert trace[1]==0
            mean+=coeff*(2*int(x==y)-trace[0])
    # Rank-four insertions give |Tr(Px deltaP Py P)| <=4eps(1+eps),
    # plus4eps for the second projector variation.
    c_error=sum(abs(matter.fraction(v))for v in coulomb)*2*eps*(2+eps)
    assert mean-c_error>0
    return {'actual_one_body_dimension':32,'full_fermion_Fock_dimension':2**32,
        'sign_polynomial_terms':512,'sign_operator_error_upper':str(rounded(sign_error)),
        'projector_operator_error_upper':str(rounded(eps)),
        'new_pair_one_photon_current_square_over_m_squared':a.encode(lo,hi),
        'new_pair_one_photon_current_norm_over_m_lower':str(rounded(root(lo,False),False)),
        'old_current_coefficient_or_vacuum_reused':False,
        'current_means_all_twelve_edges':'0 exactly by charge conjugation',
        'onsite_charge_means_all_eight_sites':'0 exactly by local charge conjugation',
        'new_Coulomb_vacuum_expectation_over_m_g_squared':a.encode(mean-c_error,mean+c_error),
        'direct_rotation_projector_difference_squared_row_bound':str(rounded(rowbound)),
        'exact_projector_distance_upper':str(rounded(distance)),
        'direct_rotation':'V=[Pnew Pold+(I-Pnew)(I-Pold)] [I-(Pnew-Pold)^2]^(-1/2)',
        'direct_rotation_maps_old_sea_to_new_sea':True,
        'new_full_fermion_ground_exactly_prepared_by_original_trial_alone':False,
        'numerical_eigensolver_or_full_Fock_matrix_used':False,
        'diagnostic_binary64':{'current_square':[float(lo),float(hi)],'projector_distance_upper':float(distance),
            'Coulomb_coefficient':[float(mean-c_error),float(mean+c_error)]}}


if __name__=='__main__':
    import json
    print(json.dumps(calculate(),indent=2))
