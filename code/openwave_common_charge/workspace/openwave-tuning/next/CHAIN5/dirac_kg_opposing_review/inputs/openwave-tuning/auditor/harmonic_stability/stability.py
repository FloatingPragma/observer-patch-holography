"""Uniform harmonic curvature bound from the actual local face columns."""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'openwave-tuning/common_dynamics/vector_field'))
import linear_response as source


def run():
    paths=[Path(__file__),Path(source.__file__),Path(source.cc.__file__),
           ROOT/'openwave-tuning/common_dynamics/vector_field/MANIFEST.json']
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    _,_,_,terms,_=source.bloch_data()
    columns=defaultdict(list)
    for edge,face,shift,coefficient in terms:
        columns[face].append((edge,shift,coefficient))
    convolution=defaultdict(int)
    for column in columns.values():
        for edge,shift,coefficient in column:
            for other,other_shift,other_coefficient in column:
                delta=tuple(a-b for a,b in zip(shift,other_shift))
                convolution[edge,other,delta]+=coefficient*other_coefficient
    majorant=np.zeros((36,36),dtype=object)
    for (edge,other,delta),coefficient in convolution.items():
        majorant[edge,other]+=abs(coefficient)
    assert np.array_equal(majorant,majorant.T)
    weights=majorant@np.ones(36,dtype=object)
    assert min(weights)>0 and max(weights)==158
    quotients=[F(int(a),int(b)) for a,b in zip(majorant@weights,weights)]
    upper=max(quotients)
    assert upper==F(1504,17)
    assert all(F(a)<=upper*F(b) for a,b in zip(majorant@weights,weights))
    # At the pi corner every coefficient is integer. Floating eigenvectors
    # only choose a rational witness; its lower bound is verified exactly.
    corner=np.zeros((36,26),dtype=object)
    for edge,face,shift,coefficient in terms:
        corner[edge,face]+=coefficient*(-1 if sum(shift)%2 else 1)
    curvature=corner@corner.T
    _,vectors=np.linalg.eigh(np.asarray(curvature,dtype=float))
    vector=np.rint(1000*vectors[:,-1]).astype(np.int64).astype(object)
    if next(int(x) for x in vector if x)!=abs(next(int(x) for x in vector if x)):
        vector=-vector
    numerator=int(vector@(curvature@vector));denominator=int(vector@vector)
    lower=F(numerator,denominator)
    assert lower>F(200,3)
    scalar_product=upper/32
    assert scalar_product==F(47,17)<4
    schedules=[]
    for steps in (5,6,8):
        alpha=F(1,2*steps);beta=F(3,steps)
        speed_squared=steps**2*F(2,3)*alpha*beta
        electric=steps*F(2,3)*alpha
        magnetic=steps*beta
        assert speed_squared==1 and electric==F(1,3) and magnetic==3
        upper_product=alpha*beta*upper
        lower_product=alpha*beta*lower
        verdict='certified_all_refinements_stable' if upper_product<4 else 'not_certified'
        if lower_product>4:verdict='pi_corner_harmonic_instability_certified'
        schedules.append({'field_substeps':steps,'alpha':str(alpha),'beta':str(beta),
                          'macro_speed_squared':str(speed_squared),
                          'canonical_electric_coefficient':str(electric),
                          'canonical_magnetic_coefficient':str(magnetic),
                          'curvature_product_upper':str(upper_product),
                          'curvature_product_lower_witness':str(lower_product),'verdict':verdict})
    assert schedules[0]['verdict']=='pi_corner_harmonic_instability_certified'
    assert schedules[1]['verdict']=='certified_all_refinements_stable'
    assert pins=={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    result={'schema':'oph.local12.uniform_harmonic_stability.v1',
            'checks':{'actual_face_convolution_reconstructed':True,
                      'integer_positive_weight_certificate':True,
                      'all_periodic_refinement_bound':True,
                      'exact_pi_corner_Rayleigh_lower':True,
                      'scalar_matched_setting_stable':True,
                      'six_substeps_same_speed_and_coupling_stable':True,
                      'five_substeps_rejected_at_actual_pi_corner':True,
                      'source_bytes_unchanged':True},
            'positive_weights':list(map(int,weights)),
            'majorant':[[int(x) for x in row] for row in majorant],
            'uniform_curvature_upper':str(upper),
            'maximum_unweighted_row_sum':158,
            'pi_corner_witness':{'vector':list(map(int,vector)),
                                'Rayleigh_lower':str(lower),'lower_decimal':float(lower),
                                'scope':'pi momentum exists on even periodic refinements; its exact symbol suffices to reject a uniform stability claim'},
            'scalar_setting':{'alpha':'1/8','beta':'1/4','product_upper':str(scalar_product),
                              'margin_below_four':str(4-scalar_product)},
            'same_Dirac_speed_coupling_schedules':schedules,
            'proof':'For every Bloch momentum, |K_ij|<=M_ij and M w<=Lambda w with positive w. The weighted maximum-norm eigenvector estimate gives every curvature eigenvalue<=Lambda. The same bound holds before periodic folding. Symplectic field kicks have cos(omega)=1-alpha beta lambda/2; nonzero-mode harmonic stability holds for0<alpha beta lambda<4.',
            'limits':['Sixsubsteps is a separately specified Floquet member, not identical to eight; accumulated time-splitting and native apparatus error must be recomputed.',
                      'This bounds harmonic linear response, not cosine-to-quadratic quantum-state error or Maxwell-vacuum convergence.',
                      'Zero-curvature modes have their separate free/Gauss/global-sector behavior; no finite oscillator gap is inferred.',
                      'Norm conservation of the finite quantum unitary is not the harmonic stability criterion.'],
            'source_sha256':pins,'official_rows_changed':[]}
    (HERE/'STABILITY.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'upper':str(upper),'corner_lower':str(lower),
                      'scalar':result['scalar_setting'],'schedules':schedules}))


if __name__=='__main__':run()
