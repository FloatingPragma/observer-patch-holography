"""Literal nearest-neighbor one-body Wilson action on a finite cubic source."""
from itertools import product
import sympy as s

from ..spatial.dirac_route import clifford
from .spectrum import exact_real,source_mass,symbol


def action(state,length,spacing):
    if type(length) is not int or length<3:raise ValueError('periodic cube of side at least three required')
    a=exact_real(spacing,True);alpha,beta=clifford();out={}
    for (position,coin,face),value in state.items():
        if len(position)!=3 or any(type(v) is not int or v not in range(length) for v in position):
            raise ValueError('source cell on the declared cube required')
        if type(coin) is not int or coin not in range(4) or type(face) is not int or face not in range(5):
            raise ValueError('Clifford and literal flux-face labels required')
        contributions=[(position,(source_mass()+3/a)*beta)]
        for axis in range(3):
            for step in (-1,1):
                target=list(position);target[axis]=(target[axis]+step)%length
                # H_(x,x+e)=-(beta+i alpha)/(2a).
                contributions.append((tuple(target),(-beta+step*s.I*alpha[axis])/(2*a)))
        for target,matrix in contributions:
            for row in range(4):
                key=tuple(target),row,face
                out[key]=s.expand(out.get(key,0)+matrix[row,coin]*value)
    return {key:s.simplify(v) for key,v in out.items() if s.simplify(v)!=0}


def certificate():
    length=3;a=s.Rational(1,2);wave=(-s.Rational(1,2)+s.I*s.sqrt(3)/2,1,1)
    p=(4*s.pi/3,0,0);target=symbol(a,p)
    for face in range(5):
        for column in range(4):
            initial={(x,column,face):s.prod(wave[j]**x[j] for j in range(3)) for x in product(range(length),repeat=3)}
            output=action(initial,length,a)
            for x in product(range(length),repeat=3):
                phase=s.prod(wave[j]**x[j] for j in range(3))
                for row in range(4):
                    assert s.simplify(s.expand_complex(output.get((x,row,face),0)-phase*target[row,column]))==0
    return {'cube_side':length,'spacing':'1/2','physical_momentum':['4*pi/3','0','0'],
            'literal_source_Bloch_columns_checked':20,
            'source_action_range':'onsite and six nearest neighbors',
            'five_flux_face_channels_unchanged_by_free_action':True}
