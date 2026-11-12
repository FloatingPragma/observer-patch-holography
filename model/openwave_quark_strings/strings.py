"""Native CAR embeddings, PW119 flux and the colour-centre phase."""
from collections import defaultdict
from itertools import product
import sympy as s
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.coherent_defect_wave.chiral_joint.model import fermion
from model.coherent_defect_wave.chiral_joint.periodic_family.model import Configuration,kinetic_terms
from model.coherent_defect_wave.chiral_joint.interacting_charged_response.model import FAMILIES
from model.coherent_defect_wave.chiral_joint.oriented_links.model import label_index,gauss_column
from model.coherent_defect_wave.chiral_joint.peter_weyl_links.model import hopping_column
from model.openwave_reference_binding.strings import shortest_paths,gauss_image
from model.coherent_defect_wave.chiral_joint.interacting_charged_transport.transport import arrival_bounds
from model.coherent_defect_wave.chiral_joint.atomic_reference_control.lattice import local_norm


def clean(d):return {k:s.simplify(v) for k,v in d.items() if s.simplify(v)!=0}


def reduced_string(family,path):
    indices=FAMILIES[family];d=len(indices);lookup=label_index();state={}
    for word in product(range(d),repeat=len(path)+1):
        links=tuple(lookup[(1,indices[word[i]],indices[word[i+1]]) if forward else
                  (-1,indices[word[i+1]],indices[word[i]])]
                    for i,(_,_,_,forward) in enumerate(path))
        state[word[0],links,word[-1]]=d**(-s.Rational(len(path)+1,2))
    assert sum(s.conjugate(v)*v for v in state.values())==1
    return state


def full_string(family,path,hole_spin=0):
    """Embed actual hole states c_i|filled> and the full quantum reference."""
    destination=path[-1][2] if path else 0;indices=FAMILIES[family];out={}
    for (reference,links,matter),amplitude in reduced_string(family,path).items():
        mask,sign=fermion((1<<240)-1,30*destination+15*hole_spin+indices[matter],False)
        full_links=[0]*24
        for step,label in zip(path,links,strict=True):full_links[step[0]]=label
        out[((mask,tuple(full_links),(0,)*8),indices[reference])]=sign*amplitude
    return out


def full_gauss(state,cell,current,omit_reference=False):
    out=defaultdict(lambda:s.Integer(0));t=representation()[2][current]
    for (key,ref),amp in state.items():
        for target,coef in gauss_column(key,cell,current).items():out[target,ref]+=amp*coef
        if cell==0 and not omit_reference:
            for row in range(15):
                if t[row,ref]:out[key,row]+=amp*t[row,ref]
    return clean(out)


def inner(left,right):return s.simplify(sum(s.conjugate(v)*right.get(k,0) for k,v in left.items()))


def hop(state):
    out=defaultdict(lambda:s.Integer(0))
    terms=kinetic_terms(Configuration())
    for (key,ref),amp in state.items():
        for target,coef in hopping_column(key,terms).items():out[target,ref]+=amp*coef
    return clean(out)


def pointer_probabilities(charge,turns=1):
    phase=s.expand_complex(s.exp(2*s.pi*s.I*turns*charge))
    real,imag=s.re(phase),s.im(phase)
    probabilities=tuple(s.simplify(x) for x in ((1+real)/4,(1-real)/4,(1+imag)/4,(1-imag)/4))
    assert sum(probabilities)==1 and all(x>=0 for x in probabilities)
    return probabilities


def certificate():
    gauges=representation()[2];q=gauges[10]+gauges[11]
    c=s.Rational(2,3)*gauges[6]+s.Rational(4,3)*gauges[7]
    rows=[];reduced_checks=0;native_checks=0;routes=[]
    for family in ('uc','dc'):
        indices=FAMILIES[family];charge=-q[indices[0],indices[0]]
        assert all(-q[i,i]==charge and gauges[10][i,i]==0 for i in indices)
        assert all((q[i,i]+c[i,i]).is_Integer for i in indices)
        for destination in range(8):
            paths=shortest_paths(destination)
            for path in paths:
                state=reduced_string(family,path)
                for node in [0]+[step[2] for step in path]:
                    for current in range(12):
                        assert not gauss_image(family,path,state,node,current)
                        reduced_checks+=1
            # Actual CAR/compiler checks, including all spectator cells, on one path per endpoint.
            actual=full_string(family,paths[0])
            for node in range(8):
                for current in range(12):
                    assert not full_gauss(actual,node,current)
                    native_checks+=1
            assert full_gauss(actual,0,11,omit_reference=True)
        initial=full_string(family,());moved=hop(initial)
        norm=inner(moved,moved);moments=[]
        for destination in (1,2,4):
            for path in shortest_paths(destination):
                amplitudes=[inner(full_string(family,path,spin),moved) for spin in (0,1)]
                moment=s.simplify(sum(s.conjugate(a)*a for a in amplitudes))
                assert moment==s.Rational(1,16)
                moments.append(moment)
                routes.append({'family':family,'cell':destination,'edge':path[0][0],
                    'source_charge':str(charge),'reference_charge':str(-charge),
                    'first_hop_spin_amplitudes':list(map(str,amplitudes)),
                    'first_hop_probability_coefficient':str(moment)})
        assert sum(moments)==norm
        rows.append({'family':family,'source_hole_charge':str(charge),
            'weak_singlet':True,'colour_dimension':len(indices),
            'electromagnetic_2pi_phase':str(s.expand_complex(s.exp(2*s.pi*s.I*charge))),
            'three_turn_phase':'1','compensating_colour_phase':str(s.expand_complex(s.exp(-2*s.pi*s.I*charge))),
            'complete_one_turn_pointer_probabilities':list(map(str,pointer_probabilities(charge))),
            'complete_three_turn_pointer_probabilities':list(map(str,pointer_probabilities(charge,3))),
            'full_first_hop_norm_squared':str(norm),
            'single_turn_phase_is_identity':False})
    b=8*local_norm()+24
    bound=arrival_bounds((s.Rational(1,16),s.Rational(1,16)),s.Rational(1,2**34),
                         s.Rational(1,2**24),15,0,360,b,0)
    assert bound['strict_scaled_amplitude_margin']>0
    return {'schema':'oph.native_fractional_charge_strings.v1','families':rows,'routes':routes,
        'all_shortest_path_Gauss_checks':reduced_checks,'actual_full_CAR_Gauss_checks':native_checks,
        'source':'Full N2 PW119 Hamiltonian with lambda=2^-34 and reference cost4; filled CAR/scalar-vacuum input is declared, not the interacting ground.',
        'source_modes':240,'link_dimension':119,'reference_modes':30,
        'full_atomic_norm_bound':str(b),'time':str(s.Rational(1,2**24)),
        'per_route_original_input_probability_bounds':list(map(str,bound['probability_interval'])),
        'per_route_probability_decimal':[float(v) for v in bound['probability_interval']],
        'complete_records':['reference-vacuum/preparation failure','native family pump failure',
                            'six orthogonal complete spin-pair string arrival effects','all remaining pump-success outcomes'],
        'selected_input_preparation_and_pump_success':'1 under the declared exact filled-state preparation; general input failures retained',
        'central_phase_identity':'exp(2pi i Q) exp(2pi i C_colour)=I on the native complete Fock carrier',
        'spatial_fractional_texture_winding_derived':False,
        'physical_quark_mass_or_hadron_identification':False,'verdict':'PASS'}
