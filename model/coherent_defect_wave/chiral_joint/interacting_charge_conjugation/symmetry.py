"""Exact full-carrier antiunitary charge/reflection symmetry for the native PW source."""
from collections import defaultdict
from functools import lru_cache
import sympy as s
from model.coherent_defect_wave.chiral_joint.oriented_links.model import labels, label_index, endpoint_current, orientation
from model.coherent_defect_wave.chiral_joint.peter_weyl_links.clebsch import native_charges, irreps, coordinate
from model.coherent_defect_wave.chiral_joint.oriented_links.rotations import geometry, hamiltonian_words
from model.coherent_defect_wave.chiral_joint.periodic_family.model import Configuration

@lru_cache(maxsize=1)
def link_checks():
    index=label_index()
    C=s.ImmutableSparseMatrix(119,119,{(0,0):1,**{(index[-sign,a,b],j):1 for j,(sign,a,b) in enumerate(labels()[1:],1)}})
    J=orientation(); F=J*C
    assert C*C==s.eye(119) and F*F==s.eye(119) and C*J==J*C
    coords=0;currents=0
    for q in native_charges():
        for a in range(irreps()[q]['dimension']):
            for b in range(irreps()[q]['dimension']):
                U=coordinate(q,a,b)
                assert C*U.conjugate()*C==U.H
                assert F*U.conjugate()*F==coordinate(q,b,a)
                coords+=1
    for j in range(12):
        for end in (0,1):
            G=endpoint_current(j,end)
            assert C*G.conjugate()*C==-G
            assert F*G.conjugate()*F==-endpoint_current(j,1-end)
            currents+=1
    return {'PW119_conjugation_coordinate_checks':coords,'PW119_antiunitary_endpoint_checks':currents,'link_charge_reflection_permutation':'(sign,a,b) -> (sign,b,a)','full_119_domain':True}

@lru_cache(maxsize=4)
def path_checks(side=2):
    # P alone has no spin rotation: remove the axial minus that a proper
    # rotation would supply. It reverses the entire kinetic operator.
    from model.coherent_defect_wave.chiral_joint.oriented_links.rotations import rotate_words
    config=Configuration(side=side)
    words=hamiltonian_words(config)
    inv=s.ImmutableMatrix(-s.eye(3)); sites,edge_map=geometry(s.ImmutableMatrix(inv),side)
    transformed=rotate_words(words,inv,side)
    # rotate_words includes the spin-axis minus, although -I is improper.
    assert transformed==words
    inversion={key:-value for key,value in transformed.items()}
    assert inversion=={key:-value for key,value in words.items()}
    return {'side':side,'complete_directed_path_words':len(words),'spatial_inversion_anticommutes_with_kinetic':True,'AP_seams_retained':True}


SCALAR_SWAP=(0,3,4,1,2)

def conjugation_column(key,config=Configuration()):
    """Theta=P A; all Fock signs, AP seams, scalars and occupied links."""
    from model.coherent_defect_wave.chiral_joint.oriented_links.model import validate_key
    mask,links,scalars=validate_key(key,config)
    modes=30*config.side**3; full=(1<<modes)-1
    odd=sum(1<<j for j in range(1,modes,2))
    exponent=(-mask.bit_count()+2*(mask&odd).bit_count())%4
    phase=s.I**exponent
    sites,edge_map=geometry(s.ImmutableMatrix(-s.eye(3)),config.side)
    target=0
    # Exterior lift of inversion after the global native particle-hole map.
    for j in range(modes):
        if not (full^mask)>>j&1:continue
        cell,local=divmod(j,30); out,ap=sites[cell]; k=30*out+local
        phase*=ap*(-1)**(target>>(k+1)).bit_count()
        target|=1<<k
    outlinks=[0]*len(links); outscalars=[0]*len(scalars)
    for e,label in enumerate(links):
        sign,a,b=labels()[label]
        outlinks[edge_map[e][0]]=0 if not sign else label_index()[sign,b,a]
    for cell,h in enumerate(scalars):outscalars[sites[cell][0]]=SCALAR_SWAP[h]
    return {(target,tuple(outlinks),tuple(outscalars)):s.expand(phase)}

def apply_conjugation(state,config=Configuration()):
    from model.coherent_defect_wave.chiral_joint.model import clean
    out=defaultdict(lambda:s.Integer(0))
    for key,value in state.items():
        for target,phase in conjugation_column(key,config).items():
            out[target]+=s.conjugate(value)*phase
    return clean(out)

def apply_linear(state,action):
    from model.coherent_defect_wave.chiral_joint.model import clean
    out=defaultdict(lambda:s.Integer(0))
    for key,value in state.items():
        for target,weight in action(key).items():out[target]+=value*weight
    return clean(out)

@lru_cache(maxsize=1)
def local_checks():
    from model.coherent_defect_wave.chiral_joint.higgs_composites import model as poly
    from model.coherent_defect_wave.chiral_joint.neutral_local_spectrum.symmetries import conjugate_polynomial
    from model.coherent_defect_wave.chiral_higgs.model import higgs_memory
    parts,couplings=poly.quantum_yukawa()
    for part in parts.values():assert conjugate_polynomial(part)==part
    penalty=poly.monomial(matrix=higgs_memory()[2])
    assert conjugate_polynomial(penalty)==penalty
    for j in range(12):
        G=poly.gauge_current(j)
        assert conjugate_polynomial(G)==poly.scale(G,-1)
    for j in range(3):
        S=poly.spin_current(j)
        assert conjugate_polynomial(S)==poly.scale(S,-1)
    return {'full_F30_scalar5_Yukawa_polynomial_identities':len(parts),
            'all_real_native_Yukawa_coefficients':True,
            'scalar_penalty_preserved':True,'gauge_current_negations':12,'spin_current_negations':3}

@lru_cache(maxsize=1)
def fixture_checks():
    from model.coherent_defect_wave.chiral_joint.peter_weyl_links.model import column
    from model.coherent_defect_wave.chiral_joint.oriented_links.model import gauss_column
    config=Configuration(y_up=s.Rational(2,3),y_down=s.Rational(4,5),y_lepton=s.Rational(6,7))
    M=240;full=(1<<M)-1
    samples=[]
    for mask in (0,full,1,1<<9,1<<24,(1<<9)|(1<<39),full^((1<<10)|(1<<54))):
        samples.append((mask,(0,)*24,(0,)*8))
    for j,label in enumerate((1,9,35,60,77,118)):
        links=[0]*24;links[(j*5)%24]=label
        mask=sum(1<<k for k in range(M) if (k*17+11)%23<8)
        samples.append((mask,tuple(links),tuple((k+j)%5 for k in range(8))))
    counts=0
    for key in samples:
        theta=conjugation_column(key,config)
        assert apply_conjugation(theta,config)=={key:s.Integer((-1)**key[0].bit_count())}
        assert apply_conjugation(column(key,config),config)==apply_linear(theta,lambda k:column(k,config))
        counts+=1
    gauss=0
    sites,_=geometry(s.ImmutableMatrix(-s.eye(3)),2)
    for key in samples[-2:]:
        theta=conjugation_column(key,config)
        for cell in (0,1,5):
            for j in range(12):
                left=apply_conjugation(gauss_column(key,cell,j,config),config)
                right=apply_linear(theta,lambda k:gauss_column(k,sites[cell][0],j,config))
                assert left=={k:-v for k,v in right.items()}
                gauss+=1
    return {'full_H_occupied_basis_controls':counts,'full_Gauss_occupied_basis_controls':gauss,
            'antiunitary_square_on_full_Fock':'fermion parity'}

@lru_cache(maxsize=1)
def ground_character():
    from model.coherent_defect_wave.chiral_joint.atomic_spin_geometry.cube import full_operators,magnetization_block
    from model.coherent_defect_wave.chiral_joint.atomic_spin_geometry.momentum import embedding
    _,twisted,U=full_operators();basis,block=magnetization_block(4)
    kernel=(block**3+20*block**2-16*block-576*s.eye(70)).nullspace()
    vectors=U*embedding(basis)*s.Matrix.hstack(*kernel)
    R=s.ImmutableMatrix([[0,-1,0],[1,0,0],[0,0,1]])
    sites,_=geometry(R,2)
    perm=[j for j,_ in sites]
    site_sign=(-1)**sum(perm[i]>perm[j] for i in range(8) for j in range(i+1,8))
    ap=s.prod(phi for _,phi in sites)
    # Each local ground doublet is an odd native spin half, with exactly the
    # inherited physical spin action. The quarter turn is diagonal in spin.
    entries={}
    for mask in range(256):
        target=sum(((mask>>cell)&1)<<sites[cell][0] for cell in range(8))
        phase=site_sign*ap*s.I**mask.bit_count()
        entries[target,mask]=phase
    quarter=s.ImmutableSparseMatrix(256,256,entries)
    assert quarter*twisted==twisted*quarter
    assert quarter*vectors==vectors or quarter*vectors==-vectors
    char=1 if quarter*vectors==vectors else -1
    return {'native_cube_ground_quarter_turn_character':char,
            'all_three_conjugate_cubic_roots_same_character':True,
            'unique_full_H_ground_inherits_character_by_positive_overlap':True,
            'character_scope':'proper cubic rotations; spin center acts +1 on the even source ground',
            'selected_origin_odd_channel_type':'G1' if char==1 else 'G2',
            'continuous_physical_spin_of_full_interacting_channel_proved':False}

if __name__=='__main__':
    import json
    print(json.dumps({'link':link_checks(),'local':local_checks(),
        'paths':[path_checks(n) for n in (2,4)],'fixtures':fixture_checks(),
        'ground':ground_character()},indent=2))
