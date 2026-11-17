#!/usr/bin/env python3
"""Exact scalar composite channels of the conditional source matter module.

No elementary Higgs module, Yukawa coupling, potential, condensate or mass is
inserted. The physical two-left-Weyl CAR attachment is a declared additional
domain; it is not the auxiliary five-mode exterior grading used to build M15.
"""
from __future__ import annotations
from evidence_contract import public_receipt
import hashlib,itertools,json,sys
from pathlib import Path
sys.dont_write_bytecode=True
import sympy as s

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE))
import representation_identification as rep

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def zero(a):return all(s.expand(v)==0 for v in a)
def simp(a):return s.SparseMatrix(a.applyfunc(s.simplify))
def load_matrix(p):return s.SparseMatrix(*p['shape'],{(i,j):s.sympify(v) for i,j,v in p['entries']})

def pair_lift(h,symmetric):
    """Second quantization on normalized two-particle occupation states."""
    n=h.rows
    pairs=list(itertools.combinations_with_replacement(range(n),2) if symmetric
               else itertools.combinations(range(n),2))
    lookup={p:i for i,p in enumerate(pairs)}
    columns={j:[] for j in range(n)}
    for (i,j),v in h.todok().items():columns[j].append((i,v))
    out={}
    for col,(a,b) in enumerate(pairs):
        if symmetric:
            occupation={a:1,b:1}
            if a==b:occupation[a]=2
            for j,nj in occupation.items():
                remaining=[a,b];remaining.remove(j)
                for i,v in columns[j]:
                    final=tuple(sorted((*remaining,i)))
                    coefficient=s.sqrt(nj*(occupation.get(i,0)+1-int(i==j)))*v
                    key=lookup[final],col;out[key]=out.get(key,0)+coefficient
        else:
            for i,v in columns[a]:
                if i!=b:
                    key=lookup[tuple(sorted((i,b)))],col
                    out[key]=out.get(key,0)+(1 if i<b else -1)*v
            for i,v in columns[b]:
                if i!=a:
                    key=lookup[tuple(sorted((a,i)))],col
                    out[key]=out.get(key,0)+(1 if a<i else -1)*v
    return s.SparseMatrix(len(pairs),len(pairs),out),pairs

def main():
    packet=json.loads((HERE/'results/representation_identification_operators.json').read_text())
    basis,dg=rep.fock_data();indices=packet['even_indices']
    assert basis==[tuple(x) for x in packet['basis_subsets']]
    def matter(h):return s.SparseMatrix(dg(s.SparseMatrix(h)).extract(indices,indices))
    e={(i,j):matter(rep.unit(5,i,j)) for i in range(3) for j in range(3)}
    nc=sum((e[i,i] for i in range(3)),s.zeros(15))
    pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)]
    epsilon=s.Matrix([[0,1],[-1,0]])
    assert all(epsilon*(-p.conjugate()/2)*epsilon.inv()==p/2 for p in pauli)
    weak=[]
    for p in pauli:
        h=s.zeros(5);h[3:,3:]=p/2;weak.append(matter(h))
    q=load_matrix(packet['even_identification_operators']['q'])
    assert q==matter(s.diag(-2,-2,-2,3,3))
    ee={k:pair_lift(v,True)[0] for k,v in e.items()}
    nn=pair_lift(s.SparseMatrix(nc),True)[0]
    cc=simp(sum((ee[i,j]*ee[j,i]/2 for i in range(3) for j in range(3)),s.zeros(120))-nn*nn/6)
    ww=[pair_lift(v,True)[0] for v in weak]
    weakc=simp(sum((v*v for v in ww),s.zeros(120)))
    charge,pairs=pair_lift(q,True)
    assert len(pairs)==120 and cc.H==cc and weakc.H==weakc
    assert zero(cc*charge-charge*cc) and zero(cc*weakc-weakc*cc)
    # First decode all color-singlet charge blocks, without field names or
    # an expected Higgs table. A color Casimir kernel is basis independent.
    singlets=[];projector=s.zeros(120);doublet_projector=s.zeros(120)
    doublet_vectors=[];vectors_by_q={}
    for value in sorted(set(charge.diagonal())):
        positions=[i for i,v in enumerate(charge.diagonal()) if v==value]
        block=cc.extract(positions,positions)
        kernel=block.nullspace()
        if not kernel:continue
        z=s.zeros(120,len(kernel))
        for j,v in enumerate(kernel):
            for i,pos in enumerate(positions):z[pos,j]=v[i]
        gram=z.H*z
        p=simp(z*gram.inv()*z.H)
        assert p*p==p and p.H==p and zero(cc*p)
        dim=int(s.trace(p));c2=s.simplify(s.trace(p*weakc)/dim)
        assert weakc*p==c2*p
        weights=[]
        for t in sorted(set(ww[2].diagonal())):
            mult=s.simplify(sum(p[i,i] for i,x in enumerate(ww[2].diagonal()) if x==t))
            if mult:weights.append({'T3':str(t),'multiplicity':int(mult),'electric_Q':str(t+value/6)})
        spin=s.simplify((s.sqrt(1+4*c2)-1)/2)
        multiplicity=s.Rational(dim)/(2*spin+1)
        assert multiplicity.is_Integer
        singlets.append({'q6Y':str(value),'Y':str(value/6),'dimension':dim,
            'weak_C2':str(c2),'weak_isospin':str(spin),
            'irreducible_multiplicity':int(multiplicity),'weights':weights,
            'complete_uniform_pair_state_probability':str(s.Rational(dim,120))})
        projector+=p
        vectors_by_q[str(value)]=z
        if c2==s.Rational(3,4):
            doublet_projector+=p;doublet_vectors.append(z)
    projector=simp(projector);doublet_projector=simp(doublet_projector)
    assert s.trace(projector)==10 and s.trace(doublet_projector)==6
    assert s.eye(120)-projector==(s.eye(120)-projector).H
    assert (s.eye(120)-projector)**2==s.eye(120)-projector
    assert [(x['q6Y'],x['dimension'],x['weak_C2']) for x in singlets]==[
        ('-6',3,'2'),('-3',2,'3/4'),('3',4,'3/4'),('12',1,'0')]
    # Explicit isometry into two fermions of the independent Spin2 x M15
    # module. In spin-major ordering, the singlet coefficients are epsilon
    # tensor symmetric-internal coefficients. This uses physical CAR, not
    # the auxiliary even exterior parity of the M15 construction.
    fermion_pairs=list(itertools.combinations(range(30),2));fi={p:i for i,p in enumerate(fermion_pairs)}
    w=s.zeros(435,120)
    for j,(a,b) in enumerate(pairs):
        if a==b:w[fi[(a,15+a)],j]=1
        else:
            w[fi[(a,15+b)],j]=1/s.sqrt(2)
            w[fi[(b,15+a)],j]=1/s.sqrt(2)
    w=s.SparseMatrix(w)
    assert w.H*w==s.eye(120)
    full_scalar_projector=w*w.H
    assert full_scalar_projector*full_scalar_projector==full_scalar_projector
    assert s.trace(s.eye(435)-full_scalar_projector)==315
    # The same normalized pair operators are not canonical elementary
    # bosons. Their vacuum creation norms and filled-state annihilation
    # norms are both one; [B,B^dagger] has expectations +1 and -1.
    holes=s.diag(*[(-1)**(a+b-1) for a,b in fermion_pairs])*w.conjugate()
    assert holes.H*holes==s.eye(120)
    spin_checks=[]
    for p in pauli:
        j30=s.SparseMatrix(s.kronecker_product(p/2,s.eye(15)))
        j435,_=pair_lift(j30,False)
        assert zero(j435*w)
        for handedness in (-1,1):assert zero(handedness*s.I*j435*w)
        spin_checks.append(True)
    # The intertwining of internal currents is checked for all twelve actual
    # lifted source generators, rather than merely charge eigenvalues.
    native_checks=0
    for encoded in packet['spin_lifted_one_body_generators']:
        g=matter(load_matrix(encoded))
        sg,_=pair_lift(g,True)
        fg,_=pair_lift(s.SparseMatrix(s.kronecker_product(s.eye(2),g)),False)
        assert zero(fg*w-w*sg)
        assert zero(sg*doublet_projector-doublet_projector*sg)
        native_checks+=1
    # Derive component origins only after the operator decoder. The unordered
    # charge-pair projector commutes with color and weak transformations.
    origin=[]
    for qa,qb in sorted(set(tuple(sorted((q[a,a],q[b,b]))) for a,b in pairs)):
        ix=[i for i,(a,b) in enumerate(pairs) if tuple(sorted((q[a,a],q[b,b])))==(qa,qb)]
        rank=s.simplify(sum(doublet_projector[i,i] for i in ix))
        if rank:
            origin.append({'constituent_q6Y':[str(qa),str(qb)],'doublet_component_dimension':int(rank),
                           'weak_doublet_multiplicity':int(rank/2)})
    assert origin==[
        {'constituent_q6Y':['-4','1'],'doublet_component_dimension':2,'weak_doublet_multiplicity':1},
        {'constituent_q6Y':['-3','6'],'doublet_component_dimension':2,'weak_doublet_multiplicity':1},
        {'constituent_q6Y':['1','2'],'doublet_component_dimension':2,'weak_doublet_multiplicity':1}]
    # A positive complete effect decomposition retains all other pair states.
    other_singlet=projector-doublet_projector;colored=s.eye(120)-projector
    assert doublet_projector+other_singlet+colored==s.eye(120)
    assert all(p*p==p and p.H==p for p in [doublet_projector,other_singlet,colored])
    assert [s.trace(p)/120 for p in [doublet_projector,other_singlet,colored]]==[
        s.Rational(1,20),s.Rational(1,30),s.Rational(11,12)]
    assert sum(map(s.Rational,['2/145','4/435','22/87','21/29']))==1
    artifact=HERE/'results/scalar_bilinear_operators.json'
    artifact.write_text(json.dumps({'symmetric_pair_basis':[list(x) for x in pairs],
        'physical_two_fermion_basis':[list(x) for x in fermion_pairs],
        'Lorentz_scalar_isometry':rep.sparse(w),'q':rep.sparse(charge),
        'color_C2':rep.sparse(cc),'weak_C2':rep.sparse(weakc),
        'color_singlet_projector':rep.sparse(projector),
        'Higgs_charge_doublet_projector':rep.sparse(doublet_projector),
        'color_singlet_bases_by_q':{k:rep.sparse(v) for k,v in vectors_by_q.items()}},indent=2)+'\n', encoding='utf-8', newline='\n')
    paths=[HERE/'representation_identification.py',HERE/'results/representation_identification_operators.json',
           HERE/'results/representation_identification.json',HERE/'spin_matter_join.py',
           ROOT/'Lean/Screen/FermionSectorAssembly.lean',
           ROOT/'code/a5_closure/super_tannakian_matter_lift_certificate.py']
    result={'status':'EXACT_CONDITIONAL_LORENTZ_SCALAR_BILINEAR_DOUBLETS_CONSTRUCTED_PASS',
        'producer_sha256':sha(__file__),'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
        'artifacts':{artifact.name:sha(artifact)},'empirical_mass_or_coupling_inputs':[],
        'construction':{'input_internal_dimension':15,'independent_left_Weyl_factor_dimension':2,
            'physical_two_fermion_dimension':435,'Lorentz_scalar_internal_dimension':120,
            'internal_symmetry':'Sym^2(M15): antisymmetric two-spinor epsilon and fermionic exchange make internal pair coefficients symmetric',
            'isometry_Wdagger_W_identity':True,'rotation_and_both_same_chirality_boost_actions_annihilate_image':True,
            'actual_source_current_intertwining_checks':native_checks},
        'operator_decoded_color_singlets':singlets,'doublet_component_origins_after_decoding':origin,
        'charge_aligned_doublet_interpolators':{'Y':'1/2','independent_channels':3,
            'construction':'Two direct positive-charge channels plus epsilon times the conjugate negative-charge channel',
            'SU2_pseudoreality_checked_exactly':True,'number_of_physical_Higgs_poles_inferred':False},
        'complete_declared_uniform_pair_state_readout':{'doublet':'1/20','other_color_singlet':'1/30','colored_complement':'11/12',
            'state_domain':'I120/120 on the declared Lorentz-scalar pair subspace',
            'all_probabilities_sum_to_one':True,'new_native_pair_preparation_or_POVM_pulse_run':False},
        'complete_declared_uniform_full_two_fermion_readout':{
            'state_domain':'I435/435 on the entire two-fermion space',
            'scalar_doublet':'2/145','other_color_singlet_scalar':'4/435',
            'colored_scalar':'22/87','non_scalar_complement':'21/29',
            'non_scalar_complement_dimension':315,
            'non_scalar_representation':'Sym^2(C2) tensor Lambda^2(M15)',
            'all_probabilities_sum_to_one':True},
        'composite_CAR_control':{'normalized_pair_columns_checked':120,
            'vacuum_B_Bdagger_expectation':1,'filled_state_Bdagger_B_expectation':1,
            'vacuum_commutator_expectation':1,'filled_state_commutator_expectation':-1,
            'composite_pair_operators_asserted_to_obey_elementary_boson_CCR':False},
        'positive_result':'The existing conditional M15 and independent same-chirality Weyl factor admit one Y=-1/2 and two Y=+1/2 color-singlet Lorentz-scalar doublets, even though End(C6) has no such module.',
        'physical_scope':{'physical_CAR_attachment_declared':True,'auxiliary_even_CAR_grading_used_as_physical_statistics':False,
            'composite_interpolating_channels_constructed':True,'normalization_uses_unit_norm_two_particle_states':True,
            'one_unique_Higgs_doublet_selected':False,'bound_state_pole_or_condensate_selected':False,
            'three_charge_aligned_interpolators_imply_three_Higgs_particles':False,
            'source_Yukawa_binding_kernel_or_scalar_potential_supplied':False,
            'mass_or_width_emitted':False,'new_elementary_Higgs_representation_inserted':False,
            'native_End_C6_obstruction_applied_to_this_larger_pair_domain':False}}
    target=HERE/'results/scalar_bilinear_route.json';target.write_text(json.dumps(public_receipt(result),indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'status':result['status'],'color_singlet_rows':singlets,'origins':origin},indent=2))

if __name__=='__main__':main()
