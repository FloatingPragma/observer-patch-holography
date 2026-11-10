"""Actual particle-hole transform and collective local Yukawa algebras."""
from collections import Counter,defaultdict
from functools import lru_cache
from itertools import combinations
import sympy as s

from ..periodic_family.model import yukawa_words
from ..higgs_composites import model as car
from ...chiral_higgs.model import higgs_memory
from ...physical_attachment.chiral import representation

SINGLETS=(0,1,4,9,12,13,14)
LEPTONS=(9,10,11,24,25,26)
QUARKS=tuple(i for i in range(30) if i not in LEPTONS)


def particle_hole(operator):
    result=defaultdict(lambda:s.Integer(0))
    for (word,row,col),coefficient in operator.items():
        transformed=tuple((i,not create if i%15 in SINGLETS else create) for i,create in word)
        for normal,sign in car.normal_order(transformed):
            result[normal,row,col]+=sign*coefficient
    return car.clean(result)


@lru_cache(maxsize=1)
def transitions():
    """H'=N_H+sum_r (dGamma(M_r)|r><0|+adjoint), on all occupations."""
    compiled=defaultdict(lambda:s.Integer(0))
    for _,row,col,word,value in yukawa_words():
        for normal,sign in car.normal_order(word):
            compiled[normal,row,col]+=value*sign
    compiled=car.add(car.clean(compiled),car.monomial(matrix=higgs_memory()[2]))
    parts,couplings=car.quantum_yukawa()
    independent=car.add(*parts.values(),car.monomial(matrix=higgs_memory()[2]))
    independent=car.clean({key:value.subs({v:1 for v in couplings.values()})
                           for key,value in independent.items()})
    assert compiled==independent
    transformed=particle_hole(compiled)
    entries={r:defaultdict(lambda:s.Integer(0)) for r in range(1,5)}
    for (word,row,col),value in transformed.items():
        if col!=0:continue
        assert row in range(1,5) and len(word)==2
        assert word[0][1] and not word[1][1]
        entries[row][word[0][0],word[1][0]]+=value
    matrices=tuple(s.ImmutableSparseMatrix(30,30,entries[r]) for r in range(1,5))
    terms=[]
    for row,matrix in enumerate(matrices,1):
        unit=s.zeros(5);unit[row,0]=1
        for (i,j),value in matrix.todok().items():
            terms.append(car.scale(car.monomial(((i,True),(j,False)),unit),value))
    positive=car.add(*terms)
    assert transformed==car.add(positive,car.adjoint(positive),car.monomial(matrix=higgs_memory()[2]))
    assert matrices[2]==matrices[0].T and matrices[3]==matrices[1].T
    assert all((i in LEPTONS)==(j in LEPTONS) for matrix in matrices for i,j in matrix.todok())
    assert particle_hole(transformed)==compiled
    return matrices


def lie_closure(generators):
    """Exact sparse elimination; every appended direction is a commutator."""
    n=generators[0].rows;basis=[];echelon=[]
    def append(matrix):
        values={i*n+j:s.Rational(value) for (i,j),value in matrix.todok().items() if value}
        for pivot,previous in echelon:
            factor=values.get(pivot,0)
            if factor:
                for key,value in previous.items():
                    result=values.get(key,0)-factor*value
                    if result:values[key]=result
                    else:values.pop(key,None)
        if not values:return
        pivot=min(values);factor=values[pivot]
        values={key:value/factor for key,value in values.items()}
        echelon.append((pivot,values))
        basis.append(s.ImmutableSparseMatrix(n,n,{divmod(key,n):value for key,value in values.items()}))
    for matrix in generators:append(matrix)
    index=0
    while index<len(basis):
        for matrix in generators:append(matrix*basis[index]-basis[index]*matrix)
        index+=1
    return tuple(basis)


def one_body_fock(matrix):
    n=matrix.rows;entries=defaultdict(lambda:s.Integer(0))
    for state in range(1<<n):
        for (i,j),coefficient in matrix.todok().items():
            if not state>>j&1:continue
            sign=(-1)**((state&((1<<j)-1)).bit_count());target=state^(1<<j)
            if target>>i&1:continue
            sign*=(-1)**((target&((1<<i)-1)).bit_count());target^=1<<i
            entries[target,state]+=sign*coefficient
    return s.ImmutableSparseMatrix(1<<n,1<<n,entries)


@lru_cache(maxsize=1)
def components():
    parent=list(range(30))
    def root(i):
        while parent[i]!=i:i=parent[i]
        return i
    for matrix in transitions():
        for i,j in matrix.todok():parent[root(i)]=root(j)
    groups=defaultdict(list)
    for i in range(30):groups[root(i)].append(i)
    result=tuple(tuple(row) for row in sorted(groups.values()))
    assert Counter(map(len,result))=={4:6,3:2}
    return result


@lru_cache(maxsize=1)
def quark_clifford():
    rows=next(row for row in components() if len(row)==4)
    f=tuple(matrix.extract(rows,rows) for matrix in transitions())
    for i,j in ((0,0),(0,1),(1,0),(1,1)):
        assert f[i]*f[j]+f[j]*f[i]==s.zeros(4)
        assert f[i]*f[j+2]+f[j+2]*f[i]==(s.eye(4) if i==j else s.zeros(4))
    gamma=tuple(matrix for i in (0,1) for matrix in
                ((f[i]+f[i+2])/2,(f[i]-f[i+2])/(2*s.I)))
    for i in range(4):
        for j in range(4):
            assert gamma[i]*gamma[j]+gamma[j]*gamma[i]==(s.eye(4)/2 if i==j else s.zeros(4))
    rotations=tuple(s.I*(a*b-b*a) for a,b in combinations(gamma,2))
    generators=gamma+rotations
    assert all(a==a.H for a in generators)
    assert s.Matrix([[s.trace(a*b) for b in generators] for a in generators])==s.eye(10)
    assert sum((a*a for a in generators),s.zeros(4))==s.Rational(5,2)*s.eye(4)
    assert len(lie_closure(f))==10
    for matrix in f:
        residual=sum((s.kronecker_product(a*matrix-matrix*a,a)+
                      s.kronecker_product(a,a*matrix-matrix*a) for a in generators),s.zeros(16))
        assert residual==s.zeros(16)
    h=(f[0]*f[2]-f[2]*f[0],f[1]*f[3]-f[3]*f[1])
    assert all(a==s.diag(*a.diagonal()) for a in h)
    weights=tuple((int(h[0][i,i]),int(h[1][i,i])) for i in range(4))
    assert set(weights)=={(-1,-1),(-1,1),(1,-1),(1,1)}
    root_counts=Counter()
    for matrix in lie_closure(f):
        roots={(weights[i][0]-weights[j][0],weights[i][1]-weights[j][1]) for i,j in matrix.todok()}
        assert len(roots)==1
        root_counts[next(iter(roots))]+=1
    assert root_counts==Counter({(0,0):2,(2,0):1,(-2,0):1,(0,2):1,(0,-2):1,
                                (2,2):1,(2,-2):1,(-2,2):1,(-2,-2):1})
    return rows,f,gamma,rotations,h,weights


def transformed_one_body(matrix):
    """Return the constant and one-body matrix after the same CAR automorphism."""
    operator=particle_hole(car.one_body(matrix));constant=operator.get(((),0,0),s.Integer(0))
    entries={}
    for (word,row,col),value in operator.items():
        if row!=0 or col!=0 or not word:continue
        assert len(word)==2 and word[0][1] and not word[1][1]
        entries[word[0][0],word[1][0]]=value
    result=s.ImmutableSparseMatrix(30,30,entries)
    assert operator==car.add(car.one_body(result),car.scale(car.monomial(),constant))
    return constant,result


@lru_cache(maxsize=1)
def lepton_data():
    matrices=tuple(one_body_fock(m.extract(LEPTONS,LEPTONS)) for m in transitions())
    squared=s.ImmutableSparseMatrix(sum((m.T*m for m in matrices),s.zeros(64)))
    rows=[]
    for n in range(7):
        ids=[j for j in range(64) if j.bit_count()==n]
        spectrum=s.Matrix(squared.extract(ids,ids)).eigenvals()
        assert all(value.is_Integer and 0<=value<=6 for value in spectrum)
        assert sum(spectrum.values())==len(ids)
        rows.append({'transformed_number':n,'raw_L':n-2,
                     'squared_transfer_spectrum':{str(k):int(v) for k,v in sorted(spectrum.items())}})
    return matrices,squared,rows


@lru_cache(maxsize=1)
def certificate():
    matrices=transitions();groups=components();rows,f,gamma,rotations,h,weights=quark_clifford()
    full=lie_closure(matrices)
    assert len(full)==18
    lepton_rows=next(row for row in groups if len(row)==3)
    lepton_generators=tuple(m.extract(lepton_rows,lepton_rows) for m in matrices)
    assert len(lie_closure(lepton_generators))==8 and all(s.trace(m)==0 for m in lepton_generators)
    # All six copies obey the same Clifford relations, including native signs.
    qf=tuple(m.extract(QUARKS,QUARKS) for m in matrices)
    assert len(lie_closure(qf))==10
    assert len(lie_closure(tuple(m.extract(LEPTONS,LEPTONS) for m in matrices)))==8
    for i in (0,1):
        for j in (0,1):
            assert qf[i]*qf[j]+qf[j]*qf[i]==s.zeros(24)
            assert qf[i]*qf[j+2]+qf[j+2]*qf[i]==(s.eye(24) if i==j else s.zeros(24))
    ell=s.diag(*( (-1 if i%15==9 else 1 if i%15 in (10,11) else 0) for i in range(30)))
    baryon=s.diag(*(1 if i%15 in (2,3,5,6,7,8) else -1 if i%15 in (0,1,4,12,13,14) else 0 for i in range(30)))
    ell_constant,ell_matrix=transformed_one_body(ell)
    baryon_constant,baryon_matrix=transformed_one_body(baryon)
    assert ell_constant==-2 and ell_matrix==s.diag(*(int(i in LEPTONS) for i in range(30)))
    assert baryon_constant==-12 and baryon_matrix==s.diag(*(int(i in QUARKS) for i in range(30)))
    # A color singlet has equal transformed occupations in the three colors.
    color_cartans=[]
    for index in (6,7):
        constant,matrix=transformed_one_body(s.kronecker_product(s.eye(2),2*representation()[2][index]))
        assert constant==0 and matrix==s.diag(*matrix.diagonal())
        color_cartans.append(matrix)
    color_weights=Counter((int(color_cartans[0][i,i]),int(color_cartans[1][i,i])) for i in QUARKS)
    assert color_weights=={(1,0):8,(-1,1):8,(0,-1):8}
    spectator_checks=0
    native_spectators=[s.kronecker_product(s.eye(2),representation()[2][i]) for i in range(8)]
    from ...chiral_higgs.model import PAULI
    native_spectators.extend(s.kronecker_product(axis/2,s.eye(15)) for axis in PAULI)
    for native in native_spectators:
        constant,matrix=transformed_one_body(native)
        assert constant==0
        restricted=matrix.extract(QUARKS,QUARKS)
        assert all(restricted*f-f*restricted==s.zeros(24) for f in qf)
        spectator_checks+=1
    return {'full_native_CAR_modes':30,'full_scalar_states':5,
        'particle_hole_modes':[i for i in range(30) if i%15 in SINGLETS],
        'exact_full_normal_ordered_H_intertwiner':True,
        'one_body_components':[list(row) for row in groups],
        'quark_complex_Lie_dimension':10,'lepton_complex_Lie_dimension':8,
        'joint_complex_Lie_dimension':len(full),
        'quark_algebra':'sp4 = so5, six spinor copies','lepton_algebra':'sl3, two fundamental copies',
        'quark_fundamental_Casimir':'5/2','lepton_fundamental_Casimir':'8/3',
        'quark_Casimir_universal_tensor_commutator_checks':4,
        'quark_collective_color_and_spin_commutators':spectator_checks,
        'transformed_charge_identities':{'quark_number':'12+B3','lepton_number':'2+L'},
        'color_singlet_requires_transformed_quark_number_divisible_by':3,
        'lepton_exact_squared_transfer':lepton_data()[2],
        'lepton_squared_transfer_operator_upper_bound':6,
        'collective_quark_squared_transfer_bound':'S_q <= 2 C_q',
        'single_particle_Clifford_identity_is_not_a_collective_CAR_identity':True}
