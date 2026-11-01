"""Exact complete-Fock centered dimer ground and charged spectral weights."""
import sympy as s


def annihilator(mode):
    result=s.zeros(16)
    for mask in range(16):
        if mask>>mode&1:
            result[mask^(1<<mode),mask]=(-1)**((mask&((1<<mode)-1)).bit_count())
    return result


def data(t=s.Rational(1,16), U=s.Integer(1)):
    t=s.sympify(t);U=s.sympify(U)
    if t.is_positive is not True or U.is_nonnegative is not True:
        raise ValueError('positive hopping and nonnegative coupling required')
    c=[annihilator(j) for j in range(4)]
    n=[a.T*a for a in c];N=sum(n,s.zeros(16))
    H=t*sum((c[j].T*c[j+2]+c[j+2].T*c[j] for j in (0,1)),s.zeros(16))+U*(n[0]*n[1]+n[2]*n[3])
    K=H-U*N/2
    discriminant=s.sqrt(U*U+16*t*t)
    ratio=4*t/(U+discriminant)
    ground=s.zeros(16,1)
    ground[9]=1;ground[6]=-1;ground[3]=-ratio;ground[12]=-ratio
    norm=s.simplify((ground.T*ground)[0]);Eg=-(U+discriminant)/2
    assert s.simplify(K*ground-Eg*ground)==s.zeros(16,1)
    density=s.simplify(ground*ground.T/norm)
    local=s.zeros(4)
    for i in range(4):
        for j in range(4):
            local[i,j]=s.simplify(sum(density[i+4*x,j+4*x] for x in range(4)))
    doublon=(1-U/discriminant)/4
    assert s.simplify(local-s.diag(doublon,s.Rational(1,2)-doublon,s.Rational(1,2)-doublon,doublon))==s.zeros(4)
    rows=[]
    for sign in (-1,1):
        energy=discriminant/2+sign*t
        projector=(K+U*s.eye(16)/2+sign*t*s.eye(16))/(2*sign*t)
        addition=c[0].T*ground;removal=c[0]*ground
        wa=s.simplify((addition.T*projector*addition)[0]/norm)
        wr=s.simplify((removal.T*projector*removal)[0]/norm)
        expected=(1-sign*4*t/discriminant)/4
        assert s.simplify(wa-expected)==0 and s.simplify(wr-expected)==0
        rows.append({'energy':energy,'addition_weight':wa,'removal_weight':wr})
    return {'H':H,'K':K,'N':N,'ground':ground,'ground_norm_squared':norm,
            'ground_energy':Eg,'local_density':local,'double_occupation_per_site':doublon,
            'charged_lines':rows,'signed_second_moment':t*t+U*U/4}


def certificate():
    x=data()
    return {'sites':2,'Fock_dimension':16,'hopping':'1/16','coupling':'1',
            'centered_ground_number':2,'centered_ground_unique':True,
            'ground_energy':str(x['ground_energy']),
            'local_density_diagonal':[str(s.simplify(v)) for v in x['local_density'].diagonal()],
            'local_full_CAR_faithful':True,
            'charged_lines':[{k:str(s.simplify(v)) for k,v in line.items()} for line in x['charged_lines']],
            'signed_second_moment':str(x['signed_second_moment']),
            'relativistic_mass_identification':False}
