"""Actual signed routes, Fourier intertwiner and controlled Dirac/KG limit."""
from connection import *
from itertools import product,islice
import sympy as s

def save(name,obj):
    (HERE/name).write_text(json.dumps(obj,indent=2,default=str)+'\n');print(name,'PASS',flush=True)

def routes():
    cases=[]
    for size in (2,3,4):
        graph=native.spatial.graph(size);edges=np.array(graph.seams,int);count=0;longest=0;currents=[]
        for cell in range(size**3):
            assert tuple(graph.cell_coordinates[cell])==(cell//size**2,(cell//size)%size,cell%size)
            for axis in range(3):
                for sign in (-1,1):
                    other,path,J=route(size,cell,axis,sign)
                    divergence=np.zeros(12*size**3,int)
                    np.add.at(divergence,edges[:,0],J);np.add.at(divergence,edges[:,1],-J)
                    expected=np.zeros_like(divergence);expected[12*cell]=1;expected[12*other]-=1
                    assert np.array_equal(divergence,expected)
                    for x,y in zip(path,path[1:]):assert y in graph.neighbours[x]
                    rev_destination,rev_path,rev_J=route(size,other,axis,-sign)
                    assert rev_destination==cell and tuple(reversed(path))==tuple(rev_path) and np.array_equal(-J,rev_J)
                    count+=1;longest=max(longest,len(path)-1)
                    currents.append({'cell':cell,'axis':axis,'sign':sign,'destination':other,'path':list(path),
                        'current':[[int(e),int(J[e])] for e in np.flatnonzero(J)]})
        # Removing an actual internal seam breaks the charge/current equation.
        J=route(size,0,0,1)[2].copy();e=next(e for e in np.flatnonzero(J) if edges[e,0]//12==edges[e,1]//12)
        wrong=J.copy();wrong[e]=0;assert not np.array_equal(wrong,J)
        delta=np.zeros(12*size**3,int);np.add.at(delta,edges[:,0],wrong-J);np.add.at(delta,edges[:,1],J-wrong);assert np.any(delta)
        cases.append({'size':size,'routes':count,'maximum_seams':longest,'all_integer_Bj':True,'all_reverse_routes':True,
                      'source_currents_sha256':hashlib.sha256(json.dumps(currents,sort_keys=True).encode()).hexdigest()})
        if size==3:save('SIZE3_ROUTES.json',currents)
    return cases

def signed_native_permutations():
    cases=[]
    for width in (1,2,3,4):
        rail=width;sc=tuple(range(width+1,width+1+max(0,width-2)));bits=width+1+len(sc)
        gates=signed_shift_gates(width,tuple(range(width)),rail,sc)
        allcols=np.arange(1<<bits);result=native.frozen.execute_reversible(gates,allcols)
        assert len(set(map(int,result)))==len(result)
        for sign in (0,1):
            for x in range(1<<width):
                source=x+(sign<<rail);expected=((x+(1 if sign==0 else -1))%(1<<width))+(sign<<rail)
                assert result[source]==expected
        # Fourier states are native Hadamard/phase product preparations.
        size=1<<width
        for j in range(size):
            wave=np.exp(2j*np.pi*j*np.arange(size)/size)/math.sqrt(size)
            for sign in (0,1):
                output=np.zeros(size,complex)
                for x in range(size):output[(x+(1 if sign==0 else -1))%size]=wave[x]
                assert np.linalg.norm(output-np.exp(-2j*np.pi*j*(1 if sign==0 else -1)/size)*wave)<2e-13
        cases.append({'coordinate_bits':width,'all_reversible_columns':len(allcols),'logical_gates':len(gates),
                      'clean_scratch_return':True,'dirty_scratch_defined_bijection':True,'all_momentum_rail_intertwiners':2*size})
    # Literal native all-column controlled signed shift at size4.
    gates=signed_shift_gates(2,(0,1),2)
    tape=events(gates_word(gates));M,fail,gram,bound=execute(tape,3,np.eye(8))
    perm=native.frozen.execute_reversible(gates,np.arange(8));target=np.eye(8)[:,perm]
    error=float(np.linalg.norm(M-target,2));assert error<2e-11
    save('SIGNED_SHIFT_NATIVE_PROGRAM.json',tape)
    return {'integer_cases':cases,'native_size4_all_columns':8,'native_pulses':len(tape),
            'native_operator_residual':error,'first_export_Gram_trace':float(np.trace(fail).real),
            'complete_Gram_residual':gram,'native_operator_error_bound_exact':str(bound)}

def factor_checks():
    # At g=1 the extension delegates the entire original iterator. Compare
    # a concrete prefix as well as the definition, without expanding it.
    prefix=list(islice(parameter_macro(1,1,1,80,coupling=1),32))
    original=list(islice(words.macro_word(1,1,1,80),32));assert prefix==original
    flavored=list(islice(flavored_parameter_macro(1,1,1,80),64))
    actual_flavored=list(islice(species.flavored_macro(1,1,1,80,MASSES),64))
    assert flavored==actual_flavored
    # g=0 deletes only zero current generators; the finite source prefix is
    # the same mass and spin-conjugated geometric cell flips.
    expected=[('Z',(0,),native.angle(0,0,F(-1,8)))]
    for axis in range(3):
        D=list(words.spin_diagonalizer(axis));expected+=D
        destinations=[matter.route(cell,axis,1)[0] for cell in range(8)]
        mask=next(iter({i^j for i,j in enumerate(destinations)}));expected+=list(words.remap(native.primitive_templates()['flip'],(2+mask.bit_length()-1,)))
        expected+=native.negate_template(D)
    actual=list(islice(parameter_macro(1,1,1,80,coupling=0),len(expected)))
    assert actual==expected
    # The normalized oscillator vacuum plus freely spreading harmonic
    # packets may be entangled with an outside reference. Factor trace holds
    # for every trace-class state; this numerical case includes correlations.
    rng=np.random.default_rng(291030);z=rng.normal(size=(24,24))+1j*rng.normal(size=(24,24));rho=z@z.conj().T;rho/=np.trace(rho)
    Uf=expm(-1j*np.diag([.2,.7,1.3]));Um=bloch((.2,.1,-.3),F(1,8),MASSES[0])
    full=np.kron(np.kron(Um,Uf),np.eye(2));small=np.kron(Um,np.eye(2))
    def reduce(x):return np.einsum('asrbst->arbt',x.reshape(4,3,2,4,3,2)).reshape(8,8)
    residual=float(np.linalg.norm(reduce(full@rho@full.conj().T)-small@reduce(rho)@small.conj().T));assert residual<1e-13
    return {'g1_original_word_delegation':True,'g1_prefix_operations_equal':len(prefix),'g1_task29_flavored_prefix_operations_equal':len(flavored),'g0_source_prefix_operations_equal':len(expected),
        'all_field_factors_retained':193,'coexact_modes':190,'harmonic_modes':3,'field_lattice_time_per_macro':1,
        'physical_time_per_macro':'h','clock_phase_per_macro':'h/8','current_scale':'g','Coulomb_scale':'g^2',
        'correlated_reference_autonomous_trace_residual':residual,
        'free_intervention':'g=0 is neutral/free; it does not prove the interacting charged source has free dispersion.',
        'spatial_family_boundary':'The g=1 native field iterator is size2. Higher-L coupled canonical sources are defined geometrically; no interacting spatial/native transfer is asserted.'}

def dispersion():
    def sym(M):return s.Matrix([[s.nsimplify(z) for z in row] for row in M])
    aa=list(map(sym,matter.ALPHA));bb=sym(matter.BETA);I=s.eye(4)
    for j,A in enumerate(aa+[bb]):
        for l,B in enumerate(aa+[bb]):assert A*B+B*A==(2*I if j==l else s.zeros(4))
    p=s.symbols('px py pz',real=True);m=s.symbols('m',real=True);H=sum((p[j]*aa[j] for j in range(3)),s.zeros(4))+m*bb
    assert (H*H-(sum(x*x for x in p)+m*m)*I).applyfunc(s.expand)==s.zeros(4)
    directions=[(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(1,2,3)];cases=[]
    for size in (64,256,1024,4096,16384):
      h=F(1,size)
      for direction in directions:
       momentum=2*math.pi*np.array(direction)
       for flavor,mass in enumerate((*MASSES,F())):
        U=bloch(momentum*float(h),h,mass);actual=np.sort(-np.angle(np.linalg.eigvals(U))/float(h))
        E=math.sqrt(float(momentum@momentum)+float(mass)**2);target=np.array([-E,-E,E,E])
        upper=[abs(mass)]+[F(44,7)*abs(n) for n in direction]
        S=sum((upper[i]*upper[j] for i in range(4) for j in range(i+1,4)),F())
        delta=h*h*S;epsilon=2*h*S
        assert delta<F(1,2) and h*sum(upper)<3
        defect=float(max(abs(actual-target)));assert defect<float(epsilon)+2e-10
        H0=sum((momentum[j]*matter.ALPHA[j] for j in range(3)),np.zeros((4,4),complex))+float(mass)*matter.BETA
        macro=float(np.linalg.norm(U-expm(-1j*float(h)*H0),2));assert macro<float(delta)+2e-13
        squared=float(max(abs(actual**2-E**2)));squared_bound=(2*sum(upper)+epsilon)*epsilon
        assert squared<float(squared_bound)+2e-8
        cases.append({'L':size,'mode_integer':direction,'flavor':flavor,'mass_exact':str(mass),'spin_eigenfrequencies':actual.tolist(),
            'positive_continuum_frequency':E,'frequency_error':defect,'frequency_error_bound_exact':str(epsilon),
            'squared_dispersion_error':squared,'squared_dispersion_bound_exact':str(squared_bound),
            'macro_operator_error':macro,'macro_operator_bound_exact':str(delta)})
    corners=[]
    for corner in product((0,1),repeat=3):
        U=bloch(np.pi*np.array(corner),F(1,8),MASSES[0]);target=(-1)**sum(corner)*expm(-1j*float(MASSES[0])/8*matter.BETA)
        err=float(np.linalg.norm(U-target));assert err<1e-13;corners.append({'corner_pi':corner,'Floquet_sign':(-1)**sum(corner),'residual':err})
    return {'symbolic_Clifford_and_KG_exact':True,'all_four_spin_components':True,'physical_family_coefficients_exact':list(map(str,MASSES)),
        'flavor_padding_mass_exact':'0','flavor_padding_kinetic_action':'same Dirac transport, not a fourth physical family',
        'controlled_cases':cases,'all_eight_Brillouin_corners':corners,
        'common_field_speed_squared':'(1/2)*3*(2/3)=1','routed_Dirac_speed_squared':'1',
        'scalar_source_identification':False,'full90_Wilson_source_identification':False}

def main():
    pins=verify();save('SOURCE_PINS.json',pins)
    save('ROUTE_GEOMETRY.json',{'cases':routes(),'native_signed_shift':signed_native_permutations()})
    save('FREE_FACTOR.json',factor_checks());save('DISPERSION.json',dispersion())
    assert verify()==pins;save('GEOMETRY_CHECKS.json',{'all_checks_pass':True,'producer_pins_unchanged':True})
if __name__=='__main__':main()
