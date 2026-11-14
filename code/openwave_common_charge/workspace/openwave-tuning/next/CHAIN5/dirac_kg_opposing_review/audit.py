"""Independent bounded-memory full-source KG review.

Reuses the frozen primitive compiler as the source definition, never the
producer's Bloch, execution, tomography or spatial comparison routines.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from collections import defaultdict, deque
from itertools import product
import hashlib, importlib.util, json, math, sys
import numpy as np
from scipy.linalg import expm

HERE=Path(__file__).resolve().parent
ROOT=HERE/'inputs'
BASE=ROOT/'openwave-tuning/next/CHAIN4/dirac_kg_connection'
I=np.eye(2,dtype=complex);X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1,-1]).astype(complex)
AA=(np.kron(X,X),np.kron(Y,X),np.kron(Z,X));BB=np.kron(I,Z)


def read(path):return json.loads(Path(path).read_text())
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def custody():
    pins=read(HERE/'INPUTS.json')['files_sha256']
    for name,sha in pins.items():assert digest(ROOT/name)==sha,name
    assert digest(BASE/'MANIFEST.json')=='398482e1a7e05fcba687cabe6056b9f870c8e364fad47abd4f8244db4a2187f8'
    return len(pins)


def load_compiler():
    sys.path.insert(0,str(BASE))
    spec=importlib.util.spec_from_file_location('review_source_definition',BASE/'connection.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def symbol(k,h,m):
    U=math.cos(float(h*m))*np.eye(4)-1j*math.sin(float(h*m))*BB
    for p,A in zip(k,AA):U=(math.cos(p)*np.eye(4)-1j*math.sin(p)*A)@U
    return U


def math_checks(masses):
    generators=(*AA,BB)
    for i,A in enumerate(generators):
        for j,B in enumerate(generators):
            assert np.array_equal(A@B+B@A,2*np.eye(4) if i==j else np.zeros((4,4)))
    # Fresh signed and off-axis modes, not the producer's six-direction set.
    directions=[(-2,1,3),(3,-2,1),(1,1,-2),(-1,-3,2),(2,2,1),(0,-2,3)]
    rows=[]
    for L,nu,m in product((128,512,2048),directions,(*masses,F())):
        h=F(1,L);p=2*np.pi*np.array(nu);H=sum((x*A for x,A in zip(p,AA)),float(m)*BB)
        E=math.sqrt(float(p@p)+float(m*m));U=symbol(p/L,h,m)
        assert np.linalg.norm(H@H-E*E*np.eye(4),2)<1e-11
        coeff=[abs(m)]+[F(44,7)*abs(x) for x in nu]
        S=sum((coeff[i]*coeff[j] for i in range(4) for j in range(i+1,4)),F())
        delta=h*h*S;frequency=2*h*S
        actual=np.sort(-np.angle(np.linalg.eigvals(U))/float(h))
        error=float(np.max(abs(actual-[-E,-E,E,E])))
        defect=float(np.linalg.norm(U-expm(-1j*float(h)*H),2))
        n=L//8;product_error=float(np.linalg.norm(np.linalg.matrix_power(U,n)-expm(-1j*float(n*h)*H),2))
        assert delta<F(1,2) and float(h)*E+2*float(delta)<math.pi
        assert defect<float(delta)+2e-12 and error<float(frequency)+2e-10
        assert product_error<float(n*delta)+2e-11
        rows.append({'L':L,'nu':nu,'mass_exact':str(m),'one_step_error':defect,
                     'one_step_bound_exact':str(delta),'frequency_error':error,
                     'frequency_bound_exact':str(frequency),'macros':n,'product_error':product_error})
    m=masses[0];k=np.array([.13,-.19,.27]);h=F(1,32)
    wrong_mass=float(np.linalg.norm(symbol(k,h,0)-symbol(k,h,m),2))
    wrong_speed=float(np.linalg.norm(symbol(k/math.sqrt(48),h,m)-symbol(k,h,m),2))
    assert wrong_mass>0.006 and wrong_speed>0.2
    return {'fresh_bounded_momentum_cases':rows,'exact_Clifford_pairs':16,
            'wrong_mass_separation':wrong_mass,'wrong_speed_separation':wrong_speed,
            'physical_clock':'one macro is h; omega=-Arg(U)/h; both field and Dirac limiting c^2=1'}


class IndependentNative:
    def __init__(self,compiler):
        self.compiler=compiler;self.pulses=0;self.certified=set()
        source=ROOT/'openwave-tuning/next/CHAIN3/native_field_connection/inputs/reverse-engineering-reality/code/source_selection_model/response.json'
        data=read(source);self.generators=[]
        for blocks in data['generators']:
            matrix=np.zeros((6,6),complex)
            for block,rows in enumerate(blocks):
                for i,row in enumerate(rows):
                    for j,a in enumerate(row):
                        matrix[3*block+i,3*block+j]=float(F(a[0]))+math.sqrt(5)*float(F(a[1]))+1j*(float(F(a[2]))+math.sqrt(5)*float(F(a[3])))
            self.generators.append(matrix)

    @lru_cache(None)
    def response(self,coefficients):
        return expm(sum((float(F(x))*G for x,G in zip(coefficients,self.generators)),np.zeros((6,6),complex)))

    @lru_cache(None)
    def groups(self,bits,targets):
        mask=sum(1<<q for q in targets);bases=[n for n in range(1<<bits) if n&mask==0]
        return np.array([[base+sum(((local>>(len(targets)-1-j))&1)<<q for j,q in enumerate(targets))
                          for local in range(1<<len(targets))] for base in bases])

    def run(self,tape,bits,start):
        M=np.asarray(start,complex).reshape(1<<bits,-1).copy();initial=M.conj().T@M
        failure=np.zeros_like(initial);bound=F();failure_stages=[]
        for event in tape:
            signature=json.dumps(event,sort_keys=True)
            if signature not in self.certified:
                self.compiler.native.verify_pulse(event);self.certified.add(signature)
            U=self.response(tuple(event['controls_exact']));groups=self.groups(bits,tuple(event['targets']))
            block=M[groups];d=groups.shape[1]
            lost=np.einsum('ij,bjk->bik',U[d:,:d],block).reshape(-1,M.shape[1])
            effect=lost.conj().T@lost;failure+=effect;failure_stages.append(float(np.trace(effect).real))
            M[groups]=np.einsum('ij,bjk->bik',U[:d,:d],block)
            bound+=F(event['full_response_operator_error_exact']);self.pulses+=1
        residual=float(np.linalg.norm(M.conj().T@M+failure-initial,2))
        assert residual<1e-10
        return M,failure,bound,residual,failure_stages


def pointer_effects(native,tape):
    M,Fp,bound,gram,_=native.run(tape,3,np.eye(8));selected=M[:,[0,4]]
    # Independent literal native endpoint repair, including spectator shells.
    endpoints=((0,0),(0,1),(1,0),(0,2))
    result={}
    for law in ('C','E'):
        groups=defaultdict(lambda:defaultdict(lambda:np.zeros((2,2),complex)))
        records={}
        for code,(u,v) in enumerate(endpoints):
            spectator={0:(0,)*10,1:(1,)+(0,)*9,2:(1,1)+(0,)*8}[2-u-v]
            K=selected[[code,code+4]]
            for coin in (0,1):
                d=u-v;j=(d-(2*coin-1)*(abs(d)%2))//2
                out=(u-j,v+j)+spectator;key=(coin,j if law=='C' else d)
                groups[key][out]+=K/math.sqrt(2)
                records[key]={'coin':coin,'label':key[1],'current':j,'bit':int(2*j+2*coin-1<0)}
        result[law]=[(records[key],sum((K.conj().T@K for K in groups[key].values()),np.zeros((2,2),complex))) for key in sorted(groups)]
        assert np.linalg.norm(sum(E for _,E in result[law])+Fp[np.ix_([0,4],[0,4])]-np.eye(2),2)<1e-12
    return result,{'pointer_all_columns':8,'Gram_residual':gram,'word_bound_exact':str(bound)}


def native_checks(c,native,masses):
    tapes=read(BASE/'NATIVE_PROGRAMS.json');records=read(BASE/'NATIVE_RECORDS.json')
    pointer,pointer_summary=pointer_effects(native,tapes['pointer'])
    for law,rows in records['fine_effects'].items():
        for recorded,(meta,E) in zip(rows,pointer[law]):
            assert recorded['record']==meta
            assert np.linalg.norm(E-np.array(recorded['effect_real'])-1j*np.array(recorded['effect_imag']),2)<1e-12
    preparations={};flavors={};analyses={};bounds={}
    for name,tape in tapes.items():
        if name.startswith('prep_') or name.startswith('flavor_'):
            M,failed,bound,gram,_=native.run(tape,2,np.eye(4));bounds[name]=bound
            if name.startswith('prep_'):preparations[name[5:]]=M[:,0]
            else:flavors[int(name[7:])]=M[:,0]
        elif name.startswith('analysis_'):
            M,failed,bound,gram,_=native.run(tape,1,np.eye(2));analyses[name[-1]]=M;bounds[name]=bound
    assert len(preparations)==28 and len(flavors)==4
    plus=analyses['X'][:,0];plusbound=bounds['analysis_X']
    sources={};source_summaries=[]
    for case in read(BASE/'NATIVE_RESULT.json')['source_cases']:
        M,Fp,bound,gram,stages=native.run(tapes[case['name']],5,np.eye(32));ideal=np.eye(32,dtype=complex)
        for flavor,m in enumerate((*masses,F())):
            k=2*np.pi*np.array(case['mode_integer'])/case['L']
            ideal[16+4*flavor:20+4*flavor,16+4*flavor:20+4*flavor]=symbol(k,F(case['h']),m)
        error=float(np.linalg.norm(M-ideal,2));assert error<3e-10
        assert bound==F(case['native_word_error_exact'])
        sources[case['name']]=M,Fp
        source_summaries.append({'source':case['name'],'all_columns':32,'independent_operator_error':error,
                                 'complete_Gram_residual':gram,'first_failure_stages':len(stages)})
    # Independent full tomography reconstruction from every archived complete
    # coarse law; no finite-source eigenvalues are used to generate these data.
    contrasts={};normalization=0.;clock_error=0.;complete_array_error=0.
    clocks={}
    timeout=F(records['reset_contract']['complete_timeout_union_bound_exact'])
    for name,n in product(sources,(1,2,4)):
        size=int(name.split('_')[0][1:]);word=[('Z',(0,),c.native.angle(0,0,-F(1,16*size)))]*n
        M,Fp,bound,gram,_=native.run(c.events(word),1,np.eye(2));clocks[name,n]=M@plus,bound
    for row in records['records']:
        source,steps,flavor,prep=row['source'],row['macros'],row['flavor'],row['preparation']
        initial=np.kron(plus,np.kron(flavors[flavor],preparations[prep]))
        v=np.linalg.matrix_power(sources[source][0],steps)@initial
        v=v.reshape(2,16);rho=analyses[row['setting']]@(v@v.conj().T)@analyses[row['setting']].conj().T
        clock,clockbound=clocks[source,steps];clock=analyses[row['clock_setting']]@clock;clockrho=np.outer(clock,clock.conj())
        sourcebound=F(next(x['native_word_error_exact'] for x in read(BASE/'NATIVE_RESULT.json')['source_cases'] if x['name']==source))
        allowance=steps*sourcebound+clockbound+bounds['prep_'+prep]+bounds['flavor_'+str(flavor)]+2*plusbound+bounds['analysis_'+row['setting']]+bounds['analysis_'+row['clock_setting']]+2*F(pointer_summary['word_bound_exact'])+timeout
        assert allowance==F(row['complete_native_TV_bound_exact'])
        for law,values in row['joint_component_laws'].items():
            values=np.asarray(values);normalization=max(normalization,abs(float(values.sum())-1))
            recomputed=np.zeros((2,2))
            for meta,E in pointer[law[0]]:
                for cm,CE in pointer[law[1]]:
                    recomputed[meta['bit'],cm['bit']]+=float(np.trace(E@rho).real)*float(np.trace(CE@clockrho).real)
            complete_array_error=max(complete_array_error,float(np.max(abs(recomputed-values))))
            assert np.max(abs(recomputed-values))<3e-10
            assert values.min()>-1e-12 and abs(float(values.sum())-1)<3e-10
            source_contrast=float(values[0,:].sum()-values[1,:].sum())
            clock_contrast=float(values[:,0].sum()-values[:,1].sum())
            size=int(row['source'].split('_')[0][1:]);phase=row['macros']/size/8
            target=(math.cos if row['clock_setting']=='X' else math.sin)(phase)
            clock_error=max(clock_error,abs(clock_contrast-target))
            key=(row['source'],row['macros'],row['flavor'],row['preparation'],row['setting'])
            if law=='CC' and row['clock_setting']=='X':contrasts[key]=source_contrast
    reconstruction=0.;reconstructed_count=0
    for name,n,flavor in product(sources,(1,2,4),range(4)):
        def z(prep):return contrasts[name,n,flavor,prep,'X']+1j*contrasts[name,n,flavor,prep,'Y']
        U=np.zeros((4,4),complex)
        for a in range(4):U[a,a]=z('basis'+str(a))
        for a in range(4):
            for b in range(a+1,4):
                ab='pair'+str(a)+str(b)+'_q';S=z(ab+'0')-z(ab+'2');D=z(ab+'1')-z(ab+'3')
                U[a,b]=(S-1j*D)/2;U[b,a]=(S+1j*D)/2
        source=sources[name][0][16+4*flavor:20+4*flavor,16+4*flavor:20+4*flavor]
        error=float(np.linalg.norm(U-np.linalg.matrix_power(source,n),2));reconstruction=max(reconstruction,error);reconstructed_count+=1
        assert error<5e-10
    # Fresh signed, three-axis native controlled source: L16, nu=(-1,2,3).
    k=(-F(1,8),F(1,4),F(3,8));h=F(1,16)
    fresh=c.events(c.controlled(c.reduced_matter_word(k,h),4))
    M,Fp,bound,gram,stages=native.run(fresh,5,np.eye(32));target=np.eye(32,dtype=complex)
    for flavor,m in enumerate((*masses,F())):
        target[16+4*flavor:20+4*flavor,16+4*flavor:20+4*flavor]=symbol(np.pi*np.array(list(map(float,k))),h,m)
    error=float(np.linalg.norm(M-target,2));assert error<3e-10
    # Deliberately erase all nontrivial source phases: rejects the target.
    erased=float(np.linalg.norm(target-np.eye(32),2));assert erased>0.5
    (HERE/'FRESH_NATIVE_PROGRAM.json').write_text(json.dumps(fresh,indent=2)+'\n')
    return {'source_replays':source_summaries,'pointer':pointer_summary,'archived_records_checked':len(records['records']),
            'all28_native_preparations_replayed':len(preparations),'all4_native_flavor_preparations_replayed':len(flavors),
            'all_joint_CE_arrays_recomputed':4*len(records['records']),'joint_array_recomputation_error':complete_array_error,
            'all_component_native_TV_budgets_recomputed_exactly':len(records['records']),
            'independent_tomography_matrices':reconstructed_count,'tomography_replay_error':reconstruction,
            'component_normalization_diagnostic':normalization,'clock_contrast_error':clock_error,
            'fresh_native_case':{'L':16,'nu':[-1,2,3],'all_columns':32,'pulses':len(fresh),
                                 'operator_error':error,'Gram_residual':gram,'bound_exact':str(bound),
                                 'phase_erasure_control_separation':erased},
            'complete_record_domain':'These are component laws on blank-reset success. Complete reset weights, timeout, preparation, source, analysis and pointer export instruments compose; components are not normalized after failure.'}


def spatial_checks(c,native,masses):
    stored=read(BASE/'SIZE3_ROUTES.json');graph=c.native.spatial.graph(3);edges=list(graph.seams)
    indexed={(r['cell'],r['axis'],r['sign']):r for r in stored}
    for row in stored:
        divergence=defaultdict(int);current=defaultdict(int)
        for x,y in zip(row['path'],row['path'][1:]):
            assert y in graph.neighbours[x]
            current[edges.index(tuple(sorted((x,y))))]+=1 if x<y else -1
        assert sorted(current.items())==list(map(tuple,row['current']))
        for edge,j in row['current']:
            x,y=edges[edge];divergence[x]+=j;divergence[y]-=j
        assert {x:j for x,j in divergence.items() if j}=={12*row['cell']:1,12*row['destination']:-1}
        reverse=indexed[row['destination'],row['axis'],-row['sign']]
        assert reverse['path']==list(reversed(row['path']))
    # Fresh non-dyadic actual carrier geometry at L5, all750 signed paths.
    graph=c.native.spatial.graph(5);coords={tuple(v):i for i,v in enumerate(graph.cell_coordinates)}
    edges=set(graph.seams);links=c.matter.spin.common.spatial.links();paths={}
    def within(start,end):
        todo=deque([(start,[start])]);seen={start}
        while todo:
            node,path=todo.popleft()
            if node==end:return path
            for nxt in graph.neighbours[node]:
                if nxt//12==start//12 and nxt not in seen:seen.add(nxt);todo.append((nxt,path+[nxt]))
        raise AssertionError('connected cell')
    for cell,axis in product(range(125),range(3)):
        destination=list(graph.cell_coordinates[cell]);destination[axis]=(destination[axis]+1)%5;other=coords[tuple(destination)]
        a,b=links[axis][0];path=within(12*cell,12*cell+a)+within(12*other+b,12*other)
        assert all(tuple(sorted(edge)) in edges for edge in zip(path,path[1:]))
        paths[cell,axis,1]=(other,path);paths[other,axis,-1]=(cell,list(reversed(path)))
    assert len(paths)==750
    # Fresh Fourier input on the actual L4 native spatial word, all16 columns.
    mode=(-1,1,1);prep=[]
    for axis in range(3):
        for j in range(2):
            wire=(2-axis)*2+j
            prep+=list(c.words.remap(c.native.primitive_templates()['hadamard'],(wire,)))
            prep.append(('Z',(wire,),c.native.angle(-F(mode[axis]*(1<<j),4))))
    preparation=c.events(prep);wave,Pf,pbound,pgram,_=native.run(preparation,6,np.eye(64)[:,0]);wave=wave[:,0]
    phase=np.exp(-1j*np.pi*sum(F(n*3,4) for n in mode))
    ideal=np.array([phase*np.exp(2j*np.pi*(mode[0]*(i//16)+mode[1]*((i//4)%4)+mode[2]*(i%4))/4)/8 for i in range(64)])
    assert np.linalg.norm(wave-ideal)<3e-12
    embedding=np.zeros((1024,16),complex)
    for flavor,spin,cell in product(range(4),range(4),range(64)):
        embedding[spin+4*cell+256*flavor,spin+4*flavor]=wave[cell]
    tape=read(BASE/'SPATIAL_NATIVE_PROGRAM.json')['whole_spatial_macro']
    out,failed,bound,gram,_=native.run(tape,10,embedding);block=np.zeros((16,16),complex)
    for flavor,m in enumerate((*masses,F())):block[4*flavor:4*flavor+4,4*flavor:4*flavor+4]=symbol(2*np.pi*np.array(mode)/4,F(1,4),m)
    defect=float(np.linalg.norm(out-embedding@block,2));assert defect<3e-10
    (HERE/'FRESH_SPATIAL_PREPARATION.json').write_text(json.dumps(preparation,indent=2)+'\n')
    # Wider clean and dirty scratch, with a separate integer gate interpreter.
    permutations=[]
    for width in (5,6):
        rail=width;scratch=tuple(range(width+1,2*width-1));bits=2*width-1
        gates=c.signed_shift_gates(width,tuple(range(width)),rail,scratch)
        values=list(range(1<<bits))
        for kind,targets in gates:
            assert kind in ('flip','cx','ccx'),kind
            controls=targets[:-1];target=targets[-1]
            values=[x^(1<<target) if all((x>>bit)&1 for bit in controls) else x for x in values]
        assert len(set(values))==len(values)
        for railvalue,x in product((0,1),range(1<<width)):
            original=x+(railvalue<<rail);expected=((x+(1 if railvalue==0 else -1))%(1<<width))+(railvalue<<rail)
            assert values[original]==expected
        permutations.append({'coordinate_bits':width,'all_dirty_columns':len(values),'clean_columns':2*(1<<width),'logical_gates':len(gates)})
    return {'stored_size3_exact_routes':len(stored),'fresh_size5_signed_routes':len(paths),
            'fresh_L4_mode':mode,'fresh_native_full_output_dimension':1024,'all_spin_flavor_columns':16,
            'full_spatial_output_error':defect,'complete_Gram_residual':gram,
            'fresh_preparation_native_pulses':len(preparation),'spatial_source_native_pulses':len(tape),
            'fresh_signed_permutation_controls':permutations}


def resource_checks(c,masses):
    contract=read(BASE/'REFINEMENT_CONTRACT.json');h=F(contract['h']);T=F(contract['T'])
    v=[max(masses)]+[F(44,7)*abs(x) for x in contract['mode_integer']]
    S=sum((v[i]*v[j] for i in range(4) for j in range(i+1,4)),F())
    assert T*h*S==F(contract['uniform_physical_flavor_source_error_exact'])
    register=read(ROOT/'openwave-tuning/next/CHAIN3/species_reuse_connection/REGISTER.json')
    radius=max((F(hi)-F(lo))/2 for lo,hi in register['effect_intervals'])
    assert T*radius==F(contract['coefficient_transfer_error_exact'])
    reset=contract['reset_contract'];timeout=reset['registers']*F(reset['base_uniform_failure_exact'])**reset['blocks_per_register']
    assert timeout==F(reset['complete_timeout_union_bound_exact'])
    pulses=contract['native_pulses_upper'];primitive=c.native.primitive_templates();info=c.spatial_layout(18)
    logical=sum(len(primitive[name]) for axis in range(3) for name,targets in c.signed_shift_gates(18,info['coordinates'][axis],1,info['scratch']))
    mass=len(list(c.species.flavor_mass_word(masses,h,info['flavor'],0)))
    rotations=2*sum(len(list(c.words.spin_diagonalizer(a))) for a in range(3))
    factor=max(len(list(c.controlled([(kind,tuple(range(2 if kind=='J' else 1)),c.native.angle(0,0,F(1,7)))],3))) for kind in ('I','X','Y','Z','H','J'))
    computed=int(T/h)*(factor*(logical+mass+rotations)+1)+18*(6+sum(abs(x) for x in contract['mode_integer']))+24+40
    assert computed==pulses
    native=pulses*c.words.universal_pulse_error(contract['precision_bits'],max(F(4),h*max(masses)))
    assert native==F(contract['native_operator_error_upper_exact'])
    total=T*h*S+T*radius+native+timeout
    assert total==F(contract['complete_observed_record_TV_bound_exact'])
    field=read(BASE/'FULL_FIELD_TERMINAL.json');fc=field['field_contract']
    complete_field=F(fc['native_export_operator_error_exact'])+F(fc['preparation']['native_reset']['complete_timeout_union_bound_exact'])
    assert complete_field==F(field['additional_reduced_record_TV_bound_exact'])
    assert not field['all_field_native_words_executed'] and not field['larger_L_field_bound_transferred']
    # Fresh entangled matter/field/reference, with a subnormalized input.
    rng=np.random.default_rng(430070);v=rng.normal(size=40)+1j*rng.normal(size=40);v/=np.linalg.norm(v);rho=.73*np.outer(v,v.conj())
    U=symbol([.11,-.08,.16],F(1,23),masses[1]);fieldU=expm(-1j*np.diag([.2,.7,1.9,2.4,3.1]));full=np.kron(np.kron(U,fieldU),np.eye(2))
    def tracefield(x):return np.einsum('afrbfs->arbs',x.reshape(4,5,2,4,5,2)).reshape(8,8)
    reduced=np.kron(U,np.eye(2));error=float(np.linalg.norm(tracefield(full@rho@full.conj().T)-reduced@tracefield(rho)@reduced.conj().T,2));assert error<1e-12
    # Rebuild the finite-size Gaussian small-g telescoping budget directly
    # from the frozen spectral and charge certificates, without g.data().
    spectral=read(ROOT/'openwave-tuning/next/CHAIN3/native_field_connection/inputs/openwave-tuning/common_dynamics/quadratic_field/FOCK.json')['spectral_certificate']
    low=min(map(F,spectral['eigenvalues_exact']))-F(spectral['actual_curvature_operator_error_upper'])
    def sqrtbox(x,bits=100):
        scale=1<<bits;v=math.isqrt(x.numerator*scale*scale//x.denominator)
        return F(v,scale),F(v+1,scale)
    omega=sqrtbox(F(3,2)*low)[0];vco=F(1,4)/omega
    sectors=read(ROOT/'openwave-tuning/next/CHAIN3/native_field_connection/COORDINATES.json')['charge_sectors']
    coulomb=max(F(row['coulomb']) for row in sectors)
    coupling_cases=read(BASE/'CONTROLS.json')['finite_size2_small_coupling_connection']
    for row in coupling_cases:
        n=row['macros'];g=F(row['coupling_exact']);kick=sqrtbox(5*(vco+F(1,2)+F(n*n,8)))[1]
        bound=min(F(2),3*n*abs(g)*kick+n*g*g*coulomb)
        assert bound==F(row['state_error_upper_exact'])
    return {'refinement_pulse_count_exact':pulses,'refinement_total_TV_bound_exact':str(total),
            'refinement_total_TV_bound':float(total),'full_field_extra_TV_bound':float(complete_field),
            'subnormalized_correlated_reference_factor_error':error,
            'independently_recomputed_small_g_bounds':len(coupling_cases),
            'huge_word_executed':False,'interacting_g1_KG_claim':False,
            'source_identity_review':'g=1 directly delegates frozen field25; flavor replacement removes precisely its first scalar mass pulse. g=0 removes all g current and g^2 Coulomb terms while retaining every autonomous field mode and clock.'}


def main():
    count=custody();c=load_compiler();masses=tuple(map(F,read(BASE/'RESULT.json')['mass_coefficients_exact']))
    native=IndependentNative(c)
    result={'schema':'oph.dirac_kg_opposing.result.v1','input_files_verified':count,
            'math':math_checks(masses),'native':native_checks(c,native,masses),
            'spatial':spatial_checks(c,native,masses),'resources_and_source':resource_checks(c,masses)}
    result['actual_native_pulses_replayed']=native.pulses;result['distinct_exact_pulse_certificates_checked']=len(native.certified)
    assert custody()==count
    result['full_literal_criterion_met']=True
    result['scope']='Neutral g=0 free intervention of native-local12-dirac-noncompact-flavor-coupling-v1; all spin and flavor code blocks; same h clock and c=1. No interacting g=1 charged mass-pole or continuum award.'
    (HERE/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('math','native','spatial','resources_and_source')},indent=2))


if __name__=='__main__':main()
