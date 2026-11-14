"""A declared periodic cubic routing of the original twelve-port cells."""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from hashlib import sha256
from math import factorial,isqrt
import json

from model.common_source import slow_geometry
from model.openwave_native_process import source as native
from model.openwave_coherent_native_architecture import source as architecture
from model.openwave_coherent_native_propagation import source as walk

@lru_cache(None)
def frame():
    c=slow_geometry.source();v=slow_geometry.raw_frame();seams=tuple(c.seams())
    for i in range(12):
        for j in range(i):
            diff=[c.q5_sub(a,b)for a,b in zip(v[i],v[j],strict=True)];d=c.q5(0)
            for x in diff:d=c.q5_add(d,c.q5_mul(x,x))
            assert(d==c.q5(4))==((j,i)in seams)
    return v

@lru_cache(None)
def links(control=False):
    c=slow_geometry.source();v=frame();zero=c.q5(0);out=[]
    for axis in range(3):
        row=[]
        for i,p in enumerate(v):
            if p[axis]!=zero:continue
            positive=c.q5_sign(c.q5_mul(p[(axis+1)%3],p[(axis+2)%3]))>0
            if control and axis==0:positive=c.q5_sign(p[1])>0
            if positive:
                q=list(p);flip=(axis+1)%3 if control and axis==0 else(axis+2)%3;q[flip]=c.q5_sub(zero,q[flip]);row.append((i,v.index(tuple(q))))
        assert len(row)==2;out.append(tuple(row))
    assert sorted(x for row in out for edge in row for x in edge)==list(range(12))
    return tuple(out)

@dataclass(frozen=True)
class CubicGraph:
    size:int
    cells:int
    ports:int
    seams:tuple
    neighbours:tuple
    inter_seams:tuple
    cell_coordinates:tuple
    @property
    def max_degree(self):return 6
    @property
    def events(self):return 2*len(self.seams)

@lru_cache(None)
def graph(size):
    if type(size)is not int or size<2:raise ValueError('periodic cubic size at least2 required')
    coords=tuple(product(range(size),repeat=3));index={c:i for i,c in enumerate(coords)};edges=set();inter=[];internal=tuple(slow_geometry.source().seams())
    for cell,c in enumerate(coords):
        for a,b in internal:edges.add((12*cell+a,12*cell+b))
        for axis,row in enumerate(links()):
            nxt=list(c);nxt[axis]=(nxt[axis]+1)%size;other=index[tuple(nxt)]
            for a,b in row:
                edge=tuple(sorted((12*cell+a,12*other+b)));assert edge not in edges;edges.add(edge);inter.append(edge)
    seams=tuple(sorted(edges));neighbours=[[]for _ in range(12*len(coords))]
    for a,b in seams:neighbours[a].append(b);neighbours[b].append(a)
    assert len(seams)==36*size**3 and len(inter)==6*size**3 and all(len(v)==6 for v in neighbours)
    g=CubicGraph(size,size**3,12*size**3,seams,tuple(map(tuple,neighbours)),tuple(sorted(inter)),coords);walk.distance(g);return g


def exact_symbol_zero(control=False):
    import sympy as s
    H=6*s.eye(12)
    for a,b in slow_geometry.source().seams():H[a,b]-=1;H[b,a]-=1
    for row in links(control):
        for a,b in row:H[a,b]-=1;H[b,a]-=1
    return H

def exact_hessian(control=False):
    import sympy as s
    H=exact_symbol_zero(control);P=s.ones(12)/12;G=(H+P).inv()-P;currents=[]
    for row in links(control):
        j=s.zeros(12,1)
        for a,b in row:j[a]-=1;j[b]+=1
        currents.append(j)
    gram=s.Matrix(3,3,lambda a,b:(currents[a].T*G*currents[b])[0]);A=s.eye(3)/6-gram/12
    if not control:
        assert all(H*j==8*j for j in currents)
        assert gram==s.eye(3)/2 and A==s.eye(3)/8
        assert H.eigenvals()=={s.Integer(0):1,s.Integer(4):3,s.Integer(6):2,s.Integer(8):6}
    encode=lambda M:[[str(v)for v in row]for row in M.tolist()]
    return {'quadratic_coefficient':encode(A),'hessian':encode(2*A),'current_pseudoinverse_gram':encode(gram),'current_vectors':[[int(v)for v in j]for j in currents],'zero_symbol_characteristic_polynomial':str(H.charpoly().as_expr()),'optical_gap_exact':'4'if not control else None,'control':control}


def symbol(k,control=False):
    import numpy as np
    H=6*np.eye(12,dtype=complex)
    for a,b in slow_geometry.source().seams():H[a,b]-=1;H[b,a]-=1
    for axis,row in enumerate(links(control)):
        p=np.exp(1j*float(k[axis]))
        for a,b in row:H[a,b]-=p;H[b,a]-=p.conjugate()
    return H

# Rational real and complex interval arithmetic for the twelve-dimensional
# spectral certificates. All arithmetic below uses Fraction endpoints.
def I(x):return(F(x),F(x))
def ia(x,y):return(x[0]+y[0],x[1]+y[1])
def ni(x):return(-x[1],-x[0])
def im(x,y):
    v=[a*b for a in x for b in y];return(min(v),max(v))
def sc(x,c):return im(x,I(c))
def ca(x,y):return(ia(x[0],y[0]),ia(x[1],y[1]))
def cn(x):return(ni(x[0]),ni(x[1]))
def cm(x,y):return(ia(im(x[0],y[0]),ni(im(x[1],y[1]))),ia(im(x[0],y[1]),im(x[1],y[0])))
def cc(x):return(x[0],ni(x[1]))
CZ=(I(0),I(0))

def pi_interval():
    def atan(q,n):
        x=F(1,q);s=sum(((-1)**j*x**(2*j+1)/F(2*j+1)for j in range(n)),F());r=x**(2*n+1)/F(2*n+1);return(s-r,s+r)
    return tuple(map(F,outward(ia(sc(atan(5,36),16),sc(atan(239,10),-4)),60)))

def trig_interval(angle):
    lo,hi=angle;x=(lo+hi)/2;spread=(hi-lo)/2;n=24
    co=sum(((-1)**j*x**(2*j)/factorial(2*j)for j in range(n)),F());si=sum(((-1)**j*x**(2*j+1)/factorial(2*j+1)for j in range(n)),F())
    cr=abs(x)**(2*n)/factorial(2*n)+spread;sr=abs(x)**(2*n+1)/factorial(2*n+1)+spread
    return tuple(tuple(map(F,outward(z,60)))for z in ((max(F(-1),co-cr),min(F(1),co+cr)),(max(F(-1),si-sr),min(F(1),si+sr))))

def symbol_interval(k):
    H=[[(I(6 if i==j else 0),I(0))for j in range(12)]for i in range(12)]
    for a,b in slow_geometry.source().seams():H[a][b]=ca(H[a][b],(I(-1),I(0)));H[b][a]=ca(H[b][a],(I(-1),I(0)))
    for axis,row in enumerate(links()):
        phase=trig_interval(k[axis])
        for a,b in row:H[a][b]=ca(H[a][b],cn(phase));H[b][a]=ca(H[b][a],cn(cc(phase)))
    return H

def outward(x,places=15):
    scale=10**places;lo=x[0].numerator*scale//x[0].denominator;hi=-((-x[1].numerator*scale)//x[1].denominator)
    return[str(F(lo,scale)),str(F(hi,scale))]

def low_branch_certificate(k,label):
    import numpy as np
    k=tuple(tuple(map(F,x))for x in k);r=max(abs(v)for x in k for v in x);assert r<=F(1,2)
    mid=[float((a+b)/2)for a,b in k];_,vec=np.linalg.eigh(symbol(mid));scale=10**12
    v=[(F(round(z.real*scale),scale),F(round(z.imag*scale),scale))for z in vec[:,0]]
    norm=sum(a*a+b*b for a,b in v);assert norm>0;H=symbol_interval(k);Hv=[]
    for row in H:
        total=CZ
        for h,(a,b)in zip(row,v,strict=True):total=ca(total,cm(h,(I(a),I(b))))
        Hv.append(total)
    q=CZ
    for(a,b),z in zip(v,Hv,strict=True):q=ca(q,cm((I(a),I(-b)),z))
    assert q[1][0]<=0<=q[1][1];q=sc(q[0],1/norm);center=sum(q)/2;residual=F()
    for(a,b),z in zip(v,Hv,strict=True):
        w=ca(z,(I(-center*a),I(-center*b)));residual+=max(abs(t)for t in w[0])**2+max(abs(t)for t in w[1])**2
    residual/=norm;beta=4-r;assert q[1]<beta
    low=max(F(),q[0]-residual/(beta-q[1]));high=q[1]
    k2=(F(),F())
    for x in k:k2=ia(k2,im(x,x)if x[0]>=0 or x[1]<=0 else(F(),max(x[0]**2,x[1]**2)))
    leading=sc(k2,F(1,8));remainder=F(3,8)*r**4
    assert high>=leading[0]-remainder and low<=leading[1]+remainder
    # Stronger check: the certified whole interval itself obeys the theorem.
    assert low>=max(F(),leading[0]-remainder) and high<=leading[1]+remainder
    weight=1-r*r/(4-r-r*r/2)**2
    return {'label':label,'k_interval_exact':[[str(x)for x in row]for row in k],'ground_eigenvalue_interval_exact':outward((low,high),20),'temple_variance_upper_exact':outward((residual,residual),50)[1],'second_eigenvalue_lower_exact':str(beta),'trial_preparation_all12_components_retained':True,'trial_norm_squared_exact':str(norm),'quadratic_term_interval_exact':outward(leading,20),'fourth_order_remainder_upper_exact':str(remainder),'uniform_internal_preparation_ground_band_weight_lower_exact':str(weight),'exact_spectral_projection_assumed_in_preparation':False,'eigenvector_generation':'Floating source eigensolution chooses a rational trial; its full-source interval certificate is rational.'}


def symmetry_certificate():
    c=slow_geometry.source();v=frame();zero=c.q5(0);perms=[]
    for label,transform,axes,sgn in [('two_coordinate_flips',lambda p:(c.q5_sub(zero,p[0]),c.q5_sub(zero,p[1]),p[2]),(0,1,2),(-1,-1,1)),('cyclic_coordinate_permutation',lambda p:(p[2],p[0],p[1]),(1,2,0),(1,1,1))]:
        perm=tuple(v.index(transform(p))for p in v);assert sorted(perm)==list(range(12));seams=set(tuple(e)for e in slow_geometry.source().seams())
        assert {tuple(sorted((perm[a],perm[b])))for a,b in seams}==seams
        for axis,row in enumerate(links()):
            for a,b in row:
                image=(perm[a],perm[b])if sgn[axis]==1 else(perm[b],perm[a]);assert image in links()[axes[axis]]
        perms.append({'transformation':label,'port_permutation':list(perm),'macro_axis_permutation':list(axes),'macro_axis_signs':list(sgn),'internal_and_oriented_intercell_edges_exact':True})
    return perms

def native_columns(g):
    """All full one-load columns, using the proved local hold complement."""
    rows=native.raw_local_branches('mixture');assert len(rows)==25
    assert sum(r.weight for r in rows if(0,0)in r.pairs)==1
    incident=[[]for _ in range(g.ports)]
    for e in g.seams:
        for x in e:incident[x].append(e)
    digest=sha256();M=len(g.seams);D=2*M
    for x in range(g.ports):
        values={x:F(M-len(incident[x]),M)}
        for e in incident[x]:
            for row in rows:
                image=native.apply_branch(1<<x,e,row)
                if image is not None:
                    assert image.bit_count()==1;y=image.bit_length()-1;values[y]=values.get(y,F())+row.weight/M
        nums=[]
        for y in range(g.ports):
            v=values.get(y,F())*D;assert v.denominator==1
            expected=D-6 if x==y else int(y in g.neighbours[x]);assert v==expected;nums.append(int(v))
        assert sum(nums)==D
        digest.update(','.join(map(str,nums)).encode()+b'\n')
    return {'complete_one_load_source_columns':g.ports,'complete_population_matrix_entries':g.ports**2,'raw_slots_per_endpoint_event':25,'actual_incident_raw_branch_evaluations':g.ports*6*25,'nonincident_events_per_column_factored_by_exact_locality':M-6,'all_zero_raw_slots_retained':True,'complete_integer_column_digest':digest.hexdigest(),'transfer_identity':'T=I-L/(2M), self numerator2M-6 and neighbour numerator1','arbitrary_reference_extension':'Inherited local recorded C/E isometry and inverse recoder, then controlled source swaps and event reflection.'}

def complete_event_control(g):
    """Every S basis column, complete R Gram identity and padded failure."""
    import numpy as np
    N,D=g.ports,g.events;dimension=N*D
    permutation=np.arange(dimension,dtype=np.int64)
    for e,(a,b)in enumerate(g.seams):
        i,j=a*D+2*e+1,b*D+2*e+1
        permutation[i],permutation[j]=permutation[j],permutation[i]
    assert np.array_equal(permutation[permutation],np.arange(dimension))
    assert np.count_nonzero(permutation!=np.arange(dimension))==D
    diagonal=F(2,D)-1;other=F(2,D)
    assert diagonal**2+(D-1)*other**2==1
    assert 2*diagonal*other+(D-2)*other**2==0
    padded=1<<(D-1).bit_length()
    return {'all_source_event_basis_columns_checked':dimension,'source_swap_involution_exact':True,'moved_arc_columns':D,'uniform_reflection_complete_Gram_diagonal_exact':'1','uniform_reflection_complete_Gram_offdiagonal_exact':'0','event_register_padded_dimension':padded,'event_padding_dimension':padded-D,'additional_padded_W_minus_one_modes':N*(padded-D),'full_padded_source_event_dimension':N*padded,'padding_rule':'S is identity and R is minus identity on event padding, so W is minus identity; terminal admission retains every padding failure without postselection.','uniform_event_preparation_and_reflection':'Declared exact coherent services; complete fine-cube native pulse compilation is not supplied.','source_move_locality':'One W step moves at most one graph edge; reflection changes only the retained event register.','source_program_order':'The sequential coherent event word is declared; disjoint local endpoint operations retain their commuting partial-order semantics.'}

def bloch_control(g):
    import numpy as np
    L=g.size;C=g.cells;coords=np.asarray(g.cell_coordinates);modes=tuple(product(range(L),repeat=3));characters=[];operators=0
    for n in modes:
        k=2*np.pi*np.asarray(n)/L;chi=np.exp(1j*(coords@k));characters.append(chi/np.sqrt(C));U=np.zeros((g.ports,12),complex)
        for cell,z in enumerate(chi):U[12*cell:12*cell+12]=z*np.eye(12)/np.sqrt(C)
        action=6*U.copy()
        for a,b in g.seams:action[a]-=U[b];action[b]-=U[a]
        H=symbol(k);assert np.max(abs(action-U@H))<2e-12
        vals=np.linalg.eigvalsh(H);assert vals.min()>-2e-12 and vals.max()<12+2e-12;operators+=12
    Fourier=np.stack(characters,axis=1);assert np.max(abs(Fourier.conjugate().T@Fourier-np.eye(C)))<2e-12
    return {'quantized_momenta':len(modes),'all12_internal_branches_per_momentum':True,'complete_source_basis_columns':operators,'dimension_check':operators==g.ports,'full_graph_vs12_by12_symbol_residual_below_2e_12':True,'Fourier_orthogonality_guard_pass':True,'floating_control_is_not_outward_spectral_certificate':True,'exact_identity_proof':'Each translated edge contributes its exact Fourier character; the L^3 characters tensor12 internal columns form a complete basis.'}

def cell_packet(g):
    """Exact all-port law from the equal twelve-port initial cell state."""
    a=[int(i<12)for i in range(g.ports)];b=[0]*g.ports;D=g.events;power=1;steps=(0,1,isqrt(D),2*isqrt(D));out=[]
    for n in range(max(steps)+1):
        if n in steps:
            nums=walk.port_numerators(g,a,b);den=12*D*power*power;assert min(nums)>=0 and sum(nums)==den
            cells=[sum(nums[12*i:12*i+12])for i in range(g.cells)];moments=[sum(v*min(c[axis],g.size-c[axis])**2 for v,c in zip(cells,g.cell_coordinates,strict=True))for axis in range(3)];assert moments[0]==moments[1]==moments[2]
            out.append({'coherent_steps':n,'complete_port_outcomes':g.ports,'total_probability_exact':'1','all_port_numerators_sha256':sha256((','.join(map(str,nums))).encode()).hexdigest(),'denominator_sha256':sha256(str(den).encode()).hexdigest(),'all_cell_position_probabilities_interval_exact':[outward((F(v,den),F(v,den)))for v in cells],'all_three_periodic_axis_second_moments_equal_exact':True,'axis_second_moment_interval_exact':outward((F(moments[0],den),F(moments[0],den)))})
        if n<max(steps):a,b=walk.cyclic_step(g,a,b);power*=D
    return {'size':g.size,'preparation':'Equal coherent amplitudes on all12 ports of cell(0,0,0), with uniform retained event register','records':out,'optical_branches_and_unprepared_fixed_event_modes_not_deleted':True}

def refinement_certificates():
    pi=pi_interval();rows=[]
    for L in(16,32,64):
        epsilon=F(1,L);k0=sc(pi,F(2,L));cert=low_branch_certificate((k0,I(0),I(0)),f'first_axis_mode_cube_{L}');lo,hi=map(F,cert['ground_eigenvalue_interval_exact']);D=72*L**3
        omega2=(2*lo/(epsilon*epsilon),2*hi/((epsilon*epsilon)*(1-hi/(2*D))))
        limit=sc(im(pi,pi),1);pmax=2*pi[1];bound=epsilon**2*pmax**4*(F(3,4)+1/(4*D*(1-epsilon**2*pmax*pmax/(4*D))))
        assert omega2[0]>=limit[0]-bound and omega2[1]<=limit[1]+bound
        rows.append({'periodic_size':L,'full_graph_enumerated':False,'complete12_by12_symbol_certified':True,'cell_spacing_exact':str(epsilon),'events_exact':D,'declared_time_per_event':'cell_spacing/sqrt(2M)','ground_band':cert,'continuum_frequency_squared_interval_exact':outward(omega2,15),'continuum_limit_frequency_squared_interval_exact':outward(limit,15),'controlled_difference_upper_exact':str(bound)})
    return rows

def run():
    h=exact_hessian();control=exact_hessian(True);assert control['quadratic_coefficient']!=h['quadratic_coefficient']
    members=[]
    for L in(2,3,4):
        g=graph(L);members.append({'size':L,'cells':g.cells,'ports':g.ports,'physical_seams':len(g.seams),'intercell_seams':len(g.inter_seams),'every_port_degree':6,'connected':True,'macro_topology':'periodic cubic three-torus, declared routing','native_source_columns':native_columns(g),'complete_bloch_decomposition':bloch_control(g),'complete_coherent_event_control':complete_event_control(g),'complete_event_spectrum':walk.full_arc_spectrum(g),'native_scalar_packet':cell_packet(g)})
    angles=(('1/10','0','0'),('1/10','1/10','0'),('1/10','1/10','1/10'),('1/7','-1/9','1/11'),('1/3','1/5','-1/7'),('1/2','0','0'))
    return {'verdict':'PASS_DECLARED_CUBIC_NATIVE_SCALAR_REFINEMENT','local_hardware':architecture.construction_descriptor(),'routing':'Original icosahedral 30 internal seams; each coordinate-zero group supplies two plus and two minus ports; flip the cyclic last nonzero coordinate across a cubic neighbour.','intercell_oriented_port_pairs':[[list(p)for p in row]for row in links()],'macro_dimension_declared':3,'original_geodesic_federation_modified':False,'exact_symmetries':symmetry_certificate(),'low_momentum_exact':h,'anisotropic_partition_control':control,'source_members':members,'held_out_low_momentum_certificates':[low_branch_certificate(tuple(I(x)for x in k),'rational_'+','.join(k))for k in angles],'refinement_family':refinement_certificates(),'controlled_low_momentum_theorem':{'domain':'max_i|k_i|<=1/2','lowest_eigenvalue':'|k|^2/8 plus remainder','absolute_remainder_upper':'(3/8)*max_i|k_i|^4','optical_second_eigenvalue_lower':'4-max_i|k_i|','declared_time_step':'epsilon/sqrt(2M)','prepared_cyclic_frequency_squared_limit':'|p|^2/4, k=epsilon*p','dimensionless_wave_speed_exact':'1/2','all_optical_and_full_event_complements_retained':True},'integration_interface':{'graph_API':'graph(L).ports,seams,neighbours,inter_seams,cell_coordinates','literal_native_endpoint_rule_and_complete25records':'unchanged','integer_Gauss_link_shift':'May be applied on each oriented seam with tagged charge/current; no automatic whole-program transfer claimed.','full_two_load_federation_enumerated':False,'existing_spherical_federation_not_relabelled_three_dimensional':True},'scope':{'three_dimensional_scalar_dispersion_limit_proved_in_stated_preparation_and_clock':True,'whole_interacting_quantum_field_theory_constructed':False,'physical_light_speed_derived':False,'physical_particle_masses_derived':False,'four_dimensional_Yang_Mills_mass_gap':False,'Einstein_gravity_derived':False,'OpenWave_criterion_awarded':False}}
