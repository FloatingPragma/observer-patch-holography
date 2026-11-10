"""Exact same-carrier kinematic coupling of reversible repair and encoded gauge transport.

Uses local20-dimensional symbolic matrices and actual60-dimensional connection
blocks; never allocates the7200-dimensional walk. A static central link is a
declared background, not an interaction, protected defect, spin or particle mass.
"""
import json, hashlib
import sympy as s
from model.common_source import slow_geometry
from .holonomy import certificate as holonomy_certificate


def require(ok,label):
    if not bool(ok):raise ArithmeticError(label)


def local(G):
    d=G.rows;I=s.eye(d)
    m=I.col_join(-G)/s.sqrt(2);p=I.col_join(G)/s.sqrt(2)
    v=s.kronecker_product(m,s.Matrix([1,0]))-s.kronecker_product(p,s.Matrix([0,1]))
    return s.simplify(s.eye(4*d)-v*v.H),m,p


def certificate():
    sx=s.Matrix([[0,1],[1,0]]);sy=s.Matrix([[0,-s.I],[s.I,0]]);sz=s.diag(1,-1)
    Q=s.diag(2,2,2,-3,-3);Z=s.diag(1,1,1,-1,-1)
    # Noncommuting exact controls test the general unitary-block identities;
    # only the central Z below is attached to the actual compiled link word.
    X=s.diag(s.diag(1,-1,-1),s.I*sz)
    Y=s.diag(s.diag(-1,1,-1),s.I*sy)
    G=X;Va=X;Vb=Y
    require(G.H*G==s.eye(5) and G.det()==1 and G*Q==Q*G,'unitary charge-preserving local link')
    U,m,p=local(G)
    require(m.H*m==s.eye(5) and p.H*p==s.eye(5) and m.H*p==s.zeros(5),'orthogonal transported parity isometries')
    require(s.simplify(U.H*U)==s.eye(20) and U==U.H,'local unitary completion')
    source_swap=s.kronecker_product(s.Matrix([[0,1],[1,0]]),s.eye(5))
    swap=s.kronecker_product(source_swap,s.eye(2))
    reversed_U=swap*local(G.H)[0]*swap
    expected=s.eye(20)-s.kronecker_product(m*m.H,s.diag(1,0))-s.kronecker_product(p*p.H,s.diag(0,1))
    require(s.simplify((U+reversed_U)/2-expected)==s.zeros(20),'orientation average cancels parity cross blocks')
    Vsource=s.diag(Va,Vb);V=s.kronecker_product(Vsource,s.eye(2));Gp=Vb*G*Va.H
    require(s.simplify(local(Gp)[0]-V*U*V.H)==s.zeros(20),'local gauge covariance')
    Qlocal=s.kronecker_product(s.eye(2),Q,s.eye(2))
    require(U*Qlocal==Qlocal*U,'exact additional internal charge conservation')
    require(s.trace(U)==10,'local20-dimensional rank-five flip trace')
    # Read the generated actual central history, not a fitted matrix table.
    gauge=holonomy_certificate();require(gauge['central_word']==[0,0],'actual compiled central word')
    central=next(row for row in gauge['compiler'] if row['word']==gauge['central_word'])
    require(central['carrier_permutation']==list(range(12)) and central['central_bit']==1,'central word has trivial port action and nontrivial cover sign')
    require(central['quaternion']==[['-1','0'],['0','0'],['0','0'],['0','0']],'central word is quaternion minus one')
    seams=list(slow_geometry.source().seams());n=12;M=len(seams);addresses=2*M;internal=5
    central_edge=seams[0]
    A=s.zeros(n);Az=s.zeros(n)
    for edge in seams:
        a,b=edge;A[a,b]=A[b,a]=1;Az[a,b]=Az[b,a]=-1 if edge==central_edge else 1
    # The full60x60 adjacency is independently assembled in port-major order.
    AG=s.zeros(n*internal);Qtotal=s.kronecker_product(s.eye(n),Q)
    for edge in seams:
        a,b=edge;link=Z if edge==central_edge else s.eye(internal)
        AG[internal*b:internal*(b+1),internal*a:internal*(a+1)]=link
        AG[internal*a:internal*(a+1),internal*b:internal*(b+1)]=link.H
    require(AG==AG.H and AG*Qtotal==Qtotal*AG,'Hermitian charge-conserving actual connection adjacency')
    require(all(AG[5*a+j,5*b+k]==(A if j<3 else Az)[a,b]*int(j==k) for a in range(n) for b in range(n) for j in range(5) for k in range(5)),'actual60x60 source operator equals its five checked blocks')
    A3=A**3;Az3=Az**3
    flat_trace=5*s.trace(A3);flux_trace=3*s.trace(A3)+2*s.trace(Az3)
    require(flat_trace==600 and flux_trace==552,'actual static-flux triangle invariant')
    require((5*s.eye(n)-A).nullspace().__len__()==1,'flat triplet constant section')
    require((5*s.eye(n)-Az).det()>0 and (5*s.eye(n)-Az).rank()==12,'frustrated doublet has no parallel zero section')
    trial=s.ones(n,1)
    rayleigh=(trial.T*(5*s.eye(n)-Az)*trial)[0]/n
    require(rayleigh==s.Rational(4,n),'single-flux variational finite-volume bound')
    # Exact embedded three-step return amplitudes: A_embed*W^3 A_embed=T3(P).
    a=s.Rational(11,12);b=s.Rational(1,60)
    P3=a**3*s.eye(n)+3*a*a*b*A+3*a*b*b*(A**2)+b**3*A3
    Pz3=a**3*s.eye(n)+3*a*a*b*Az+3*a*b*b*(Az**2)+b**3*Az3
    P=a*s.eye(n)+b*A;Pz=a*s.eye(n)+b*Az
    p,q=central_edge
    flat_amp=(4*P3-3*P)[p,p];flux_amp=(4*Pz3-3*Pz)[p,p]
    require(flux_amp-flat_amp==-s.Rational(1,6750),'actual endpoint embedded return amplitude changes')
    require((flux_trace-flat_trace)/s.Integer(60)**3==-s.Rational(1,4500),'exact trace cubic kernel change')
    triangles=[(i,j,k) for i in range(n) for j in range(i+1,n) for k in range(j+1,n) if A[i,j]*A[j,k]*A[k,i]]
    flux_faces=[f for f in triangles if set(central_edge).issubset(f)]
    require(len(triangles)==20 and len(flux_faces)==2,'exact two nontrivial triangular Wilson loops')
    source_parity_dim=2*n*internal;dimension=source_parity_dim*addresses
    primitive_trace=source_parity_dim-2*internal
    traceW=(2-addresses)*primitive_trace
    normal=s.Rational(5,4)*dimension+traceW
    success=s.simplify(normal/(dimension*s.Rational(9,4)))
    require((dimension,primitive_trace,traceW,normal,success)==(7200,110,-6380,2620,s.Rational(131,810)),'same coupled-source stationary filter numbers')
    return {'status':'EXACT finite SAME-carrier reversible-repair and encoded-transport coupling',
      'declared_architecture_changes':['five-level encoded fibre attached to the moving one-load source','transported parity seam rule','fixed internal link background','owned coherent address reflection and complete heralded preparation'],
      'internal_charge':'Q5=diag(2,2,2,-3,-3)','internal_charge_conserved':True,
      'general_local_control_dimension':20,'general_local_unitarity_and_covariance_checked':True,
      'carrier_ports':n,'seams':M,'coupled_walk_dimension':dimension,
      'ordered_seam_reflector_negative_rank':internal,'ordered_seam_reflector_trace':primitive_trace,
      'discriminant_blocks':['I-L_G/(2*M)','I-(Degree+A_G)/(2*M)'],
      'same_coupled_W_trace':traceW,
      'same_coupled_reference':{'rho':'((5/4)I+(W_G+W_G*)/2)/2620','eigenvalue_lower':'1/10480','success_probability':str(success),'failure_tail':'(679/810)^k','normalization':str(normal),'uncoupled_reference_inherited':False},
      'central_link':{'edge':list(central_edge),'matrix':'diag(I3,-I2)','actual_compiled_word':gauge['central_word'],'quaternion':central['quaternion'],'source_receipt_sha256':hashlib.sha256(json.dumps(gauge,sort_keys=True,separators=(',',':')).encode()).hexdigest()},
      'static_flux_witness':{'nontrivial_triangular_faces':flux_faces,'Wilson_loop_trace':str(s.trace(Z)),'flat_Wilson_loop_trace':'5','trace_A_cubed_flat':str(flat_trace),'trace_A_cubed_flux':str(flux_trace),'trace_P_cubed_difference':'-1/4500','trace_embedded_W_cubed_difference':'-1/1125','doublet_endpoint_three_step_return_amplitude_flat':str(flat_amp),'doublet_endpoint_three_step_return_amplitude_flux':str(flux_amp),'amplitude_difference':str(flux_amp-flat_amp),'probability_difference':str(s.factor(flux_amp**2-flat_amp**2))},
      'parallel_sections':{'triplet_kernel_dimension':3,'doublet_kernel_dimension':0,'finite_doublet_L_determinant':str((5*s.eye(n)-Az).det()),'determinant_scope':'one12-dimensional doublet component','full24_doublet_L_determinant':str((5*s.eye(n)-Az).det()**2),'physical_rest_mass_claim':False,
          'single_flux_Rayleigh_quotient':str(rayleigh),
          'general_single_negative_edge_bound':'lambda_min(L_signed) <= 4/N',
          'walk_phase_upper':'arccos(1-2/(M*N)) in the complete-attempt walk clock',
          'uniform_physical_mass_gap_inferred':False},
      'gauge_scope':'Covariance holds for arbitrary block-unitary passive frames with transformed links. Local frames in the finite120-element cover preserve its source link alphabet. A common global continuous block center also preserves every link. Arbitrary independent continuous local frames require an enlarged link orbit.',
      'readout_scope':'Prepared source-return probabilities change. Full7200-dimensional tomography does not follow from the1440-dimensional support proof; conserved internal charges split the support and require the code/referral apparatus contract.',
      'scope':{'kinematic_gauge_coupling_on_same_source_state_clock':True,'dynamical_gauge_field_or_backreaction':False,'protected_flux_defect':False,'interacting_many_particle_scattering':False,'physical_spin_or_statistics':False,'continuum_or_GeV_mass':False}}
