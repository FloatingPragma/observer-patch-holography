"""A reversible parity-pump dilation and exact native-kernel finite walk.

Declared owned coherent record + address reflection. Not the committed pump,
not a protected particle, and not a continuum or stationary-state selection.
"""
from fractions import Fraction as F
import sympy as s
from .common import require, local_pump_kraus
from model.common_source import slow_geometry


def certificate():
    # Complete the actual local pure-pump isometry, including double loads.
    # Source order00,10,01,11,20,02 and record order0,1,2,3.
    e=[s.eye(6)[:,k] for k in range(6)];r=[s.eye(4)[:,k] for k in range(4)]
    plus=(e[1]+e[2])/s.sqrt(2);minus=(e[1]-e[2])/s.sqrt(2)
    pairs=[(s.kronecker_product(minus,r[0]),s.kronecker_product(plus,r[1])),
           (s.kronecker_product(e[4],r[0]),s.kronecker_product(e[3],r[2])),
           (s.kronecker_product(e[5],r[0]),s.kronecker_product(e[3],r[3]))]
    U=s.eye(24)
    for a,b in pairs:
        v=a-b;U-=v*v.T
    require(s.simplify(U.T*U)==s.eye(24) and U==U.T,'complete involutory local dilation')
    expected=local_pump_kraus()
    for k,K in enumerate(expected):
        block=s.Matrix(6,6,lambda i,j:U[4*i+k,4*j])
        require(s.simplify(block)==K,'literal pump Kraus rows recovered on blank memory')
    load=s.diag(*[a+b for a,b in [(0,0),(1,0),(0,1),(1,1),(2,0),(0,2)]])
    require(U*s.kronecker_product(load,s.eye(4))==s.kronecker_product(load,s.eye(4))*U,'total load conservation')
    # Actual-carrier coherent address walk. Q_e,p swaps oriented minus0 andplus1.
    seams=list(slow_geometry.source().seams());n=12;M=len(seams);D=2*M
    addresses=[(a,b) for edge in seams for a,b in (edge,edge[::-1])]
    adj=s.zeros(n)
    for a,b in seams:adj[a,b]=adj[b,a]=1
    degree=s.diag(*[sum(adj[i,j] for j in range(n)) for i in range(n)])
    disc=s.zeros(2*n);traces=[];first_U=None
    for a,b in addresses:
        mm=s.zeros(2*n,1);pp=s.zeros(2*n,1)
        mm[2*a]=1/s.sqrt(2);mm[2*b]=-1/s.sqrt(2)
        pp[2*a+1]=pp[2*b+1]=1/s.sqrt(2)
        v=mm-pp;require((v.T*v)[0]==2,'literal unitary reflector normalization')
        Uk=s.eye(2*n)-v*v.T;disc+=Uk;traces.append(s.trace(Uk))
        if first_U is None:first_U=Uk
    disc=s.simplify(disc/D)
    P=s.eye(n)-(degree-adj)/(2*M);Psign=s.eye(n)-(degree+adj)/(2*M)
    target=s.Matrix(2*n,2*n,lambda i,j:(P if i%2==0 else Psign)[i//2,j//2] if i%2==j%2 else 0)
    require(disc==target,'exact reversible pump discriminant')
    # Rational amplitudes normalized by1/sqrt(D); no floating spectral test.
    def idx(a,r,k):return (2*a+r)*D+k
    def reflect(x):
        y=x[:]
        for i in range(2*n):
            mean=sum(x[i*D:(i+1)*D],F(0))/D
            for k in range(D):y[i*D+k]=2*mean-x[i*D+k]
        return y
    def controlled(x):
        y=x[:]
        for k,(a,b) in enumerate(addresses):
            # v=(a0-b0-a1-b1)/sqrt2; U=I-vv*.
            sites=[idx(a,0,k),idx(b,0,k),idx(a,1,k),idx(b,1,k)]
            signs=[1,-1,-1,-1]
            c=sum(sign*x[j] for sign,j in zip(signs,sites))/2
            for sign,j in zip(signs,sites):y[j]-=sign*c
        return y
    def populations(x):return [sum(z*z for z in x[2*a*D:2*(a+1)*D])/D for a in range(n)]
    initial=[F(0)]*(2*n*D)
    for k in range(D):initial[idx(0,0,k)]=F(1)
    first=controlled(reflect(initial));second=controlled(reflect(first))
    p1=populations(first);p2=populations(second)
    require(p1==[F(P[i,0]) for i in range(n)],'first blank-record event has native population shadow')
    require(sum(p2)==1,'two-step exact normalization')
    classical=P*P;difference=F(p2[0])-F(classical[0,0])
    require(difference!=0,'coherent reuse is not original all-time Markov shadow')
    # Source embedding recurrence follows Q²=I and R=2AA*−I:
    # W A=Q A; W Q A=2 Q A disc−A.
    x=s.symbols('x');char=s.factor(P.charpoly(x).as_expr())
    native_expected=(x-1)*(x-s.Rational(9,10))**5*((x-s.Rational(11,12))**2-s.Rational(1,720))**3
    require(s.expand(char-native_expected)==0,'exact carrier native branch spectrum')
    sign_expected=(x-s.Rational(5,6))*(x-s.Rational(14,15))**5*((x-s.Rational(11,12))**2-s.Rational(1,720))**3
    require(s.expand(Psign.charpoly(x).as_expr()-sign_expected)==0,'exact parity-one discriminant spectrum')
    # Rank/intersection count for two reflections:23 interior discriminant modes
    # each produce a reciprocal pair; the one endpoint lambda=1 is unpaired.
    cyclic_dimension=2*(2*n)-1
    full_plus_one=1+(D-(2*n-1))
    full_minus_one=(2*n*D-D)-2*n
    require(full_plus_one+full_minus_one+2*(2*n-1)==2*n*D,'complete unitary spectrum count')
    require(set(traces)=={22},'literal controlled repair traces')
    traceW=(s.Rational(2,D)-1)*sum(traces)
    require(traceW==-1276,'exact coherent walk trace')
    alpha=s.Rational(1,2);dimension=2*n*D
    Z=s.simplify((1+alpha**2)*dimension+2*alpha*traceW)
    require(Z==524,'filtered stationary normalization')
    probability=s.simplify(Z/(dimension*(1+alpha)**2))
    require(probability==s.Rational(131,810),'complete heralded filter success probability')
    diagonal_inside=s.Rational(1,2)*(s.Rational(2,D)-1)
    diagonal_outside=s.Rational(2,D)-1
    a,b=addresses[0]
    require(first_U[2*a,2*a]==s.Rational(1,2),'literal endpoint keep diagonal')
    outside=next(i for i in range(n) if i not in (a,b))
    require(first_U[2*outside,2*outside]==1,'literal outside keep diagonal')
    require(first_U[2*a,2*b]==s.Rational(1,2),'literal offdiagonal readout witness')
    rho_offdiagonal=s.simplify(alpha*2*first_U[2*a,2*b]*(s.Rational(2,D)-1)/Z)
    require(rho_offdiagonal==-s.Rational(29,31440),'available record projector noncommutes with prepared state')
    return {'status':'EXACT local total-two unitary completion and actual one-load finite walk',
      'local_unitary_dimension':24,'local_record_dimension':4,'all_four_blank_record_kraus_blocks_match':True,
      'source_port_dimension':n,'one_load_parity_record_dimension':2,'coherent_ordered_seam_addresses':D,
      'full_walk_dimension':2*n*D,'walk':'W=Q(2AA*-I), A=I_source,parity tensor uniform_ordered_seam',
      'discriminant_blocks':['I-L/(2*M)','I-(Degree+Adjacency)/(2*M)'],
      'blank_parity_discriminant':'exact native one-load kernel',
      'native_kernel_eigenvalues':['1 multiplicity1','11/12+sqrt(5)/60 multiplicity3','11/12-sqrt(5)/60 multiplicity3','9/10 multiplicity5'],
      'walk_phase_relation':'cos(omega)=lambda; reciprocal pairs away from endpoints',
      'full_walk_eigenvalue_plus_one_multiplicity':full_plus_one,
      'full_walk_eigenvalue_minus_one_multiplicity':full_minus_one,
      'full_walk_oscillatory_dimension':2*(2*n-1),
      'full_A_cyclic_subspace_dimension':cyclic_dimension,
      'blank_parity_A_cyclic_subspace_dimension':2*n-1,
      'auxiliary_record_modes_are_not_particle_species':True,
      'first_event_source_populations_exact_native':True,
      'second_event_source_return_probability':str(p2[0]),
      'classical_two_event_source_return_probability':str(classical[0,0]),
      'second_event_return_difference':str(difference),
      'all_time_native_population_shadow':False,
      'trace_W':str(traceW),
      'W_diagonal_endpoint':str(diagonal_inside),'W_diagonal_outside':str(diagonal_outside),
      'stationary_filter':{'K_success':'(2I+W)/3','K_failure':'sqrt(2)*(I-W)/3',
          'rho':'((5/4)I+(W+W*)/2)/524',
          'normalization':str(Z),'eigenvalue_lower':str((1-alpha)**2/Z),
          'success_probability_from_maximally_mixed':str(probability),
          'independent_fresh_restart_failure_tail':'(679/810)^k',
          'rho_a0k_b0k':str(rho_offdiagonal),
          'scope':'Declared finite stationary preparation for the same W, with explicit outcome and restart cost; not inherited from the committed pump or a modular/physical-time identification.'},
      'added_primitives':['owned coherent outcome memory','coherent seam/endpoint address','uniform-address reflection','reversible record reuse'],
      'faithful_phase_pump_stationary_state_inherited':False,
      'protected_particle_or_continuum_claim':False,
      'scope':'A coherent architectural alternative, exact finite wave recurrence and explicit faithful stationary preparation. Full apparatus control, physical clock interpretation, federation locality, protected sectors and particle response remain unproved.'}
