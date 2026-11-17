"""Finite coefficient implementation; every matching phase is exactly rational.

Hadamard, controlled addressing and phase memory operations are declared
coherent services. A blank-record census is not an execution of Hadamard.
"""
import math
from decimal import Decimal, localcontext, ROUND_HALF_EVEN
from fractions import Fraction as F
import numpy as np
from native_source import native
from amplified_native import require,routing,hadamard_choices
from selected_source import hadamard,contrast_dilation,EDGES

BITS=40


def qadd(x,y):return x[0]+y[0],x[1]+y[1]
def qmul(x,y,d):return x[0]*y[0]+x[1]*y[1]*d,x[0]*y[1]+x[1]*y[0]
def qscale(x,c):return x[0]*c,x[1]*c
def qsign(x,d):
    a,b=x
    if not b or not d:return (a>0)-(a<0)
    if not a:return (b>0)-(b<0)
    if (a>0)==(b>0):return (a>0)-(a<0)
    difference=a*a-b*b*d
    return ((difference>0)-(difference<0))*((a>0)-(a<0))
def qinverse(x,d):
    denominator=x[0]*x[0]-x[1]*x[1]*d
    require(denominator!=0,'Quadratic inverse needs rational simplification')
    return x[0]/denominator,-x[1]/denominator
def qrecord(x,d):return {'rational':str(x[0]),'sqrt_coefficient':str(x[1]),'radicand':str(d)}


def decimal(value):
    return Decimal(value.numerator)/Decimal(value.denominator)


def phase_from_tangent(tangent,exact_square,radicand=F(0)):
    with localcontext()as ctx:
        ctx.prec=90
        numerator=int((tangent*(1<<BITS)).to_integral_value(rounding=ROUND_HALF_EVEN))
    t=F(numerator,1<<BITS)
    lower=F(2*numerator-1,1<<(BITS+1));upper=F(2*numerator+1,1<<(BITS+1))
    lo=qadd(exact_square,(-lower*lower,F(0)))
    hi=qadd((upper*upper,F(0)),qscale(exact_square,-1))
    require(lower>=0 and qsign(lo,radicand)>=0 and qsign(hi,radicand)>=0,'Chosen tangent is outside exact nearest cell')
    require((qsign(lo,radicand)!=0 and qsign(hi,radicand)!=0)or numerator%2==0,'Nearest tie is not even')
    re=(1-t*t)/(1+t*t);im=2*t/(1+t*t)
    require(re*re+im*im==1,'Rational phase unitarity')
    return complex(float(re),float(im)),{'bits':BITS,'numerator':numerator,
        'tangent':str(t),'real':str(re),'imag':str(im),'nearest_rounding':'half_even',
        'exact_tangent_square':qrecord(exact_square,radicand),
        'nearest_cell_lower_square_margin':qrecord(lo,radicand),
        'nearest_cell_upper_square_margin':qrecord(hi,radicand),
        'certified_phase_operator_error_from_ideal_upper':str(F(1,1<<BITS))},(re,im)


def matching_phase(counts,fixed=None):
    constant=sum(counts[0]);trace=sum(counts[i][i]for i in range(3))-constant
    square=sum(counts[i][j]*counts[j][i]for i in range(3)for j in range(3))-constant**2
    discriminant=2*square-trace**2
    require(discriminant>=0,'Optical characteristic discriminant')
    record={'trace':str(trace),'square':str(square),'discriminant':str(discriminant)}
    rn,rd=math.isqrt(discriminant.numerator),math.isqrt(discriminant.denominator)
    if rn*rn==discriminant.numerator and rd*rd==discriminant.denominator:
        d=F(0);exact_s=((abs(trace)+F(rn,rd))/2,F(0))
    else:d=discriminant;exact_s=(abs(trace)/2,F(1,2))
    require(qsign(qadd(exact_s,(-F(1,5),F(0))),d)>=0 and qsign(qadd(exact_s,(-F(1),F(0))),d)<=0,
        'Source norm outside declared [1/5,1] compiler domain')
    preboost=qsign(qadd(exact_s,(-F(1,2),F(0))),d)<0
    signal=qadd(qscale(qmul(qmul(exact_s,exact_s,d),exact_s,d),4),qscale(exact_s,-3))if preboost else exact_s
    denominator=qadd(qscale(qmul(signal,signal,d),4),(-F(1),F(0)))
    require(qsign(denominator,d)>=0,'Exact amplification norm is below half')
    with localcontext()as ctx:
        ctx.prec=90
        s=(abs(decimal(trace))+decimal(discriminant).sqrt())/2
        if fixed is not None:
            return fixed['phase_complex'],fixed['phase_record'],fixed['preboost'],fixed['calibration_s']
        signal=abs(4*s**3-3*s)if preboost else s
        if qsign(denominator,d)==0:
            z=complex(-1,0);phase_record={'bits':BITS,'endpoint':'minus_one','real':'-1','imag':'0',
                'exact_endpoint_certificate':qrecord(denominator,d),
                'certified_phase_operator_error_from_ideal_upper':'0'}
        else:
            require(signal>Decimal('.5'),'Decimal proposal contradicts exact nonboundary admission')
            z,phase_record,_=phase_from_tangent(1/(4*signal**2-1).sqrt(),qinverse(denominator,d),d)
        phase_record={**phase_record,'source_characteristic':record}
    return z,phase_record,preboost,str(s)


def current_reversal(edges,law,width,columns=None,sign=1):
    """V, coin feedback, V, current-add, V†; all workspaces finish blank.

V† is declared coherent inverse access to the recorded source isometry.
Two literal forward events alone return the load to its initial position.
"""
    require(sign in(-1,1),'Source charge sign')
    if columns is not None:
        expected={(p,e,c)for p in range(3)for e in((0,1),(1,2),(0,2))for c in(0,1)}
        keys=[(r['input_position'],tuple(r['edge']),r['coin'])for r in columns]
        require(len(keys)==18 and len(set(keys))==18 and set(keys)==expected,'Calibration columns have missing, duplicate or extra coverage')
        require(all(r.get('sign',sign)==sign and r.get('law',law)==law for r in columns),'Calibration source metadata mismatch')
    observed=None if columns is None else {(r['input_position'],tuple(r['edge']),r['coin']):r for r in columns}
    permutation=np.arange(3*width);records=[]
    counts=[[F(0)for _ in range(3)]for _ in range(3)]
    for x in range(3):
        for ticket,edge in enumerate(edges):
            for c in (0,1):
                q=tuple(sign*int(p==x)for p in range(3))
                if observed is None:
                    audit,event=native(q,edge,c,law)
                else:
                    row=observed[x,tuple(edge),c];event=[edge[0],edge[1],c,row['current']]
                    audit=tuple(sign*int(p==row['decoded_output'])for p in range(3))
                a,b,_,j=event;decoded=list(q);decoded[a]-=j;decoded[b]+=j
                require(tuple(decoded)==audit,'Current reconstruction')
                y=decoded.index(sign);cc=c^int(j!=0)
                if observed is None:
                    reverse_audit,reverse=native(tuple(decoded),edge,cc,law)
                else:
                    row=observed[y,tuple(edge),cc];reverse=[edge[0],edge[1],cc,row['current']]
                    reverse_audit=tuple(sign*int(p==row['decoded_output'])for p in range(3))
                require(reverse[3]==-j and reverse_audit==q,'Current workspace not erased')
                inverse_restored=list(reverse_audit)
                inverse_restored[a]+=reverse[3];inverse_restored[b]-=reverse[3]
                require(tuple(inverse_restored)==tuple(decoded),'Inverse-source final load')
                permutation[x*width+2*ticket+c]=y*width+2*ticket+cc
                counts[y][x]+=F(1,2*len(edges))
                records.append({'input':x,'ticket':ticket,'forward':event,'output':y,
                    'feedback_coin':cc,'reverse':reverse,'final_current_workspace':j+reverse[3],
                    'charge_sign':sign,'after_two_forward_position':reverse_audit.index(sign),
                    'after_inverse_position':inverse_restored.index(sign),
                    'source_accesses':{'forward':2,'inverse':1},'inverse_scratch_final':0})
    require(np.array_equal(np.sort(permutation),np.arange(3*width)),'Source permutation')
    require(np.array_equal(permutation[permutation],np.arange(3*width)),'Source involution')
    return permutation,counts,records


def weighted(law,columns=None,sign=1,blank_evidence=None):
    if blank_evidence is not None:
        require({r[2]for r in blank_evidence}=={0,1}and all(r[:2]==[29,30]and r[3]==0 for r in blank_evidence),
            'Measured blank current coupon coverage')
    with localcontext()as ctx:
        ctx.prec=90
        phase,phase_record,(re,im)=phase_from_tangent((Decimal(2)/3).sqrt(),(F(2,3),F(0)))
    h=hadamard(8);zero=np.eye(8,dtype=complex);zero[0,0]=phase
    good=np.diag([phase]*5+[1.]*3)
    choice=np.kron(h@zero@h@good@h,hadamard(2))
    # Exactly derived finite-choice weights; inactive labels are not discarded.
    z2re,z2im=re*re-im*im,2*re*im
    ar,ai=5*z2re+6*re-3,5*z2im+6*im
    active=5*(ar*ar+ai*ai)/512
    br,bi=5*z2re-2*re+5,5*z2im-2*im
    inactive=3*(br*br+bi*bi)/512
    require(active+inactive==1 and 0<active<=1,'Exact selector weights')
    width=16;permutation,base,records=current_reversal(EDGES,law,width,columns,sign)
    blank_records=[]
    for position in range(3):
        for label in range(5,8):
            for coin in (0,1):
                first,event=native((0,0),(0,1),coin,law)
                second,reverse=native(first,(0,1),coin,law)
                require(first==second==(0,0) and event[3]==reverse[3]==0,'Inactive-label blank source')
                blank_records.append({'target_spectator':position,'label':label,'coin':coin,
                    'physical_blank_edge':[29,30],'forward':event,'second_forward':reverse,
                    'inverse_on_image_final':[0,0],'all_current_workspaces':0,
                    'meaning':'Declared-vacuum source-column audit; separate actual blank records are attached, not replaced by this call.'})
    reversal=np.eye(3*width)[permutation]
    preparation=np.kron(np.eye(3),choice)
    unit=preparation.conj().T@reversal@preparation
    counts=[[active*base[i][j]+inactive*int(i==j)for j in range(3)]for i in range(3)]
    p=np.array(counts,dtype=float);indices=np.arange(3)*width
    require(np.linalg.norm(unit[np.ix_(indices,indices)]-p)<1e-12,'Finite weighted source block')
    return {'U':unit,'P':p,'counts':counts,'width':width,'source_records':records,'inactive_source_records':blank_records,
        'selector_phase':phase_record,'selector_active_weight':str(active),
        'selector_inactive_weight':str(inactive),'native_reversals_per_U':1,
        'source_law':law,'source_kind':'finite_weighted','charge_sign':sign,
        'operational_blank_records':blank_evidence,'vacuum_assumption':'Supplied and admission-preserved; not derived from zero current.'}


def instrument(source,fixed=None):
    width=source['width'];mean=source['P'];counts=source['counts']
    phase,phase_record,preboost,s=matching_phase(counts,fixed)
    initial=np.zeros((3*width,3),complex);initial[np.arange(3)*width,np.arange(3)]=1
    if 'U'in source:
        def primitive(x):return source['U']@x
    else:
        def primitive(x):
            values=hadamard_choices(x,width);out=np.empty_like(values)
            out[source['permutation']]=values
            return hadamard_choices(out,width)
    def mark(x,z):
        out=x.copy();out[np.arange(3)*width]*=z;return out
    def signal(x):
        return primitive(mark(primitive(mark(primitive(x),-1)),-1))if preboost else primitive(x)
    def full(z):return signal(mark(signal(mark(signal(initial),z)),z))
    forward=full(phase).reshape(3,width,3).transpose(1,0,2)
    inverse=full(phase.conjugate()).reshape(3,width,3).transpose(1,0,2)
    for ops in (forward,inverse):
        require(np.linalg.norm(sum(op.conj().T@op for op in ops)-np.eye(3))<1e-10,'Complete finite instrument')
    block=inverse[0]@forward[0]
    effective=4*mean@mean@mean-3*mean if preboost else mean
    alpha=2*(1-phase.real)
    predicted=np.eye(3)-(np.eye(3)-effective@effective)@np.linalg.matrix_power(np.eye(3)-alpha*effective@effective,2)
    require(np.linalg.norm(block-predicted)<1e-10,'Finite compiler polynomial')
    return {**source,'K':block,'forward':forward,'inverse':inverse,'phase_complex':phase,
        'phase_record':phase_record,'preboost':preboost,'calibration_s':s,
        'native_reversals_per_filter':18 if preboost else 6,
        'source_forward_calls_per_filter':27 if preboost else 9,
        'source_inverse_calls_per_filter':27 if preboost else 9}


def star(slot,m,law,fixed=None,columns=None,sign=1):
    spin_to_energy=(1,2,0)
    edges=[tuple(sorted((spin_to_energy[a],spin_to_energy[b])))for a,b in routing(slot,m)]
    width=2*len(edges)
    permutation,counts,records=current_reversal(edges,law,width,columns,sign)
    return instrument({'P':np.array(counts,dtype=float),'counts':counts,'width':width,
        'permutation':permutation,'source_records':records,'source_law':law,
        'source_kind':'star','slot':slot,'M':m,'charge_sign':sign,
        'coordinate_system':'energy_target_order','physical_target_ports':[0,24,25],
        'spin_marker_to_energy_vertex':list(spin_to_energy)},fixed)


def scalar_control(compiled):
    """Exact positive scalar for swap on the complete optical sector."""
    p=compiled['counts'];trace=sum(p[i][i]for i in range(3))-1
    require(trace==0,'Control is not traceless on optical plane')
    u=(sum(p[i][j]*p[j][i]for i in range(3)for j in range(3))-1)/2
    if compiled['preboost']:u=u*(4*u-3)**2
    re=F(compiled['phase_record']['real']);alpha=2*(1-re)
    scalar=1-(1-u)*(1-alpha*u)**2
    require(0<scalar<=1,'Nonpositive scalar control')
    return scalar
