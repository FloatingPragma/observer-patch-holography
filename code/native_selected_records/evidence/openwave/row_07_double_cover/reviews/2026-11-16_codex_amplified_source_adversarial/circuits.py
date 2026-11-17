"""Independent circuit reconstruction for the finite source review.

No producer import, state-specific ray or desired rotation sign is an input.
The source law, signed endpoint currents, apparatus coefficients and the
public choice/address program completely determine the matrices below.
"""
from decimal import Decimal, localcontext, ROUND_HALF_EVEN
from fractions import Fraction as F
import math
import numpy as np

EDGES=((0,1),(0,1),(0,2),(0,2),(1,2))
TRI=((0,1),(1,2),(0,2))


def need(ok,message):
    if not ok:raise AssertionError(message)


def array(value):
    return np.asarray(value['real'])+1j*np.asarray(value['imag']) if isinstance(value,dict) else np.asarray(value)


def near(a,b,message,tolerance=2e-10):
    error=float(np.max(abs(np.asarray(a)-np.asarray(b))))
    need(error<tolerance,message+': '+str(error))
    return error


def hadamard(width):
    need(width>0 and not width&(width-1),'Choice width is not dyadic')
    return np.array([[(-1)**((i&j).bit_count()) for j in range(width)]for i in range(width)])/math.sqrt(width)


def native_column(position,edge,coin,law,sign):
    q=[sign*int(i==position)for i in range(3)]
    a,b=edge;difference=q[a]-q[b]
    j=(difference-(2*coin-1)*(abs(difference)%2))//2 if law=='native' else difference
    q[a]-=j;q[b]+=j
    return q.index(sign),j,coin^int(j!=0)


def reversal(edges,width,law,sign):
    permutation=list(range(3*width));counts=[[F(0)for _ in range(3)]for _ in range(3)]
    columns=[]
    for position in range(3):
        for ticket,edge in enumerate(edges):
            for coin in (0,1):
                out,j,feedback=native_column(position,edge,coin,law,sign)
                permutation[position*width+2*ticket+coin]=out*width+2*ticket+feedback
                counts[out][position]+=F(1,2*len(edges))
                columns.append((position,ticket,edge,coin,j,out,feedback))
    need(sorted(permutation)==list(range(3*width)),'Literal source permutation')
    need(all(permutation[permutation[i]]==i for i in range(3*width)),'Literal source involution')
    return np.array(permutation),counts,columns


def phase(counts,bits=40):
    constant=sum(counts[0]);trace=sum(counts[i][i]for i in range(3))-constant
    square=sum(counts[i][j]*counts[j][i]for i in range(3)for j in range(3))-constant**2
    disc=2*square-trace**2
    need(disc>=0,'Source discriminant')
    with localcontext()as context:
        context.prec=110
        convert=lambda x:Decimal(x.numerator)/Decimal(x.denominator)
        norm=(abs(convert(trace))+convert(disc).sqrt())/2
        preboost=norm<Decimal('0.5')
        boosted=abs(4*norm**3-3*norm)if preboost else norm
        need(boosted>=Decimal('0.5'),'Norm outside one-boost domain')
        if boosted==Decimal('0.5'):return F(-1),F(0),preboost
        tangent=(1/(4*boosted**2-1)).sqrt()
        integer=int((tangent*(2**bits)).to_integral_value(rounding=ROUND_HALF_EVEN))
    t=F(integer,2**bits)
    return (1-t*t)/(1+t*t),2*t/(1+t*t),preboost


def weighted_source(law,sign,selector_real,selector_imag):
    re,im=F(selector_real),F(selector_imag);need(re*re+im*im==1,'Finite selector phase norm')
    z=complex(float(re),float(im));h=hadamard(8)
    zero=np.eye(8,dtype=complex);zero[0,0]=z
    good=np.diag([z]*5+[1]*3)
    choice=np.kron(h@zero@h@good@h,hadamard(2));width=16
    permutation,base,columns=reversal(EDGES,width,law,sign)
    r=np.eye(3*width)[permutation]
    prepared=np.kron(np.eye(3),choice)
    u=prepared.conj().T@r@prepared
    ar=5*(re*re-im*im)+6*re-3;ai=10*re*im+6*im
    br=5*(re*re-im*im)-2*re+5;bi=10*re*im-2*im
    active=5*(ar*ar+ai*ai)/512;inactive=3*(br*br+bi*bi)/512
    need(active+inactive==1,'Finite selector completeness')
    counts=[[active*base[i][j]+inactive*int(i==j)for j in range(3)]for i in range(3)]
    return u,counts,columns,active,inactive


def contrast(u,counts):
    n=len(u);eye=np.eye(n);a=(eye+u)/2;b=(eye-u)/2
    full=np.block([[b,a],[a,b]])
    width=n//3
    order=[anc*3*width+p*width+j for p in range(3)for anc in range(2)for j in range(width)]
    return full[np.ix_(order,order)], [[(F(int(i==j))-counts[i][j])/2 for j in range(3)]for i in range(3)]


def star_source(slot,m,law,sign):
    vertex,k=divmod(slot%(3*m),m)
    edges=[e for e in TRI if vertex in e]*(m-k)+[e for e in TRI if (vertex+1)%3 in e]*k
    width=4*m;perm,counts,columns=reversal(edges,width,law,sign)
    h=hadamard(width)
    def operation(x):
        blocks=np.einsum('ij,pjk->pik',h,x.reshape(3,width,-1),optimize=True).reshape(3*width,-1)
        blocks=blocks[perm].reshape(3,width,-1)
        return np.einsum('ij,pjk->pik',h,blocks,optimize=True).reshape(3*width,-1)
    return operation,width,counts,columns


def amplified(primitive,width,counts,re,im,preboost):
    u=primitive if callable(primitive) else lambda x:primitive@x
    zero=np.arange(3)*width
    start=np.zeros((3*width,3),complex);start[zero,np.arange(3)]=1
    z=complex(float(re),float(im))
    def marked(x,phase_value):
        y=x.copy();y[zero]*=phase_value;return y
    def signal(x):
        return u(marked(u(marked(u(x),-1)),-1))if preboost else u(x)
    def sequence(value):return signal(marked(signal(marked(signal(start),value)),value))
    f=sequence(z).reshape(3,width,3).transpose(1,0,2)
    g=sequence(z.conjugate()).reshape(3,width,3).transpose(1,0,2)
    near(sum(x.conj().T@x for x in f),np.eye(3),'Independent forward completeness')
    near(sum(x.conj().T@x for x in g),np.eye(3),'Independent inverse completeness')
    near(g[0],f[0].conj().T,'Accepted adjoint')
    p=np.array(counts,dtype=float);t=4*p@p@p-3*p if preboost else p
    alpha=2*(1-float(re));residual=np.eye(3)-alpha*t@t
    polynomial=np.eye(3)-(np.eye(3)-t@t)@residual@residual
    near(g[0]@f[0],polynomial,'Independent finite polynomial')
    return f,g,g[0]@f[0]
