"""Complete source-derived weighted-triangle contrast/cooling circuits.

Exploratory successor; the source permutation comes only from current columns.
The compiler accepts a code-block norm, never a preferred ray or target phase.
"""
import math
from fractions import Fraction
import numpy as np
from native_source import native
from amplified_native import require, OPTICAL

EDGES = ((0,1),(0,1),(0,2),(0,2),(1,2))


def hadamard(n):
    out = np.ones((1,1))
    while len(out) < n:
        out = np.block([[out,out],[out,-out]])/math.sqrt(2)
    require(len(out)==n,'Hadamard dimension')
    return out


def uniform_five():
    """Exact finite amplitude amplification of five native blank words in eight."""
    h = hadamard(8)
    phase = complex(1/5,2*math.sqrt(6)/5)
    zero = np.eye(8,dtype=complex); zero[0,0] = phase
    good = np.diag([phase]*5+[1.]*3)
    unit = h@zero@h@good@h
    initial = unit[:,0]
    require(np.max(abs(initial[5:]))<1e-13,'Five-ticket selector leakage')
    require(np.max(abs(initial[:5]-initial[0]))<1e-13,'Five-ticket amplitudes unequal')
    require(np.linalg.norm(unit.conj().T@unit-np.eye(8))<1e-13,'Selector not unitary')
    histories = []
    for word in range(8):
        events = []
        for bit in range(3):
            q,record = native((0,0),(0,1),(word>>bit)&1)
            require(q==(0,0) and record[3]==0,'Blank selector source')
            events.append(record)
        histories.append({'word':word,'events':events,'accepted_label':word if word<5 else None})
    return np.kron(unit,hadamard(2)),histories


def weighted_source(law):
    choice,histories = uniform_five()
    width = len(choice)
    reversal = np.eye(3*width,dtype=complex)
    counts = [[Fraction(0)for _ in range(3)]for _ in range(3)]
    records = []
    for x in range(3):
        for ticket,edge in enumerate(EDGES):
            for coin in (0,1):
                q = tuple(int(i==x)for i in range(3))
                audit,event = native(q,edge,coin,law)
                a,b,c,j = event
                decoded = list(q); decoded[a]-=j; decoded[b]+=j
                require(tuple(decoded)==audit,'Current-only output decoder')
                y = decoded.index(1); cc = c^int(j!=0)
                source,target = x*width+2*ticket+c,y*width+2*ticket+cc
                reversal[:,source] = 0; reversal[target,source] = 1
                counts[y][x] += Fraction(1,10)
                records.append([x,ticket,a,b,c,j,y,cc])
    require(np.linalg.norm(reversal@reversal-np.eye(3*width))==0,'Current reversal involution')
    preparation = np.kron(np.eye(3),choice)
    unitary = preparation.conj().T@reversal@preparation
    indices = np.arange(3)*width
    mean = np.array(counts,dtype=float)
    require(np.linalg.norm(unitary[np.ix_(indices,indices)]-mean)<1e-13,'Source block mismatch')
    return {'U':unitary,'P':mean,'counts':counts,'width':width,
            'source_records':records,'selector_histories':histories,
            'native_reversals_per_U':1,'source_law':law}


def contrast_dilation(source):
    """X H diag(I,U_P) H; success code is ancilla0/choice0, all flags retained."""
    u = source['U']; n = len(u)
    controlled = np.block([[np.eye(n),np.zeros((n,n))],[np.zeros((n,n)),u]])
    h = np.kron(hadamard(2),np.eye(n))
    x = np.kron(np.array([[0,1],[1,0]]),np.eye(n))
    unitary = x@h@controlled@h
    width = source['width']
    # Reorder ancilla/position/choice into position/ancilla/choice.
    order = np.arange(2*3*width).reshape(2,3,width).transpose(1,0,2).reshape(-1)
    unitary = unitary[np.ix_(order,order)]
    counts = [[(Fraction(int(i==j))-source['counts'][i][j])/2 for j in range(3)]for i in range(3)]
    return {**source,'U':unitary,'P':np.array(counts,dtype=float),'counts':counts,
            'width':2*width,'source_kind':'contrast_lcu'}


def amplify(source,fixed_norm=None):
    u,width = source['U'],source['width']
    counts = source['counts']
    constant = sum(counts[0])
    trace = sum(counts[i][i]for i in range(3))-constant
    square = sum(counts[i][j]*counts[j][i]for i in range(3)for j in range(3))-constant**2
    disc = 2*square-trace**2
    require(disc>=0,'Invalid record characteristic polynomial')
    s = (abs(float(trace))+math.sqrt(float(disc)))/2 if fixed_norm is None else fixed_norm
    require(.199999999999<s<=1.000000000001,'Unsupported source norm')
    indices = np.arange(3)*width
    zpi = np.eye(len(u),dtype=complex);zpi[indices,indices]=-1
    preboost = s<.5-1e-12
    signal = u@zpi@u@zpi@u if preboost else u
    norm = abs(4*s**3-3*s)if preboost else s
    require(norm>=.5-1e-12,'Insufficient preliminary boost')
    cosine = 1-1/(2*norm**2)
    require(1-cosine*cosine>=-1e-12,'Invalid source matching phase')
    phase = complex(cosine,math.sqrt(max(0.,1-cosine*cosine)))
    mark = np.eye(len(u),dtype=complex);mark[indices,indices]=phase
    complete = signal@mark@signal@mark@signal
    forward = complete[:,indices].reshape(3,width,3).transpose(1,0,2)
    inverse = complete.conj().T[:,indices].reshape(3,width,3).transpose(1,0,2)
    require(np.linalg.norm(complete.conj().T@complete-np.eye(len(u)))<1e-11,'Amplifier unitary')
    for ops in (forward,inverse):
        require(np.linalg.norm(sum(a.conj().T@a for a in ops)-np.eye(3))<1e-11,'Full complement missing')
    block = inverse[0]@forward[0]
    p = source['P']; effective = 4*p@p@p-3*p if preboost else p
    polynomial = np.eye(3)-(np.eye(3)-effective@effective)@np.linalg.matrix_power(np.eye(3)-effective@effective/norm**2,2)
    require(np.linalg.norm(block-polynomial)<1e-10,'Full circuit/polynomial disagreement')
    return {**source,'forward':forward,'inverse':inverse,'K':block,'s':s,
            'norm_record':{'trace':str(trace),'square':str(square),'discriminant':str(disc)},
            'phase':[phase.real,phase.imag],'preboost':preboost,
            'native_reversals_per_filter':18 if preboost else 6,
            'optical_filter':OPTICAL.T@block@OPTICAL}
