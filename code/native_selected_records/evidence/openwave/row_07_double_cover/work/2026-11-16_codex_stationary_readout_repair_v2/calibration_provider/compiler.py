"""Derive classical means and conditional coherent maps from measured columns only."""
from fractions import Fraction as F
EDGES=((0,1),(1,2),(0,2))


def require(ok,message):
    if not ok:raise ValueError(message)


def table(columns):
    require(len(columns)==18,'exact18 observed columns')
    indexed={}
    for row in columns:
        p,c,y,j=row['input_position'],row['coin'],row['decoded_output'],row['current']
        edge=tuple(row['edge'])
        require(type(p) is int and p in range(3) and type(c) is int and c in (0,1),'column input')
        require(edge in EDGES and type(y) is int and y in range(3) and type(j) is int and j in (-1,0,1),'column output')
        key=(p,EDGES.index(edge),c)
        require(key not in indexed,'duplicate column');indexed[key]=(y,j)
    require(set(indexed)=={(p,e,c) for p in range(3) for e in range(3) for c in (0,1)},'complete column domain')
    return indexed


def compile_menu(columns,tickets,weights=None,blank_records=None):
    indexed=table(columns);tickets=list(tickets);n=len(tickets)
    require(n>0 and all(type(e) is int and e in range(3) for e in tickets),'declared tickets')
    weights=[F(1,n)]*n if weights is None else [F(w) for w in weights]
    require(len(weights)==n and min(weights)>=0 and sum(weights)<=1,'declared ticket probabilities')
    inactive=1-sum(weights)
    if inactive:
        require(blank_records is not None,'inactive reference must be checked')
        require({r[2] for r in blank_records}=={0,1} and all(r[3]==0 for r in blank_records),'observed blank coins/current')
    P=[[F(inactive*int(i==j)) for j in range(3)] for i in range(3)]
    gram=[[F(inactive*int(i==j)) for j in range(3)] for i in range(3)]
    primitive=[]
    for ticket,e in enumerate(tickets):
        for coin in (0,1):
            rows=[indexed[p,e,coin] for p in range(3)]
            for p,(y,j) in enumerate(rows):P[y][p]+=weights[ticket]/2
            for p in range(3):
                for q in range(3):
                    if rows[p]==rows[q]:gram[p][q]+=weights[ticket]/2
            primitive.append({'ticket':ticket,'edge':list(EDGES[e]),'coin':coin,
                              'outputs':[r[0] for r in rows],'currents':[r[1] for r in rows]})
    require(gram==[[F(int(i==j)) for j in range(3)] for i in range(3)],'observed current instrument is not isometric')
    width=2*n+1;permutation=[]
    for p in range(3):
        for h in range(width):
            if h==2*n:y,hh=p,h
            else:
                e,coin=tickets[h//2],h%2;y,current=indexed[p,e,coin]
                hh=2*(h//2)+(coin^int(current!=0))
            permutation.append(width*y+hh)
    require(sorted(permutation)==list(range(3*width)),'observed current reversal is not a permutation')
    require(all(permutation[permutation[i]]==i for i in range(3*width)),'observed current reversal is not involutive')
    return {'tickets':tickets,'ticket_probabilities':[str(w) for w in weights],
            'inactive_probability':str(inactive),'P':[[str(x) for x in row] for row in P],
            'fine_isometry_gram':[[str(x) for x in row] for row in gram],
            'R_folded_permutation':permutation,'R_dimension':3*width,'primitive_columns':primitive,
            'folding':'Inactive fair-coin identity is restricted to its declared coherent coin vector; orthogonal padding is retained by the full instrument.',
            'scope':'Measured classical endpoint/current columns determine these maps. Coherent control, erasure, inverse and reflection remain declared services.'}


def catalog(columns,blank_records):
    cooling=[0,0,2,2,1]
    result={'cooling5':compile_menu(columns,cooling),
            'weighted_cooling_diagnostic':compile_menu(columns,cooling,[F(3,20)]*5,blank_records),
            'stars':[]}
    for m in (4,8,32):
        for slot in range(3*m):
            v,r=divmod(slot,m);w=(v+1)%3
            tickets=[e for e,edge in enumerate(EDGES) if v in edge]*(m-r)
            tickets += [e for e,edge in enumerate(EDGES) if w in edge]*r
            result['stars'].append({'M':m,'slot':slot,**compile_menu(columns,tickets)})
    return result
