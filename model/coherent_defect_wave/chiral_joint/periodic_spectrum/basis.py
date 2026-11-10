"""Exact full Gauss kernels on the actual two-oriented-link carrier."""
from collections import defaultdict
from functools import lru_cache
import sympy as s
from ..two_cell.basis import boundary_charges,current_data,_local_weight_states,_integer
from ..two_cell.symmetries import periodic_current_column

DIMENSION=3686400


@lru_cache(maxsize=None)
def weight_blocks(boundary,left_two_T3=0,right_two_T3=0):
    ql,qr=boundary_charges(boundary);local=_local_weight_states();_,_,edge=current_data()
    blocks=defaultdict(list)
    for l01 in range(6):
        for l10 in range(6):
            ly=_integer(6*(edge[0][11][l01,l01]+edge[1][11][l10,l10]))
            ry=_integer(6*(edge[1][11][l01,l01]+edge[0][11][l10,l10]))
            lz=_integer(2*(edge[0][10][l01,l01]+edge[1][10][l10,l10]))
            rz=_integer(2*(edge[1][10][l01,l01]+edge[0][10][l10,l10]))
            for ml,hl,ll,sl in local.get((-6*ql-ly,left_two_T3-lz),()):
                for mr,hr,lr,sr in local.get((-6*qr-ry,right_two_T3-rz),()):
                    blocks[ll+lr,sl+sr].append((ml|(mr<<6),l01,l10,hl,hr))
    return {label:tuple(sorted(keys)) for label,keys in sorted(blocks.items())}


def block_labels(boundary=0):return tuple(weight_blocks(boundary))


def diagonal_basis(boundary,L,twoSz):return weight_blocks(boundary).get((L,twoSz),())


@lru_cache(maxsize=None)
def character_multiplicities(boundary=0):
    weights={(a,b):weight_blocks(boundary,a,b) for a in (0,2) for b in (0,2)}
    labels=set().union(*(set(v) for v in weights.values()));result={}
    for label in sorted(labels):
        n=len(weights[0,0].get(label,()))-len(weights[2,0].get(label,()))-len(weights[0,2].get(label,()))+len(weights[2,2].get(label,()))
        assert n>=0
        if n:result[label]=n
    return result


def raising_column(key,site):
    a=periodic_current_column(key,site,8);b=periodic_current_column(key,site,9)
    return {target:s.expand(a.get(target,0)+s.I*b.get(target,0))
            for target in set(a)|set(b) if s.expand(a.get(target,0)+s.I*b.get(target,0))}


@lru_cache(maxsize=None)
def raising_matrix(boundary,L,twoSz):
    keys=diagonal_basis(boundary,L,twoSz);columns=[];rows=set()
    for key in keys:
        col={(site,target):value for site in (0,1) for target,value in raising_column(key,site).items()}
        columns.append(col);rows.update(col)
    rows=tuple(sorted(rows));positions={row:i for i,row in enumerate(rows)}
    return s.ImmutableSparseMatrix(len(rows),len(keys),{(positions[row],j):v
        for j,col in enumerate(columns) for row,v in col.items()}),rows


@lru_cache(maxsize=None)
def physical_block(boundary,L,twoSz):
    m,rows=raising_matrix(boundary,L,twoSz);n=m.cols;parent=list(range(n))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    row_columns=defaultdict(list)
    for row,col in m.todok():row_columns[row].append(col)
    for columns in row_columns.values():
        for col in columns[1:]:parent[find(col)]=find(columns[0])
    components=defaultdict(list)
    for col in range(n):components[find(col)].append(col)
    component_rows=defaultdict(list)
    for row,cols in row_columns.items():component_rows[find(cols[0])].append(row)
    entries={};coordinates=[];ranks=[];rank=0
    for columns in sorted(components.values(),key=lambda c:c[0]):
        small=s.Matrix(m.extract(component_rows[find(columns[0])],columns))
        rref,pivots=small.rref();rank+=len(pivots)
        free=[j for j in range(len(columns)) if j not in pivots]
        for freecol in free:
            outcol=len(coordinates);coordinates.append(columns[freecol]);entries[columns[freecol],outcol]=1
            for r,pivot in enumerate(pivots):
                value=-rref[r,freecol]
                if value:entries[columns[pivot],outcol]=value
        ranks.append((len(columns),small.rows,len(pivots)))
    b=s.ImmutableSparseMatrix(n,len(coordinates),entries)
    assert m*b==s.zeros(m.rows,b.cols)
    assert b.extract(coordinates,range(b.cols))==s.eye(b.cols)
    assert b.cols+rank==n and b.cols==character_multiplicities(boundary).get((L,twoSz),0)
    return b,s.ImmutableSparseMatrix(b.H*b),tuple(coordinates),{
        'boundary':boundary_charges(boundary),'L':L,'twoSz':twoSz,
        'diagonal_dimension':n,'physical_dimension':b.cols,'exact_constraint_rank':rank,
        'raising_image_rows':m.rows,'largest_constraint_component':max((r[0] for r in ranks),default=0),
        'constraint_components':[{'domain':a,'image':b,'rank':c} for a,b,c in ranks],
        'independent_character_count_matches':True}


def physical_basis(boundary,L,twoSz):return physical_block(boundary,L,twoSz)[:2]


def coordinate_rows(boundary,L,twoSz):return physical_block(boundary,L,twoSz)[2]


def metadata(boundary=0,execute_kernel=False):
    sizes=character_multiplicities(boundary)
    return {'boundary':boundary_charges(boundary),'full_active_dimension':DIMENSION,
        'diagonal_dimension':sum(map(len,weight_blocks(boundary).values())),
        'physical_dimension':sum(sizes.values()),'largest_physical_block':max(sizes.values(),default=0),
        'blocks':[physical_block(boundary,*label)[3] if execute_kernel else
                  {'L':label[0],'twoSz':label[1],'physical_dimension':n,
                   'diagonal_dimension':len(diagonal_basis(boundary,*label))}
                  for label,n in sizes.items()]}
