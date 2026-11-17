"""Continuing-source interface; input preparation belongs only to the harness."""
from admission import admit
from controller import calibrate, check_blank
from source import source

CELL=((0,1),(0,5),(0,7),(0,10),(0,11),(1,5),(1,7),(1,8),(1,9),
      (2,3),(2,4),(2,6),(2,10),(2,11),(3,4),(3,6),(3,8),(3,9),
      (4,5),(4,9),(4,11),(5,9),(5,11),(6,7),(6,8),(6,10),(7,8),(7,10),(8,9),(10,11))
TRI=((0,1),(1,2),(0,2))


def graph(bound,m):
    if bound<2 or m<1:raise ValueError('public graph parameters')
    n=12*bound
    core=[(12*c+a,12*c+b) for c in range(bound) for a,b in CELL]
    core += [(12*c+v,12*(c+1)+v) for c in range(bound-1) for v in range(12)]
    target=[0,n,n+1];copy=[n+2,n+3,n+4];blank=[n+5,n+6]
    rings=[list(range(n+7,n+7+3*m)),list(range(n+7+3*m,n+7+6*m))]
    edges=list(core)
    for vs in (target,copy):edges += [(vs[a],vs[b]) for a,b in TRI]
    edges += list(zip(target,copy))+[tuple(blank)]
    for ring in rings:edges += list(zip(ring,ring[1:]+ring[:1]))
    bridges=[(0,blank[0]),(0,rings[0][0]),(0,rings[1][0])];edges+=bridges
    return {'N':n+7+6*m,'core_N':n,'core':core,'edges':edges,'target':target,
            'copy':copy,'blank':blank,'rings':rings,'idle_bridges':bridges}


def prepare(first_cell,g):
    """One declared hidden initial input, before the source closure is created."""
    q=list(first_cell)+[0]*(g['N']-12)
    for ring in g['rings']:q[ring[2*(len(ring)//3)]]=1
    return q


def start(initial,g,admission_seed,calibration_seed,probe_law='native',cap=64,forced=None):
    """Return live request closure, audit-only inspector and observer transcript.

    The observer receives neither initial state nor either public replay seed.
    No state assignment occurs after the source is created.
    """
    request,audit=source(initial,g['edges'],admission_seed,calibration_seed,probe_law,forced)
    record={}
    try:
        entry=admit(request,g,cap)
        record['admission']=entry
        if entry['outcome']!='signed_target_unit':
            record['outcome']='not_admitted';return request,audit,record
        calibration=calibrate(request,g['target'],entry['sign'],1,cap)
        record['calibration']=calibration
        if calibration['outcome']!='calibrated':
            record['outcome']=calibration['outcome'];return request,audit,record
        blank=check_blank(request,g['blank'],cap);record['blank']=blank
        record['outcome']='ready' if blank['outcome']=='blank_checked' else blank['outcome']
    except Exception as error:
        record.update(outcome='controller_error',error=repr(error))
    return request,audit,record
