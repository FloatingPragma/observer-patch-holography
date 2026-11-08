"""Record-only local masks, with explicit coverage and constant-sector erasures."""
from fractions import Fraction as F
from math import comb
from model.native_program import interpreter as source
from model.observer import native as endpoint
from model.native_sector_readout import readout as inherited
from model.weighted_source_family import family as native
from model.native_rooted_refinement import routing


def integer(x,low=0,high=None):
    if type(x) is not int or x<low or (high is not None and x>high):
        raise ValueError('an exact integer in the declared domain is required')
    return x


def plan(ports,incident):
    """A graph-local theorem interface; source attachment is checked separately."""
    if (type(ports) is not tuple or len(ports)<2 or any(type(x) is not int or x<0 for x in ports)
            or tuple(sorted(set(ports)))!=ports or type(incident) is not tuple):
        raise ValueError('canonical distinct port and incident-seam tuples required')
    keep=set(ports)
    if (not incident or any(
        type(e) is not tuple or len(e)!=2 or any(type(x) is not int or x<0 for x in e)
        or e[0]>=e[1] or not keep.intersection(e) for e in incident)
        or tuple(sorted(set(incident)))!=incident):
        raise ValueError('ordered distinct incident seams required')
    internal=tuple(e for e in incident if set(e)<=keep)
    adjacent={v:[] for v in ports}
    for a,b in internal:adjacent[a].append(b);adjacent[b].append(a)
    seen={ports[0]};queue=[ports[0]];tree=[]
    for a in queue:
        for b in sorted(adjacent[a]):
            if b not in seen:seen.add(b);queue.append(b);tree.append(tuple(sorted((a,b))))
    if seen!=keep:raise ValueError('the retained induced graph must be connected')
    return {'ports':ports,'incident_seams':incident,'tree':tuple(tree)}


def authenticate(p):
    if type(p) is not dict or set(p)!={'ports','incident_seams','tree'}:
        raise ValueError('exact local readout plan fields required')
    expected=plan(p['ports'],p['incident_seams'])
    if type(p['tree']) is not tuple or any(type(e) is not tuple or any(type(i) is not int for i in e) for e in p['tree']):
        raise ValueError('typed tree addresses required')
    if source.canonical(p)!=source.canonical(expected):raise ValueError('readout tree changed')
    return p


def source_plan(level,ports):
    integer(level,0,source.MAX_EVALUATED_LEVEL)
    family=native.federation_family(level,(1,)*12)
    if type(ports) is not tuple or any(type(p) is not int or not 0<=p<len(family.weights) for p in ports):
        raise ValueError('literal source ports required')
    return plan(ports,tuple(e for e in family.seams if set(e).intersection(ports)))


def rooted_plan(window,inner_ports):
    """All incident edges must be certified inside the supplied outer window.

    The caller supplies a union of complete interior cells, never an inferred
    untested routing limit. Every neighbour of each retained cell is checked.
    """
    routing.authenticate_window(window)
    if type(inner_ports) is not tuple or any(type(i) is not int or not 0<=i<len(window['source_ports']) for i in inner_ports):
        raise ValueError('canonical inner-window port indices required')
    names=window['canonical_ports'];selected={names[i][0] for i in inner_ports}
    expected=tuple(i for i,name in enumerate(names) if name[0] in selected)
    if inner_ports!=expected:raise ValueError('complete canonical cells required')
    allcells={name[0] for name in names}
    for cell,neighbours in window['routing']:
        if cell in selected and not {n for n,_ in neighbours}<=allcells:
            raise ValueError('the outer window omits an incident boundary seam')
    keep=set(inner_ports);addresses=window['source_ports']
    incident=tuple((addresses[a],addresses[b]) for a,b in window['internal_seams'] if keep.intersection((a,b)))
    return plan(tuple(sorted(addresses[i] for i in inner_ports)),tuple(sorted(tuple(sorted(e)) for e in incident)))


def record(edge,coin,label,selector='uniform'):
    if type(edge) is not tuple or len(edge)!=2 or any(type(x) is not int or x<0 for x in edge) or edge[0]>=edge[1]:
        raise ValueError('ordered source edge required')
    integer(coin,0,1);integer(label,-1,1)
    if selector not in ('uniform','fixed'):raise ValueError('declared selector tag required')
    rows=source.circuit.record_values('E',coin)
    raw=next(k for k,v in rows.items() if v==label)
    current=(label-(2*coin-1)*(abs(label)%2))//2
    return source.Record('repair',selector,edge,'E',coin,raw,label,current)


def valid_record(r,p,selector):
    if type(r) is not source.Record or r.seam not in p['incident_seams']:
        raise ValueError('an incident original source record required')
    if any(type(x) is not int for x in r.seam) or type(r.coin) is not int or r.coin not in (0,1):
        raise ValueError('typed source labels required')
    if r.law not in ('C','E') or type(r.raw) is not int or not 0<=r.raw<(1<<source.circuit.RECORD_BITS[r.law]):
        raise ValueError('complete native raw slot required')
    label=source.circuit.record_values(r.law,r.coin).get(r.raw)
    current=None if label is None else label if r.law=='C' else (label-(2*r.coin-1)*(abs(label)%2))//2
    wanted=source.Record('repair',selector,r.seam,r.law,r.coin,r.raw,label,current)
    if not inherited._same_record(r,wanted):raise ValueError('native record fields disagree')
    return r.law=='E' and label in (-1,0,1)


def decode(p,events,attempts,*,selector='uniform',population=None,total=None):
    """Only records, addresses and a public total-attempt count are inputs.

    Each event is (global tick, raw record); every incident event is required.
    A caller authenticates that coverage against the actual complete history.
    Returned masks use the compact bit order of the canonical ports tuple.
    """
    authenticate(p);integer(attempts)
    known = population is not None or total is not None
    if known:
        bounds=feasibility(population,total,len(p['ports']))
        if p['ports'][-1]>=population:raise ValueError('port lies outside the declared global carrier')
    if type(events) is not tuple or selector not in ('uniform','fixed'):
        raise ValueError('literal chronological events and selector required')
    previous=0;shift={i:0 for i in p['ports']};differences={}
    for event in events:
        if type(event) is not tuple or len(event)!=2:raise ValueError('tick/record tuple required')
        tick,r=event;integer(tick,previous+1,attempts);previous=tick
        if not valid_record(r,p,selector):
            return {'status':'ZERO_RAW_WORD','initial_mask':None,'terminal_mask':None,'missing_tree_edges':()}
        a,b=r.seam
        if known and bounds['global_constant'] and r.label!=0:
            return {'status':'ZERO_RAW_WORD','initial_mask':None,'terminal_mask':None,'missing_tree_edges':()}
        if r.seam in p['tree'] and r.seam not in differences:
            differences[r.seam]=r.label-shift[a]+shift[b]
        if a in shift:shift[a]-=r.current
        if b in shift:shift[b]+=r.current
    missing=tuple(e for e in p['tree'] if e not in differences)
    if known and bounds['global_constant']:
        mask=0 if total==0 else (1<<len(p['ports']))-1
        return {'status':'DECODED','initial_mask':mask,'terminal_mask':mask,'missing_tree_edges':missing}
    if missing:return {'status':'COVERAGE_ERASURE','initial_mask':None,'terminal_mask':None,'missing_tree_edges':missing}
    offset={p['ports'][0]:0}
    for a,b in p['tree']:
        if a in offset:offset[b]=offset[a]-differences[(a,b)]
        else:offset[a]=offset[b]+differences[(a,b)]
    if not any(offset.values()) and not (known and not bounds['ambiguous_constant']):
        return {'status':'CONSTANT_ERASURE','initial_mask':None,'terminal_mask':None,'missing_tree_edges':()}
    candidates=[]
    for base in (0,1):
        initial={i:base+offset[i] for i in p['ports']};now=dict(initial)
        if any(v not in (0,1) for v in initial.values()):continue
        if known and not max(0,total-(population-len(initial)))<=sum(initial.values())<=min(total,len(initial)):
            continue
        for _,r in events:
            a,b=r.seam
            na=now[a] if a in now else now[b]+r.label
            nb=now[b] if b in now else now[a]-r.label
            if na not in (0,1) or nb not in (0,1) or na-nb!=r.label:break
            output,(_,j)=endpoint.native((na,nb),(0,1),r.coin)
            if j!=r.current:break
            if a in now:now[a]=output[0]
            if b in now:now[b]=output[1]
        else:candidates.append((sum(initial[i]<<k for k,i in enumerate(p['ports'])),sum(now[i]<<k for k,i in enumerate(p['ports']))))
    if not candidates:
        return {'status':'ZERO_RAW_WORD','initial_mask':None,'terminal_mask':None,'missing_tree_edges':()}
    if len(candidates)!=1:raise ValueError('records do not support one feasible local input')
    return {'status':'DECODED','initial_mask':candidates[0][0],'terminal_mask':candidates[0][1],'missing_tree_edges':()}


def uniform_program(level,total,attempts):
    integer(level);integer(total);integer(attempts,1,source.MAX_ATTEMPTS)
    return source.program(federation={'level':level,'total':total},biased=False,
        commands=[{'op':'repair','selector':{'kind':'uniform'},'law':'E','repeat':attempts}],
        capacity=1,overflow='export_oldest')


def uniform_ceiling(level,total):
    """Actual one-command attempt and retained-record admission ceiling."""
    seed=uniform_program(level,total,1)
    bits=seed['declared_resources']['record_bits_per_slot_upper']
    ceiling=min(source.MAX_ATTEMPTS,(2**40)//bits-1)
    accepted=uniform_program(level,total,ceiling)
    try:uniform_program(level,total,ceiling+1)
    except ValueError:pass
    else:raise ArithmeticError('the owning admission ceiling changed')
    if accepted['declared_resources']['record_storage_bits_upper']>2**40:
        raise ArithmeticError('retained record budget exceeded')
    return ceiling


def fixed_program(level,total,p):
    authenticate(p);integer(level);integer(total)
    if len(p['tree'])>source.MAX_NODES:raise ValueError('literal tree exceeds source node cap')
    return source.program(federation={'level':level,'total':total},biased=False,
        commands=[{'op':'repair','selector':{'kind':'fixed','seam':list(e)},'law':'E','repeat':1} for e in p['tree']],
        capacity=1,overflow='export_oldest')


def observed_events(p,records):
    authenticate(p)
    if type(records) is not tuple or any(type(r) is not source.Record for r in records):
        raise ValueError('complete source record tuple required')
    return tuple((t,r) for t,r in enumerate(records,1) if r.kind=='repair' and r.seam in p['incident_seams'])


def feasibility(population,total,retained):
    integer(population,2);integer(total,0,population);integer(retained,2,min(population,4096))
    empty=total<=population-retained;full=total>=retained
    return {'empty_feasible':empty,'full_feasible':full,'ambiguous_constant':empty and full,
            'global_constant':total in (0,population)}


def canonical_erasure(population,total,retained,*,use_global_total=False):
    bounds=feasibility(population,total,retained)
    if type(use_global_total) is not bool:raise ValueError('an explicit global-specification choice required')
    if use_global_total and not bounds['ambiguous_constant']:
        return {'exact':F(),'product_upper':F()}
    # Falling-factorial ratios avoid a binomial with population-sized indices.
    empty=full=F(1)
    for j in range(retained):
        empty*=F(max(0,population-total-j),population-j)
        full*=F(max(0,total-j),population-j)
    p=F(total,population)
    bound=(1-p)**retained+p**retained
    if empty+full>bound:raise ArithmeticError('canonical erasure upper bound failed')
    return {'exact':empty+full,'product_upper':bound}


def coverage(seams,tree_edges,attempts,*,exact=False):
    integer(seams,1);integer(tree_edges,1,seams);integer(attempts)
    blocks=attempts//seams
    if blocks>4096:raise ValueError('bounded displayed rational coverage certificate required')
    upper=min(F(1),F(tree_edges,2**blocks))
    lower=max(F(),1-F(attempts,seams))
    result={'upper':upper,'lower':lower,'per_seam_attempt_time':F(attempts,seams)}
    if exact:
        if tree_edges>64 or attempts*seams.bit_length()>200000:
            raise ValueError('exact inclusion-exclusion output exceeds its finite arithmetic budget')
        value=1-sum(((-1)**j*comb(tree_edges,j)*F(seams-j,seams)**attempts for j in range(tree_edges+1)),F())
        if not lower<=value<=upper:raise ArithmeticError('coverage bounds failed')
        result['exact']=value
    return result


def error_bound(coverage_error,constant_weight,*,global_constant=False):
    if type(coverage_error) is not F or type(constant_weight) is not F or not 0<=coverage_error<=1 or not 0<=constant_weight<=1:
        raise ValueError('exact probability bounds required')
    if type(global_constant) is not bool:raise ValueError('a typed global-constant flag required')
    return F() if global_constant else coverage_error+(1-coverage_error)*constant_weight


def raw_query(value,p,slots,density,*,use_global_total=False):
    """Execute one unnormalized complete native word, including zero outcomes.

    This preserves the original continuation kernel and FIFO. It is not an
    expanded-history execution, normalized outcome sampler or host-cap lift.
    """
    value=source.admit(value);authenticate(p)
    if type(use_global_total) is not bool:raise ValueError('an explicit global-specification choice required')
    c=value['carrier']
    if c['kind']!='federation' or c['level']>source.MAX_EVALUATED_LEVEL:
        raise ValueError('literal evaluated federation query required')
    if source.canonical(p)!=source.canonical(source_plan(c['level'],p['ports'])):
        raise ValueError('incomplete or foreign incident-edge attachment')
    commands=value['commands']
    if type(slots) is not tuple or any(type(slot) is not int for slot in slots):
        raise ValueError('typed raw slots required')
    uniform=(len(commands)==1 and commands[0].get('selector')=={'kind':'uniform'})
    expected=uniform_program(c['level'],c['total'],len(slots)) if uniform else fixed_program(c['level'],c['total'],p)
    if source.canonical(value)!=source.canonical(expected):raise ValueError('exact declared readout protocol required')
    attempts=value['declared_resources']['attempts_upper']
    if type(slots) is not tuple or len(slots)!=attempts or not 1<=attempts<=64:
        raise ValueError('one bounded raw slot per attempted repair required')
    actual=source._validate_density(density,c);memory=archive=()
    for t,slot in enumerate(slots):
        command=commands[0] if uniform else commands[t]
        outcomes=source._outcomes(value,command,actual,memory,archive)
        integer(slot,0,len(outcomes)-1)
        out=outcomes[slot];actual,memory,archive=out.density,out.memory,out.archive
    records=archive+memory
    if len(records)!=attempts:raise ArithmeticError('source archive omitted an attempt')
    events=observed_events(p,records)
    specification={'population':240*4**c['level'],'total':c['total']} if use_global_total else {}
    decoded=decode(p,events,attempts,selector='uniform' if uniform else 'fixed',**specification)
    return {'program_digest':source.digest(value),'status':'COMPLETE_RAW_WORD_QUERY',
            'uses_declared_global_total':use_global_total,'slots':slots,'density':actual,'archive':archive,'memory':memory,
            'observed_events':events,'decoded':decoded}
