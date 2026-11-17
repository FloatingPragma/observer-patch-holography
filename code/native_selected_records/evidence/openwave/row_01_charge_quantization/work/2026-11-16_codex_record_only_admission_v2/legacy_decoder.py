"""Exact b6df current-only functions, copied for a metadata isolation test."""

def tree(n, edges):
    seen={0}; queue=[0]; result=[]
    for a in queue:
        for u,v in edges:
            b=v if u==a else u if v==a else None
            if b is not None and b not in seen:
                seen.add(b);queue.append(b);result.append((a,b))
    assert len(seen)==n
    return result

def observe(n, edges, body, probe, cap):
    change=[0]*n; relative=[None]*n; relative[0]=0; cursor=0
    for a,b in tree(n,edges):
        correction=change[a]-change[b]; first=None
        for _ in range(cap):
            if cursor==len(probe): return {'outcome':'probe_timeout'}
            aa,bb,c,j=probe[cursor];cursor+=1
            assert (aa,bb)==(a,b) and c in (0,1) and type(j) is int
            change[a]-=j;change[b]+=j
            if first is None: first=(c,j)
            elif c!=first[0]:
                relative[b]=relative[a]-(2*first[1]+j-correction);break
        else: return {'outcome':'probe_timeout'}
    assert cursor==len(probe)
    initial=relative.copy()
    for a,b,c,j in reversed(body): initial[a]+=j;initial[b]-=j
    offset=-initial[12];initial=[x+offset for x in initial]
    live=initial.copy();first=None
    for t,record in enumerate(body+[None]):
        if first is None and (set(live)<={0,1} or set(live)<={-1,0}): first=t
        if record:
            a,b,c,j=record;live[a]-=j;live[b]+=j
    before=live.copy()
    for a,b,c,j in probe:live[a]-=j;live[b]+=j
    return {'outcome':'decoded','initial':initial,'before_probe':before,
            'after_probe':live,'first_unit_event':first,'total':sum(initial)}
