"""Continuing-unit calibration: only event records and public coordinates enter."""
EDGES=((0,1),(1,2),(0,2))


def require(ok, message):
    if not ok:raise ValueError(message)


def calibrate(request, target, sign, position, cap=64):
    require(len(target)==3 and len(set(target))==3,'public triangle')
    require(sign in (-1,1) and type(position) is int and 0<=position<3,'observed admission')
    require(type(cap) is int and cap>=2,'coupon cap')
    target=tuple(target);entry=position;events=[];columns=[];routing_phases=0

    def take(edge,stage,group):
        nonlocal position
        before=position;record=request(tuple(target[i] for i in edge),stage)
        require(isinstance(record,(list,tuple)) and len(record)==4,'event record shape')
        a,b,coin,current=record
        require((a,b)==tuple(target[i] for i in edge) and type(coin) is int and coin in (0,1),'address/coin record')
        require(type(current) is int,'integer current record')
        local=[sign*int(i==position) for i in range(3)]
        local[edge[0]]-=current;local[edge[1]]+=current
        require(local.count(sign)==1 and local.count(0)==2,'record contradicts admitted isolated unit')
        position=local.index(sign)
        events.append({'stage':stage,'group':list(group),'before_position':before,
                       'record':list(record),'after_position':position})
        return coin,current

    def route(destination,group):
        nonlocal routing_phases
        routing_phases+=1
        if position==destination:return True
        edge=tuple(sorted((position,destination)))
        for _ in range(cap):
            take(edge,'calibration_route',group)
            if position==destination:return True
        return False

    def finish(outcome,group=None):
        return {'outcome':outcome,'target':list(target),'sign':sign,'entry_position':entry,
                'return_position':position,'cap':cap,'columns':columns,'events':events,
                'routing_phases':routing_phases,'endpoint_calls':len(events),'unfinished_group':group}

    for input_position in range(3):
        for e,edge in enumerate(EDGES):
            seen={};group=(input_position,e)
            for attempt in range(cap):
                if not route(input_position,group):return finish('route_timeout',list(group))
                coin,current=take(edge,'calibration_probe',group)
                row={'input_position':input_position,'edge':list(edge),'coin':coin,
                     'current':current,'decoded_output':position}
                if coin in seen:
                    if seen[coin]!=row:return finish('inconsistent_column',list(group))
                else:
                    seen[coin]=row;columns.append(row)
                if len(seen)==2:break
            else:return finish('coupon_timeout',list(group))
    if not route(1,('return',)):return finish('return_timeout',['return'])
    columns.sort(key=lambda r:(r['input_position'],EDGES.index(tuple(r['edge'])),r['coin']))
    return finish('calibrated')


def check_blank(request, blank, cap=64):
    """Check actual zero currents; absolute vacuum remains the admission premise."""
    require(len(blank)==2 and blank[0]!=blank[1] and cap>=2,'public blank/cap')
    records=[];coins=set()
    for _ in range(cap):
        record=request(tuple(blank),'calibration_blank');records.append(list(record))
        require(len(record)==4 and tuple(record[:2])==tuple(blank),'blank record address')
        require(type(record[2]) is int and record[2] in (0,1) and type(record[3]) is int,'blank coin/current')
        if record[3]!=0:return {'outcome':'reference_rejected','records':records,'endpoint_calls':len(records)}
        coins.add(record[2])
        if len(coins)==2:return {'outcome':'blank_checked','records':records,'endpoint_calls':len(records)}
    return {'outcome':'blank_timeout','records':records,'endpoint_calls':len(records)}
