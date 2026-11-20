"""Station-owned clock exchanges and route delivery in the supplied addressed apparatus."""
import numpy as np

from ..kernel import require


def addressed(start, legs, size, payload):
    """The sealed payload travels with the token through successive neighbor services."""
    where = list(start)
    ticks = 0
    for axis,distance in legs:
        sign = 1 if distance>=0 else -1
        for _ in range(abs(distance)):
            where[axis] = (where[axis]+sign)%size
            ticks += 1
    return tuple(where),ticks,payload


def stations(size, axis, times):
    home = (0,0,0)
    order = (5*np.arange(size)+3+axis)%size
    exchanges, stamps = [],[]
    for key,j in enumerate(order.tolist()):
        offset = 17+(3*key%23)
        samples = []
        for sent in (0,11,39):
            at,elapsed,tag = addressed(home,[(axis,j)],size, (key,sent))
            remote = 2*(sent+elapsed)+offset
            end,back,received = addressed(at,[(axis,-j)],size,(tag,remote))
            require(end==home and received[0]==(key,sent),'failed station receipt')
            samples.append((sent,received[1],sent+elapsed+back))
        exchanges.append(samples)
        stamps.append(2*np.asarray(times)+offset)
    _,lap,_ = addressed(home,[(axis,size),(axis,-size)],size,'lap')
    return order,np.asarray(exchanges),np.asarray(stamps).T,lap


def deliver(size,axis,order,local):
    columns = []
    for key,j in enumerate(order.tolist()):
        at = [0,0,0]
        at[axis] = j
        end,_,message = addressed(at,[(axis,-j)],size,(key,local[:,j].copy()))
        require(end==(0,0,0) and key==message[0],'lost station report')
        columns.append(message[1])
    return np.asarray(columns).T


def reconstruct(exchanges,stamps,lap):
    require(exchanges.ndim==3 and exchanges.shape[1:]==(3,3) and len(exchanges)>0 and
            exchanges.dtype.kind=='i' and stamps.dtype.kind=='i' and stamps.ndim==2 and
            stamps.shape[1]==len(exchanges) and len(stamps)>=2,'incomplete clock menu')
    sent,remote,returned = exchanges.transpose(2,0,1)
    require(np.all(sent[:,1:]>sent[:,:-1]) and np.all(returned>=sent),'invalid clock exchange order')
    ranges = (returned-sent)/2
    require(np.all(ranges==ranges[:,:1]) and lap>0,'inconsistent range records')
    coords = ranges[:,0]
    require(np.array_equal(np.sort(coords),np.arange(lap/2)),'missing or repeated station')
    midpoint = sent+ranges
    scale = (remote[:,1]-remote[:,0])/(midpoint[:,1]-midpoint[:,0])
    require(np.all(np.isfinite(scale)) and np.all(scale>0),'invalid station clock rate')
    offset = remote[:,0]-scale*midpoint[:,0]
    require(np.array_equal(remote[:,2],scale*midpoint[:,2]+offset),'third clock exchange failed')
    home = (stamps-offset)/scale
    require(np.all(home==home[:,:1]) and np.all(np.diff(home[:,0])>0),'asynchronous clock records')
    return coords,home[:,0],lap/2
