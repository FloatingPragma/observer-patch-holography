"""Station-owned clock exchanges and route delivery in the supplied addressed apparatus."""
from fractions import Fraction

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


def station_times(exchanges, stamps):
    """Recover one station's range and home timestamps from three reciprocal trips."""
    require(exchanges.shape == (3, 3) and exchanges.dtype.kind == "i" and
            stamps.ndim == 1 and stamps.dtype.kind == "i", "incomplete station clock records")
    # Python integers preserve offsets and differences across the signed-counter range.
    sent, remote, returned = zip(*exchanges.tolist())
    elapsed = returned[0] - sent[0]
    require(sent[0] < sent[1] < sent[2] and elapsed >= 0 and elapsed % 2 == 0 and
            all(end - start == elapsed for start, end in zip(sent, returned)), "inconsistent range records")
    home_step, remote_step = sent[1] - sent[0], remote[1] - remote[0]
    require(remote_step > 0 and
            (remote[2] - remote[0]) * home_step == (sent[2] - sent[0]) * remote_step,
            "third clock exchange failed")
    distance = elapsed // 2
    midpoint = sent[0] + distance
    return distance, [midpoint + Fraction((tick - remote[0]) * home_step, remote_step)
                      for tick in stamps.tolist()]


def reconstruct(exchanges, stamps, lap):
    require(exchanges.ndim == 3 and exchanges.shape[1:] == (3, 3) and len(exchanges) > 0 and
            stamps.ndim == 2 and stamps.shape[1] == len(exchanges) and len(stamps) >= 2,
            "incomplete clock menu")
    lap = np.asarray(lap)
    require(lap.shape == () and lap.dtype.kind in "iu" and lap.item() == 2 * len(exchanges),
            "lap differs from twice the station count")
    coords = []
    home = None
    for column, station in enumerate(exchanges):
        distance, times = station_times(station, stamps[:, column])
        coords.append(distance)
        if home is None:
            home = times
        require(times == home, "asynchronous clock records")
    require(sorted(coords) == list(range(len(exchanges))), "missing or repeated station")
    observed = np.asarray(home, dtype=float)
    require(np.all(np.diff(observed) > 0), "unresolved observation times")
    return np.asarray(coords), observed, lap.item() / 2
