"""Install the historical arrow order and acquire a signed spatial difference."""

from collections import Counter
from fractions import Fraction as Q
from math import expm1, log, log1p

from .local import require


def edge(a,b):
    return tuple(sorted((a,b)))


def triangles(board):
    root = 0
    return sorted((m,f) for m in board.neighbours[root] for f in board.neighbours[root]
                  if f in board.neighbours[m] and m//12 == f//12 == 0)


def install(board, triangle, rank=None):
    """Actual addresses on the accepted branch. Blank swaps and the head read are sharp."""
    rank = range(board.n) if rank is None else rank
    r = 0
    m,f = triangle
    require(triangle in triangles(board), "Missing native parking triangle")
    rm,rf,mf = edge(r,m), edge(r,f), edge(m,f)
    marks = {}
    port = r
    events = []

    def write(e, head=None):
        require(e not in marks and port not in e, "Collar converter needs two vacant endpoints")
        marks[e] = head if head is not None else min(e,key=rank.__getitem__)
        events.append(("write",e,marks[e]))

    def move(e, marked):
        nonlocal port
        require(port in e and ((e in marks) == marked), "Wrong transfer contact")
        to = next(p for p in e if p != port)
        if marked:
            require(marks[e] == port, "Sharp move forbidden by the written head")
            marks[e] = to
        events.append(("head_hop" if marked else "blank_swap",e,to))
        port = to

    for e in board.seams:
        if r not in e and e != mf:
            write(e)
    move(rm,False)
    for e in board.seams:
        if r in e and e != rm:
            if e == rf and min(e,key=rank.__getitem__) != r:
                raise ValueError("Parking schedule requires the donor root earlier than the far port")
            write(e, f if e == rf else None)
    move(mf,False)
    write(rm)
    move(rf,True)
    write(mf)
    require(port == r and set(marks) == set(board.seams), "Incomplete material")
    require(all(h == min(e,key=rank.__getitem__) for e,h in marks.items()), "Installed a different donor")
    return marks, events


def query(v, kind, yes):
    """Unnormalized selective flag on (parent, left, right); compatible with exact rationals."""
    if kind=="occupied":
        hit = [v[0], 0, 0]
    else:
        tip = {"left":1, "right":2}[kind]
        a = (v[0]+v[tip])/2
        hit = [a,0,0]
        hit[tip] = a
    return hit if yes else [x-y for x,y in zip(v,hit,strict=True)]


def attempt(sign):
    require(sign in (-1,1), "Use one of the two real analyzer phases")
    v = [Q(1),Q(0),Q(0)]
    first_wrong = []
    for kind,want in (("left",sign==-1),("right",False),("occupied",False),
                      ("left",False),("occupied",False)):
        wrong = query(v,kind,not want)
        first_wrong.append([q*q for q in wrong])
        v = query(v,kind,want)
    require(v[0]==0 and v[1]!=0 and v[2]==sign*v[1], "Success has the wrong signed ray")
    prob = sum(q*q for q in v)
    fail_at = [sum(t) for t in zip(*first_wrong,strict=True)]
    require(prob == Q(1,8) and prob+sum(fail_at)==1, "Incomplete source attempt")
    return prob,fail_at,first_wrong


def branch_law(edges, depth, n, sign, p):
    """All failure categories on the original trial, including the final exhausted attempt."""
    if type(n) is not int or n<1 or not 0<=p<=1 or not 1 <= depth:
        raise ValueError("Invalid finite source")
    # With 0 < p <= 1 the failure after n fuel trials is (1-p)^n.
    wf = (1-p)**n
    hf = 2.0**(-n)
    s,fail_at,_ = attempt(sign)
    leaf = float(sum(fail_at[1:]))
    retry = 1-float(s)-leaf*hf
    def geom(k):
        return -expm1(k*log(retry))/(1-retry) if k else 0.
    left = 1.
    failure = {}
    def spend(key, q):
        nonlocal left
        failure[key] = left*q
        left *= 1-q
    # Before the commissioning hop only the final opposite edge remains blank.
    spend("writer_before_return", -expm1((edges-1)*log1p(-wf)) if wf<1 else 1.)
    spend("commissioning_hop", hf)
    spend("last_write", wf)
    spend("route_to_parent",-expm1(depth*log1p(-hf)))
    failure["failed_reset"] = left*leaf*hf*geom(n-1)
    failure["no_accepted_word"] = left*(1-float(s))*retry**(n-1)
    accepted = left*float(s)*geom(n)
    require(abs(accepted+sum(failure.values())-1)<5e-15, "Source lost a complementary outcome")
    return accepted,failure


def exact_floor(edges,depth,n):
    # At most n-1 returns after selective failures, and 1 installation hop.
    # This union bound holds on all measurement histories, without rounding.
    return max(Q(0), 1-Q(7,8)**n-Q(edges+depth+n,2**n))


def bare_floor(edges,depth,n):
    return max(Q(0),1-Q(7,8)**n-edges*Q(3,4)**n-Q(depth+n,2**n))
