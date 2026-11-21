"""A 26-service loading schedule for the signed four-corner circulation.

Age zero is the compulsory bright birth update. The age-one local gates
cancel the six unwanted arms of that first hop. On the four surviving
corners each later service uses both horizontal and vertical edges. All
gates act before a complete native service or at the terminal appointment.
Each rail has its own inverse contact coin. The source phases and the local
gates are supplied inputs. No two-site control or rail measurement is used.

In the symmetric four-site sector let a and c be the nonnegative root
and far-corner amplitudes, and b the norm at the two middle sites together.
An on-site plate divides each middle site's amplitude between the edge to
the root and the edge to the far corner using cos(beta) and sin(beta).
The two independent angles below specify the service-dependent splitting.
They solve the two terminal conditions a=c=1/2. Their rational enclosure
was chosen by an exploratory root search on the native pulse angles.
The interval calculation proves existence and uniqueness within that box.
The executable gate matrices use the rounded center of the enclosure.
"""
from __future__ import annotations

import numpy as np

from . import clock_loader as base
from .. import law

END = 26
CENTER = ("0.1878769362445742693320310536844898", "1.1136838729257317506159348841672850")
RADIUS = "1e-30"


def _setup(cos0, sin0, cos1, sin1, sqrt):
    """Refocus six root-neighbor edges and keep the other two."""
    d = sin0 / sqrt(8)
    m = d*cos1/sin1
    g = sqrt((cos0*cos0 - 6*m*m)/2)
    b = sqrt(2)*(sin1*g + cos1*d)
    return d, m, g, b


def _endpoint(ctx, angles):
    """Terminal amplitudes and their Jacobian, using outward interval operations."""
    cos, sin, sqrt = ctx.cos, ctx.sin, ctx.sqrt
    one = ctx.mpf(1)
    theta = [ctx.pi*(1-cos(ctx.pi*j/128))/4 for j in range(1, END+1)]
    _, _, _, b = _setup(cos(theta[0]), sin(theta[0]), cos(theta[1]), sin(theta[1]), sqrt)
    a, c = sqrt(1-b*b), ctx.mpf(0)
    derivatives = [(0, 0, 0), (0, 0, 0)]
    for j in range(2, END):
        u = (one*(j-2)/(END-3))**2
        beta = (1-u)*angles[0] + u*angles[1]
        cb, sb = cos(beta), sin(beta)
        co, si = cos(theta[j]), sin(theta[j])
        x, y = si*a + co*b*cb, co*b*sb - si*c
        next_b = sqrt(x*x+y*y)
        out = []
        for (da, db, dc), dbeta in zip(derivatives, (1-u, u)):
            du = db*cb - b*sb*dbeta
            dv = db*sb + b*cb*dbeta
            dx, dy = si*da + co*du, co*dv - si*dc
            out.append((co*da - si*du, (x*dx+y*dy)/next_b, co*dc + si*dv))
        derivatives = out
        a, b, c = co*a-si*b*cb, next_b, co*c+si*b*sb
        if not (a > 0 and b > 0 and c > 0):
            raise ArithmeticError("shell amplitudes leave the positive orientation chart")
    return (a-one/2, c-one/2), [[d[i] for d in derivatives] for i in (0, 2)]


def certify_angles():
    """An outward interval contraction proof using a rational preconditioner.

    For F=(a-1/2,c-1/2) and the enclosed box X, T(x)=x-MF(x) is
    contracting and maps X into its interior. M is nonsingular, so the
    unique fixed point in X solves both terminal conditions exactly.
    """
    from mpmath.ctx_iv import MPIntervalContext
    ctx = MPIntervalContext()
    ctx.dps = 65
    radius = ctx.mpf(RADIUS)
    center = [ctx.mpf(s) for s in CENTER]
    box = [v + radius*ctx.mpf([-1, 1]) for v in center]
    f, _ = _endpoint(ctx, center)
    _, jac = _endpoint(ctx, box)
    matrix = [[ctx.mpf(s) for s in row] for row in
              (("5.262717", "-17.769744"), ("-.409285", "8.358533"))]
    det = matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0]
    if not det > 0:
        raise ArithmeticError("singular preconditioner")
    norms, shifts = [], []
    for i in range(2):
        shifts.append(abs(sum(matrix[i][k]*f[k] for k in range(2))))
        norms.append(sum(abs(int(i == j)-sum(matrix[i][k]*jac[k][j] for k in range(2)))
                         for j in range(2)))
    if not all(v < ctx.mpf(".00001") for v in norms):
        raise ArithmeticError("interval map is not contracting at the declared bound")
    if not all(delta + lipschitz*radius < radius for delta, lipschitz in zip(shifts, norms)):
        raise ArithmeticError("interval map fails to preserve the angle box")
    return {"center": list(CENTER), "radius": RADIUS, "contraction_upper": "0.00001",
            "fixed_point_interior": True}


def _ray_gate(start, desired):
    """Map one nonzero complex ray to another with its declared phase."""
    u = start/np.linalg.norm(start)
    v = desired/np.linalg.norm(desired)
    overlap = np.vdot(v, u)
    phase = overlap/abs(overlap) if overlap != 0 else 1.
    difference = u-phase*v
    length = np.linalg.norm(difference)
    if length < 2e-15:
        return phase.conjugate()*base.IDENTITY
    w = difference/length
    return phase.conjugate()*(base.IDENTITY-2*np.outer(w,w.conj()))


def schedule(carrier, root):
    """Compile the local gate table, preserving the native first-service rail phase."""
    carrier, root = base.chart(carrier, root)
    sites = base.corners(carrier, root)
    slots = base.IDENTITY
    phi = base.SOURCE.rail_contact_phases
    undo = np.stack([law.contact_coin(slots, -p).T for p in phi])
    theta = np.array([law.hop_angle(j, base.SOURCE.maturation) for j in range(END)])
    co, si = np.cos(theta), np.sin(theta)
    d, m, g, _ = _setup(co[0], si[0], co[1], si[1], np.sqrt)
    root_input = np.full(8, -m, complex)
    root_input[[0,2]] = g
    layer = {root: undo @ _ray_gate(law.EQUAL, root_input)}
    for slot in range(8):
        layer[base.neighbour(root, slot, carrier)] = undo
    gates = {1: layer}

    # The first two services leave zero amplitude at all six unwanted neighbors.
    middle = 1j*(si[1]*g+co[1]*d)
    incoming = [co[1]*root_input - si[1]*d, middle*slots[1], middle*slots[3],
                np.zeros(8, complex)]
    root_ray = (slots[0]+slots[2])/np.sqrt(2)
    far_ray = -(slots[1]+slots[3])/np.sqrt(2)
    left_angle, right_angle = map(float, CENTER)
    for j in range(2, END):
        a = np.linalg.norm(incoming[0])
        b = np.sqrt(2)*np.linalg.norm(incoming[1])
        c = np.linalg.norm(incoming[3])
        u = ((j-2)/(END-3))**2
        beta = left_angle*(1-u) + right_angle*u
        cb, sb = np.cos(beta), np.sin(beta)
        toward_root, toward_far = 1j*b*cb/np.sqrt(2), 1j*b*sb/np.sqrt(2)
        desired = [a*root_ray, toward_root*slots[1]+toward_far*slots[2],
                   toward_root*slots[3]+toward_far*slots[0], c*far_ray]
        gates[j] = {site: undo @ _ray_gate(before, after)
                    for site, before, after in zip(sites, incoming, desired) if np.any(before)}
        x = 1j*(si[j]*a+co[j]*b*cb)/np.sqrt(2)
        y = 1j*(co[j]*b*sb-si[j]*c)/np.sqrt(2)
        incoming = [(co[j]*a-si[j]*b*cb)*root_ray, x*slots[1]+y*slots[2],
                    x*slots[3]+y*slots[0], (co[j]*c+si[j]*b*sb)*far_ray]
    target = base.circulation(carrier, root)
    gates[END] = {}
    for site, before in zip(sites, incoming):
        terminal = _ray_gate(before, target[site])
        gates[END][site] = np.stack([terminal, terminal])
    return gates
