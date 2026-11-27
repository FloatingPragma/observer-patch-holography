import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, expm, schur

from .apparatus import Guide, aperture
from .emitter import absorber, launch, selected_parent, three_exits

HERE = Path(__file__).resolve().parent
IDENTITY = np.eye(2)
X = np.array([[0,1],[1,0]])
Y = np.array([[0,-1j],[1j,0]])
Z = np.diag([1,-1])


def tensor(values):
    out = np.array([[1.]])
    for item in values:
        out = np.kron(out, item)
    return out


def spins(n, fraction):
    xs = [tensor([IDENTITY]*i+[X]+[IDENTITY]*(n-i-1)) for i in range(n)]
    zs = [tensor([IDENTITY]*i+[Z]+[IDENTITY]*(n-i-1)) for i in range(n)]
    gamma = np.array([tensor([X]*i+[v]+[IDENTITY]*(n-i-1)) for i in range(n) for v in (Z,Y)])
    link = sum((fraction if i == 0 else 1)*zs[i] @ zs[i+1] for i in range(n-1))
    u = expm(-1j*link/16) @ expm(1j*sum(xs)/16)
    return xs, zs, gamma, u


def vacuum(g, gamma):
    result = np.eye(len(gamma[0]), dtype=complex)
    for mode in g.modes.T:
        destroy = np.einsum('i,ijk->jk', mode.conj()/2, gamma)
        result = result @ (np.eye(len(result))-destroy.conj().T @ destroy)
    values, vectors = eigh((result+result.conj().T)/2)
    np.testing.assert_allclose(values[:-1], 0, atol=5e-12)
    np.testing.assert_allclose(values[-1], 1, atol=5e-12)
    return vectors[:, -1]


def test_direct_material_spins_whole_aperture_and_both_parent_outcomes():
    n, tick = 4, 7
    g = Guide(n,'native',.5)
    cold = Guide(n,'native',0)
    xs, zs, gamma, u = spins(n,.5)
    # Every generator row is checked against the true full-spin two-layer product.
    for i, gm in enumerate(gamma):
        np.testing.assert_allclose(u.conj().T @ gm @ u, np.einsum('j,jab->ab',g.r[i],gamma), atol=2e-13)
    vac = vacuum(g,gamma)
    old = vacuum(cold,gamma)
    tail = zs[-1]
    tri, v = schur(u, output='complex')
    angles = np.angle(tri.diagonal())
    target = g.phi[0]/g.period
    signed = v.conj().T @ tail @ v
    tap = v @ (signed*aperture(target*g.period+angles[None,:]-angles[:,None],2048)) @ v.conj().T
    r = np.linalg.matrix_power(g.r,tick)
    for kind, source in [('parent',zs[0] @ old),('packet',zs[1] @ old),('sham',old)]:
        state = np.linalg.matrix_power(u,tick) @ source
        c = r @ launch(cold,kind) @ r.T
        for sign in (1,-1):
            branch = (state + sign*xs[0] @ state)/2
            prob = np.vdot(branch,branch).real
            observed, collapsed = selected_parent(c,sign)
            np.testing.assert_allclose(observed,prob,atol=3e-12)
            tagged = tap @ branch
            yes = abs(np.vdot(vac,tagged))**2
            success = np.vdot(tagged,tagged).real
            actual = observed*three_exits(collapsed,absorber(g,target),kind=='sham')
            np.testing.assert_allclose(actual,[yes,success-yes,prob-success],atol=3e-11)


def test_two_stage_finite_residual_filter_and_start_delay():
    n, tick = 4, 9
    g = Guide(n,'native')
    _, zs, gamma, u = spins(n,1.)
    vac = vacuum(g,gamma)
    tri, v = schur(u,output='complex')
    angles = np.angle(tri.diagonal())
    ref = np.angle(np.vdot(vac,u @ vac))
    theta = g.phi[0]
    tap = v @ ((v.conj().T @ zs[-1] @ v)*aperture(theta+angles[None,:]-angles[:,None],2048)) @ v.conj().T
    filt = v @ np.diag(aperture(angles-ref,4096)) @ v.conj().T
    for sham, state in [(False,zs[0] @ vac),(True,vac)]:
        post = tap @ (np.linalg.matrix_power(u,tick) @ state)
        bright = np.linalg.norm(filt @ post)**2
        passed = np.vdot(post,post).real
        expected, exact = g.absorb(theta/g.period,tick,sham)
        np.testing.assert_allclose(expected,[bright,passed-bright,1-passed],atol=3e-11)
        np.testing.assert_allclose(exact,abs(np.vdot(vac,post))**2,atol=3e-11)


def test_joint_launch_denominators_bright_bound_and_charge_coherence():
    data = json.loads((HERE/'emitter.json').read_text())
    for panel in data['panels']:
        for bank in panel['sources']:
            assert sum(bank['joint_counts']) == sum(bank['uninspected_counts']) == data['shots']
        parent, sham, packet = panel['sources']
        if panel['source_link_fraction'] == .5:
            assert parent['joint_lower_99'] > panel['bright_only_tag_upper']
            assert parent['uninspected_counts'][0]/data['shots'] > panel['bright_only_tag_upper']
            assert sham['joint_counts'][0] == 0
        if panel['source_link_fraction'] == 0:
            assert parent['joint_counts'][0] == 0 < packet['joint_counts'][0]
    charge = json.loads((HERE/'charge.json').read_text())['settings']
    for setting in charge.values():
        assert max(setting['coherent_test_errors']) < 3e-12
    np.testing.assert_allclose(np.array(charge['parity']['head_20_probabilities'])[[1,2,5,6]],[0,1,0,1],atol=1e-12)
    np.testing.assert_allclose(np.array(charge['sign']['head_20_probabilities'])[[1,2,5,6]],[0,.5,1,.5],atol=1e-12)
