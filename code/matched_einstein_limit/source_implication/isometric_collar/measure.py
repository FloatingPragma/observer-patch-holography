"""Read the head-space effects of a charge-pooled four-station capped radar menu."""

from itertools import permutations
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse.linalg import expm_multiply

from ..geometry_gain import pairing, pooled_ranges, ranging


def basis(fixed_marginals):
    pauli = [np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
    return np.array([np.kron(x, y)/2 for a, x in enumerate(pauli) for b, y in enumerate(pauli)
                     if (a and b) or (not fixed_marginals and (a or b))])


def effects(u, step):
    menus = []
    for horizon in (2, 4, 6):
        stop = round(horizon/step)
        retards = []
        for a, b in permutations(range(4), 2):
            left = ranging.first_hits(u[:stop+1], a, b)
            right = ranging.first_hits(u[:stop+1], b, a)
            weighted = np.zeros((4, 4), complex)
            for n in range(2, stop):
                for j in range(1, n):
                    amp = right[n-j] @ left[j]
                    weighted += ((n*step-horizon)/2)*(amp.conj().T @ amp)
            retards.append(weighted)
        menus.append(retards)
    return np.array(menus)


def radar():
    board = pairing.formation.federation(1)
    route = pairing.formation.region(board)
    joint = pairing.prepared_states(board, route, 1)[0]
    rho = joint/np.trace(joint).real
    banks = []
    for sign in (-1, 1):
        h = ranging.generator(board, route, sign, 8, False)
        u = ranging.kernels(h, route, .25, 24)
        banks.append(np.array([effects(u, .25), effects(u[::2], .5)]))
    return rho, np.mean(banks, axis=0)


def components(effects, operators):
    return np.einsum("...ij,kji->...k", effects, operators).real


def source_columns():
    form = pairing.formation
    board = form.federation(1)
    route = form.region(board)
    h = form.hamiltonian(board, route, 1, 8)
    left, right = form.starting_coins(route, 1)
    src = np.zeros(h.shape[0], complex)
    src[4*route[0]+2*left+right] = 1
    psi = expm_multiply(-4j*h, src, traceA=0).reshape(-1,4)
    psi /= np.linalg.norm(psi)
    return np.einsum("xi,xj->xij", psi, psi.conj())


def operational(rho, all_menus, log):
    """After the mobile position is read, keep the original heads or use an independent formed pair."""
    branches = source_columns()
    np.testing.assert_allclose(branches.sum(axis=0), rho, atol=1e-14)
    weights = np.trace(branches, axis1=1, axis2=2).real
    selected = np.flatnonzero(weights>1e-7)
    columns = (branches - weights[:,None,None]*rho)[selected]
    effect = all_menus[1,2]
    records = components(effect, columns)
    entropy = components(log, columns)
    _, s, v = np.linalg.svd(records, full_matrices=False)
    # Remove the constrained coordinates twice before optimizing the tiny residual score.
    target = entropy.copy()
    for _ in range(2):
        target -= (target @ v.T) @ v
    fit = linprog(-target/max(abs(target)), A_eq=v, b_eq=np.zeros(len(v)),
                  bounds=(-.9,.9), method="highs")
    if not fit.success:
        raise ValueError(fit.message)
    factors = fit.x
    factors -= (factors @ v.T) @ v
    D = np.einsum("x,xij->ij", factors/2, columns)
    shifts = components(effect, D[None])[...,0]
    def S(state):
        w = np.linalg.eigvalsh(state)
        return -w @ np.log(w)
    held = components(all_menus,D[None])[...,0]
    return dict(addresses=selected.tolist(), biases=factors.tolist(), singular=s.tolist(),
                norm=float(np.linalg.norm(D)), entropy_derivative=float(components(log,D[None])[0]),
                entropy_centered=(S(rho+D)-S(rho-D))/2, max_range_shift=float(max(abs(shifts))),
                cap_four_max_range_shift=float(np.max(abs(held[1,1]))),
                marginal_norms=[float(np.linalg.norm(x)) for x in pairing.formation.marginals(D)],
                min_output_eigenvalue=float(min(np.linalg.eigvalsh(rho+D)[0],np.linalg.eigvalsh(rho-D)[0])),
                tangent_real=D.real.tolist(),tangent_imag=D.imag.tolist())


def run():
    rho, banks = radar()
    ranging.ordinary_write()
    archive = json.loads(Path(ranging.__file__).with_name("results.json").read_text())
    native = [arm for arm in archive["runs"] if arm["contrast"]==8 and not arm["left_odd"]]
    states = {"formed":rho, "independent":np.kron(*pairing.formation.marginals(rho))}
    for i in range(2):
        for j,t in enumerate((2,4,6)):
            for name, state in states.items():
                archived = pooled_ranges([arm["panels"][i] for arm in native],name,t)
                clock = np.array([archived[a,b] for a,b in permutations(range(4),2)])
                np.testing.assert_allclose(t/2+components(banks[i,j],state[None])[:,0],clock,
                                           rtol=0,atol=8e-14)
    log = -pairing.log_density(rho)
    output = []
    for fixed in (False, True):
        ops = basis(fixed)
        entropy = components(log, ops)
        for name, menu in [((tick, horizon), banks[i,j]) for i,tick in enumerate((.25, .5))
                           for j,horizon in enumerate((2, 4, 6))] + [("all", banks.reshape(-1,4,4))]:
            rows = components(menu, ops)
            scales = np.linalg.norm(rows, axis=1)
            normalized = rows/scales[:, None]
            _, s, vt = np.linalg.svd(normalized, full_matrices=True)
            rank = sum(s>1e-11)
            kernel = vt[rank:]
            tangent = kernel.T @ (kernel @ entropy)
            output.append(dict(menu=name, fixed=bool(fixed), rank=int(rank), singular=s.tolist(),
                               entropy_hidden=float(np.linalg.norm(tangent)), minimum_row=float(min(scales)),
                               sv_raw=np.linalg.svd(rows,compute_uv=False).tolist()))
    practical = operational(rho,banks,log)
    D = np.array(practical["tangent_real"])+1j*np.array(practical["tangent_imag"])
    c, odd = -.004650970878318474, .25
    s_even, s_odd = .5+c/(1-odd), .5-c/odd
    tau = states["independent"]+s_even*(rho-states["independent"])
    ops = basis(False)
    rows = components(banks[1,2],ops)
    kernel = np.linalg.svd(rows,full_matrices=True)[2][12:]
    direction = kernel.T @ (kernel @ components(log,ops))
    direction /= np.linalg.norm(direction)
    held_log = -pairing.log_density(tau)
    held = dict(odd=odd,c=c,s_even=s_even,s_odd=s_odd,even_formed_fraction=(1-odd)*s_even,
                never_dephased_formed_fraction=(1-odd)*s_even-odd*s_odd,
                same_direction_even_slope=float(components(held_log,ops)@direction),
                practical_even_slope=float(np.trace(held_log@D).real),
                min_even_eigenvalue=float(np.linalg.eigvalsh(tau)[0]))
    return dict(spectrum=np.linalg.eigvalsh(rho).tolist(), menus=output,
                operational_half_tick_six=practical, held_selected_source=held)


if __name__ == "__main__":
    print(json.dumps(run()))
