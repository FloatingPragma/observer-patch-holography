"""ROW_10_FABLE_B, second script: the repository's neutral family coefficient rebuilt from its stated architecture
with the reviewer's port vectors, its levels, and its invariance under the carrier's rotations.

Stated architecture (the candidate's PROOF.md section 6, quoting model/openwave_neutrino_oscillation): Y = sum of
port projectors on the orbit [1, 5, 11, 10, 7] of a fifth-turn with weights 20736, 1728, 144, 12, 1 over 22621;
M_nu = V* Y V^+ with V = I + (-i - 1) P_2.  The reviewer's embedding differs from the repository's chart by an
orthogonal map at most, which changes neither the levels nor the invariance group.  audit.py is the reviewer's code.

    python audit_neutral_coefficient.py     writes RESULT_NEUTRAL_COEFFICIENT.json; stdout is RUN_NEUTRAL_COEFFICIENT.log
"""
import json, math, os, sys
sys.dont_write_bytecode = True
import numpy as np
import audit as A

HERE = os.path.dirname(os.path.abspath(__file__))
S5 = math.sqrt(5.0)
ORBIT, W = [1, 5, 11, 10, 7], [20736, 1728, 144, 12, 1]
adjacent_in_order = all(ORBIT[(i + 1) % 5] in A.ADJ[ORBIT[i]] for i in range(5)) and all(p in A.ADJ[0] for p in ORBIT)


def conj_vec(v):
    out = []
    for x in v:
        p, q = A.ident(x / 2)                      # x/2 = (p + q sqrt5)/4
        out.append((p - q * S5) / 2)
    return np.array(out)


def build(vecs):
    Y = sum(w * np.outer(vecs[p], vecs[p]) / (vecs[p] @ vecs[p]) for p, w in zip(ORBIT, W)) / sum(W)
    P2 = np.outer(vecs[2], vecs[2]) / (vecs[2] @ vecs[2])
    V = np.eye(3) + (-1j - 1) * P2
    return Y, V.conj() @ Y @ V.conj().T, P2


res = {'orbit_is_the_ring_of_neighbours_of_port_0_in_cyclic_order': bool(adjacent_in_order), 'sum_of_weights': sum(W)}
for name, vecs, reps in (('3', A.VEC, A.R3), ("3'", np.array([conj_vec(v) for v in A.VEC]), A.R3P)):
    Y, M, P2 = build(vecs)
    lv = np.linalg.eigvalsh(Y)
    r = (lv[1] ** 2 - lv[0] ** 2) / (lv[2] ** 2 - lv[0] ** 2)
    comY = sorted(float(np.linalg.norm(R @ Y - Y @ R)) for i, R in enumerate(reps) if i != A.E_ID)
    invM = sorted(float(np.linalg.norm(R @ M @ R.T - M)) for i, R in enumerate(reps) if i != A.E_ID)
    port_action = max(float(np.abs(reps[i] @ vecs[k] - vecs[A.G[i][k]]).max()) for i in range(60) for k in range(12))
    res[name] = {'levels_of_Y': lv.tolist(), 'trace': float(np.trace(Y)), 'ratio_r': float(r), 'ratio_over_reference_0.029805': float(r / 0.029805),
                 'M_nu_is_symmetric': float(np.abs(M - M.T).max()), 'smallest_commutator_of_a_rotation_with_Y_frobenius': comY[0],
                 'smallest_change_of_M_nu_under_a_rotation_frobenius': invM[0], 'rotations_act_on_these_port_vectors': port_action,
                 'port_1_and_port_2_are_antipodes': float(np.abs(vecs[1] + vecs[2]).max())}
res['candidate_quotes'] = {'levels_on_3': [0.0020642, 0.0636598, 0.9342760], 'levels_on_the_conjugate': [0.0047330, 0.0608109, 0.9344561], 'r': [0.004638, 0.004209],
                           'smallest_commutator_with_Y': 0.0044, 'smallest_change_of_M_nu': 0.0035}
res['levels_agree_with_the_quotes'] = bool(max(abs(a - b) for a, b in zip(res['3']['levels_of_Y'], res['candidate_quotes']['levels_on_3'])) < 5e-7
                                           and max(abs(a - b) for a, b in zip(res["3'"]['levels_of_Y'], res['candidate_quotes']['levels_on_the_conjugate'])) < 5e-7)
json.dump(res, open(os.path.join(HERE, 'RESULT_NEUTRAL_COEFFICIENT.json'), 'w'), indent=1)
print(json.dumps(res, indent=1))
