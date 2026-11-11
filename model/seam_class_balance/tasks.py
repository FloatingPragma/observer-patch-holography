"""task lists for staged.py and staged_face_step.py. Run as a module, it writes tasks_gates.json, tasks_prod.json and
tasks_fm.json into the working directory. Reads no measured value."""
import json
import math

X_U = 24.3166
BETA_U = X_U / math.pi                    # Wilson coupling at the chain value, beta = 4/g^2 = x/pi
BETAS = [2.5, 3.0, 4.0, 5.0, 6.0, BETA_U, 9.0]


def par_of(kind, beta):
    return beta if kind == 0 else 4.0 / beta


def tag(kind):
    return {0: "w", 1: "hk", 2: "ab"}[kind]


seed = [1000]


def task(name, cx, kind, beta, N, ncyc, nloc=20, ntherm=40, untw=400, **kw):
    seed[0] += 1
    t = {"name": name, "cx": cx, "kind": kind, "par": par_of(kind, beta), "beta_equiv": beta, "N": N, "ncyc": ncyc,
         "nloc": nloc, "ntherm": ntherm, "nhit": 4, "seed": seed[0], "untw_cyc": untw, "nshell": 1}
    t.update(kw)
    return t


gates = [
    task("g1_ico2_hk_a", "ico2", 1, 4.0 / 0.5168, 16, 3000, nloc=10, untw=2000, nshell=2),
    task("g1_ico2_hk_b", "ico2", 1, 4.0 / 0.1, 32, 3000, nloc=10, untw=2000, nshell=2),
    task("g1_ico2_w_a", "ico2", 0, BETA_U, 16, 3000, nloc=10, untw=2000, nshell=2),
    task("g1_ico2_w_b", "ico2", 0, 30.0, 32, 3000, nloc=10, untw=2000, nshell=2),
    task("g2_cub3_ab", "cub3_L4", 2, 4.0, 16, 1500, untw=400),
    task("g2_cub4_ab", "cub4_L4", 2, 4.0, 16, 800, untw=400),
    task("g2_icoT_ab", "icoT_NT4", 2, 4.0, 16, 1500, untw=400),
    task("g2_cub4_ab_U", "cub4_L4", 2, BETA_U, 24, 600, untw=400),
    task("g3_cub4_w23", "cub4_L6", 0, 2.3, 0, 0, untw=800),
    task("g3_cub4_w25", "cub4_L6", 0, 2.5, 0, 0, untw=800),
]
seed[0] += 3        # the production list held three further engine gates here; the seeds below are those of the runs
gates.append(task("g5_cub4_w25_direct", "cub4_L4", 0, 2.5, 16, 1500, untw=6000))

prod = []
# main scan: three complexes, two SU(2) actions, seven couplings. Interleaved so that every parallel process sees every kind.
for ib, beta in enumerate(BETAS):
    N = 32 if beta > 8.5 else 24
    for kind in (1, 0):
        prod.append(task("cub4_L6_%s_b%d" % (tag(kind), ib), "cub4_L6", kind, beta, N, 400))
        prod.append(task("cub3_L8_%s_b%d" % (tag(kind), ib), "cub3_L8", kind, beta, N, 600))
        prod.append(task("icoT_NT8_%s_b%d" % (tag(kind), ib), "icoT_NT8", kind, beta, N, 1000))
# with the doublet: beta_G at the chain value, beta_H = 2 beta_G (PROOF.md section 9), both sides of the line
for kind in (1, 0):
    for h, start in ((0.0, "pole"), (0.5, "pole"), (1.5, "pole"), (1.5, "equator"), (3.0, "equator"), (4.0, "equator")):
        prod.append(task("higgs_%s_h%02d_%s" % (tag(kind), int(round(10 * h)), start), "cub4_L4", kind, BETA_U, 24, 500,
                         untw=1000, beta_h=2 * BETA_U, h=h, start=start))
    # stiffer doublet on the Higgs side: the frozen-link limit
    for ratio in (8.0, 32.0):
        prod.append(task("higgs_%s_stiff%d" % (tag(kind), int(ratio)), "cub4_L4", kind, BETA_U, 40, 500, untw=1000,
                         beta_h=ratio * BETA_U, h=4.0, start="equator"))
# finite-size and tick-number checks at the chain coupling
for kind in (1, 0):
    prod.append(task("cub4_L4_%s_b5" % tag(kind), "cub4_L4", kind, BETA_U, 24, 600))
    prod.append(task("cub3_L6_%s_b5" % tag(kind), "cub3_L6", kind, BETA_U, 24, 600))
    for NT in (2, 4, 16):
        prod.append(task("icoT_NT%d_%s_b5" % (NT, tag(kind)), "icoT_NT%d" % NT, kind, BETA_U, 24, 1000))
# the Cartan circle through the same engine at the chain coupling (cross-check of the exact gate)
prod.append(task("cub4_L6_ab_b5", "cub4_L6", 2, BETA_U, 24, 400))
prod.append(task("cub3_L8_ab_b5", "cub3_L8", 2, BETA_U, 24, 600))
prod.append(task("icoT_NT8_ab_b5", "icoT_NT8", 2, BETA_U, 24, 1000))
# second reading (A at the midpoint of the face step), engine engine_face_step through staged_face_step.py
fm = [task("fmg_el_cub4_L4_ab_b5", "cub4_L4", 2, BETA_U, 24, 600, mode="el"),
      task("fmg_cub4_L4_ab_b5", "cub4_L4", 2, BETA_U, 40, 500, mode="fm"),
      task("fmg_cub3_L6_ab_b5", "cub3_L6", 2, BETA_U, 32, 600, mode="fm"),
      task("fmg_icoT_NT8_ab_b5", "icoT_NT8", 2, BETA_U, 32, 1000, mode="fm")]
for ib in (5, 2, 4, 6):
    beta = BETAS[ib]
    for kind in (1, 0):
        fm.append(task("fm_cub4_L4_%s_b%d" % (tag(kind), ib), "cub4_L4", kind, beta, 40 if kind == 1 else 32, 400, mode="fm"))
        fm.append(task("fm_cub3_L6_%s_b%d" % (tag(kind), ib), "cub3_L6", kind, beta, 32, 600, mode="fm"))
        fm.append(task("fm_icoT_NT8_%s_b%d" % (tag(kind), ib), "icoT_NT8", kind, beta, 32, 1000, mode="fm"))

if __name__ == "__main__":
    json.dump(gates, open("tasks_gates.json", "w"), indent=1)
    json.dump(prod, open("tasks_prod.json", "w"), indent=1)
    json.dump(fm, open("tasks_fm.json", "w"), indent=1)
    print(len(gates), "gate tasks,", len(prod), "production tasks,", len(fm), "tasks of the second reading")

