"""Every admissible quantum record law for the native seam repair.

The native classical rule fixes which load configurations a seam repair can
produce and with what probability. A quantum record law must also say which
coherences between configurations survive the event. This module classifies
all such laws on the 78-state total-two carrier and checks the classification
exactly against the repository's executed C and E instruments.

Requirements on a law Phi for one seam e=(a,b):

1. Seam locality: Phi acts on the endpoint loads and as the identity on the
   other ten ports.
2. Exact classical law: Phi(|x><x|) is the native repair distribution.
3. Local number-phase covariance: Phi commutes with exp(i sum_p t_p N_p) for
   every real vector t. Equivalently every surviving matrix unit keeps its
   signed load difference.
4. Endpoint-swap covariance: the rotation exchanging a and b, with the
   current reversed, maps the law to itself.

Then Phi(|x><x'|) = sum_j G_j(key x, key x') |x+j D><x'+j D|, where D moves
one unit from a to b and key x=(x_a,x_b). The Choi matrix is the direct sum
of the G_j, so Phi is completely positive exactly when every G_j is positive
semidefinite. Trace preservation fixes the diagonal G_j(k,k)=P(j|k). On the
total-two carrier five multipliers remain:

    k = G_0((0,0),(1,1))   empty seam versus full seam,
    g = G_0((0,0),(0,1))   empty seam versus a seam whose single load waits,
    h = G_0((1,1),(0,1))   full seam versus a waiting seam,
    m = G_0((1,0),(0,1))   the two waiting orientations (real by symmetry),
    n = G_-1((0,1),(0,2))  a moving single load versus a moving collision.

Named members are built from explicit Kraus operators: C0 records only the
current, C records coin and current, E records coin and pre-repair contrast,
T records coin and both endpoint loads, equivalently coin, current and seam
total. The endpoint-prevalue recorder is a
control that keeps the classical law and fails requirement 3.

Arithmetic is exact in Q(sqrt 2)(i). No floating-point value enters.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations, product

from model.observer import native
from model.predictions import observer_instrument as oi

SEAM = native.E
STATES = native.states()
INDEX = {state: i for i, state in enumerate(STATES)}
KEYS = ((0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2))
CURRENTS = (-1, 0, 1)
SETTINGS = oi.SETTINGS


class S:
    """Exact element a + b sqrt2 + i(c + d sqrt2) of Q(sqrt 2)(i)."""

    __slots__ = ("a", "b", "c", "d")

    def __init__(self, a=0, b=0, c=0, d=0):
        self.a, self.b, self.c, self.d = F(a), F(b), F(c), F(d)

    @staticmethod
    def of(value):
        return value if isinstance(value, S) else S(value)

    def parts(self):
        return self.a, self.b, self.c, self.d

    def __add__(self, other):
        o = S.of(other)
        return S(self.a+o.a, self.b+o.b, self.c+o.c, self.d+o.d)

    __radd__ = __add__

    def __neg__(self):
        return S(-self.a, -self.b, -self.c, -self.d)

    def __sub__(self, other):
        return self + (-S.of(other))

    def __rsub__(self, other):
        return S.of(other) - self

    def __mul__(self, other):
        o = S.of(other)

        def mul(x1, y1, x2, y2):
            return x1*x2 + 2*y1*y2, x1*y2 + y1*x2

        pr, qs = mul(self.a, self.b, o.a, o.b), mul(self.c, self.d, o.c, o.d)
        ps, qr = mul(self.a, self.b, o.c, o.d), mul(self.c, self.d, o.a, o.b)
        return S(pr[0]-qs[0], pr[1]-qs[1], ps[0]+qr[0], ps[1]+qr[1])

    __rmul__ = __mul__

    def __pow__(self, exponent):
        if type(exponent) is not int or exponent < 0:
            raise ValueError("nonnegative integer exponents only")
        result = S(1)
        for _ in range(exponent):
            result = result * self
        return result

    def conj(self):
        return S(self.a, self.b, -self.c, -self.d)

    def __eq__(self, other):
        return self.parts() == S.of(other).parts()

    def __hash__(self):
        return hash(self.parts())

    def __bool__(self):
        return any(self.parts())

    def real(self):
        return S(self.a, self.b)

    def is_real(self):
        return not self.c and not self.d

    def sign(self):
        """Exact sign of a real element."""
        if not self.is_real():
            raise ValueError("sign of a nonreal number")
        a, b = self.a, self.b
        if a >= 0 and b >= 0:
            return 0 if not (a or b) else 1
        if a <= 0 and b <= 0:
            return -1
        if a > 0:
            return 1 if a*a > 2*b*b else -1
        return 1 if 2*b*b > a*a else -1

    def __str__(self):
        def surd(x, y):
            if not y:
                return str(x)
            root = "sqrt2" if y == 1 else f"{y}*sqrt2"
            return root if not x else f"{x}+{root}".replace("+-", "-")
        real, imag = surd(self.a, self.b), surd(self.c, self.d)
        if not (self.c or self.d):
            return real
        imag = f"i*({imag})"
        return imag if not (self.a or self.b) else f"{real}+{imag}"

    __repr__ = __str__


ROOT_HALF = S(0, F(1, 2))          # 1/sqrt2
PHASE = {"phase_0": S(1), "phase_pi": S(-1), "phase_pi_over_2": S(0, 0, 1)}


def key(state):
    return state[SEAM[0]], state[SEAM[1]]


def current(seam_key, coin):
    d = seam_key[0] - seam_key[1]
    return (d - (2*coin - 1)*(abs(d) % 2)) // 2


def move(state, j):
    out = list(state)
    out[SEAM[0]] -= j
    out[SEAM[1]] += j
    return tuple(out)


PROBABILITY = {(k, j): F(sum(current(k, c) == j for c in (0, 1)), 2)
               for k in KEYS for j in CURRENTS}
FIBER = {j: tuple(k for k in KEYS if PROBABILITY[k, j]) for j in CURRENTS}


def check_native_rule():
    """The local formula reproduces the repository rule on every state and coin."""
    for state in STATES:
        for coin in (0, 1):
            output, (recorded_coin, j) = native.native(state, SEAM, coin)
            if (recorded_coin, j) != (coin, current(key(state), coin)) or output != move(state, j):
                raise AssertionError("seam formula differs from the native rule")
    return {"states": len(STATES), "coins": 2, "agrees_with_native_rule": True}


# ----------------------------------------------------------------- tables

PARAMETERS = ("k", "g", "h", "m", "n")


def table(k, g, h, m, n):
    values = dict(zip(PARAMETERS, map(S.of, (k, g, h, m, n))))
    if not values["m"].is_real():
        raise ValueError("endpoint-swap covariance makes m real")
    return values


def multipliers(t):
    """The coherence multipliers G_j(k,k') of a parameter table."""
    k, g, h, m, n = (t[p] for p in PARAMETERS)
    half = S(F(1, 2))
    G = {}

    def put(j, x, y, value):
        G[j, x, y] = value
        G[j, y, x] = value.conj()

    for j in CURRENTS:
        for x in FIBER[j]:
            G[j, x, x] = S(PROBABILITY[x, j])
    put(0, (0, 0), (1, 1), k)
    put(0, (0, 0), (1, 0), g)
    put(0, (0, 0), (0, 1), g)
    put(0, (1, 1), (1, 0), h)
    put(0, (1, 1), (0, 1), h)
    put(0, (1, 0), (0, 1), m)
    put(-1, (0, 1), (0, 2), n)
    put(1, (1, 0), (2, 0), n)
    if G[0, (1, 0), (1, 0)] != half:
        raise AssertionError("unexpected classical probability")
    return G


def apply_table(t, rho):
    """Forgotten-record channel of a table on a sparse operator."""
    G = multipliers(t)
    out = {}
    for (x, y), value in rho.items():
        sx, sy = STATES[x], STATES[y]
        for j in CURRENTS:
            coefficient = G.get((j, key(sx), key(sy)))
            if coefficient is not None and coefficient:
                target = INDEX[move(sx, j)], INDEX[move(sy, j)]
                out[target] = out.get(target, S()) + coefficient*value
    return {unit: value for unit, value in out.items() if value}


def compose(first, second):
    """Table of 'second after first'. Derived in README; checked exactly below."""
    return table(first["k"]*second["k"], first["g"]*second["g"],
                 first["h"]*second["h"], first["m"]*second["m"],
                 first["n"]*second["h"].conj())


def power(t, r):
    if type(r) is not int or r < 1:
        raise ValueError("positive repetition counts only")
    result = t
    for _ in range(r-1):
        result = compose(result, t)
    return result


def determinant(matrix):
    size = len(matrix)
    if size == 1:
        return matrix[0][0]
    total = S()
    for column in range(size):
        minor = [row[:column]+row[column+1:] for row in matrix[1:]]
        term = matrix[0][column]*determinant(minor)
        total = total + term if column % 2 == 0 else total - term
    return total


def positive_semidefinite(t):
    """All principal minors of every G_j are nonnegative (exact)."""
    G = multipliers(t)
    for j in CURRENTS:
        keys = FIBER[j]
        for size in range(1, len(keys)+1):
            for chosen in combinations(keys, size):
                minor = determinant([[G[j, x, y] for y in chosen] for x in chosen])
                if minor.sign() < 0:
                    return False
    return True


def text(t):
    return {p: str(t[p]) for p in PARAMETERS}


# ------------------------------------------------------------ instruments

def coin_instrument(record_of, sign_of=lambda seam_key, coin: 1):
    """Kraus maps {record: {input: (output, amplitude)}} from a coin rule."""
    kraus = {}
    for i, state in enumerate(STATES):
        for coin in (0, 1):
            j = current(key(state), coin)
            record = record_of(key(state), coin, j)
            branch = kraus.setdefault(record, {})
            if i in branch:
                raise AssertionError("two coins share one record branch for one input")
            branch[i] = INDEX[move(state, j)], ROOT_HALF*sign_of(key(state), coin)
    return kraus


def minimal_instrument():
    """C0: one Kraus map per current, amplitude sqrt P(j|key)."""
    kraus = {}
    for i, state in enumerate(STATES):
        for j in CURRENTS:
            probability = PROBABILITY[key(state), j]
            if probability:
                amplitude = S(1) if probability == 1 else ROOT_HALF
                kraus.setdefault(("current", j), {})[i] = INDEX[move(state, j)], amplitude
    return kraus


INSTRUMENTS = {
    "C0": minimal_instrument,
    "C": lambda: coin_instrument(lambda k, c, j: ("coin_current", c, j)),
    "E": lambda: coin_instrument(lambda k, c, j: ("coin_contrast", c, k[0]-k[1])),
    "T": lambda: coin_instrument(lambda k, c, j: ("coin_occupation", c, k)),
    "C_signed": lambda: coin_instrument(
        lambda k, c, j: ("coin_current", c, j),
        lambda k, c: -1 if (k, c) in {((1, 0), 1), ((0, 1), 0)} else 1),
}
CONTROL = lambda: coin_instrument(lambda k, c, j: ("coin_prevalue", c, k[0]))


def completeness(kraus):
    """sum_r K_r^dagger K_r = I on the 78-state carrier, exactly."""
    gram = {}
    for branch in kraus.values():
        by_output = {}
        for i, (y, amplitude) in branch.items():
            by_output.setdefault(y, []).append((i, amplitude))
        for entries in by_output.values():
            for (i, a), (l, b) in product(entries, repeat=2):
                gram[i, l] = gram.get((i, l), S()) + a.conj()*b
    for (i, l), value in gram.items():
        if value != (1 if i == l else 0):
            return False
    return all(gram.get((i, i)) == 1 for i in range(len(STATES)))


def apply_instrument(kraus, rho):
    out = {}
    for branch in kraus.values():
        for (x, y), value in rho.items():
            if x in branch and y in branch:
                (u, a), (v, b) = branch[x], branch[y]
                out[u, v] = out.get((u, v), S()) + a*b.conj()*value
    return {unit: value for unit, value in out.items() if value}


def definite_currents(kraus):
    """Each Kraus branch translates every input by a single current."""
    for branch in kraus.values():
        moved = {key(STATES[y])[0] - key(STATES[x])[0] for x, (y, _) in branch.items()}
        if len(moved) != 1:
            return False
    return True


def difference_preserving(kraus):
    """Every output matrix unit keeps the input's signed load difference."""
    for x, y in product(range(len(STATES)), repeat=2):
        for u, v in apply_instrument(kraus, {(x, y): S(1)}):
            gap_in = [p-q for p, q in zip(STATES[x], STATES[y])]
            gap_out = [p-q for p, q in zip(STATES[u], STATES[v])]
            if gap_in != gap_out:
                return False, {"input": [list(STATES[x]), list(STATES[y])],
                               "output": [list(STATES[u]), list(STATES[v])]}
    return True, None


def extract_table(kraus):
    """Read k,g,h,m,n from Kraus maps; demand that every carrier pair agrees."""
    seen = {}
    for branch in kraus.values():
        for x, y in product(branch, repeat=2):
            (u, a), (v, b) = branch[x], branch[y]
            j = key(STATES[x])[0] - key(STATES[u])[0]
            if key(STATES[y])[0] - key(STATES[v])[0] != j:
                raise AssertionError("branch mixes currents")
            label = j, key(STATES[x]), key(STATES[y]), x, y
            seen[label] = seen.get(label, S()) + a*b.conj()
    local = {}
    for (j, kx, ky, x, y), value in seen.items():
        if local.setdefault((j, kx, ky), value) != value:
            raise AssertionError("multiplier depends on spectator ports")
    get = lambda j, x, y: local.get((j, x, y), S())
    return table(get(0, (0, 0), (1, 1)), get(0, (0, 0), (0, 1)), get(0, (1, 1), (0, 1)),
                 get(0, (1, 0), (0, 1)), get(-1, (0, 1), (0, 2)))


def same_channel_everywhere(kraus, t):
    """The instrument's channel equals the table's on all 6,084 matrix units."""
    for x, y in product(range(len(STATES)), repeat=2):
        unit = {(x, y): S(1)}
        if apply_instrument(kraus, unit) != apply_table(t, unit):
            return False
    return True


def same_as_repository(name, kraus):
    """C and E coincide branch by branch with the executed predictions module."""
    support = oi.kraus_support(name)
    mine = {}
    for record, branch in kraus.items():
        label = (record[1], record[2])
        mine[label] = tuple(sorted((i, y) for i, (y, a) in branch.items()))
        if any(a != ROOT_HALF for _, a in branch.values()):
            return False
    return mine == {record: tuple(sorted(entries)) for record, entries in support.items()}


# ----------------------------------------------------------- the ladder

def diagonal_after(kraus, refine):
    """Record one more classical function of the post-repair seam occupation.

    refine(record, output_key) returns (new_record, amplitude). The amplitudes
    must form a complete diagonal instrument on each branch's range.
    """
    result = {}
    for record, branch in kraus.items():
        for i, (y, a) in branch.items():
            new_record, factor = refine(record, key(STATES[y]))
            if factor:
                result.setdefault(new_record, {})[i] = y, a*factor
    return result


def same_instrument(left, right, relabel):
    mapped = {}
    for record, branch in left.items():
        target = relabel(record)
        if target in mapped:
            return False
        mapped[target] = branch
    return mapped == right


def ladder():
    """C0 -> C (coin), C -> E (seam parity), C -> T (seam total), exactly."""
    C0, C, E, T = (INSTRUMENTS[name]() for name in ("C0", "C", "E", "T"))

    def coin_refinement(record, out_key):
        # The coin is a diagonal readout on C0's output: weight 1/2 on a
        # balanced output, a sharp label on an odd one.
        j = record[1]
        results = []
        source = (out_key[0]+j, out_key[1]-j)
        for coin in (0, 1):
            weight = F(int(current(source, coin) == j), 2) / PROBABILITY[source, j]
            amplitude = S(1) if weight == 1 else (ROOT_HALF if weight == F(1, 2) else S())
            results.append((("coin_current", coin, j), amplitude))
        return results

    coin_split = {}
    for record, branch in C0.items():
        for i, (y, a) in branch.items():
            for new_record, factor in coin_refinement(record, key(STATES[y])):
                if factor:
                    coin_split.setdefault(new_record, {})[i] = y, a*factor
    parity = diagonal_after(C, lambda r, k: ((r[1], r[2], sum(k) % 2), S(1)))
    total = diagonal_after(C, lambda r, k: ((r[1], r[2], sum(k)), S(1)))

    def contrast_label(r):
        coin, j, p = r
        return "coin_contrast", coin, 2*j + (2*coin-1)*p

    def occupation_label(r):
        coin, j, s = r
        d = 2*j + (2*coin-1)*(s % 2)
        return "coin_occupation", coin, ((s+d)//2, (s-d)//2)

    return {"C_is_C0_followed_by_a_coin_readout": coin_split == C,
            "E_is_C_followed_by_seam_parity": same_instrument(parity, E, contrast_label),
            "T_is_C_followed_by_seam_total": same_instrument(total, T, occupation_label)}


def coin_partition_laws():
    """Coin-recording laws whose records refine (coin, current), with swap symmetry."""
    fibers = {}
    for k in KEYS:
        fibers.setdefault(current(k, 1), []).append(k)

    def partitions(items):
        if not items:
            yield []
            return
        first, rest = items[0], items[1:]
        for part in partitions(rest):
            for i in range(len(part)):
                yield part[:i] + [[first]+part[i]] + part[i+1:]
            yield [[first]] + part

    laws = []
    for choice in product(*(list(partitions(sorted(v))) for _, v in sorted(fibers.items()))):
        blocks = [frozenset(b) for part in choice for b in part]
        label = {k: i for i, block in enumerate(blocks) for k in block}

        def record_of(k, c, j, label=label):
            block = label[k] if c == 1 else label[(k[1], k[0])]
            return ("partition", c, block)
        laws.append(coin_instrument(record_of))
    return laws


# ----------------------------------------------------------- observables

def prepared(setting, second):
    """(|X> + e^{i theta}|second>)/sqrt2 or the equal mixture; X = loads at 2,3."""
    x = oi.X
    rho = {(x, x): S(F(1, 2)), (second, second): S(F(1, 2))}
    if setting != "diagonal_control":
        phase = PHASE[setting]
        rho[x, second] = phase.conj()*F(1, 2)
        rho[second, x] = phase*F(1, 2)
    return rho


def plus_probability(output, second):
    """<phi_0|rho|phi_0> with phi_0=(|X>+|second>)/sqrt2."""
    pair = (oi.X, second)
    return sum((output.get((a, b), S()) for a in pair for b in pair), S())*F(1, 2)


def wait_coherence_prediction(t, setting):
    """3/8 + Re(g e^{-i theta})/2, the closed form derived in README."""
    if setting == "diagonal_control":
        return S(F(3, 8))
    return S(F(3, 8)) + (t["g"]*PHASE[setting].conj()).real()*F(1, 2)


FULL = INDEX[tuple(1 if p in (0, 1) else 0 for p in range(12))]


def full_coherence_prediction(t, setting):
    """1/2 + Re(k e^{-i theta})/2 for an empty versus a full seam."""
    if setting == "diagonal_control":
        return S(F(1, 2))
    return S(F(1, 2)) + (t["k"]*PHASE[setting].conj()).real()*F(1, 2)


def classical_kernel(kraus):
    rows = []
    for x in range(len(STATES)):
        out = apply_instrument(kraus, {(x, x): S(1)})
        rows.append(tuple(sorted((u, str(v)) for (u, w), v in out.items() if u == w)))
    return tuple(rows)


# ------------------------------------------------------------------ run

def generic_tables():
    """Exact test points: two admissible laws off the ladder and one outside.

    C0_quarter_turn gives the waiting configurations amplitude i/sqrt2 instead
    of 1/sqrt2; it is rank one per current and therefore admissible. The
    mixture is a convex combination of C and T. The outside point violates
    positivity; the composition law is algebraic and is checked there too.
    """
    quarter = S(0, 0, 0, F(1, 2))            # i/sqrt2
    return {"C0_quarter_turn": table(1, quarter.conj(), quarter.conj(), F(1, 2), quarter),
            "C_T_mixture": table(F(1, 2), F(1, 4), F(1, 4), 0, F(1, 4)),
            "outside": table(S(0, 0, 1), F(3, 4), 1, F(1, 2), F(5, 7))}


REPRESENTATIVE = {(0, 0): (2, 3), (1, 0): (0, 2), (0, 1): (1, 2),
                  (1, 1): (0, 1), (2, 0): (0, 0), (0, 2): (1, 1)}


def representative(seam_key):
    loads = [0]*12
    for port in REPRESENTATIVE[seam_key]:
        loads[port] += 1
    state = tuple(loads)
    if key(state) != seam_key:
        raise AssertionError("bad representative")
    return INDEX[state]


def run():
    native_check = check_native_rule()
    members, tables = {}, {}
    for name, build in INSTRUMENTS.items():
        kraus = build()
        t = extract_table(kraus)
        tables[name] = t
        covariant, _ = difference_preserving(kraus)
        members[name] = {
            "record_values_on_seam": len(kraus),
            "complete": completeness(kraus),
            "definite_current_per_branch": definite_currents(kraus),
            "number_phase_covariant": covariant,
            "table": text(t),
            "channel_equals_table_on_all_matrix_units": same_channel_everywhere(kraus, t),
            "positive_semidefinite": positive_semidefinite(t),
            "idempotent": compose(t, t) == t,
            "wait_coherence_plus": {s: str(plus_probability(apply_table(t, prepared(s, oi.Z)), oi.Z))
                                    for s in SETTINGS},
            "full_coherence_plus": {s: str(plus_probability(apply_table(t, prepared(s, FULL)), FULL))
                                    for s in ("phase_0", "phase_pi")},
        }
    for name in ("C", "E"):
        members[name]["equals_executed_instrument"] = same_as_repository(name, INSTRUMENTS[name]())
        members[name]["equals_executed_prediction"] = all(
            S(oi.prediction(name, s)) == plus_probability(apply_table(tables[name], prepared(s, oi.Z)), oi.Z)
            for s in SETTINGS)

    control = CONTROL()
    control_covariant, control_witness = difference_preserving(control)
    kernels = {name: classical_kernel(build()) for name, build in INSTRUMENTS.items()}
    kernels["prevalue_control"] = classical_kernel(control)

    # apply_table depends on spectators only through the moved state, so one
    # representative per seam-key pair covers the law; two pairs also run on
    # every matrix unit of the carrier.
    composition_checks = 0
    candidates = dict(tables, **generic_tables())
    units = [(representative(a), representative(b)) for a in KEYS for b in KEYS]
    for (_, first), (_, second) in product(candidates.items(), repeat=2):
        for x, y in units:
            unit = {(x, y): S(1)}
            if apply_table(second, apply_table(first, unit)) != apply_table(compose(first, second), unit):
                raise AssertionError("composition law fails")
        composition_checks += 1
    full_carrier_pairs = (("C", "C"), ("C0", "E"))
    for left, right in full_carrier_pairs:
        for x, y in product(range(len(STATES)), repeat=2):
            unit = {(x, y): S(1)}
            if (apply_table(tables[right], apply_table(tables[left], unit))
                    != apply_table(compose(tables[left], tables[right]), unit)):
                raise AssertionError("composition law fails on the full carrier")

    for name, t in candidates.items():
        for s in SETTINGS:
            if plus_probability(apply_table(t, prepared(s, oi.Z)), oi.Z) != wait_coherence_prediction(t, s):
                raise AssertionError(f"wait-coherence formula fails for {name}")
            if plus_probability(apply_table(t, prepared(s, FULL)), FULL) != full_coherence_prediction(t, s):
                raise AssertionError(f"full-coherence formula fails for {name}")

    partition = [extract_table(law) for law in coin_partition_laws()]
    distinct = {tuple(t[p] for p in PARAMETERS) for t in partition}
    idempotent_partitions = [t for t in partition if compose(t, t) == t]

    completion = {}
    for name in ("C0", "C", "C_signed"):
        t = tables[name]
        completion[name] = {"after_8_repairs": text(power(t, 8)),
                            "limit": "E" if t["k"] == 1 else ("T" if t["k"] == 0 else "none")}
    C, E = tables["C"], tables["E"]
    halving = all(
        power(C, r) == table(*(E[p] + (C[p]-E[p])*F(2)**(1-r) for p in PARAMETERS))
        for r in range(1, 9))

    bounds = {
        "wait_coherence_bound_squared": "1/2",
        "C0_attains_the_bound": tables["C0"]["g"]*tables["C0"]["g"].conj() == S(F(1, 2)),
        "every_member_within_the_bound": all(
            (t["g"]*t["g"].conj() - F(1, 2)).sign() <= 0 for t in tables.values()),
        "coin_recording_laws_bound_by_one_half": all(
            (t["g"]*t["g"].conj() - F(1, 4)).sign() <= 0 for t in partition),
        "outside_point_rejected": not positive_semidefinite(generic_tables()["outside"]),
        "generic_admissible_points_accepted": all(
            positive_semidefinite(generic_tables()[n]) for n in ("C0_quarter_turn", "C_T_mixture")),
        "exact_full_coherence_forces_g_equal_h": all(
            t["g"] == t["h"] for t in tables.values() if t["k"] == 1),
    }

    result = {
        "schema": "oph.record_laws.v1",
        "carrier": {"states": len(STATES), "binary_states": len(native.states(binary=True)),
                    "seam": list(SEAM), "matrix_units": len(STATES)**2},
        "native_rule": native_check,
        "classical_probabilities": {f"{k}|{j}": str(PROBABILITY[k, j])
                                    for k in KEYS for j in CURRENTS if PROBABILITY[k, j]},
        "members": members,
        "prevalue_control": {
            "complete": completeness(control),
            "definite_current_per_branch": definite_currents(control),
            "number_phase_covariant": control_covariant,
            "witness": control_witness,
            "same_classical_law_as_C": kernels["prevalue_control"] == kernels["C"]},
        "same_classical_law_for_every_member": len({kernels[n] for n in INSTRUMENTS}) == 1,
        "ladder": ladder(),
        "composition": {"law": "(k,g,h,m,n) then (k',g',h',m',n') = (kk',gg',hh',mm',n conj(h'))",
                        "ordered_pairs_checked_on_every_seam_key_pair": composition_checks,
                        "ordered_pairs_checked_on_all_matrix_units": [list(p) for p in full_carrier_pairs]},
        "idempotent_members": sorted(n for n in members if members[n]["idempotent"]),
        "completion": completion,
        "C_power_identity_C^r=E+2^(1-r)(C-E)_for_r_1_to_8": halving,
        "coin_partition_laws": {
            "count": len(partition), "distinct_channels": len(distinct),
            "idempotent": len(idempotent_partitions),
            "wait_coherence_values": sorted({str(t["g"]) for t in partition}),
            "with_C_interference": sum(t["g"] == F(1, 2) for t in partition)},
        "bounds": bounds,
        "interference_closed_forms": {
            "wait": "plus(theta) = 3/8 + Re(g exp(-i theta))/2",
            "full": "plus(theta) = 1/2 + Re(k exp(-i theta))/2"},
        "scope": {"single_seam_classification": True,
                  "carrier_covariance_by_identical_law_on_every_seam": True,
                  "network_coherence_rates_beyond_C_and_E_computed": False,
                  "laboratory_data_used": False},
    }
    for name, row in members.items():
        if not (row["complete"] and row["definite_current_per_branch"] and row["number_phase_covariant"]
                and row["channel_equals_table_on_all_matrix_units"] and row["positive_semidefinite"]):
            raise AssertionError(f"member {name} failed an admissibility check")
    if control_covariant or not result["prevalue_control"]["same_classical_law_as_C"]:
        raise AssertionError("the prevalue control no longer separates covariance from the classical law")
    if not (all(result["ladder"].values()) and halving and result["same_classical_law_for_every_member"]):
        raise AssertionError("ladder, power identity or classical equality failed")
    if result["idempotent_members"] != ["E", "T"]:
        raise AssertionError("unexpected idempotent members")
    if not all(v for k, v in bounds.items() if k != "wait_coherence_bound_squared"):
        raise AssertionError("bound checks failed")
    return result
