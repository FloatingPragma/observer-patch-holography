"""Bounded integer half-seam records joined to the port/transit generator.

J02 Form B, with integer rather than modular records. The finite regulator
rejects hops beyond |F_h| <= cutoff. Every rejection is retained. Native commit
measures position AND records, applies the literal one-load C/E event on ports,
and leaves transit states unchanged. This explicit monitored instrument is a
supplied architecture, not the seam-decay channel or its fast-commit theorem.
"""
from collections import deque
from dataclasses import dataclass
import math
import operator
import numpy as np
from .construction import Construction


@dataclass
class JointRecords:
    construction: Construction
    initial_port: int = 0
    cutoff: int = 1
    max_states: int = 4096

    def __post_init__(self):
        c = self.construction
        if type(self.cutoff) is not int or self.cutoff < 1:
            raise ValueError("positive integer record cutoff required")
        try:
            self.initial_port = operator.index(self.initial_port)
            self.max_states = operator.index(self.max_states)
        except TypeError as exc:
            raise ValueError("integer initial port and max_states required") from exc
        if self.max_states < 1:
            raise ValueError("positive max_states required")
        if not 0 <= self.initial_port < c.ports:
            raise ValueError("invalid initial port")
        start = (self.initial_port, (0,) * (2 * len(c.seams)))
        states, seen, queue = [start], {start}, deque([start])
        while queue:
            state = queue.popleft()
            for target, _ in self.half_hops(state):
                if target not in seen:
                    if len(states) >= self.max_states:
                        raise ValueError("joint sector exceeds max_states")
                    seen.add(target)
                    states.append(target)
                    queue.append(target)
        self.states = tuple(states)
        self.index = {state: i for i, state in enumerate(states)}

    def half_hops(self, state):
        pos, records = state
        c = self.construction
        for edge, endpoints in enumerate(c.seams):
            for side, port in enumerate(endpoints):
                if pos == port:
                    destination, delta = c.ports + edge, 1
                elif pos == c.ports + edge:
                    destination, delta = port, -1
                else:
                    continue
                h = 2 * edge + side
                f = list(records)
                f[h] += delta
                if abs(f[h]) <= self.cutoff:
                    yield (destination, tuple(f)), (1 if side == 0 else -1)

    def gauss(self, state):
        pos, f = state
        c = self.construction
        g = [int(pos == x) for x in range(c.ports + len(c.seams))]
        for edge, (a, b) in enumerate(c.seams):
            g[a] += f[2*edge]
            g[b] += f[2*edge+1]
            g[c.ports+edge] -= f[2*edge] + f[2*edge+1]
        return tuple(g)

    @property
    def hamiltonian(self):
        c = self.construction
        h = np.zeros((len(self.states), len(self.states)))
        for col, state in enumerate(self.states):
            h[col,col] = c.mass * (1 if state[0] < c.ports else -1)
            for target, sign in self.half_hops(state):
                h[self.index[target], col] += sign
        return h

    def initial_density(self):
        rho = np.zeros((len(self.states), len(self.states)), complex)
        rho[0,0] = 1
        return rho

    def native_target(self, state, seam, coin):
        """One-load native C-current event in the half-seam register basis."""
        pos, records = state
        c = self.construction
        try:
            seam, coin = operator.index(seam), operator.index(coin)
        except TypeError as exc:
            raise ValueError("integer event indices required") from exc
        if not 0 <= seam < len(c.seams) or coin not in (0,1):
            raise ValueError("invalid event")
        if pos >= c.ports:
            return state, "transit_noop", 0
        n = tuple(int(pos == x) for x in range(c.ports))
        n2, _, current = c.committed_event(n, (0,)*len(c.seams), seam, coin)
        f = list(records)
        f[2*seam] += current
        f[2*seam+1] -= current
        target = (n2.index(1), tuple(f))
        if target not in self.index:
            return state, "boundary_rejected", current
        return target, "native", current

    def commit(self, rho, seam=None, coin=None, law=None):
        """CPTP measured native instrument, including noops and rejections.

        Each outcome has K=sqrt(p)|target><source|. Source record, seam and
        coin are retained in the returned ledger. Thus sum K* K = identity;
        collisions of target states do not silently identify distinct outcomes.
        """
        rho = self._density(rho)
        for value in (seam,coin):
            if value is not None:
                try:
                    operator.index(value)
                except TypeError as exc:
                    raise ValueError("integer event indices required") from exc
        if law not in (None, "C", "E"):
            raise ValueError("record law must be C or E")
        events = [(e,q,l) for e in range(len(self.construction.seams)) for q in (0,1)
                  for l in (("C", "E") if law is None else (law,))
                  if (seam is None or seam == e) and (coin is None or coin == q)]
        if not events:
            raise ValueError("empty event distribution")
        out = np.zeros_like(rho)
        ledger = []
        for source, state in enumerate(self.states):
            for e,q,record_law in events:
                target, status, current = self.native_target(state,e,q)
                probability = float(rho[source,source].real) / len(events)
                out[self.index[target],self.index[target]] += probability
                ledger.append({"source": source, "target": self.index[target],
                               "seam": e, "coin": q,
                               "proposed_current": current,
                               "current": current if status == "native" else 0,
                               "native_event_executed": status == "native",
                               "applied_current": current if status == "native" else 0,
                               "law": record_law,
                               "native_record": (q, current if record_law == "C" else
                                   int(state[0] == self.construction.seams[e][0])-int(state[0] == self.construction.seams[e][1]))
                                   if status == "native" else None,
                               "status": status, "probability": probability})
        return out, ledger

    def _density(self, rho):
        rho = np.asarray(rho, complex)
        n = len(self.states)
        if rho.shape != (n,n) or not np.all(np.isfinite(rho)):
            raise ValueError("invalid density dimensions or values")
        if not np.allclose(rho,rho.conj().T) or not np.isclose(np.trace(rho),1):
            raise ValueError("density must be Hermitian with unit trace")
        if np.linalg.eigvalsh(rho).min() < -1e-10:
            raise ValueError("density is not positive")
        return rho

    def step(self, rho, time, commit_probability):
        """Declared discrete unitary-then-monitored-event protocol.

        No continuous-time or fast-commit convergence is asserted. At p=0
        this is joint coherent evolution, at time=0,p=1 the native port law
        where the boundary is not touched. Records of no commit are retained.
        """
        if not math.isfinite(time) or not math.isfinite(commit_probability) or not 0 <= commit_probability <= 1:
            raise ValueError("invalid protocol parameters")
        rho = self._density(rho)
        eigenvalues, vectors = np.linalg.eigh(self.hamiltonian)
        unitary = (vectors * np.exp(-1j*time*eigenvalues)) @ vectors.T
        evolved = unitary @ rho @ unitary.conj().T
        committed, ledger = self.commit(evolved)
        for row in ledger:
            row["probability"] *= commit_probability
        ledger.append({"status": "no_commit", "probability": 1-commit_probability})
        return (1-commit_probability)*evolved + commit_probability*committed, ledger


def cutoff_experiment():
    """Executed numerical regulator comparison, not a cutoff-removal theorem."""
    rows = []
    c = Construction(3, ((0,1),(1,2),(2,0)), .23)
    for cutoff in (1,2):
        joint = JointRecords(c,cutoff=cutoff)
        rho = joint.initial_density()
        rejected = 0.
        for _ in range(8):
            rho, ledger = joint.step(rho,.37,.29)
            rejected += sum(r['probability'] for r in ledger if r['status']=='boundary_rejected')
        populations = np.zeros(c.ports+len(c.seams))
        for i,(pos,_) in enumerate(joint.states):
            populations[pos] += rho[i,i].real
        rows.append({"cutoff":cutoff,"dimension":len(joint.states),
                     "position_probabilities":populations.tolist(),
                     "trace":float(np.trace(rho).real),
                     "sum_expected_rejections_over_steps":rejected})
    return {"graph":"oriented triangle", "mass":.23, "steps":8,
            "time_per_step":.37,"commit_probability":.29,"rows":rows,
            "maximum_position_probability_drift":float(np.max(np.abs(
                np.array(rows[0]['position_probabilities'])-rows[1]['position_probabilities']))),
            "scope":"Two finite cutoffs only; cumulative expected rejection count is not an error bound."}


def run():
    joint = JointRecords(Construction(3, ((0,1),(1,2),(2,0)), .23))
    rho, ledger = joint.step(joint.initial_density(), .17, .37)
    gauss = joint.gauss(joint.states[0])
    native, _ = joint.commit(joint.initial_density())
    populations = np.zeros(3)
    for i,(pos,_) in enumerate(joint.states):
        if pos < 3:
            populations[pos] += native[i,i].real
    expected = (np.eye(3)-joint.construction.incidence@joint.construction.incidence.T/6)[:,0]
    return {"scope": "Bounded integer Gauss sector, monitored one-load native C-current instrument, discrete coherent/commit protocol; no continuum or fast-commit claim.",
            "checks": {"joint_hermitian": bool(np.allclose(joint.hamiltonian,joint.hamiltonian.T)),
                       "exact_gauss": all(joint.gauss(s)==gauss for s in joint.states),
                       "native_population": bool(np.allclose(populations,expected)),
                       "trace": bool(np.isclose(np.trace(rho),1)),
                       "positive": bool(np.linalg.eigvalsh(rho).min()>-1e-12),
                       "ledger_complete": bool(np.isclose(sum(r['probability'] for r in ledger),1))},
            "mutations": {"missing_half_record_detected": joint.gauss((1,(0,)*6)) != gauss},
            "full_criterion_met": False,
            "cutoff_comparison": cutoff_experiment(),
            "supplied": {"integer_cutoff": joint.cutoff,"states":len(joint.states),
                         "commit_monitor": "position and all half-seam records",
                         "transit_event": "retained noop", "boundary": "coherent forbidden half-hops omitted (Hermitian hard wall); monitored native attempts retained as rejections",
                         "commit_probability": .37, "time": .17,
                         "record_law": "equal C (coin,current) and E (coin,pre-repair contrast); separate seam-offset dynamics not implemented"}}


if __name__ == '__main__':
    import json
    print(json.dumps(run(),indent=2))
