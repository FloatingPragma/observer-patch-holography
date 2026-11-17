"""Executed finite chart agreement for the coherent gauge observer.

The frame wrapper is a declared apparatus service. Its matrices are compiled
from the same finite carrier-cover words as the internal transport. A charted
detector is executed by undoing that frame, using the owned register detector,
then restoring the frame. No arbitrary positive-effect service is inferred.
"""
from functools import lru_cache
from pathlib import Path
import argparse
import hashlib
import json
import sympy as s
from model.particle_common_contract.provenance import import_closure
from .common import require
from .gauge_walk import local
from .holonomy import certificate as holonomy_certificate
from .program import Observer, State

ROOT = Path(__file__).resolve().parents[2]


def input_hashes():
    """Pin the executable import closure and dynamically loaded carrier bytes."""
    seeds = {
        'model/coherent_seam_repair/agreement.py',
        'model/particle_common_contract/provenance.py',
        'model/requirements.txt',
        'Lean/Screen/PrimitivePortFrameQuotient.lean',
        'Lean/Screen/PortFrameGram.lean',
        'evidence/native_geometric_source_20260925/spec.json',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/__init__.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/__init__.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/__init__.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py',
    }
    paths = {ROOT/name for name in import_closure(seeds)}
    paths = {path for path in paths if path.name in {p.name for p in path.parent.iterdir()}}
    return {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(paths)}


def _expand(value):
    return s.expand(value)


def _matrix_zero(matrix):
    return all(_expand(x) == 0 for x in matrix)


def frame_from_lift(row):
    """Actual five-level lift, reconstructed from its exact quaternion record."""
    q = [s.Rational(a) + s.Rational(b)*s.sqrt(5) for a, b in row['quaternion']]
    w, x, y, z = q
    rotation = s.Matrix([
        [1-2*(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y)],
        [2*(x*y+w*z), 1-2*(x*x+z*z), 2*(y*z-w*x)],
        [2*(x*z-w*y), 2*(y*z+w*x), 1-2*(x*x+y*y)],
    ]).applyfunc(_expand)
    w, x, y, z = [_expand(a.xreplace({s.sqrt(5): -s.sqrt(5)})) for a in q]
    spin = s.Matrix([[w+s.I*x, y+s.I*z], [-y+s.I*z, w-s.I*x]])
    matrix = s.diag(rotation, spin).applyfunc(_expand)
    require(_matrix_zero(matrix.H*matrix-s.eye(5)), 'cover frame is unitary')
    require(_expand(matrix.det()) == 1, 'cover frame has determinant one')
    return matrix


class ChartedObserver:
    """The same observer in a source-bound, locally controlled code frame.

    ``words`` maps a source port to a compiled finite-cover word. The links
    transform as G_ba -> V_b G_ba V_a*. The raw pointer label is interpreted
    after the matching frame wrapper; bare pointer readout is a different query.
    """

    def __init__(self, observer, words):
        if not isinstance(observer, Observer) or observer.fibre != 5:
            raise ValueError('charted gauge observer requires the five-level source')
        rows = {tuple(row['word']): row for row in holonomy_certificate()['compiler']}
        self.observer = observer
        self.words = {}
        self.frames = [s.eye(5) for _ in range(12)]
        for port, word in words.items():
            if type(port) is not int or not 0 <= port < 12:
                raise ValueError('frame port outside the carrier')
            if any(type(label) is not int for label in word) or tuple(word) not in rows:
                raise ValueError('frame is not a compiled cover word')
            self.words[port] = tuple(word)
            self.frames[port] = frame_from_lift(rows[tuple(word)])
        self.links = {}
        for a, b in observer.addresses:
            G = s.diag(1, 1, 1, -1, -1) if tuple(sorted((a,b))) == observer.central_edge else s.eye(5)
            self.links[(a,b)] = (self.frames[b]*G*self.frames[a].H).applyfunc(_expand)

    def reframe(self, state, inverse=False):
        self.observer._validate(state)
        values = list(state.amplitudes)
        for port in self.words:
            frame = self.frames[port].H if inverse else self.frames[port]
            for parity in (0,1):
                for address in range(self.observer.address_count):
                    indices = [self.observer.index(port,parity,address,j) for j in range(5)]
                    vector = [state.amplitudes[i] for i in indices]
                    for j, index in enumerate(indices):
                        values[index] = _expand(sum(frame[j,k]*vector[k] for k in range(5)))
        return State(tuple(values), state.norm_denominator)

    def repair(self, state):
        """Execute the transformed link law directly, without sandwiching Q."""
        self.observer._validate(state)
        x, out = state.amplitudes, list(state.amplitudes)
        for address, (a,b) in enumerate(self.observer.addresses):
            G = self.links[(a,b)]
            a0 = [self.observer.index(a,0,address,j) for j in range(5)]
            a1 = [self.observer.index(a,1,address,j) for j in range(5)]
            b0 = [self.observer.index(b,0,address,j) for j in range(5)]
            b1 = [self.observer.index(b,1,address,j) for j in range(5)]
            u = [_expand(x[a0[j]]-x[a1[j]]-
                  sum(G[k,j].conjugate()*(x[b0[k]]+x[b1[k]]) for k in range(5)))
                 for j in range(5)]
            gu = [_expand(sum(G[j,k]*u[k] for k in range(5))) for j in range(5)]
            for j in range(5):
                out[a0[j]], out[a1[j]] = _expand(x[a0[j]]-u[j]/2), _expand(x[a1[j]]+u[j]/2)
                out[b0[j]], out[b1[j]] = _expand(x[b0[j]]+gu[j]/2), _expand(x[b1[j]]+gu[j]/2)
        return State(tuple(out), state.norm_denominator)

    def step(self, state):
        return self.repair(self.observer.reflect(state))

    def preparation_filter(self, state):
        self.observer._validate(state)
        evolved = self.step(state).amplitudes
        yes = tuple(_expand((2*x+y)/3) for x,y in zip(state.amplitudes,evolved))
        no = tuple(_expand(s.sqrt(2)*(x-y)/3) for x,y in zip(state.amplitudes,evolved))
        return (({'kind':'preparation','outcome':'success'},State(yes,state.norm_denominator)),
                ({'kind':'preparation','outcome':'failure'},State(no,state.norm_denominator)))

    def quarter_turn(self, state, address):
        unframed = self.reframe(state,inverse=True)
        acted, record = self.observer.quarter_turn(unframed,address)
        return self.reframe(acted), record

    def register_measurement(self, state, address):
        unframed = self.reframe(state,inverse=True)
        return tuple((record,self.reframe(branch))
                     for record,branch in self.observer.register_measurement(unframed,address))


def _states_equal(left, right):
    return left.norm_denominator == right.norm_denominator and all(
        _expand(a-b) == 0 for a,b in zip(left.amplitudes,right.amplitudes))


@lru_cache(maxsize=1)
def certificate():
    controller = holonomy_certificate()
    # This rational Hurwitz quaternion belongs to the generated source cover.
    row = next(r for r in controller['compiler'] if r['quaternion'] == [['1/2','0']]*4)
    native = Observer(fibre=5,central_edge=(0,1))
    charted = ChartedObserver(native,{0:row['word']})
    Q5 = s.diag(2,2,2,-3,-3)
    require(all(_matrix_zero(V*Q5-Q5*V) for V in charted.frames), 'frames preserve charge blocks')

    local_checks = 0
    for a,b in native.addresses:
        G = s.diag(1,1,1,-1,-1) if tuple(sorted((a,b))) == native.central_edge else s.eye(5)
        V = s.kronecker_product(s.diag(charted.frames[a],charted.frames[b]),s.eye(2))
        require(_matrix_zero(local(charted.links[(a,b)])[0]-V*local(G)[0]*V.H),
                'every oriented seam has the same rechart covariance')
        local_checks += 1

    initial = native.embedded_port(0,internal=3)
    reframed = charted.reframe(initial)
    require(_states_equal(charted.reframe(reframed,inverse=True),initial), 'frame inverse execution')
    require(reframed.norm_squared() == initial.norm_squared() == 1, 'frame preserves probability')
    require(_states_equal(charted.step(reframed),charted.reframe(native.step(initial))),
            'direct transformed links execute the conjugated carrier walk')
    address = native.index(0,0,0,internal=3)
    wrong_probability = native.register_measurement(reframed,address)[0][1].norm_squared()
    correct_probability = native.register_measurement(initial,address)[0][1].norm_squared()
    require(wrong_probability == s.Rational(1,120) and correct_probability == s.Rational(1,60),
            'untransported detector is a strict countercontrol')

    probabilities = []
    old_filtered, new_filtered = native.preparation_filter(initial), charted.preparation_filter(reframed)
    for (old_label,old),(new_label,new) in zip(old_filtered,new_filtered):
        require(old_label == new_label and _states_equal(new,charted.reframe(old)), 'both preparation outcomes intertwine')
        # Record-adaptive setting, with failure retained instead of discarded.
        probe = address if old_label['outcome']=='success' else native.index(1,0,0,internal=3)
        old, ro = native.quarter_turn(old,probe)
        new, rn = charted.quarter_turn(new,probe)
        require(ro == rn and _states_equal(new,charted.reframe(old)), 'optional framed phase')
        old, new = native.step(old), charted.step(new)
        for (lo,so),(ln,sn) in zip(native.register_measurement(old,address),charted.register_measurement(new,address)):
            require(lo == ln and _states_equal(sn,charted.reframe(so)), 'complete adaptive history intertwining')
            po,pn = so.norm_squared(),sn.norm_squared()
            require(po == pn, 'exact source-generated history probabilities agree')
            probabilities.append({'preparation':old_label['outcome'],'detector':lo['outcome'],'probability':str(po)})
    require(sum(s.Rational(r['probability']) for r in probabilities) == 1, 'all four recorded outcomes sum to one')

    return {
        'schema':'oph.coherent_seam_agreement.v1',
        'input_sha256':input_hashes(),
        'status':'EXACT finite accepted-chart agreement with executed apparatus frame wrapper',
        'source_program':'coherent-seam-repair-v1',
        'source_dimension':native.dimension,
        'frame_word_at_port0':row['word'],
        'frame_quaternion':row['quaternion'],
        'oriented_local_operator_identities':local_checks,
        'direct_recharted_link_execution_checked':True,
        'complete_adaptive_outcomes':probabilities,
        'wrong_untransported_detector_probability':str(wrong_probability),
        'correct_transported_detector_probability':str(correct_probability),
        'heralded_reference':'rho_Gprime=S rho_G S* on the reduced walk algebra',
        'reference_extension':'operator identities tensor with identity on every untouched reference',
        'retained_archives':'preparation and frame labels are transformed as typed records; no joint stationary-archive claim',
        'declared_service':'controlled finite-cover frame on owned source/code registers, its inverse, and the existing complete detector',
        'physical_clock_or_growing_federation_claim':False,
        'A1_complete_response_certified':False,
        'A3_state_selection_derived':False,
        'Born_rule_derived_from_classical_repair':False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--write',action='store_true')
    action.add_argument('--check',nargs='?',const=str(Path(__file__).with_name('agreement_results.json')))
    args = parser.parse_args()
    content = json.dumps(certificate(),indent=2,sort_keys=True,allow_nan=False)+'\n'
    output = Path(__file__).with_name('agreement_results.json')
    if args.write:
        output.write_text(content)
    elif not Path(args.check).exists() or Path(args.check).read_text()!=content:
        raise SystemExit('finite chart-agreement certificate differs')
    print('PASS complete finite chart agreement with a transported apparatus')


if __name__=='__main__':
    main()
