"""Signed message grammar and the honest observer's local signing guards."""

from functools import lru_cache

from cryptography.exceptions import InvalidSignature
from cryptography.hazmat.primitives.asymmetric.ed25519 import (
    Ed25519PrivateKey,
    Ed25519PublicKey,
)

from .transport import LockMemory, encode, identifier

EMPTY = ""


class Signers:
    """Only the acquisition owns private keys. A replay uses the public deposit."""

    def __init__(self, enrollment=None):
        if enrollment is None:
            self.private = [Ed25519PrivateKey.generate() for _ in range(4)]
            self.keys = [s.public_key().public_bytes_raw().hex() for s in self.private]
            self.signatures = {}
        else:
            self.private = None
            self.keys = enrollment["keys"]
            self.signatures = dict(enrollment["signatures"])
        self.verifiers = [
            Ed25519PublicKey.from_public_bytes(bytes.fromhex(k)) for k in self.keys
        ]
        self.used = set()

    def sign(self, stage, i, view, value, ref=EMPTY):
        body = [stage, i, view, value, ref]
        k = identifier(body)
        if self.private is not None:
            self.signatures[k] = self.private[i].sign(encode(body)).hex()
        msg = [*body, self.signatures[k]]
        if not self.valid(msg, stage, view, value, ref):
            raise ValueError("Invalid signature in signing bank")
        self.used.add(k)
        return msg

    @lru_cache(maxsize=1024)
    def verifies(self, i, body, signature):
        try:
            self.verifiers[i].verify(bytes.fromhex(signature), body)
        except (ValueError, InvalidSignature):
            return False
        return True

    def valid(self, msg, stage, view, value, ref):
        return (
            isinstance(msg, list)
            and len(msg) == 6
            and type(msg[1]) is int
            and 0 <= msg[1] < len(self.keys)
            and encode(msg[:5]) == encode([stage, msg[1], view, value, ref])
            and isinstance(msg[5], str)
            and self.verifies(msg[1], encode(msg[:5]), msg[5])
        )

    def archive(self):
        if self.used != self.signatures.keys():
            raise ValueError("Unused signed objects in deposit")
        return {"keys": self.keys, "signatures": self.signatures}


class Rules:
    def __init__(self, signers, threshold):
        self.auth = signers
        self.q = threshold

    def quorum(self, votes, phase, view, value, ref):
        return (
            isinstance(votes, list)
            and len(votes) == self.q
            and all(self.auth.valid(m, phase, view, value, ref) for m in votes)
            and len({m[1] for m in votes}) == self.q
        )

    def cert(self, c, kind):
        if not isinstance(c, list) or len(c) != 5:
            return False
        tag, v, x, parent, votes = c
        if (
            tag != kind
            or type(v) is not int
            or v < 0
            or type(x) is not int
            or x not in (0, 1)
        ):
            return False
        if kind == "PC":
            return parent is None and self.quorum(votes, "prepare", v, x, EMPTY)
        previous, phase = {"LC": ("PC", "accept"), "DC": ("LC", "commit")}[kind]
        return (
            self.cert(parent, previous)
            and parent[1:3] == [v, x]
            and self.quorum(votes, phase, v, x, identifier(parent))
        )

    def proposal(self, packet):
        if not isinstance(packet, list) or len(packet) != 2:
            return False
        leader, reports = packet
        if (
            not isinstance(leader, list)
            or len(leader) != 6
            or not isinstance(reports, list)
        ):
            return False
        _, i, view, x, _, _ = leader
        if type(view) is not int or view < 0 or type(x) is not int or x not in (0, 1):
            return False
        if i != view % 4 or not self.auth.valid(
            leader, "propose", view, x, identifier(reports)
        ):
            return False
        valid, forced = self.selection(view, reports)
        return valid and forced in (None, x)

    def report_signer(self, view, report):
        if not isinstance(report, list) or len(report) != 2:
            return None
        mark, lock = report
        if lock is None:
            value, ref = -1, EMPTY
        else:
            if not self.cert(lock, "LC") or lock[1] >= view:
                return None
            value, ref = lock[2], identifier(lock)
        if not self.auth.valid(mark, "change", view, value, ref):
            return None
        return mark[1]

    def first_quorum(self, messages, signer):
        accepted = {}
        for message in messages:
            identity = signer(message)
            if identity is not None and identity not in accepted:
                accepted[identity] = message
                if len(accepted) == self.q:
                    break
        return list(accepted.values())

    def selection(self, view, reports):
        if len(reports) != self.q:
            return False, None
        signers = {self.report_signer(view, report) for report in reports}
        if None in signers or len(signers) != self.q:
            return False, None
        locks = [report[1] for report in reports if report[1] is not None]
        if not locks:
            return True, None
        top = max(lock[1] for lock in locks)
        values = {lock[2] for lock in locks if lock[1] == top}
        return (True, values.pop()) if len(values) == 1 else (False, None)


class Observer:
    def __init__(self, i, rules, board):
        self.i, self.rules = i, rules
        self.view = -1
        self.value = None
        self.memory = LockMemory(i, board)
        self.decision = None
        self.decision_certificate = None

    @property
    def lock(self):
        return self.memory.read(self.i)[1]

    def change(self, view):
        if view < self.view:
            return None
        if view > self.view:
            self.view, self.value = view, None
        lock = self.lock
        x, ref = (-1, EMPTY) if lock is None else (lock[2], identifier(lock))
        return [self.rules.auth.sign("change", self.i, view, x, ref), lock]

    def respond(self, stage, request):
        if stage == "prepare":
            if not self.rules.proposal(request):
                return None
            view, x = request[0][2:4]
            ref = EMPTY
        else:
            kind = {"accept": "PC", "commit": "LC"}[stage]
            if not self.rules.cert(request, kind):
                return None
            view, x = request[1:3]
            ref = identifier(request)
        if self.view != view or self.value not in (None, x):
            return None
        self.value = x
        if stage == "commit":
            self.memory.install(request)
        return self.rules.auth.sign(stage, self.i, view, x, ref)

    def finalize(self, certificate):
        if not self.rules.cert(certificate, "DC"):
            return False
        x = certificate[2]
        if self.decision not in (None, x):
            return False
        self.decision = x
        self.decision_certificate = certificate
        return True
