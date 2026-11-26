"""Execute the prepared-lock protocol against the received material packet bytes."""

import argparse
from copy import deepcopy
import gzip
from itertools import combinations
import json
from pathlib import Path

from .certificates import EMPTY, Observer, Rules, Signers
from .transport import Wire, encode, identifier, layout


def must(condition, reason):
    if not condition:
        raise ValueError(reason)


def lock_at(node):
    lock = node.lock
    return None if lock is None else [*lock[1:3], identifier(lock)]


class Trial:
    def __init__(self, auth, q=3, faulty=()):
        self.rules = Rules(auth, q)
        board = layout().fresh([0] * 6)
        self.nodes = [Observer(i, self.rules, board) for i in range(6)]
        self.wire = Wire()
        self.log = []
        self.faulty = set(faulty)

    def reports(self, view, members, enter=None):
        leader = view % 4
        bank = {}
        for i in members if enter is None else enter:
            report = self.nodes[i].change(view)
            if i in members:
                bank[i] = self.wire.send(i, leader, report)
        return [bank[i] for i in members]

    def propose(self, view, fresh, reports):
        arrivals = len(reports)
        reports = self.rules.first_quorum(
            reports, lambda m: self.rules.report_signer(view, m)
        )
        valid, forced = self.rules.selection(view, reports)
        self.log.append(["proposal_selection", view, arrivals, valid, forced, fresh])
        if not valid:
            return None
        return self.sign_offer(view, fresh if forced is None else forced, reports)

    def sign_offer(self, view, x, reports):
        m = self.rules.auth.sign("propose", view % 4, view, x, identifier(reports))
        return [m, reports]

    def vote(self, stage, request, i, leader):
        observed = self.wire.send(leader, i, request)
        node = self.nodes[i]
        before = lock_at(node)
        if i in self.faulty:
            if stage == "prepare":
                view, x = observed[0][2:4]
                ref = EMPTY
            else:
                view, x, ref = *observed[1:3], identifier(observed)
            vote = self.rules.auth.sign(stage, i, view, x, ref)
        else:
            vote = node.respond(stage, observed)
        self.log.append([stage, i, node.view, before, lock_at(node), vote is not None])
        return self.wire.send(i, leader, vote) if vote is not None else None

    def collect(self, stage, packet, members, leader):
        votes = [self.vote(stage, packet, i, leader) for i in members]
        if stage == "prepare":
            kind, parent, (v, x) = "PC", None, packet[0][2:4]
        else:
            kind, parent, (v, x) = (
                {"accept": "LC", "commit": "DC"}[stage],
                packet,
                packet[1:3],
            )
        ref = EMPTY if parent is None else identifier(parent)

        def signer(m):
            return m[1] if self.rules.auth.valid(m, stage, v, x, ref) else None

        bundle = [kind, v, x, parent, self.rules.first_quorum(votes, signer)]
        return bundle if self.rules.cert(bundle, kind) else None

    def prepare_lock(self, offer, prepares, accepts):
        leader = offer[0][1]
        pc = self.collect("prepare", offer, prepares, leader)
        must(pc is not None, "Expected a raw prepare certificate")
        lc = self.collect("accept", pc, accepts, leader)
        must(lc is not None, "Expected an assembled lock certificate")
        return pc, lc

    def decision(self, sender, receiver, cert, parked=False):
        msg = self.wire.send(sender, receiver, cert, parked)
        outcome = self.nodes[receiver].finalize(msg)
        self.log.append(["decision", receiver, outcome, self.nodes[receiver].decision])
        return outcome

    def issue(self, view, x, members):
        offer = self.propose(view, x, self.reports(view, members))
        _, lock = self.prepare_lock(offer, members, members)
        return self.collect("commit", lock, members, view % 4)

    def archive(self):
        for n in self.nodes[:4]:
            must(n.memory.read(5)[1] == n.lock, "Independent lock reader disagreed")
        return {
            "packets": self.wire.packets,
            "links": self.wire.links,
            "wire": self.wire.events,
            "observer_events": self.log,
            "final": [n.decision for n in self.nodes],
            "locks": [lock_at(n) for n in self.nodes],
            "lock_reads": [n.memory.events for n in self.nodes[:4]],
        }


def orphan_and_decision(auth):
    t = Trial(auth)
    everyone = [0, 1, 2, 3]
    offer = t.propose(0, 0, t.reports(0, [0, 1, 2], everyone))
    pc, first = t.prepare_lock(offer, [0, 1, 2], [1, 2, 3])
    must(
        all(n.lock is None for n in t.nodes), "Provisional acceptance installed a lock"
    )
    must(t.vote("commit", pc, 1, 0) is None, "Raw prepare admitted for commit")
    must(t.collect("commit", first, [1], 0) is None, "One commit became a decision")
    held = lock_at(t.nodes[1])
    must(held[:2] == [0, 0], "Orphan lock missing")
    reports = t.reports(1, [0, 2, 3], everyone)
    must(all(r[1] is None for r in reports), "Held quorum contained a lock")
    fresh = t.propose(1, 1, reports)
    pc1, lock = t.prepare_lock(fresh, [1, 2, 3], [0, 1, 3])
    must(
        lock_at(t.nodes[1]) == held,
        "Prepare or accept prematurely replaced orphan lock",
    )
    for stage, old_request in [("prepare", offer), ("accept", pc), ("commit", first)]:
        must(
            t.vote(stage, old_request, 1, 0) is None,
            "Old request regained signing authority",
        )
    must(t.nodes[1].change(0) is None, "Local view decreased on delayed change")
    dc = t.collect("commit", lock, [1, 2, 3], 1)
    must(dc is not None, "Missing view-one decision")
    must(
        t.decision(1, 4, dc) and t.decision(4, 5, t.nodes[4].decision_certificate),
        "Independent decision or relay failed",
    )
    must(not t.decision(1, 5, dc, parked=True), "Parked channel unexpectedly decoded")
    next_reports = t.reports(2, [0, 1, 2], everyone)
    must(
        t.vote("prepare", t.sign_offer(2, 0, next_reports), 3, 2) is None,
        "Decided lock omitted from selection",
    )
    must(
        t.vote("prepare", t.propose(2, 0, next_reports), 3, 2) is not None,
        "Forced continuation refused",
    )
    duplicate_report = [next_reports[1], next_reports[1], next_reports[2]]
    must(
        t.vote("prepare", t.sign_offer(2, 1, duplicate_report), 3, 2) is None,
        "Duplicate justifier accepted",
    )
    forged_report = deepcopy(next_reports)
    forged_report[1][1] = pc1
    must(
        t.vote("prepare", t.sign_offer(2, 1, forged_report), 3, 2) is None,
        "Raw PC justified a view",
    )
    alternatives = {}
    for name in (
        "duplicate",
        "short",
        "bad_signature",
        "different_lock",
        "wrong_value",
    ):
        wrong = deepcopy(dc)
        if name == "duplicate":
            wrong[4][1] = wrong[4][0]
        elif name == "short":
            wrong[4].pop()
        elif name == "bad_signature":
            sig = wrong[4][0][5]
            wrong[4][0][5] = ("0" if sig[0] != "0" else "1") + sig[1:]
        elif name == "different_lock":
            wrong[3][4].reverse()
            must(t.rules.cert(wrong[3], "LC"), "Changed parent should itself be valid")
        else:
            wrong[2] ^= 1
        alternatives[name] = [t.decision(1, j, wrong) for j in (4, 5)]
    must(
        not any(any(a) for a in alternatives.values()),
        "Invalid decision certificate accepted",
    )
    return t.archive()


def mixed_locks(auth):
    t = Trial(auth)
    members = [0, 1, 2, 3]
    p = t.propose(0, 0, t.reports(0, [1, 2, 3], members))
    _, lc = t.prepare_lock(p, [0, 1, 3], [0, 2, 3])
    t.vote("commit", lc, 0, 0)
    p = t.propose(1, 1, t.reports(1, [1, 2, 3], members))
    _, lc = t.prepare_lock(p, [1, 2, 3], [0, 1, 2])
    t.vote("commit", lc, 2, 1)
    reports = t.reports(2, [0, 2, 3], members)
    must(
        [r[1][1] if r[1] else None for r in reports] == [0, 1, None],
        "Highest-lock fixture changed",
    )
    must(
        t.vote("prepare", t.sign_offer(2, 0, reports), 3, 2) is None,
        "Oldest reported lock authorized",
    )
    p = t.propose(2, 0, reports)
    _, lc = t.prepare_lock(p, [0, 1, 2], [0, 2, 3])
    dc = t.collect("commit", lc, [0, 1, 3], 2)
    must(t.decision(2, 4, dc) and t.decision(2, 5, dc), "Highest-lock decision failed")
    return t.archive()


def silent(auth, q):
    t = Trial(auth, q)
    # All nonfaulty owners enter the silent leader's view and advance on timeout.
    t.reports(0, [1, 2, 3])
    available = [1, 2, 3]
    changes = t.reports(1, available)
    offer = t.propose(1, 1, changes)
    if q == 4:
        must(
            offer is None,
            "Silent fourth signer gave progress",
        )
    else:
        _, lc = t.prepare_lock(offer, available, available)
        dc = t.collect("commit", lc, available, 1)
        must(
            all(t.decision(1, i, dc) for i in [1, 2, 3, 4, 5]), "Decision relay failed"
        )
    return t.archive()


def equivocation(auth):
    t = Trial(auth, 3, faulty=[0])
    reports = t.reports(0, [0, 1, 2], [0, 1, 2, 3])
    x = t.sign_offer(0, 0, reports)
    y = t.sign_offer(0, 1, reports)
    pc, lc = t.prepare_lock(x, [0, 1, 2], [0, 1, 2])
    must(t.vote("prepare", y, 1, 0) is None, "Nonfaulty owner cast a second vote")
    must(
        t.collect("prepare", y, [0, 1, 3], 0) is None, "Two conflicting prepare quorums"
    )
    dc = t.collect("commit", lc, [0, 1, 2], 0)
    must(t.decision(0, 4, dc) and t.decision(0, 5, dc), "Nonfaulty decision failed")
    return t.archive()


def conflicting(auth):
    t = Trial(auth, 2, faulty=[0])
    for x, nonfaulty, final in [(0, 1, 4), (1, 2, 5)]:
        dc = t.issue(0, x, [0, nonfaulty])
        must(
            t.decision(0, final, dc),
            "Low threshold did not expose conflicting certificate",
        )
    must([n.decision for n in t.nodes[4:]] == [0, 1], "Expected incompatible decisions")
    return t.archive()


def intersection_counts():
    result = {}
    for q in [2, 3, 4]:
        qs = list(map(set, combinations(range(4), q)))
        sizes = [len(a & b - {bad}) for a in qs for b in qs for bad in range(4)]
        result[str(q)] = {
            "cases": len(sizes),
            "without_honest_witness": sizes.count(0),
            "minimum_honest_overlap": min(sizes),
        }
    return result


def run(enrollment=None):
    studies = {}
    auth_bank = {}
    for name, exercise in [
        ("orphan", orphan_and_decision),
        ("equivocation", equivocation),
        ("mixed", mixed_locks),
        ("silent", lambda a: silent(a, 3)),
        ("high", lambda a: silent(a, 4)),
        ("low", conflicting),
    ]:
        auth = Signers(enrollment[name] if enrollment is not None else None)
        studies[name] = exercise(auth)
        auth_bank[name] = auth.archive()
    return {
        "schema": 1,
        "encoding": "Canonical ASCII JSON, most significant bit first; each wire row records both pointers via its two-entry link table.",
        "studies": studies,
        "intersections": intersection_counts(),
        "authentication": auth_bank,
    }


def main():
    path = Path(__file__).with_name("records.json.gz")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", action="store_true")
    args = parser.parse_args()
    if args.acquire:
        result = run()
        with path.open("xb") as f:
            f.write(gzip.compress(encode(result), mtime=0))
    else:
        saved = json.loads(gzip.decompress(path.read_bytes()))
        must(
            encode(run(saved["authentication"])) == encode(saved),
            "Public replay differed from deposit",
        )
    print(
        "Verified six prepared-lock traces, signed records and both independent decision readers."
    )


if __name__ == "__main__":
    main()
