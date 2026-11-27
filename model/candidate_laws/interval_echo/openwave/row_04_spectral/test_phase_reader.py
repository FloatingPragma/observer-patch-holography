"""Invert the terminal reader independently using the raw forward-source columns."""

import numpy as np
import pytest

from . import raw_verify as raw


def test_replay_rejects_impossible_click_and_altered_refusal():
    from .experiment import compare

    receipt = raw.read_json(raw.HERE / "phase.json")
    row = next(row for row in receipt["boards"][0]["records"] if row["outcomes"][-1] == 0)
    outcomes = [*row["outcomes"][:-1], -1e-12]
    with pytest.raises(ValueError, match="Out-of-range probability"):
        compare(row, {**row, "outcomes": outcomes})
    with pytest.raises(ValueError, match="Numerical replay differs"):
        compare(row, {**row, "source_refused": row["source_refused"] + 1e-12})


def test_every_phase_acceptance_and_original_denominator():
    receipt=raw.read_json(raw.HERE/"phase.json")
    assert len(receipt["record_order"])==8 and receipt["record_order"][-1]=="signed_ray_click"
    for level,board in enumerate(receipt["boards"]):
        tree=raw.component(level)
        terms=raw.spectroscopy()
        U=raw.propagation(level,8,37)[0]
        seen=set()
        targets=dict(terms)
        for row in board["records"]:
            reader,rs,source,arm=row["reader"],row["read_sign"],row["source"],row["arm"]
            seen.add((reader,rs,source,arm))
            loc=[tree.words.index(w) for w in targets[reader]]
            inp=[tree.words.index(w) for w in targets[source]]
            # The unique observed parent supplies the retrodictive forward ket.
            parent=targets[reader][0][:-1]
            flag=raw.board(level)[2][0][1]
            assert [w for (p,flips),w in zip(tree.states,tree.words,strict=True)
                    if p==parent[-1] and (flag in flips)==(1 in parent)]==[parent]
            ray,_=raw.joint_attempt(rs,(1,1))
            backward=np.zeros(len(U))
            backward[loc]=[float(ray[1]),float(ray[2])]
            # This is the forward amplitude including its 1/8 source success.
            response=backward @ U[:,inp]
            if arm=="dephased":
                conditional=float(np.vdot(response,response).real/2)
            else:
                conditional=float(abs(response[0]+(1 if arm=="plus" else -1)*response[1])**2/2)
            mass,fail=raw.ledger(len(raw.board(level)[0]),len(targets[source][0])-2,96,1 if arm=="plus" else -1,.5)
            probs=row["outcomes"]
            assert len(probs)==8 and min(probs)>=0 and abs(sum(probs)+row["source_refused"]-1)<5e-12
            assert probs[-1]==pytest.approx(float(mass)*conditional,abs=2e-12,rel=0)
            assert row["source_refused"]==pytest.approx(float(sum(fail.values())),rel=1e-12,abs=0)
            # Also propagate every complementary branch, using the complete
            # physical (position, flipped coins) states rather than effect matrices.
            if arm=="dephased":
                state=U[:,inp]/np.sqrt(2)
            else:
                state=(U[:,inp[0]]+(1 if arm=="plus" else -1)*U[:,inp[1]])[:,None]/np.sqrt(2)
            keep=np.array([(flag in flips)==(1 in parent) for _,flips in tree.states])
            census=[float(np.linalg.norm(state[~keep])**2)]
            state=state*keep[:,None]
            inc=raw.board(level)[2]
            ids={s:i for i,s in enumerate(tree.states)}
            p=parent[-1]
            for tip,yes in ((None,False),(targets[reader][0][-1],False),(None,False),
                            (targets[reader][1][-1],False),(targets[reader][0][-1],rs==-1)):
                if tip is None:
                    hit=state*(tree.sites==p)[:,None]
                else:
                    hit=np.zeros_like(state)
                    e=inc[p][tip]
                    for i,(where,coins) in enumerate(tree.states):
                        if where not in (p,tip):
                            continue
                        rest=tip if where==p else p
                        j=ids.get((rest,coins^{e}))
                        if j is not None:
                            hit[i]=(state[i]+state[j])/2
                other=state-hit
                state,rejected=(hit,other) if yes else (other,hit)
                census.append(float(np.linalg.norm(rejected)**2))
            census.extend(float(np.linalg.norm(state[v])**2) for v in (tree.sites!=p, tree.sites==p))
            np.testing.assert_allclose(probs,float(mass)*np.array(census),rtol=0,atol=3e-12)
        assert len(seen)==len(board["records"])==24
