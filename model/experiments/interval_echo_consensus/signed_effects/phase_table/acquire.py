"""Replay native complex instruments, their coarse records and state readback."""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import json

import numpy as np
import sympy as s

from . import measure as native
from .. import fixed_pointer

HERE = native.HERE
ROOT = native.ROOT
RECEIPT = HERE/"results.json"


@lru_cache(maxsize=None)
def raw(context,mode,sign):
    return native.raw_branch_maps(context,mode,sign)


def encode_matrix(matrix):
    return {"real":matrix.real.tolist(),"imag":matrix.imag.tolist()}


def preparation_word(name):
    p = native.pulse
    return {"port0":[],"port1":[p("contact",0,1,1)],
            "plus":[p("contact",0,1,Fraction(1,2)),p("reference",0,0,Fraction(1,16))],
            "plus_i":[p("contact",0,1,Fraction(1,2))],
            "minus_i":[p("contact",0,1,Fraction(3,2))]}[name]


def input_vector(name,sign):
    return native.source_unitary(tuple(preparation_word(name)),sign)[:2,0]


def exact_input(name):
    return {"port0":s.Matrix([1,0]),"port1":s.Matrix([0,1]),
            "plus":s.Matrix([1,1])/s.sqrt(2),"plus_i":s.Matrix([1,s.I])/s.sqrt(2),
            "minus_i":s.Matrix([1,-s.I])/s.sqrt(2)}[name]


@lru_cache(maxsize=None)
def exact_projector(context):
    row = next(row for row in native.context_table() if row["context"]==context)
    return s.Matrix([[s.Rational(v[0])+s.Rational(v[1])*s.sqrt(3)
                      +s.I*(s.Rational(v[2])+s.Rational(v[3])*s.sqrt(3)) for v in line]
                     for line in row["effect"]])


@lru_cache(maxsize=None)
def exact_weight(context,preparation):
    vec = exact_input(preparation)
    return s.simplify((vec.conjugate().T*exact_projector(context)*vec)[0])


def ideal_branches(context,mode):
    """Independent analytic maps, including physical outcome phases."""
    p = native.projector(context)
    q = np.eye(2)-p
    h = native.address.COMMON_PHASE
    if mode=="sharp":
        return {0:h**14*p,2:h**14*q}
    if mode=="coin":
        phase = h**14*np.exp(-17j*np.pi/6)
        return {0:phase*np.sqrt(3)/2*np.eye(2),2:1j*phase/2*np.eye(2)}
    phase = h**(20 if mode=="fine" else 14)*np.exp(1j*np.pi/3)
    if mode=="coarse":
        return {0:phase*(np.sqrt(3)/2*p+q/2),
                2:1j*phase*(p/2+np.sqrt(3)/2*q)}
    return {0:phase*np.sqrt(3)/2*p,1:phase/2*q,
            2:1j*phase/2*p,3:1j*phase*np.sqrt(3)/2*q}


def reset_word(raw_record):
    return [native.pulse("contact",site["park"],site["cross"],1)
            for bit,site in enumerate(native.address.SITES) if (raw_record>>bit)&1]


def x_analyzer_word():
    # Common real rotation with first row (1,1)/sqrt(2).
    return [native.pulse("reference",0,0,Fraction(1,16)),
            native.pulse("contact",0,1,Fraction(1,2)),
            native.pulse("reference",0,0,Fraction(3,16))]


def native_readback_word(raw_record):
    return reset_word(raw_record)+x_analyzer_word()+native.address.pulse_word(native.address.schedule())


def readback_probabilities(context,mode,preparation,sign):
    source = input_vector(preparation,sign)
    analyze = native.source_unitary(tuple(x_analyzer_word()),sign)
    probabilities = np.zeros((17,17))
    for label,k in raw(context,mode,sign).items():
        after = analyze@k@source
        probabilities[label,:12] = abs(after)**2
    return probabilities


def cleared_record_check(mode):
    """Execute the unread-bit inverse on every physical routed code word."""
    calls = native.address.schedule()
    unread = [call for call in calls if not (native.keep_mask(mode)>>call["bit"])&1]
    max_phase_error = 0.0
    for sign in (-1,1):
        for address in range(12):
            for prior in range(16):
                loads,phase = native.address.execute_basis(calls,sign,address,prior)
                for call in reversed(unread):
                    edges = list(zip(call["route"],call["route"][1:]))
                    for a,b in edges:
                        loads[a],loads[b] = loads[b],loads[a]
                    site = native.address.SITES[call["site"]]
                    if loads[site["source"]]!=1:
                        a,b = site["park"],site["cross"]
                        loads[a],loads[b] = loads[b],loads[a]
                    phase *= native.address.COMMON_PHASE.conjugate()
                    for a,b in reversed(edges):
                        loads[a],loads[b] = loads[b],loads[a]
                expected = native.address.initial_word(sign,address,prior^(address&native.keep_mask(mode)))
                native.sequential.require(np.array_equal(loads,expected),"unread record failed to reset")
                max_phase_error = max(max_phase_error,abs(phase-native.address.COMMON_PHASE**(20-len(unread))))
    return {"physical_code_inputs":384,"cleared_bits":[bit for bit in range(4)
             if not (native.keep_mask(mode)>>bit)&1],"retained_native_copy_calls":20-len(unread),
            "phase_comparison_error":float(max_phase_error)}


def source_hashes():
    names = ["code/phase_operation_producer/PHASE_OPERATION_RECEIPT.v1.json",
             "model/candidate_laws/contrast_memory_spherical_v1/law.py",
             "model/candidate_laws/interval_echo/candidate.py",
             "model/executable_model/carriers/FEDERATION_LEVEL0.npz",
             "model/experiments/interval_echo_consensus/signed_effects/full_address/measure.py",
             "model/experiments/interval_echo_consensus/signed_effects/phase_table/measure.py",
             "model/experiments/interval_echo_consensus/signed_effects/phase_table/acquire.py"]
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}


def serialized(word):
    return [[kind,a,b,str(t)] for kind,a,b,t in word]


def sample_records(context,mode,preparation,shots,rng):
    p = float(exact_weight(context,preparation))
    counts = np.zeros(17,dtype=int)
    if mode=="fine":
        j0 = int(rng.binomial(shots,p))
        y0j0 = int(rng.binomial(j0,0.75))
        y0j1 = int(rng.binomial(shots-j0,0.25))
        counts[:4] = [y0j0,y0j1,j0-y0j0,shots-j0-y0j1]
    else:
        yes = p if mode=="sharp" else (0.75 if mode=="coin" else 0.25+p/2)
        counts[0] = int(rng.binomial(shots,yes))
        counts[2] = shots-counts[0]
    return counts.tolist()


def ideal_record_probabilities(context,mode,preparation):
    p = float(exact_weight(context,preparation))
    result = np.zeros(17)
    if mode=="fine":
        result[:4] = [0.75*p,0.25*(1-p),0.25*p,0.75*(1-p)]
    else:
        result[0] = p if mode=="sharp" else (0.75 if mode=="coin" else 0.25+p/2)
        result[2] = 1-result[0]
    return result


def readback_banks(rng):
    outputs = []
    shots = 32768
    for sign in (-1,1):
        for mode in ("coarse","fine","coin"):
            law = readback_probabilities("phase",mode,"plus",sign)
            counts = np.zeros((17,17),dtype=int)
            if mode=="fine":
                first = rng.multinomial(shots,[0.375,0.125,0.125,0.375])
                final0 = 0.5
                rows = range(4)
            else:
                yes = 0.75 if mode=="coin" else 0.5
                first0 = int(rng.binomial(shots,yes))
                first = [first0,shots-first0]
                final0 = 1.0 if mode=="coin" else (2+np.sqrt(3))/4
                rows = (0,2)
            for row,total in zip(rows,first,strict=True):
                hit = int(rng.binomial(int(total),final0))
                counts[row,:2] = [hit,int(total)-hit]
            ideal = np.zeros((17,17))
            first_law = ([0.375,0.125,0.125,0.375] if mode=="fine"
                         else ([0.75,0.25] if mode=="coin" else [0.5,0.5]))
            for row,probability in zip(rows,first_law,strict=True):
                ideal[row,:2] = [probability*final0,probability*(1-final0)]
            native.sequential.require(np.max(abs(law-ideal))<4e-12,
                                      "conditional state-readback branches differ")
            expected = (2+np.sqrt(3))/4 if mode=="coarse" else (0.5 if mode=="fine" else 1)
            native.sequential.require(abs(law[:,0].sum()-expected)<4e-12,
                                      "native state-readback contrast failed")
            outputs.append({"sign":sign,"mode":mode,"input":"plus","context":"phase",
                            "all_17_by_17_probability_centers":law.tolist(),
                            "all_17_by_17_counts":counts.tolist(),"shots":shots})
    return outputs


def report():
    rng = np.random.default_rng(625011)
    shots = 4096
    banks,tables = [],[]
    max_map_error = max_sign_error = max_union_error = max_input_error = 0.0
    for row in native.context_table():
        context = row["context"]
        entry = {"context":context,"frame_word":serialized(native.context_word(context)),
                 "sharp_target":row["effect"],"branch_matrices":{}}
        modes = ("sharp","coarse","fine","coin") if context=="phase" else ("sharp","coarse","fine")
        for mode in modes:
            expected = ideal_branches(context,mode)
            minus,plus = raw(context,mode,-1),raw(context,mode,1)
            entry["branch_matrices"][mode] = {str(label):encode_matrix(k[:2])
                                                for label,k in minus.items() if label in expected}
            for sign,ops in ((-1,minus),(1,plus)):
                native.sequential.require(np.linalg.norm(sum(k.conj().T@k for k in ops.values())-np.eye(2),2)<4e-12,
                                          "incomplete complex instrument")
                for label,k in ops.items():
                    target = expected.get(label,np.zeros((2,2)))
                    max_map_error = max(max_map_error,float(np.linalg.norm(k[:2]-target,2)),float(np.linalg.norm(k[2:],2)))
                for preparation in ("port0","port1","plus","plus_i","minus_i"):
                    state = input_vector(preparation,sign)
                    wanted = np.array(exact_input(preparation).evalf(),complex).ravel()
                    max_input_error = max(max_input_error,float(np.linalg.norm(
                        np.outer(state,state.conj())-np.outer(wanted,wanted.conj()),2)))
                    centers = [float(np.linalg.norm(k@state)**2) for k in ops.values()]
                    native.sequential.require(np.max(abs(np.array(centers)-ideal_record_probabilities(
                        context,mode,preparation)))<4e-12,"native record law differs from sampled ideal law")
                    counts = sample_records(context,mode,preparation,shots,rng)
                    banks.append({"context":context,"mode":mode,"sign":sign,"input":preparation,
                                  "all_17_probability_centers":centers,"all_17_counts":counts,"shots":shots})
            max_sign_error = max(max_sign_error,max(float(np.linalg.norm(minus[j]-plus[j],2)) for j in range(17)))
        fine = raw(context,"fine",-1)
        effects = {j:k.conj().T@k for j,k in fine.items()}
        coarse = raw(context,"coarse",-1)
        max_union_error = max(max_union_error,float(np.linalg.norm(effects[0]+effects[1]-coarse[0].conj().T@coarse[0],2)),
                              float(np.linalg.norm(effects[0]+effects[3]-3*np.eye(2)/4,2)))
        tables.append(entry)
    native.sequential.require(max_map_error<4e-12 and max_sign_error<4e-12 and max_union_error<4e-12
                              and max_input_error<4e-12,
                              "native maps, sign coherence or coarse effects failed")
    budgets = []
    for context in native.context_table():
        for mode in ("sharp","coarse","fine","coin"):
            prefix = native.context_word(context["context"])+native.dilation_word(mode)+native.record_word(mode)
            accepted = range(4) if mode=="fine" else (0,2)
            for label in accepted:
                word = prefix+native.post_word(context["context"],mode,label)
                budgets.append({"context":context["context"],"mode":mode,"record":label,
                                "waits":len(word),"generator_norm_sum":sum(33 if k=="contact" else 32 for k,_,_,_ in word)})
    sharp_words = {context["context"]:{str(raw):serialized(native.post_word(context["context"],"sharp",raw))
                                      for raw in (0,2)} for context in native.context_table()}
    return {
        "candidate":"Interval Echo protected native control menu",
        "input":"one unknown-sign defect on source rails0,1; fixed neutral pointer architecture",
        "contexts":"eight committed slots, four distinct sharp projectors; plus unsharp E=I/4+P/2",
        "basis_fixture":"declared table is a target; this construction and its chosen weakness do not derive PR04 or Born probabilities",
        "context_table_sha256":hashlib.sha256(native.TABLE.read_bytes()).hexdigest(),
        "source_sha256":source_hashes(),"tables":tables,
        "native_dilation_words":{mode:serialized(native.dilation_word(mode)) for mode in ("sharp","coarse","fine","coin")},
        "shared_record_words":{mode:serialized(native.record_word(mode)) for mode in ("coarse","fine")},
        "weak_common_phase_correction":serialized(fixed_pointer.diagonal_word(0,7,Fraction(5,51))),
        "outcome1_decode_word":serialized(native.transposition(0,2)+native.transposition(1,3)),
        "sharp_post_words":sharp_words,
        "readback_analyzer_word":serialized(x_analyzer_word()),
        "readback_reset_words":{str(r):serialized(reset_word(r)) for r in range(4)},
        "record_alphabet":"all17 raw words; coarse accepted0/2, fine0..3, others rejection/unused",
        "record_reset":{mode:cleared_record_check(mode) for mode in ("coarse","fine")},
        "maximum_native_vs_analytic_branch_operator_error":max_map_error,
        "maximum_same_word_sign_block_difference":max_sign_error,
        "maximum_alternative_union_effect_error":max_union_error,
        "maximum_prepared_ray_density_error":max_input_error,
        "complete_native_history_budgets":budgets,
        "duration_error_contract":"complete cq half-trace error <= min(1,epsilon_prep+B*delta) "
            "with B the maximum recorded generator norm sum for the protocol and absolute wait error delta; "
            "final sharp read assumed",
        "banks":banks,"readback_banks":readback_banks(rng),"seed":625011,
        "count_semantics":"synthetic exact-ideal-law tree draws checked against unnormalized native propagation centers; no literal native event count",
        "postselection":"none; rejection/unused raw outcomes and all second-read outcomes retained",
        "custody":"producer-side software binding; final sharp occupation rule and clock/control calibration supplied",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true")
    args = parser.parse_args()
    result = report()
    if args.write:
        with RECEIPT.open("x") as output:
            output.write(json.dumps(result,indent=2,allow_nan=False)+"\n")
    else:
        native.sequential.compare(json.loads(RECEIPT.read_text()),result)
        print("native phase table: sharp/unsharp maps, alternative records and state readback replay")


if __name__=="__main__":
    main()
