"""Seal ideal whole-grid reads, then run the finite emitter and complete local pair observer."""
import argparse
import hashlib
from pathlib import Path

import numpy as np

from tools.strict_json import load
from ..contrast_memory_spherical_v1 import carrier, reader, sectors
from .gate import ROOT, archive, counters, real_step, receipt, require
from .compiler import Compiler

HERE = Path(__file__).resolve().parent


def pins():
    native = ROOT/"model/contrast_memory_spherical_v1"
    files = [HERE/n for n in ("__init__.py","__main__.py","compiler.py","gate.py","protocol.json")]
    files += [native/n for n in ("carrier.py","sectors.py","reader.py","law.py","carriers/FEDERATION_LEVEL0.npz")]
    def digest(p):
        data = p.read_bytes()
        if p.suffix!=".npz":
            data = data.replace(b"\r\n",b"\n")
        return hashlib.sha256(data).hexdigest()
    return {p.relative_to(ROOT).as_posix():digest(p) for p in files}


def ideal_law(grid, compiler, sector, aperture):
    state = np.zeros(sector.n,complex)
    state[compiler.slots] = grid.ravel()
    a,b = compiler.slots[aperture],compiler.slots[0]
    state[[a,b]] = state[[b,a]]
    outcomes = np.zeros(1+len(sector.car.seams))
    for i,(a,b) in enumerate(sector.car.seams,1):
        if compiler.marker in (a,b):
            outcomes[i] = abs(state[a+b-compiler.marker])**2
    outcomes[0] = max(0.,1-outcomes[1:].sum())
    return outcomes


def run(stage,root,record):
    require(pins()==load(HERE/"sources.json"),"changed seam program or gate")
    spec = load(HERE/"protocol.json")
    car = carrier.federation(0)
    sector = sectors.PairSector(car)
    comp = Compiler(car,spec["marker"],spec["shape"])
    trials = spec["waits"]
    queries = spec["query_slots"]
    shots = spec["shots"]
    grid = np.zeros((*comp.shape,8),complex)
    grid.flat[0] = 1
    forecasts = {}
    probabilities = {}
    for tick in range(max(trials)+1):
        if tick in trials:
            for q in queries:
                key = f"{tick}_{q}"
                law = ideal_law(grid,comp,sector,q)
                probabilities[key] = law
                forecasts[key] = counters(law,shots)
        grid = real_step(grid)
    if record and stage=="seal":
        root.mkdir(exist_ok=False,parents=True)
    archive(root/"ideal.npz",forecasts,record and stage=="seal")
    if stage=="seal":
        return
    checked = {}
    tape = {}
    unitary = np.column_stack([real_step(col.reshape(grid.shape)).ravel()
                               for col in np.eye(grid.size,dtype=complex)])
    event = sectors.SeamEvents(sector,spec["kappa"])
    refport,source = carrier.EMITTER_SEAM[0]
    require((refport,source)==(comp.marker,comp.slots[0]), "emitter outside the declared embedding")
    entry = sector.index(refport,source)
    emitted,_ = sectors.emit(sector.hamiltonian(spec["kappa"]),entry,spec["fuel"])
    bound = sectors.emitter_error_bound(*car.degree[[refport,source]],spec["fuel"])
    ideal_input = np.zeros_like(emitted)
    ideal_input[1+entry] = -1j
    require(np.linalg.norm(emitted-ideal_input)<=bound, "finite emitter exceeded its norm bound")
    marker_seam = 1+car.seams.index((refport,source))
    control_words = {label:[int(reader.gate(flag,0,**reader.CONTROLS[label])) for flag in (0,1)]
                     for label in spec["reader_controls"]}
    query_words = {q:comp.exchange(comp.slots[q],source,4) for q in queries}
    def apply(state,program):
        for edge,angle in program:
            event.apply(state,edge,angle*np.pi/4)
    resources = {}
    for member in spec["members"]:
        program = comp.build(member)
        resources[member] = comp.resources(program)
        columns = comp.columns(program)
        target = np.zeros_like(columns)
        target[comp.slots] = unitary
        difference = float(np.linalg.norm(columns-target,ord="fro"))
        finite = emitted[1:].copy()
        max_tv = 0.
        for tick in range(max(trials)+1):
            if tick in trials:
                for q in queries:
                    state = finite.copy()
                    apply(state,query_words[q])
                    require(abs(np.vdot(state,state).real+abs(emitted[0])**2-1)<1e-10,"lost fuel or pair")
                    law = sectors.pair_outcomes(sector,state)
                    require(law.min()>=-1e-12,"negative observer probability")
                    law = np.maximum(law,0)
                    label = f"{tick}_{q}"
                    counts = counters(law,shots)
                    tape[f"{member}_{label}"] = counts
                    for ctl,words in control_words.items():
                        # The same fine seam counts feed both endpoint observers.
                        yes = int(counts[marker_seam])
                        grouped = np.zeros(reader.APP,dtype=np.int32)
                        grouped[words[0]] += shots-yes
                        grouped[words[1]] += yes
                        tape[f"{member}_{label}_{ctl}"] = grouped
                    max_tv = max(max_tv,float(abs(law-probabilities[label]).sum()/2))
            if tick<max(trials):
                apply(finite,program)
        if member=="nominal":
            require(difference<1e-10 and max_tv<=bound,"gate or complete observer law outside forecast")
        else:
            require(difference>1e-3 and max_tv>bound,"changed scheduler failed to distinguish the source")
        # Only material differences are serialized. The nominal threshold is checked above.
        checked[member] = dict(meets_operator_gate=difference<1e-10,max_seam_tv=max_tv)
    archive(root/"opened.npz",tape,record)
    result = dict(programs=resources,readout_programs={str(q):comp.resources(p) for q,p in query_words.items()},
                  comparisons=checked,emitter_norm_error=float(np.linalg.norm(emitted-ideal_input)),
                  emitter_bound=float(bound),unconverted_fuel_probability=float(abs(emitted[0])**2),
                  unmarked_pair_probability=float(np.sum(abs(emitted[1:][(sector.aa!=refport)&(sector.bb!=refport)])**2)),
                  observer_words=control_words,ports=car.n,seams=len(car.seams),pair_dimension=sector.P)
    receipt(root/"result.json",result,record)
    for name in checked:
        print(name,checked[name],resources[name]["events"],resources[name]["depth"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage",choices=("seal","open"),default="open")
    parser.add_argument("--record",type=Path)
    args = parser.parse_args()
    run(args.stage,args.record or HERE/"records",args.record is not None)


if __name__=="__main__":
    main()
