"""Release the clocked preparation, then collect every native marker-occupation outcome."""
import argparse
from hashlib import sha256
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from tools.strict_json import load
from model.candidate_laws.contrast_memory_spherical_v1 import sectors
from model.executable_model import carrier
from .clocked_apparatus import clock_laws, patches
from .clocked_pins import verify_sources
from .clocked_run import finite_source
from .compiler import Compiler
from .gate import archive, counters, receipt, require
from .preparation import HERE, execute_source, preparation_pins, rotation


def leakage(car, marker):
    """Assemble all transitions that leave the marked subspace, reducing only empty rows."""
    partners = [v for v in range(car.n) if v != marker]
    columns = {v:i for i,v in enumerate(partners)}
    rows, row_ids, col_ids = {}, [], []
    nb = car.neighbours[marker]
    for end in nb:
        for other in partners:
            if other != end:
                row_ids.append(rows.setdefault(tuple(sorted((end,other))),len(rows)))
                col_ids.append(columns[other])
    b = sparse.coo_matrix((np.full(len(row_ids),-.5),(row_ids,col_ids)),shape=(len(rows),len(partners))).tocsr()
    d = len(nb)
    diagonal = np.full(len(partners),d/4)
    near = np.array([columns[v] for v in nb])
    diagonal[near] = (d-2)/4
    j = sparse.coo_matrix((np.full(d*d,.25),(np.repeat(near,d),np.tile(near,d))),shape=(len(partners),)*2)
    residual = (b.T@b-sparse.diags(diagonal)-j).tocsr()
    residual.eliminate_zeros()
    require(d>=3 and residual.nnz==0, "leakage identity failed")
    return dict(n=car.n,d=d,nonzero_output_rows=len(rows),
                squared_singular_values=[[(d-2)/4,d-1],[d/4,len(partners)-d],[(d-1)/2,1]])


def hamiltonian(data, host, ports):
    sector = sectors.PairSector(data)
    lost = host.degree[ports]-data.degree
    return sector.hamiltonian(8)-sparse.diags(8*(lost[sector.aa]+lost[sector.bb]),dtype=float)


def terminal_words(comp):
    return [rotation(comp,p,comp.slots[0],np.pi) for p in range(comp.carrier.n) if p!=comp.marker]


def routed_read(compiler, events, initial):
    """Literal native first-yes detector, used to check the diagonal instrument reduction."""
    state = initial.copy()
    test = events.sector.index(compiler.marker,compiler.slots[0])
    yes = []
    for route in terminal_words(compiler):
        for edge,angle in route:
            events.apply(state,edge,angle)
        yes.append(abs(state[test])**2)
        state[test] = 0
        for edge,angle in reversed(route):
            events.apply(state,edge,angle)
    return np.asarray(yes), state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record",type=Path)
    args = parser.parse_args()
    verify_sources()
    require(preparation_pins()==load(HERE/"preparation_sources.json"),"source changed")
    manifest = load(HERE/"free_marker_sources.json")
    files = ("free_marker.py","free_marker_protocol.json")
    require(manifest == {n:sha256((HERE/n).read_bytes().replace(b"\r\n",b"\n")).hexdigest() for n in files},
            "free release differs from frozen sources")
    spec = load(HERE/"clocked_protocol.json")["fixture"]
    protocol = load(HERE/"free_marker_protocol.json")
    require(protocol["exposure_per_clock_opportunity"]=="1/8 in native Hamiltonian units",
            "unknown exposure unit")
    cap = spec["alarm_after_periods"]
    host,data,ports,front,bank = patches(spec["level"])
    clock,_ = clock_laws(host,front,bank,spec)
    sec,emitted,emitter = finite_source(data,host,ports,spec["fuel"])
    comp = Compiler(data,0,[1,1,1,spec["shape"][3]])
    neg = list(range(data.n-len(comp.slots),data.n))
    ev = sectors.SeamEvents(sec,8)
    marked = np.array([sec.index(0,p) for p in range(1,data.n)])
    off = sec.aa!=0
    cut = carrier.Carrier("marker held",data.n,[pair for pair in data.seams if 0 not in pair])
    operators = {"all_data_contacts":hamiltonian(data,host,ports),
                 "marker_contacts_closed":hamiltonian(cut,host,ports)}
    require(operators["marker_contacts_closed"][marked][:,off].nnz==0, "marker contact isolation failed")
    graphs = [carrier.cell(),*(carrier.federation(level) for level in (0,1,2))]
    summary = dict(leakage=[leakage(car,0) for car in [*graphs,data]],
                   terminal_queries=len(marked),terminal_worst_pi_contacts=2*sum(map(len,terminal_words(comp))),
                   all_no_excess_over_hold={})
    opened = {}
    shots = spec["shots"]
    for order in spec["private_orders"]:
        wave = np.repeat(np.exp(2j*np.pi*order*np.arange(comp.shape[3])/comp.shape[3]),8)/np.sqrt(len(comp.slots))
        launch,early,_ = execute_source(comp,ev,neg,wave,"native",emitted[1:])
        live = abs(emitted[0])**2 + np.vdot(launch,launch).real
        quiet = {}
        for member in protocol["members"]:
            if member == "hold":
                states = np.broadcast_to(launch,(cap+1,sec.P))
            else:
                states = expm_multiply(-1j*operators[member],launch,start=0,stop=cap/8,num=cap+1)
            found = abs(states[:,marked])**2
            none = abs(emitted[0])**2 + (abs(states[:,off])**2).sum(axis=1)
            require(np.max(abs(none+found.sum(axis=1)-live))<1e-11, "free terminal alphabet incomplete")
            for r in protocol["thresholds"]:
                w = clock[r]
                law = np.r_[early,live*(1-w.sum()),w@found,w@none]
                counts = counters(law,shots)
                opened[f"{order}_{member}_{r}"] = counts
                quiet[member,r] = counts[-1]/shots
        excess = {m:[float(quiet[m,r]-quiet["hold",r]) for r in protocol["thresholds"]]
                  for m in operators}
        require(max(map(abs,excess["marker_contacts_closed"]))<=2/shots, "closed marker leaked")
        summary["all_no_excess_over_hold"][str(order)] = excess
    summary["two_emitter_and_counter_bound"] = 2*(emitter["bound"]+1/shots)
    root = args.record or HERE/"free_marker_records"
    if args.record:
        root.mkdir(parents=True,exist_ok=False)
    archive(root/"counts.npz",opened,bool(args.record))
    receipt(root/"result.json",summary,bool(args.record))
    print(f"{len(opened)} complete release distributions")


if __name__=="__main__":
    main()
