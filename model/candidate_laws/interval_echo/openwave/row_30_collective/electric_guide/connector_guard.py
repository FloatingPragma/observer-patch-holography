"""Exact action of one passive connector on the two double/vacancy dimers."""

import json
from itertools import product

import sympy as sp

from model.candidate_laws.contrast_memory_spherical_v1 import law


def exact_columns():
    loads = law.N_SRC
    written = loads * law.N_LAB
    write = sp.SparseMatrix(law._write_pattern()) * sp.sqrt(2) / 2
    r, d2 = (sp.SparseMatrix(x) for x in law._label_ops())
    swap, contrast = (sp.SparseMatrix(x) for x in law._load_ops())
    echo = sp.SparseMatrix(sp.BlockMatrix([
        [sp.zeros(loads), write.T],
        [write, sp.eye(written)-write*write.T],
    ]).as_explicit())
    source_ids = [law.SRC.index(word) for word in product((0, 2), repeat=2)]
    enrolled = [s+loads*m for m in range(1+law.N_LAB) for s in source_ids]
    injection = sp.eye(law.DIM)[:, enrolled]
    for memory, blank in (((sp.eye(law.N_LAB)-r)/2, (sp.eye(loads)-swap)/2),
                          (d2, contrast)):
        memory_generator = sp.kronecker_product(memory, sp.eye(loads))
        native = sp.diag(sp.zeros(loads), memory_generator)
        claimed = sp.diag(blank, memory_generator)*injection
        delta = echo*(native*(echo*injection))-claimed
        if any(sp.expand(x)!=0 for x in delta):
            raise ArithmeticError("The passive connector changes a claimed column")
    return dict(local_dimension=law.DIM, checked_inputs=len(enrolled),
                operator_coefficients=["exchange", "contrast"],
                coefficient_entries_checked=2*law.DIM*len(enrolled),
                written_memories=law.N_LAB, unlike_blank_exchange="-1/2",
                orientation_coupling="p * (2*kappa + 1/4)",
                isolated_unlike_leakage="p * sin(t/2)^2")


if __name__ == "__main__":
    print(json.dumps(exact_columns(), indent=2))
