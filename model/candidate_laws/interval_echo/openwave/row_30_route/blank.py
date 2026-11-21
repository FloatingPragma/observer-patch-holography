"""Blank memory closure for arbitrary loads and the exact free one-unit dispersion."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation

HERE = Path(__file__).resolve().parent


def run():
    law = candidate.native_law
    kappa = candidate.KAPPA
    symbolic = law.exact_checks()["checks"]
    required = ("h_kappa_J_equals_J_seam_term", "blank_code_reducing",
                "reference_pair_D2_reduces_to_contrast_square", "reference_stays_at_zero")
    if not all(symbolic[k] for k in required):
        raise ArithmeticError("Local exact verification failed")
    exchange, contrast = law.local_term_parts()
    J = np.eye(law.DIM)[:, :law.N_SRC]
    native = (exchange + kappa * contrast) @ J
    contrast_blank = contrast @ J
    ref = np.array([law.SRC.index((n, 0)) for n in range(3)])
    one = np.array([law.SRC.index((1, 0)), law.SRC.index((0, 1))])
    block = native[np.ix_(one, one)]
    ref_one = kappa * contrast_blank[one[0], one[0]]
    ref_error = np.max(abs(contrast_blank[:, ref] - J[:, ref] * np.arange(3)**2))
    seams_error = np.max(abs(native - J @ law.seam_term(kappa)))
    if max(ref_error, seams_error) > 1e-12 or np.max(abs(native[:,law.SRC.index((0,0))])) > 1e-12:
        raise ArithmeticError("Full native seam or blank zero-reference columns differ")
    boards = []
    for level in range(3):
        car = federation(level)
        a,b = np.array(car.seams).T
        data = np.repeat(block.ravel(),len(a))
        H = sparse.coo_matrix((data, (np.r_[a,a,b,b],np.r_[a,b,a,b])),shape=(car.n, car.n)).tocsr()
        H += sparse.diags((6 - car.degree) * ref_one)
        lap = car.laplacian()
        offset = 6*kappa
        expected = offset*sparse.eye(car.n) + lap/2
        residue = abs(H-expected).max()
        implementation_error = abs(candidate.generator(car)-expected).max()
        # For a connected graph, Tr(L)^2/(n-1) equals Tr(L^2) iff all positive levels agree.
        trace = int(car.degree.sum())
        trace2 = int(car.degree @ car.degree + trace)
        separation = int((car.n-1)*trace2 - trace**2)
        if connected_components(car.adjacency(), directed=False)[0] != 1 or separation <= 0:
            raise ArithmeticError("No certified pair of distinct positive wave numbers")
        if max(residue,implementation_error) > 1e-12:
            raise ArithmeticError("Whole-board identity fails")
        boards.append(dict(cells=car.cells, ports=car.n, seams=len(a), one_unit_offset=offset,
                           all_nonzero_modes_with_multiplicity=car.n-1,
                           positive_levels_variance_numerator=separation, exact_graph_traces=[trace,trace2],
                           max_native_matrix_residual=float(f"{residue:.2g}"),
                           max_registered_matrix_residual=float(f"{implementation_error:.2g}")))
    return dict(candidate=candidate.NAME,kappa=kappa, exact_local_identities={k:symbolic[k] for k in required},
                native_dimension=law.DIM, checked_blank_columns=law.N_SRC,
                reference_zero_input_loads=[0,1,2],square_law_residual_quadratic_coefficient=.25,
                boards=boards)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result=run()
    if args.output:
        with args.output.open("x") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")
    elif result != json.loads((HERE/"blank.json").read_text()):
        raise SystemExit("Blank-collar verification diverged from archive")
    for board in result["boards"]:
        print(board)
