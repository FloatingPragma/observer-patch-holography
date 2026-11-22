"""The private-record graph in its prepared one-load sector."""
import numpy as np
from scipy import sparse
from model.executable_model import CELL_SEAMS
from model.experiments.interval_echo_consensus.sign_graph import questions

NAME = "private_record_walk"
TITLE = "Private record walk"
ROLE = "experimental"
SUMMARY = "Each load exchange and each enabled private sign reversal is one edge of the same unweighted graph."
DOCUMENTS = {"readme": "README.md", "specification": "protocol.json"}
LABELLING = CELL_SEAMS
FEDERATION_LEVELS = ("cell", 0, 1, 2)
PARAMETERS = {
    "moves": {"value": "primitive", "provenance": "designed",
              "quantity": "private binary comparison at each portable address, with one common edge rate"},
    "reference": {"value": "flatten single-load costs", "provenance": "designed",
                  "quantity": "local unary signs up to the largest single-port comparison count"},
    "preparation": {"value": "protocol.json", "provenance": "designed",
                    "quantity": "responding signs, sharp launch and tagged sign preparations"},
}
PREMISES = {}
JOINT = {}


def basis(car):
    return [("port", v) for v in range(car.n)]


def represent(rotation, car):
    return np.eye(car.n)[:, rotation]


def generator(car):
    """One odd parity per sign. Each incident comparison contributes its unit splitting."""
    offers, _ = questions(car, "primitive")
    counts = [0]*car.n
    for a, b, _ in offers:
        counts[a] += 1
        if b is not None:
            counts[b] += 1
    return car.laplacian()/2 + sparse.diags(counts, format="csr", dtype=float)


def response(car):
    return generator(car)


def run(argv=None):
    from .release import main
    return main(argv)
