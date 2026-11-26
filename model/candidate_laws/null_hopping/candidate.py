"""The contrast-memory generator with both agreement terms switched off."""
from model.candidate_laws.contrast_memory_spherical_v1 import candidate as reference

NAME = "null_hopping"
TITLE = "Bare hopping"
ROLE = "control"
SUMMARY = "Exchange at zero agreement energy isolates the part of the response fixed by the port graph."
DOCUMENTS = {"readme": "README.md", "specification": "README.md"}
LABELLING = reference.LABELLING
FEDERATION_LEVELS = reference.FEDERATION_LEVELS
KAPPA = 0.0
PARAMETERS = {
    "kappa": {"value": KAPPA, "provenance": "designed",
              "quantity": "agreement energy in units of the seam exchange"},
}
PREMISES = {}
JOINT = {}

basis = reference.basis
represent = reference.represent


def generator(car):
    return reference.generator(car, kappa=KAPPA)


def response(car):
    return generator(car).toarray()


def run(argv=None):
    from .replay import main
    return main(argv)
