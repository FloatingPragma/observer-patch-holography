"""Frozen apparatus-error diagnostic; exact architecture is the main theorem."""
import cmath
import json
from model import I, EDGES, CALLS, add, apply, density, fine, mm, scaled, star, trace, transpose


def adjoint(a):
    return tuple(tuple(complex(x).conjugate() for x in r) for r in transpose(a))


def cp(a, rho):
    return mm(mm(a, rho), adjoint(a))


rho = density((1, 1, 1))
warm = scaled(I, 0)
for m in fine(EDGES, "native"):
    warm = add(warm, scaled(apply(m, rho), 1 / 6))
pedge = ((.5, .5, 0), (.5, .5, 0), (0, 0, 1))
q = density((1, -1, 0))
out = []
for delta in (0, 1e-9, -1e-9, 1e-7, -1e-7, 1e-5, -1e-5, 1e-3, -1e-3):
    l = scaled(add(I, pedge, -cmath.exp(1j * delta)), 12 ** -.5)
    prepared = cp(l, warm)
    pprep = trace(prepared).real
    for turns in (1, 2):
        a = b = I
        v = 2
        for hop in range(3 * turns):
            v = (v + 1) % 3
            for event in range(4):
                a = mm(star(v, "native")[0], a)
                b = mm(star(2, "native")[0], b)
        b = scaled(b, 1 / 8 ** turns)
        aa, bb = cp(a, prepared), cp(b, prepared)
        cross = mm(mm(a, prepared), adjoint(b))
        na, nb = trace(aa).real, trace(bb).real
        out.append({"preparation_phase_error_radians": delta, "turns": turns,
                    "preparation_success": pprep,
                    "prepared_contrast_fidelity": trace(mm(q, prepared)).real / pprep,
                    "rotated_contrast_fidelity": trace(mm(q, aa)).real / na,
                    "visibility": 2 * trace(cross).real / (na + nb),
                    "unconditional_success": (na + nb) / (2 * 3 ** (12 * turns))
                                             * (255 / 256) ** (3 * turns)})
print(json.dumps({"kind": "frozen apparatus-error diagnostic, no gate selection",
                  "cases": out, "literal_endpoint_calls": CALLS}, indent=2))
