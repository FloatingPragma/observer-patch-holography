"""Independent-field phase reading, conditional on a declared regularity class.

The current record sets resolution only. This reader never makes its field.
No finite set of samples certifies an arbitrary callback's derivative bound.
"""
import cmath
import math

REGISTERED_SAMPLE_BUDGET = 100000
TOL = 1e-9


def declared_field(geometry, field_current, *, deformed=True):
    """A separate producer; field_current need not equal the reader's record."""
    if len(field_current) != 30 or any(type(j) is not int for j in field_current):
        raise ValueError('thirty integer field windings required')
    winding = tuple(field_current)
    edges = tuple(tuple(e) for e in geometry['dual_edges'])

    def field(e, t):
        a, b = edges[e]
        ga, gb = (.31 * math.sin(v + 1) if deformed else 0.0 for v in (a, b))
        amp = .23 * math.cos(a + b + 1) if deformed else 0.0
        return cmath.exp(1j * (ga + (gb-ga)*t + 2*math.pi*winding[e]*t
                             + amp*math.sin(2*math.pi*t)))
    return field


def read_phase(geometry, resolution_record, region, *, field,
               sample_budget=REGISTERED_SAMPLE_BUDGET):
    """Read supplied samples; every successful degree is regularity-conditional.

    Registered class: a continuous phase lift on each edge with
    sup|phi'_e| < pi*n_e/2, and a glued U(1) field on the oriented contours.
    The supplied declared_field family with |F'_e| <= 2|F_e|+1 satisfies it.
    A callback is not authenticated by this API. Its off-grid class membership
    is a source obligation, not an inference from the sampled guard.
    """
    base = {'phase_degree': None, 'sample_budget': sample_budget,
            'regularity_verified_from_samples': False,
            'scope': 'conditional on supplied field regularity; no arbitrary-field authentication'}

    def finish(status, reason, **details):
        return {**base, 'status': status, 'reason': reason, **details}

    try:
        if (not isinstance(resolution_record, (list, tuple)) or len(resolution_record) != 30
                or any(type(j) is not int for j in resolution_record)):
            return finish('rejected_input', 'thirty actual integer resolution currents required')
        if type(sample_budget) is not int or not 0 <= sample_budget <= REGISTERED_SAMPLE_BUDGET:
            return finish('rejected_input', 'budget must be an integer between zero and the registered cap')
        if (not isinstance(region, (list, tuple)) or any(type(v) is not int or not 0 <= v < 12 for v in region)
                or len(set(region)) != len(region) or not callable(field)):
            return finish('rejected_input', 'a region of distinct ports and a supplied field are required')
        loops = geometry['dual_port_boundaries']
        edges = geometry['dual_edges']
        if len(loops) != 12 or len(edges) != 30:
            return finish('rejected_geometry', 'wrong carrier size')
        for loop in loops:
            oriented = []
            for e, sign in loop:
                if type(e) is not int or not 0 <= e < 30 or type(sign) is not int or sign not in (-1, 1):
                    return finish('rejected_geometry', 'malformed oriented boundary')
                a, b = edges[e]
                oriented.append((a, b) if sign == 1 else (b, a))
            if not oriented or any(oriented[k][1] != oriented[(k+1) % len(oriented)][0]
                                   for k in range(len(oriented))):
                return finish('rejected_geometry', 'an oriented cell boundary is not closed')
        edge_ids = sorted({e for v in region for e, _ in loops[v]})
        resolutions = {e: 8*abs(resolution_record[e])+8 for e in edge_ids}
        required = sum(n+1 for n in resolutions.values())
        base.update(required_samples=required, resolutions=resolutions)
        # This integer-only check precedes evaluating field or converting a huge current.
        if required > sample_budget:
            return finish('resource_stop', 'complete readout exceeds registered sample budget')
        samples = {}
        for e, n in resolutions.items():
            values = [complex(field(e, k/n)) for k in range(n+1)]
            if any(not math.isfinite(z.real) or not math.isfinite(z.imag)
                   or abs(abs(z)-1) >= TOL for z in values):
                return finish('invalid_field', 'samples must be finite unit complex numbers')
            samples[e] = values
        total = 0
        float_total = 0.0
        max_increment = 0.0
        max_defect = 0.0
        for v in region:
            contour = []
            for e, sign in loops[v]:
                values = samples[e] if sign == 1 else list(reversed(samples[e]))
                if contour:
                    max_defect = max(max_defect, abs(contour[-1]-values[0]))
                    contour.extend(values[1:])
                else:
                    contour.extend(values)
            max_defect = max(max_defect, abs(contour[-1]-contour[0]))
            if max_defect >= TOL:
                return finish('open_curve', 'field does not glue around the cell', gluing_defect=max_defect)
            increments = [cmath.phase(b*a.conjugate()) for a, b in zip(contour, contour[1:])]
            max_increment = max(max_increment, max(map(abs, increments), default=0.0))
            if max_increment >= math.pi/2:
                return finish('increment_guard', 'sampled phase increment reaches the pi/2 guard',
                              max_principal_increment=max_increment)
            value = math.fsum(increments)/(2*math.pi)
            degree = round(value)
            if abs(value-degree) >= TOL:
                return finish('non_integer', 'closed-contour sum is not numerically integral', float_degree=value)
            total += degree
            float_total += value
        return finish('read', 'supplied field read under the stated regularity contract',
                      phase_degree=total, float_degree=float_total,
                      max_principal_increment=max_increment, gluing_defect=max_defect)
    except (ValueError, TypeError, KeyError, IndexError, OverflowError) as error:
        return finish('rejected_input', str(error))


def compare_degree(reading, claimed_degree):
    """A real comparison: the reader has not used claimed_degree to make its field."""
    return (reading.get('status') == 'read' and type(claimed_degree) is int
            and reading['phase_degree'] == claimed_degree)
