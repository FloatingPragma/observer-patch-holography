"""Passive dark storage and recurring control on the finite age register of the torus law."""
import json
import math

import numpy as np

from .. import formation, law
from . import clock_factory as factory
from . import clock_loader as loader


L = formation.AGE_PERIOD


def _register(value, label):
    if not isinstance(value, (int, np.integer)) or isinstance(value, bool) or value < 0:
        raise ValueError(f"{label} must be a nonnegative integer")
    return int(value)


def _phase(start, length):
    """Relative routing-odd phase for a contiguous interval inside one age cycle."""
    opening = loader.SOURCE.maturation[2]
    partial = max(0,min(start+length,opening)-start)
    angles = math.fsum(law.hop_angle(j,loader.SOURCE.maturation) for j in range(start,start+partial))
    return (-1)**((length-partial)%2) * np.exp(-2j*angles)


def transfer(waves, age, steps, *, saturated=False):
    """Exact passive transfer on both dark routing parities, including age wraps.

    The saturated control keeps the terminal full hop forever. The two tests
    of darkness refer to both wave and F-wave at every site, before powering.
    """
    age, steps = _register(age,"age"), _register(steps,"steps")
    waves = np.asarray(waves)
    if (age >= L or waves.ndim != 6 or waves.shape[0] != 2 or waves.shape[-1] != 8
            or 0 in waves.shape or not np.isfinite(waves).all()):
        raise ValueError("two finite four-axis rail fields and a valid register age required")
    routed = np.stack([law.route(wave) for wave in waves])
    for v in waves:
        # A rail may have arbitrarily small branch amplitude. Check its
        # direction before squaring components to form a norm.
        scale = max(np.max(abs(v.real)),np.max(abs(v.imag)))
        if scale == 0:
            continue
        unit = v.real/scale + 1j*(v.imag/scale)
        if max(np.linalg.norm(unit@law.EQUAL),
               np.linalg.norm(law.route(unit)@law.EQUAL)) > 2e-11*np.linalg.norm(unit):
            raise ValueError("passive dark transfer requires both pointwise brightnesses to vanish")
    if saturated:
        relative = _phase(age,steps)
    else:
        wraps, tail = divmod(steps,L)
        head = min(tail,L-age)
        relative = (-1j)**(wraps%4) * _phase(age,head) * _phase(0,tail-head)
    result = (waves+routed)/2 + relative*(waves-routed)/2
    for r, phi in enumerate(loader.SOURCE.rail_contact_phases):
        # The two source phases make exactly 128 and 129 dark-rail turns
        # per register period. Reducing integers avoids large-angle loss.
        turns = round(L*phi/(4*np.pi))
        if abs(turns*4*np.pi/L-phi)>1e-15:
            raise ValueError("the enrolled finite rail period has changed")
        result[r] *= (-1)**(steps%2)*np.exp(2j*np.pi*((turns*steps)%L)/L)
    return result


def capture_vertex(ordinary, matched, *, inverse=False):
    """One eligible native vertex, retaining every program, rail and reference index.

    The last axis is the matched slot. Ordinary input has an empty wire head,
    and its partner has an occupied head and age zero. The program digit
    belongs to constituent A in either sector.
    """
    ordinary, matched = np.asarray(ordinary),np.asarray(matched)
    if (matched.shape!=(*ordinary.shape,8) or not np.isfinite(ordinary).all()
            or not np.isfinite(matched).all() or type(inverse) is not bool):
        raise ValueError("matching finite amplitudes and a boolean direction required")
    angle = -formation.ANGLE if inverse else formation.ANGLE
    c,s = np.cos(angle),np.sin(angle)
    inlet = matched @ law.EQUAL
    return (c*ordinary-s*inlet,
            matched+(s*ordinary+(c-1)*inlet)[...,None]*law.EQUAL)


def return_probe(waves, head):
    """The occupied, age-zero reverse capture with an empty unmarked input.

    The pair released at the first separation is counted in the ordinary
    output. The retained field includes the attenuation in the same slot ray.
    """
    kept = np.array(waves, dtype=complex, copy=True)
    if kept.ndim != 6 or kept.shape[0] != 2 or kept.shape[-1] != 8 or not np.isfinite(kept).all():
        raise ValueError("reverse capture requires two finite four-axis fields")
    head = tuple(_register(x,"head coordinate") for x in head)
    if len(head) != 4 or any(x >= size for x,size in zip(head,kept.shape[1:-1])):
        raise ValueError("head outside the carrier")
    site = (slice(None), *head)
    released, kept[site] = capture_vertex(np.zeros(2),kept[site])
    return kept, abs(released)**2


def program_step(waves, age, lap, capacity, gates):
    """One service on a basis block of a supplied co-moving program digit.

    The same map on every digit block extends linearly to coherent inputs.
    Capture and terminal readout act on their own registers.
    """
    age, lap = _register(age,"age"), _register(lap,"program digit")
    capacity = _register(capacity,"program capacity")
    if age >= L or lap >= capacity:
        raise ValueError("program address outside the finite counter")
    new_age = (age+1) % L
    new_lap = (lap+int(new_age == 0)) % capacity
    plate = gates.get(new_age) if new_lap == 0 else None
    return factory.material_service(waves,age,plate),new_age,new_lap


def program_trial(born, gates):
    """Compare five program states with the one-state repeating control."""
    rows = []
    for capacity in (1,5):
        waves, age, lap = born.copy(),0,0
        for _ in range(factory.END):
            waves,age,lap = program_step(waves,age,lap,capacity,gates)
        for wrap in range(1,capacity+1):
            # Skip only verified passive darkness. The boundary transition
            # itself, the inverse capture and every active plate are explicit.
            waves = transfer(waves,age,L-1-age)
            waves,age,lap = program_step(waves,L-1,lap,capacity,gates)
            # At age zero this first vertex sees a dark pair for either digit.
            waves,release = return_probe(waves,(2,2,1,1))
            for _ in range(factory.END):
                waves,age,lap = program_step(waves,age,lap,capacity,gates)
            rows.append({"capacity":capacity,"wrap":wrap,"post_wrap_digit":lap,
                         "age":age,"release":float(release.sum()),
                         "bright_weight":float(np.linalg.norm(waves@law.EQUAL)**2)})
    return rows


def experiment():
    shape, root = (5,5,3,3),(2,2,1,1)
    gates = factory.uniform_schedule()
    born = np.zeros((2,*shape,8),complex)
    born[(slice(None),*root)] = law.EQUAL/np.sqrt(2)
    state = born.copy()
    check = []
    def gate(condition, name):
        if not condition:
            check.append(name)
    for age in range(128):
        state = factory.advance(state,age,gates)
    book = []
    widths = (64,4096,65)
    releases = []
    transported = np.stack([law.route(w) for w in state])
    for cycle in range(5):
        fields = transfer(state,128,cycle*L)
        expected = (state+transported)/2 + (-1j)**cycle*(state-transported)/2
        gate(np.max(abs(fields-expected))<2e-11,f"both dark parities, cycle {cycle}")
        if cycle:
            elapsed = cycle*L
            eligible = [width for width in widths if elapsed % math.lcm(width,L) == 0]
            zero_age = transfer(state,128,elapsed-128)
            retained, released = return_probe(zero_age,root)
            error = float(np.sum(released))
            everywhere = float(np.max(abs(zero_age @ law.EQUAL)))
            releases.append({"elapsed":elapsed,"age":0,"wire_at_occupied_head":eligible,
                             "release_at_natal_head":error,"largest_bright_amplitude":everywhere})
            gate(eligible==[64,4096] and error<2e-22 and everywhere<2e-11
                 and np.max(abs(retained-zero_age))<2e-11,"eligible occupied-head return")
        power = np.sum(abs(fields)**2,axis=(0,-1))
        ports = fields[(slice(None),*root)]
        local = [float(np.sum(abs(ports[0]+phase*ports[1])**2)/4) for phase in (1,1j,-1,-1j)]
        book.append({"age":128,"elapsed":128+cycle*L,"root_ports":local,
                     "born_root_miss":float(1-power[root])})
        gate(abs(sum(local)-power[root])<2e-11,"complete readout census")
    before_wrap = transfer(state,128,L-128)
    kept, loss = return_probe(before_wrap,root)
    gate(loss.sum()<2e-22 and np.max(abs(kept-before_wrap))<2e-11,"true age-zero return")
    _, positive = return_probe(born,root)
    gate(abs(positive.sum()-.5)<2e-13,"occupied bright reverse vertex")
    restarted = kept.copy()
    for j in range(factory.END):
        restarted = factory.advance(restarted,j,gates)
    restarted_bright = float(np.sum(abs(restarted@law.EQUAL)**2))
    saturated = transfer(state,128,L,saturated=True)
    native = transfer(state,128,L)
    overlap = abs(np.vdot(native,saturated))**2
    gate(abs(overlap-.5)<2e-11,"one-wrap saturated control fidelity")
    gate(restarted_bright>1e-6,"static controller restarts on stored material")
    program = program_trial(born,gates)
    for row in program:
        rearmed = row["post_wrap_digit"] == 0
        gate(row["release"] < 2e-22 and
             (row["bright_weight"]>1e-6 if rearmed else row["bright_weight"]<2e-22),
             f"finite commissioning digit {row['capacity']}, wrap {row['wrap']}")
    return {"age_period":L,"wire_periods":list(widths),"eligible_returns":releases,
            "records_conditional_on_first_write":book,
            "sharp_original_first_write":float(np.sin(formation.ANGLE)**2),
            "saturated_control_fidelity":float(overlap),"bright_return_release":float(positive.sum()),
            "restarted_local_brightness":restarted_bright,"finite_program":program,"failures":check}


if __name__ == "__main__":
    result = experiment()
    print(json.dumps(result,indent=2))
    raise SystemExit(int(bool(result["failures"])))
