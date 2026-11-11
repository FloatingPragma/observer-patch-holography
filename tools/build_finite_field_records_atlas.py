"""Render the separately scoped finite-record guide from its reviewed receipts."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
from html import escape
import json
from pathlib import Path
import subprocess
import textwrap

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'docs/research/OPH_FINITE_FIELD_RECORDS.svg'
BASE = ROOT / 'model/field_record_addenda'
COLORS = ['#ffbd63', '#43e0ce', '#ae9aff', '#89b8ff', '#f7a6c6']


def short(value, *, upper=False):
    """Eight-significant-digit enclosure, without binary rounding a bound inward."""
    q = Fraction(value)
    if q <= 0:
        raise ValueError('positive bound required')
    exponent = len(str(q.numerator)) - len(str(q.denominator))
    if q < Fraction(10) ** exponent:
        exponent -= 1
    scaled = q * Fraction(10) ** (7 - exponent)
    n = -(-scaled.numerator // scaled.denominator) if upper else scaled.numerator // scaled.denominator
    if -3 <= exponent < 0:
        digits = str(n).zfill(8-exponent)
        return (digits[:exponent-7]+'.'+digits[exponent-7:]).rstrip('0').rstrip('.')
    mantissa = f'{n // 10**7}.{n % 10**7:07d}'.rstrip('0').rstrip('.')
    return mantissa if exponent == 0 else f'{mantissa} × 10^{exponent}'


def render():
    from check_field_record_addenda_claims import check
    check()
    names = ('native_quark_spatial_spectrum', 'native_code_ground_preparation',
             'native_orthogonal_wilson_spatial_join', 'native_orthogonal_wilson_charged_band',
             'native_orthogonal_wilson_charged_spatial_records')
    receipts = {name: json.loads((BASE / name / 'results.json').read_text()) for name in names}
    q, prep, spatial, charged, charged_spatial = (receipts[name] for name in names)
    quark = q['families']['uc']
    energy = quark['all_atomic_and_full_H_bins'][quark['best_certified_bin_index']]['source_excitation_relative_actual_neutral_ground']
    ground = prep['bounds_exact']
    witness = spatial['bounds_exact']
    cb = charged['bounds_exact']
    hashes = {name: hashlib.sha256((BASE / name / 'results.json').read_bytes()).hexdigest() for name in names}
    height = 815+515*len(names)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="3600" height="{height}" viewBox="0 0 3600 {height}" role="img" aria-labelledby="title desc">',
        '<title id="title">Finite records: from an OPH source to a detector</title>',
        '<desc id="desc">Five distinct finite constructions: charged transfer and energy, interacting cube preparation, revised spatial field records, charged finite clusters, and charged spatial records. Complete outcomes and exact error bounds are retained. Physical particle identities and GeV masses remain unproved.</desc>',
        '<metadata>'+escape(json.dumps({'receipt_sha256': hashes, 'physical_particle_identifications': [], 'source_derived_masses_GeV': []}, sort_keys=True))+'</metadata>',
        '<defs><linearGradient id="bg"><stop stop-color="#102139"/><stop offset="1" stop-color="#070d18"/></linearGradient><marker id="arrow" markerWidth="9" markerHeight="9" refX="7" refY="4" orient="auto"><path d="M0 0 L8 4 L0 8" fill="#7296ba"/></marker></defs>',
        f'<rect width="3600" height="{height}" fill="url(#bg)"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#e5efff}.body{font-size:27px;fill:#bacbe2}.small{font-size:23px;fill:#a2b9d7}.title{font-size:36px;font-weight:bold}.metric{font-size:35px;font-weight:bold}.label{font-size:25px;letter-spacing:2px;font-weight:bold}a text{text-decoration:underline}</style>']

    def text(x, y, value, cls='body', color=None):
        attr = '' if color is None else f' style="fill:{color}"'
        parts.append(f'<text x="{x}" y="{y}" class="{cls}"{attr}>{escape(value)}</text>')

    def lines(x, y, value, width=49, cls='body', step=39):
        for j, line in enumerate(textwrap.wrap(value, width=width)):
            text(x, y+j*step, line, cls)

    def link(x, y, path, label):
        parts.append(f'<a href="{escape(path, quote=True)}" target="_blank">')
        text(x, y, label, 'small', '#71cbff')
        parts.append('</a>')

    def card(x, y, color, number, title, body, metric, foot):
        parts.append(f'<rect x="{x}" y="{y}" width="770" height="370" rx="22" fill="#102035" stroke="#294461" stroke-width="2"/>')
        parts.append(f'<circle cx="{x+47}" cy="{y+48}" r="22" fill="{color}"/>')
        text(x+39, y+56, str(number), 'small', '#07111e')
        text(x+84, y+59, title, 'title')
        lines(x+32, y+112, body)
        text(x+32, y+269, metric, 'metric', color)
        lines(x+32, y+316, foot, width=57, cls='small', step=30)

    text(80, 80, 'THE FINITE FIELD LABORATORY', 'label', '#43e0ce')
    text(80, 165, 'From source dynamics to complete detector records', 'title')
    lines(80, 226, 'Read each lane left to right. Each lane declares its own Hamiltonian, preparation, references and clock. A change of source is shown explicitly; results from different lanes are not silently combined into a particle claim.', width=180, step=39)
    link(80, 325, 'OPH_FINITE_FIELD_RECORDS.md', 'Open the complete guide and dependency map')
    link(1050, 325, 'OPH_PARTICLE_ATLAS_V3.svg', 'Open the full particle candidate atlas')

    lanes = [
        ('A · NATIVE CHARGED SOURCE', 'Same interacting Hamiltonian; ordered arrival and energy experiment.', names[0], [
            ('Prepare a charge', 'The actual uᶜ, dᶜ or Q operator changes a native block. An explicit reference carries the compensating charge; the preparation partner remains.', 'uᶜ  +2/3   dᶜ  −1/3   Q  −1/6', 'These are source Y changes, not confirmed quark flavours.'),
            ('Let it transfer', 'A gauge link couples two complete blocks. The source evolves autonomously for one time unit before the arrival detector acts.', 'P₀ B† T Pₑ T B P₀ = P₀ / 128', 'Positive transfer residue; every nonarrival outcome stays.'),
            ('Measure the energy', 'After arrival, measure the active source plus its fixed reference. The energy operator is identity on the retained partner.', f'{short(energy[0])}–{short(energy[1],upper=True)}', 'uᶜ/dᶜ source excitation interval; reference cost 4 removed.'),
            ('Keep the whole law', 'The original record includes failed preparation, nonarrival, every spectral bin and all dark residues. Changing measurement order is a separate experiment.', '88 / 88 / 160 record labels', 'Arrival backaction is included. This is not a rest mass.')]),
        ('B · COMPLETED INTERACTING CUBE', 'Same completed cube and raw energy convention; explicit finite apparatus services.', names[1], [
            ('Retain every sector', 'The finite source contains native code, occupation complement, charge references and photon cap. A positive completion selects its code ground.', 'Code + complement + references', 'No discarded states are hidden by a projection.'),
            ('Prepare, then herald', 'Four independently prepared copies undergo recorded code admission. Three finite energy-agreement filters select a common low-energy state.', '4 copies · 4096 pointer labels', 'All rejections remain; six service slots use the same clock.'),
            ('Bound the error', 'The successful prepared state approximates the coupled ground. Exact finite control services and background-cancellation work are declared.', '< '+short(ground['conditional_ground_purification_error_upper'],upper=True), 'Conditional purification error, with total photon cap 2^64.'),
            ('Read the interaction', 'A complete selected matter–field record has positive original probability after all preparation factors. The small cube is a finite quantum system.', '> '+short(ground['original_preparation_and_T111_probability_lower']), 'Original selected-record probability; not a scattering cross section.')]),
        ('C · REVISED SPATIAL FIELD FAMILY', 'Explicit orthogonal regulator; fresh prepared sea; actual coupled evolution on growing boxes.', names[2], [
            ('Repair the free slope', 'An orthogonal regulator Γ separates the lattice correction from the mass matrix β. The free low-momentum squared coefficient becomes one.', 'E² − m² ≈ |p|²', 'Free comparison only. Vectorlike; finite-spacing P and T break.'),
            ('Prepare in one frame', 'A finite algebraic map transports the free sea while retaining the same field, charge and Gauss operators. Every bank and code flag is retained.', 'h(θ)² ≥ m² I', 'The open-boundary proof does not assume [K,W] = 0.'),
            ('Measure coupled records', 'A recorded electric kick or identity is followed by the actual coupled wait and a complete bounded-current readout. Local laws converge at fixed spacing.', '> '+short(witness['complete_finite_law_TV_separation_lower']), 'Complete finite-law separation from the zero-coupling control.'),
            ('Fix sector selection', 'A fixed cube penalty fails on large boxes. The separate source choice c(V)=18V+1 places the global ground inside the code at every finite cap.', 'Complement separation ≥ m', 'A variational sector result, not the interacting ground vector.')]),
        ('D · CHARGED FINITE CUBE', 'Explicit completion 23/2; same charge frame and raw offset; photon cap three.', names[3], [
            ('Find a stable cluster', 'This instrument uses a proved eight-state cluster in each charge sector. Its uniform interaction estimate is inconclusive for the first two states alone.', '> '+short(cb['charged_cluster_gap_over_m_lower'])+' m', 'Gap within Q = ±1. Eight states do not mean eight species.'),
            ('Prepare and compensate', 'Recorded admission and energy agreement prepare two sources. An even fermionic pump adds or removes charge; R = −Q and its cost four remain.', 'Joint parity and Gauss charge kept', 'The opposite-charge spectra are paired by conjugation.'),
            ('Record two energies', 'A charged-source clock and retained neutral-partner clock resolve the desired joint subspace. Every failure and every pointer pattern remains.', '> '+short(cb['original_accepted_energy_pair_probability_lower']), 'Original acceptance probability; five complete driven slots.'),
            ('Bound what was read', 'Accepted records lie in the charged cluster and neutral ground with a rigorous fidelity bound. The reference-corrected interval contains binding effects.', '> '+short(cb['accepted_joint_cluster_ground_fidelity_lower']), 'Joint fidelity. No continuous spin, rest mass or GeV assignment.')]),
        ('E · CHARGED PREPARED SPATIAL FAMILY', 'Same growing field and charge frame; nonstationary prepared input; static reference in each member.', names[4], [
            ('Prepare a local charge', 'A literal local creation or removal acts on the prepared free sea. An even pointer pump retains every failure and supplies the Gauss-compensating reference.', 'Pump success ≥ 11/20', 'This state is not the interacting ground.'),
            ('Send the boundary away', 'Enlarge the box around the observer. Keep the reference fixed within each experiment, then compare sources with that reference farther away.', 'Complete local records converge', 'Fixed spacing and bounded time; remote fine flags remain.'),
            ('Account for its energy', 'The reference contributes cost four and a geometry-dependent Coulomb self energy. Subtract those displayed terms when comparing charged and neutral preparations.', 'Eref = 4 + ½ m g² G(root,root)', 'Corner, face and centre have different finite self energies.'),
            ('Read a positive excess', 'Exact quantum exchange and finite photon bounds control the prepared energy above the matching neutral prepared sea. This is an expectation, not a spectral mass.', '≥ '+short(charged_spatial['prepared_energy']['reference_subtracted_limit_over_m'][0])+' m', 'A local charge signal is proved at preparation; the mass remains open.')])]

    for i, (label, subtitle, name, cards) in enumerate(lanes):
        y = 390 + i*515
        color = COLORS[i]
        text(80, y, label, 'label', color)
        text(80, y+44, subtitle, 'body')
        for j, content in enumerate(cards):
            x = 80 + j*885
            card(x, y+75, color, j+1, *content)
            if j < 3:
                parts.append(f'<path d="M{x+785} {y+255} H{x+865}" stroke="#7296ba" stroke-width="3" marker-end="url(#arrow)"/>')
        link(80, y+483, '../../model/field_record_addenda/'+name+'/PROOF.md', 'Proof, exact source definitions and all retained outcomes')

    bottom = 410+515*len(lanes)
    parts.append(f'<rect x="80" y="{bottom}" width="3425" height="290" rx="22" fill="#263024" stroke="#64794f"/>')
    text(115, bottom+54, 'WHAT STILL SEPARATES THESE RECORDS FROM A MEASURED PARTICLE', 'label', '#d1e99f')
    lines(115, bottom+105, 'Prepare the same interacting state across size → extend the finite charged energies to spatial translations → control interacting dispersion and infrared resolution → match physical spin, statistics and currents → establish dimensionless ratios → attach an energy convention and compare with experiment.', width=180, step=41)
    text(115, bottom+248, 'No Standard Model identity, source-derived GeV mass, interacting relativistic ultraviolet limit or early-universe abundance is asserted by this poster.', 'small', '#f3dbb0')
    text(80, height-45, 'Exact receipt hashes are embedded in this SVG. Every number is an outward-rounded bound from the linked finite construction.', 'small')
    parts.append('</svg>')
    return '\n'.join(parts)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--png', action='store_true')
    args = parser.parse_args()
    content = render()
    if args.check:
        if OUT.read_text() != content:
            raise SystemExit('Finite field-record atlas is stale')
    else:
        OUT.write_text(content)
    if args.png:
        subprocess.run(['rsvg-convert', '--width', '2400', '--output', str(OUT.with_suffix('.png')), str(OUT)], check=True)
    print('PASS: finite field-record atlas matches exact sources and explicit physical scope')


if __name__ == '__main__':
    main()
