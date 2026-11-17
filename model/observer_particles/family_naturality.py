#!/usr/bin/env python3
"""Apply the actual source family-band theorem to invariant mass operators.

This tests the conditional family attachment, not an unconditional physical
family assignment. A chosen port is a symmetry-breaking control, not a law
selected by the source. All arithmetic is exact in Q(sqrt(5)).
"""
from __future__ import annotations
from evidence_contract import public_receipt
import hashlib
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True
import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PUBLIC = ROOT
sys.path.insert(0, str(PUBLIC / 'code/a5_closure'))
import family_band_attachment_certificate as bands


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def exact(x):
    return s.Rational(str(x.a)) + s.Rational(str(x.b)) * s.sqrt(5)


def clean(matrix):
    return matrix.applyfunc(s.simplify)


def serialize(matrix):
    return [[str(x) for x in row] for row in matrix.tolist()]


def main():
    carrier, rotations, pin = bands.load_carrier()
    adjacency = bands.lift(bands.adjacency_int(carrier))
    raw = bands.spectral_projectors(adjacency)
    p = s.Matrix([[exact(x) for x in row] for row in raw['3']])
    assert p*p == p and p.T == p and s.trace(p) == 3
    basis = s.Matrix.hstack(*p.columnspace())
    gram = clean(basis.T*basis)
    left = clean(gram.inv()*basis.T)
    assert clean(left*basis) == s.eye(3)
    representations = []
    for g in rotations:
        u = s.zeros(12)
        for j, i in enumerate(g):
            u[i,j] = 1
        r = clean(left*u*basis)
        assert clean(u*basis-basis*r) == s.zeros(12,3)
        assert clean(r.T*gram*r-gram) == s.zeros(3)
        representations.append(r)
    # Schur's lemma is checked on this concrete representation: the average
    # squared character is one. The finite source rotation table is the group
    # used by the owning exact certificate, not a chosen continuous SU(3).
    characters = [s.trace(r) for r in representations]
    inner = s.simplify(sum(c*c for c in characters)/60)
    assert inner == 1
    # Build the actual commutant equations. Retain two independent source
    # rotations sufficient to leave one unknown scalar, then verify that
    # scalar commutes with the entire group.
    units = [s.eye(3)[:,i]*s.eye(3)[j,:] for i in range(3) for j in range(3)]
    constraints = s.zeros(0,9)
    witness_rotations = []
    for index, r in enumerate(representations):
        block = s.Matrix.hstack(*[s.Matrix(9,1,list(clean(r*e-e*r))) for e in units])
        candidate = constraints.col_join(block)
        reduced, pivots = candidate.T.rref(simplify=s.simplify)
        if len(pivots) > constraints.rows:
            constraints = candidate[list(pivots),:]
            witness_rotations.append(index)
        if constraints.rows == 8:
            break
    kernel = constraints.nullspace()
    assert len(kernel) == 1
    only = s.Matrix(3,3,list(kernel[0]))
    assert only == s.eye(3)*only[0,0] and only[0,0] != 0
    assert all(clean(r*only-only*r) == s.zeros(3) for r in representations)
    # The full matrix deformation representation is 1 + 3 + 5. This says
    # which additional nonscalar data can split a family; it chooses none.
    chi_wedge = [s.simplify((s.trace(r)**2-s.trace(r*r))/2) for r in representations]
    chi_sym0 = [s.simplify((s.trace(r)**2+s.trace(r*r))/2-1) for r in representations]
    chars = [[s.Integer(1)]*60,chi_wedge,chi_sym0]
    character_gram = [[s.simplify(sum(x*y for x,y in zip(a,b))/60) for b in chars] for a in chars]
    assert character_gram == [[1,0,0],[0,1,0],[0,0,1]]
    assert chi_wedge == characters
    # A scalar-valued Higgs under this source family action forces a Yukawa
    # intertwiner between the same triplet representation to be y I. A
    # bilinear invariant has matrix y G in this nonorthonormal basis.
    assert all(clean(r.T*gram*r) == gram for r in representations)
    normalized_mass_squared = []
    for y in (s.Integer(1),s.Integer(2)):
        ybilinear = y*gram
        mass2 = clean(gram.inv()*ybilinear.T*gram.inv()*ybilinear)
        assert mass2 == y*y*s.eye(3)
        normalized_mass_squared.append({'chosen_coefficient':str(y),'squared_singular_values':[str(y*y)]*3})
    # Source-locality alone permits nonscalar attachments. An independently
    # chosen port gives a rank-one form on the same band. Its exact twirl
    # removes all splitting; no preferred port is inferred from the source.
    port = s.zeros(12);port[0,0] = 1
    perturbation = clean(left*port*basis)
    assert clean(perturbation*perturbation) == perturbation/4
    assert s.trace(perturbation) == s.Rational(1,4)
    twirl = clean(sum((r*perturbation*r.inv() for r in representations),s.zeros(3))/60)
    assert twirl == s.eye(3)/12
    stabilizer = sum(clean(r*perturbation-perturbation*r) == s.zeros(3) for r in representations)
    assert stabilizer < 60
    # A fixed conditional control, Y=I+4S, has singular values (1,1,2).
    # It is G-self-adjoint; this is a demonstration of broken symmetry, not
    # a fit to a particle spectrum or an empirical-valued source candidate.
    chosen = s.eye(3)+4*perturbation
    assert clean(chosen.T*gram-gram*chosen) == s.zeros(3)
    assert clean(chosen*chosen).eigenvals() == {s.Integer(1):2,s.Integer(4):1}
    artifact = HERE/'results/family_naturality_operators.json'
    artifact.write_text(json.dumps({'source_band_projector':serialize(p),'basis':serialize(basis),
        'kinetic_metric':serialize(gram),'source_rotations':list(map(serialize,representations)),
        'commutant_constraint_witness_indices':witness_rotations,
        'commutant_constraint_matrix':serialize(constraints),'chosen_port_form':serialize(perturbation)},indent=2)+'\n', encoding='utf-8', newline='\n')
    paths = [Path(bands.__file__),PUBLIC/'Lean/Screen/A5FamilyBand.lean',
        PUBLIC/'code/a5_closure/manifests/echosahedral_federation_reference.json',
        PUBLIC/'code/a5_closure/manifests/family_band_attachment_reference.json',
        PUBLIC/'paper/tex_fragments/UNIFIED_OBSERVER_PHYSICS_SPINE.tex',
        PUBLIC/'Lean/Screen/WeylYukawaConventions.lean']
    result = {
        'status':'PASS_EXACT_SOURCE_FAMILY_NATURALITY_AND_SPLITTING_CONTROL',
        'producer_sha256':sha(__file__),'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
        'artifacts':{artifact.name:sha(artifact)},'empirical_mass_or_coupling_inputs':[],
        'source_carrier_pin':pin,'source_rotations_checked':60,'source_band_dimension':3,
        'character_norm':str(inner),'commutant_complex_dimension':len(kernel),
        'endomorphism_deformation_irreducible_dimensions':[1,3,5],
        'endomorphism_character_gram':[[str(x) for x in row] for row in character_gram],
        'conditional_mass_conclusion':{
            'premises':['family attached to the existing selected source rank3 band',
                'the same unbroken source A5 action on both family factors',
                'family-invariant Higgs/vacuum','invariant kinetic metrics and Yukawa intertwiner'],
            'Yukawa_form':'y I3 in common endomorphism coordinates; y G in bilinear coordinates',
            'squared_singular_values':'|y|^2/kL/kR, repeated three times',
            'relative_family_masses_if_identified':'1:1:1 for nonzero y',
            'overall_y_or_kinetic_normalizations_selected':False,
            'same_constraints_coefficient_controls':normalized_mass_squared},
        'chosen_local_symmetry_breaking_control':{
            'source_port_index_chosen':0,'form_rank':1,'form_nonzero_eigenvalue':'1/4',
            'full_group_twirl':'I3/12','form_stabilizer_size':stabilizer,
            'chosen_Y':'I3+4S_port','chosen_relative_singular_values':['1','1','2'],
            'source_selected_control':False,
            'meaning':'covariance transports an anisotropic field; invariance of one state is a stronger condition. A selected symmetry-breaking background could split the band, but none is supplied here.'},
        'scope':{'actual_source_family_band_used':True,'new_abstract_Schur_theorem_claimed':False,
            'family_attachment_unconditional':False,'excludes_symmetry_breaking_routes':False,
            'all_OPH_mass_realizations_excluded':False,'physical_SM_masses_derived':False}}
    (HERE/'results/family_naturality.json').write_text(json.dumps(public_receipt(result),indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'status':result['status'],'commutant_dimension':len(kernel),
        'conditional_family_mass_ratios':'1:1:1','symmetry_breaking_control_stabilizer':stabilizer},indent=2))


if __name__ == '__main__':
    main()
