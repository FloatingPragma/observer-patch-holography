#!/usr/bin/env python3
"""Read internal multiplets of the existing geometric response processor.

No mass, Standard Model matter table, kinetic action or laboratory datum is
an input. The supplied geometric response is distinguished from a source
reconstruction of the strengthened A1/A2 response from primitive repair.
"""
from __future__ import annotations
from evidence_contract import public_receipt
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import sys
sys.dont_write_bytecode=True
import numpy as np
import sympy as sp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PUBLIC=ROOT
sys.path.insert(0,str(PUBLIC))
from model.observer.source_scalar import compiler


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def scalar(cell):
    a,b,c,d=map(sp.Rational,cell)
    return a+b*sp.sqrt(5)+sp.I*(c+d*sp.sqrt(5))


def full_blocks(blocks):
    return sp.diag(*(sp.Matrix([[scalar(x) for x in row] for row in block]) for block in blocks))


def clean(matrix):return matrix.applyfunc(sp.simplify)


@lru_cache(maxsize=1)
def decoded_generators():
    """Twelve exact skew generators recovered from 35 response probes each.

    Probe inputs are matrix units. Their commutator outputs determine the
    generator modulo scalar iI; fix zero total trace. No stored generator,
    named Lie type, matter representation or gauge coupling enters decoding.
    """
    packet=json.loads((PUBLIC/'code/source_selection_model/response.json').read_text())
    answer=[]
    for records in packet['matrix_unit_response_probes']:
        assert len(records)==35
        values=iter(map(scalar,records));d=sp.zeros(6)
        for i in range(6):
            for j in range(6):
                if i!=j:d[i,j]=next(values)
        for i in range(5):d[i,i]=next(values)
        d=clean(d-sp.trace(d)*sp.eye(6)/6)
        assert clean(d+d.conjugate().T)==sp.zeros(6)
        answer.append(sp.ImmutableMatrix(d))
    assert len(answer)==12
    return tuple(answer)


def adjoint(matrix):
    """Commutator on all36 operator coordinates, row-major vectorization."""
    return sp.kronecker_product(matrix,sp.eye(6))-sp.kronecker_product(sp.eye(6),matrix.T)


def matrix_unit(i,j):
    a=sp.zeros(6);a[i,j]=1;return a


def canonical_coordinates():
    # These convenient coordinates are only admitted AFTER exact equality
    # of their span to the reconstructed source span is proved in main().
    color=[]
    for i,j in ((0,1),(0,2),(1,2)):
        color.extend(((matrix_unit(i,j)+matrix_unit(j,i))/2,
                      -sp.I*(matrix_unit(i,j)-matrix_unit(j,i))/2))
    color.extend((sp.diag(1,-1,0,0,0,0)/2,sp.diag(1,1,-2,0,0,0)/2))
    weak=[]
    for i,j in ((4,5),(5,3),(3,4)):
        weak.append(-sp.I*(matrix_unit(i,j)-matrix_unit(j,i)))
    relative=sp.diag(1,1,1,-1,-1,-1)/2
    return tuple(color),tuple(weak),relative


def projector(operator,eigenvalue,spectrum):
    result=sp.eye(operator.rows)
    for other in spectrum:
        if other!=eigenvalue:result=result*(operator-other*sp.eye(operator.rows))/(eigenvalue-other)
    return clean(result)


def geometric_check(generators):
    path=PUBLIC/'code/source_selection_model/response.py'
    spec=importlib.util.spec_from_file_location('gauge_identity_geometric_source',path)
    source=importlib.util.module_from_spec(spec);spec.loader.exec_module(source)
    verts,_,_,antipode=source.incidence.port_model()
    frame=source.algebra.FrameRealization(SimpleNamespace(antipode=antipode),tuple(range(12)),verts)
    packet=json.loads((PUBLIC/'code/source_selection_model/response.json').read_text())
    for p,field in enumerate(source.algebra.BASIS_FIELDS):
        regenerated=full_blocks(source.encode_blocks(source.geometric_response(field,frame)))
        assert regenerated==full_blocks(packet['generators'][p])
        assert clean(regenerated-sp.trace(regenerated)*sp.eye(6)/6-generators[p])==sp.zeros(6)
    for row in packet['ordered_mixed_responses']:
        p,q=row['ports'];actual=generators[p]*generators[q]-generators[q]*generators[p]
        assert clean(actual-full_blocks(row['mixed_coefficient']))==sp.zeros(6)
    return len(packet['ordered_mixed_responses'])


def main():
    out=HERE/'results';out.mkdir(exist_ok=True)
    generators=decoded_generators()
    jets=geometric_check(generators)
    flat=sp.Matrix.hstack(*(sp.Matrix(d).reshape(36,1) for d in generators))
    color,weak,relative=canonical_coordinates()
    canonical=tuple(sp.I*t for t in (*color,*weak,relative))
    basis=sp.Matrix.hstack(*(d.reshape(36,1) for d in canonical))
    assert flat.rank()==basis.rank()==12
    assert flat.row_join(basis).rank()==12
    coefficients=clean((basis.conjugate().T*basis).inv()*basis.conjugate().T*flat)
    assert clean(basis*coefficients-flat)==sp.zeros(36,12)
    assert all(sp.im(c)==0 for c in coefficients)
    # Canonical coordinates are now verified source combinations, not an
    # inserted gauge model. Read their complete operator representation.
    cc=sum((adjoint(t)**2 for t in color[:-1]),sp.zeros(36))+adjoint(color[-1])**2/3
    cw=sum((adjoint(t)**2 for t in weak),sp.zeros(36))
    charge=adjoint(relative)
    assert cc*cw==cw*cc and charge*cc==cc*charge and charge*cw==cw*charge
    z=sp.Symbol('lambda')
    assert sp.expand(cc.charpoly().as_expr()-z**10*(z-3)**8*(z-sp.Rational(4,3))**18)==0
    assert sp.expand(cw.charpoly().as_expr()-z**10*(z-2)**21*(z-6)**5)==0
    specs=((0,0,0,2),(3,0,0,8),(0,2,0,3),(0,6,0,5),
           (sp.Rational(4,3),2,1,9),(sp.Rational(4,3),2,-1,9))
    rows=[];projectors=[]
    for c,w,y,dimension in specs:
        p=projector(cc,c,(0,sp.Rational(4,3),3))*projector(cw,w,(0,2,6))*projector(charge,y,(-1,0,1))
        assert p*p==p and p.conjugate().T==p and sp.trace(p)==dimension
        for d in generators:assert clean(p*adjoint(d)-adjoint(d)*p)==sp.zeros(36)
        projectors.append(p)
        rows.append({'color_Casimir':str(c),'weak_Casimir':str(w),
                     'relative_central_charge':str(y),'complex_operator_dimension':dimension})
    assert sum(projectors,sp.zeros(36))==sp.eye(36)
    assert all(a*b==sp.zeros(36) for i,a in enumerate(projectors) for b in projectors[:i])
    # The source-current subspace occupies adjoints8+3 and ONE of the two
    # singlets. The remaining operator directions form a24dim complement.
    current_projection=clean(basis*(basis.conjugate().T*basis).inv()*basis.conjugate().T)
    assert current_projection**2==current_projection and sp.trace(current_projection)==12
    current_dimensions=[int(sp.trace(p*current_projection)) for p in projectors]
    assert current_dimensions==[1,8,3,0,0,0]
    # Every color-singlet operator is neutral under the source center. This
    # is an exact domain restriction, stronger than a dimension count.
    color_singlets=projector(cc,0,(0,sp.Rational(4,3),3))
    assert sp.trace(color_singlets)==10 and charge*color_singlets==sp.zeros(36)
    # A declared dimensionless test order parameter, not an emitted VEV:
    # measure its gauge-orbit tangent Gram without kinetic coefficients.
    order_parameter=weak[2]
    orbit=[t*order_parameter-order_parameter*t for t in (*color,*weak,relative)]
    orbit_gram=sp.Matrix(12,12,lambda i,j:sp.trace(orbit[i].conjugate().T*orbit[j]))
    assert orbit_gram==sp.diag(*([0]*8+[2,2,0,0]))
    assert cw.charpoly().as_expr().subs(sp.Symbol('lambda'),sp.Rational(3,4))!=0
    # The native weak sector is the spin1 representation. Its2pi action is
    # identity on the whole source fiber and hence on every End/tensor word.
    assert weak[2].eigenvals()=={-1:1,1:1,0:4}
    full_turn=(2*sp.pi*sp.I*weak[2]).exp()
    assert full_turn==sp.eye(6)
    # Existing source pulses genuinely preserve the recovered sectors. This
    # checks the full6level execution rather than a supplied logical gate.
    errors=[]
    for kind,angle in (('X',.23),('Y',-.31),('Z',.19),('J',.27),('H',-.13)):
        u=compiler.native_pulse(kind,angle)
        a=np.kron(u,u.conj())
        errors.extend(float(np.linalg.norm(a@np.array(p,complex)-np.array(p,complex)@a,2)) for p in projectors)
    assert max(errors)<2e-12
    artifact=out/'gauge_operator_multiplets.npz'
    np.savez_compressed(artifact,
        generators=np.array([np.array(g,complex) for g in generators]),
        color_Casimir=np.array(cc,complex),weak_Casimir=np.array(cw,complex),
        relative_charge=np.array(charge,complex),
        projectors=np.array([np.array(p,complex) for p in projectors]),
        current_projection=np.array(current_projection,complex))
    paths=[PUBLIC/'code/source_selection_model/response.json',PUBLIC/'code/source_selection_model/response.py',
           PUBLIC/'code/a5_closure/port_current_inner_certificate.py',
           PUBLIC/'code/a5_closure/response_grammar_completeness_certificate.py',
           PUBLIC/'model/observer/source_scalar/compiler.py',PUBLIC/'code/m1_source_realization/receipt.json',
           PUBLIC/'claims/axiom_registry.yaml',PUBLIC/'docs/AXIOM_REFERENCE.md',
           PUBLIC/'code/a5_closure/source_current_tomography_stage1.md',
           PUBLIC/'code/a5_closure/source_current_order_sensitive_inventory.py',
           PUBLIC/'code/particles/flavor/exterior_higgs_response_spec.json',
           PUBLIC/'code/sm_fermion_current/README.md',PUBLIC/'model/common_source/gauge_matter.py',
           PUBLIC/'code/yang_mills/manifests/physical_compact_gauge_uninstantiated.json']
    result={'status':'EXACT_SOURCE_RESPONSE_INTERNAL_GAUGE_MULTIPLETS_IDENTIFIED',
        'producer_sha256':sha(__file__),'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
        'artifacts':{artifact.name:sha(artifact)},'empirical_or_mass_inputs':[],
        'public_head':'bb469bd80432fa24672bbb605afc1d171e46fa9d',
        'source_reconstruction':{'probe_count_per_port':35,'ports':12,'work_dimension':6,
            'response_rank':12,'ordered_mixed_coefficients_checked':jets,
            'independent_geometric_formula_matches':True,'overall_scalar_phase_removed':True,
            'raw_primitive_repair_histories_generated_these_response_probes':False,
            'source_status':'The existing declared geometric response law and its retained algebraic response probes, also used by the six-level source compiler.'},
        'internal_current_identity':{'compact_Lie_algebra':'u(1)+su(3)+so(3), with so(3) locally isomorphic to su(2)',
            'simple_and_central_current_dimensions':[8,3,1],
            'full_operator_dimension_complex':36,'current_dimension_real':12,
            'Hermitian_operator_complement_dimension_real':24,
            'source_basis_coordinates_in_verified_canonical_basis':[[str(x) for x in row] for row in coefficients.tolist()],
            'coupling_or_hypercharge_normalization_selected':False,
            'central_charge_convention':'eigenvalue separation of the two source blocks set to1; a coordinate convention, not measured hypercharge'},
        'operator_multiplets':rows,'current_dimensions_in_same_order':current_dimensions,
        'color_preserving_order_parameter_test':{
            'color_singlet_operator_dimension':10,
            'relative_central_U1_action_on_every_color_singlet_operator_is_zero':True,
            'identity':'ad(Y) P_color_singlet = 0 on all36 operator coordinates',
            'conditional_consequence':'For order parameters O in this color-singlet operator domain, G_ab=Tr([T_a,O]^dagger[T_b,O]) has a zero central-U1 row and column. These operators cannot realize hypercharge-charged weak-doublet mixing through that orbit metric.',
            'dimensionless_triplet_control':'O=Jz in the verified weak current ideal',
            'triplet_orbit_Gram_diagonal':list(map(str,orbit_gram.diagonal())),
            'triplet_orbit_Gram_rank':int(orbit_gram.rank()),
            'source_selected_condensate_or_Higgs_action':False,
            'gauge_particle_mass_matrix_claimed':False},
        'weak_doublet_test':{'Casimir_three_quarters_eigenspace_dimension':0,
            'native_2pi_internal_weak_rotation_is_identity':True,
            'scope':'End(C6) and finite ordinary tensor/subrepresentation constructions from this specific integer-isospin fiber; not arbitrary charged, projective, nonlocal or separately admitted matter sectors.',
            'all_compiled_codes_or_architectural_transfer_services_exhausted':False,
            'excludes_Standard_Model_matter_in_full_OPH':False,
            'existing_conditional_exterior_doublet_producer_rejected':False},
        'source_execution':{'native_full6_pulses_tested':5,'operator_projectors_tested_per_pulse':6,
            'maximum_multiplet_invariance_error':max(errors),'new_physical_detector_or_propagation_Hamiltonian_compiled':False},
        'particle_attachment':{'internal_adjoint_multiplets_identified':True,
            'transverse_vector_propagating_poles_or_helicities_identified':False,
            'physical_photon_gluon_W_Z_identified':False,'Higgs_particle_identified':False,
            'existing_3plus1_field_description_discarded':False,
            'same_source_physical_current_Ward_and_pole_map_completed':False,
            'reason':'Internal adjoint response directions are candidate gauge-current quantum numbers, not spacetime vector polarizations. The source carrier response alone supplies neither the physical propagating gauge field nor its same-current spectral kernel. The declared exterior Higgs module is a separate attachment.'}}
    path=out/'gauge_particle_identification.json';path.write_text(json.dumps(public_receipt(result),indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'status':result['status'],'multiplets':rows,'native_invariance_error':max(errors)},indent=2))


if __name__=='__main__':main()
