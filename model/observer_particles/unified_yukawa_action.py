#!/usr/bin/env python3
"""One declared Higgs doublet coupled to the actual conditional M15 currents.

Derive invariant tensors, rather than entering a table of SM Yukawa vertices.
This is a conditional local action construction, not a selected vacuum, a
spatial chiral regulator, or an empirical mass calculation.
"""
from __future__ import annotations
from evidence_contract import public_receipt
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True
import numpy as np
import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRIOR = HERE
sys.path.insert(0, str(PRIOR))
import representation_identification as rep
import scalar_bilinear_route as pairs


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def mat(p): return s.SparseMatrix(*p['shape'], {(i,j):s.sympify(v) for i,j,v in p['entries']})
def clean(a): return s.SparseMatrix(a.applyfunc(s.simplify))
def zero(a): return all(v == 0 or s.expand(v) == 0 or s.simplify(v) == 0 for v in a.todok().values())
def eye(n): return s.SparseMatrix.eye(n)
def kron(a,b):
    return s.SparseMatrix(a.rows*b.rows,a.cols*b.cols,
        {(i*b.rows+k,j*b.cols+l):v*w for (i,j),v in a.todok().items()
         for (k,l),w in b.todok().items()})


def invariant_vertices():
    raw = json.loads((PRIOR / 'results/representation_identification_operators.json').read_text())
    scalar = json.loads((PRIOR / 'results/scalar_bilinear_operators.json').read_text())
    basis, dg = rep.fock_data()
    ix = raw['even_indices']
    g5 = [mat(x) for x in raw['spin_lifted_one_body_generators']]
    currents = [clean(dg(g).extract(ix, ix)) for g in g5]
    hcurrents = [g[3:,3:] for g in g5]
    lifted = [pairs.pair_lift(g, True)[0] for g in currents]
    pair_basis = scalar['symmetric_pair_basis']
    projector = mat(scalar['Higgs_charge_doublet_projector'])
    q = mat(raw['even_identification_operators']['q'])
    t3 = mat(raw['even_identification_operators']['T3'])
    result = []
    # Read all source doublet origins from the projector and charge operator.
    origins = sorted(set(tuple(sorted((q[a,a],q[b,b]))) for a,b in pair_basis))
    for origin in origins:
        selected = [i for i,(a,b) in enumerate(pair_basis)
                    if tuple(sorted((q[a,a],q[b,b]))) == origin]
        diagonal = s.diag(*[int(i in selected) for i in range(120)])
        projected = clean(diagonal * projector)
        if s.trace(projected) == 0:
            continue
        z = s.SparseMatrix.hstack(*projected.columnspace())
        assert z.cols == 2
        left = (z.H*z).inv()*z.H
        native_pair = [clean(left*g*z) for g in lifted]
        assert all(zero(z*small-g*z) for small,g in zip(native_pair,lifted))
        pair_charge = sum(origin)
        assert abs(pair_charge) == 3
        conjugate_higgs = bool(pair_charge > 0)
        higgs = [-g.T for g in hcurrents] if conjugate_higgs else hcurrents
        equations = s.Matrix.vstack(*[
            kron(g,eye(2))+kron(eye(2),h)
            for g,h in zip(native_pair,higgs)])
        kernel = equations.nullspace()
        assert len(kernel) == 1
        # Row-vectorization identity preserves all240 tensor equations without
        # materializing a dense symbolic240-by240 Kronecker operator.
        tensor_matrix = clean(z*s.Matrix(2,2,list(kernel[0])))
        tensor = s.SparseMatrix(240,1,lambda i,j:tensor_matrix[i//2,i%2])
        assert all(zero(g*tensor_matrix+tensor_matrix*h.T)
                   for g,h in zip(lifted,higgs))
        matrices = []
        for component in range(2):
            c = s.zeros(15)
            for p,(a,b) in enumerate(pair_basis):
                value = s.conjugate(tensor[2*p+component])
                c[a,b] = value if a == b else value/s.sqrt(2)
                c[b,a] = c[a,b]
            matrices.append(clean(c))
        # A neutral H=(0,v) fixes only a convention for the free y coefficient.
        vacuum = matrices[1]
        nonzero = [v for v in vacuum if v != 0]
        assert nonzero
        normalization = s.sqrt(s.conjugate(nonzero[0])*nonzero[0])
        matrices = [clean(c/normalization) for c in matrices]
        assert all(abs(complex(s.N(v))) == 1 for v in matrices[1] if v != 0)
        charge = t3+q/6
        assert zero(charge.T*matrices[1]+matrices[1]*charge)
        # Infinitesimal invariance of the actual covector in field variables.
        for g,h in zip(currents,higgs):
            for k in range(2):
                derivative = g.T*matrices[k]+matrices[k]*g
                derivative += sum((h[j,k]*matrices[j] for j in range(2)), s.zeros(15))
                assert zero(derivative)
        result.append({'origin':origin,'conjugate_higgs':conjugate_higgs,
                       'matrices':matrices,'tensor':tensor,
                       'higgs_currents':higgs,
                       'left_indices':[i for i in range(15) if q[i,i] == origin[0]],
                       'right_indices':[i for i in range(15) if q[i,i] == origin[1]],
                       'neutral_rank':matrices[1].rank()})
    assert len(result) == 3
    assert sorted(x['neutral_rank'] for x in result) == [2,6,6]
    # The two-component invariant epsilon makes every symmetric internal
    # coefficient into an antisymmetric physical CAR pairing kernel.
    eps = s.Matrix([[0,1],[-1,0]])
    pauli = (s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1))
    for row in result:
        for c in row['matrices']:
            b = kron(eps,c)
            assert b.T == -b
            for p in pauli:
                rotation = -s.I*p/2
                boost = p/2
                for j in (rotation,boost):
                    j30 = kron(s.SparseMatrix(j),eye(15))
                    assert zero(j30.T*b+b*j30)
    return result, currents, q, t3


def family_kernel(vertices, matrices, component=1):
    """Arbitrary complex family matrices, with symmetric all-left kernel."""
    k = s.zeros(45)
    for vertex,y in zip(vertices,matrices):
        c = vertex['matrices'][component]
        for a in vertex['left_indices']:
            for b in vertex['right_indices']:
                if c[a,b]:
                    for i in range(3):
                        for j in range(3):
                            k[3*a+i,3*b+j] += c[a,b]*y[i,j]
                            k[3*b+j,3*a+i] += c[a,b]*y[i,j]
    assert k.T == k
    return s.SparseMatrix(k)


def main():
    vertices,currents,q,t3 = invariant_vertices()
    print('Exact one-Higgs vertices and Lorentz/CAR identities: PASS',flush=True)
    # Two dimensionless controls with identical representation and Higgs norm.
    # No measured masses, CKM entries, physical time or empirical scales enter.
    a = [s.diag(1,2,3),s.diag(2,3,4),s.diag(3,4,5)]
    b = [s.diag(1,4,9),s.diag(2,6,12),s.diag(3,8,15)]
    # Check nonsymmetric complex family couplings, including all family units.
    # Gauge acts trivially on the declared family factor. Charge-oriented
    # blocks prevent a spurious color-dependent transpose of Y.
    family_checks = 0
    family_currents = [kron(g,eye(3)) for g in currents]
    for vertex in vertices:
        for i,j in itertools.product(range(3),repeat=2):
            y = s.zeros(3);y[i,j] = 1+s.I
            cs = [family_kernel([vertex],[y],k) for k in range(2)]
            for gf,h in zip(family_currents,vertex['higgs_currents']):
                for k in range(2):
                    derivative = gf.T*cs[k]+cs[k]*gf
                    derivative += sum((h[l,k]*cs[l] for l in range(2)),s.SparseMatrix.zeros(45))
                    assert zero(derivative)
                    family_checks += 1
    print(f'Exact family-current-component identities: {family_checks} PASS',flush=True)
    controls = []
    for name,ys in [('A',a),('B',b)]:
        k = family_kernel(vertices,ys)
        em = kron(t3+q/6,eye(3))
        assert zero(em.T*k+k*em)
        numerical = np.array(k,dtype=complex)
        singular = np.linalg.svd(numerical,compute_uv=False)
        expected = sorted([0.]*3 + [float(y[i,i]) for v,y in zip(vertices,ys)
                           for i in range(3) for _ in range(v['neutral_rank'])],reverse=True)
        assert len(expected) == 45
        assert np.max(abs(singular-expected)) < 1e-12
        controls.append({'name':name,'Yukawa_matrices':[[[str(x) for x in r] for r in y.tolist()] for y in ys],
                         'singular_values_in_supplied_unit_Higgs_norm':list(map(float,singular)),
                         'massless_neutrino_count_in_this_minimal_attachment':3,
                         'electromagnetic_pairing_identity_exact':True})
    assert controls[0]['singular_values_in_supplied_unit_Higgs_norm'] != controls[1]['singular_values_in_supplied_unit_Higgs_norm']
    # A change of Yukawa entries is a change of the Hamiltonian/reference,
    # not a proposed freedom at a fixed full modular clock.
    artifact = HERE / 'results/unified_yukawa_operators.json'
    artifact.write_text(json.dumps({'vertices':[{
        'constituent_q6Y':list(map(str,v['origin'])),
        'uses_conjugate_of_same_Higgs_doublet':v['conjugate_higgs'],
        'coefficient_matrices':[rep.sparse(c) for c in v['matrices']]} for v in vertices],
        'three_family_control_kernels':[rep.sparse(family_kernel(vertices,ys)) for ys in (a,b)]},indent=2)+'\n', encoding='utf-8', newline='\n')
    source_paths = [PRIOR/'representation_identification.py',PRIOR/'scalar_bilinear_route.py',
                   PRIOR/'results/representation_identification_operators.json',
                   PRIOR/'results/scalar_bilinear_operators.json',
                   ROOT/'Lean/Screen/WeylYukawaConventions.lean',
                   ROOT/'Lean/Screen/FermionSectorAssembly.lean',
                   HERE/'results/family_naturality.json']
    result = {'status':'PASS_EXACT_ONE_HIGGS_SOURCE_CURRENT_VERTICES_AND_SCOPED_FAMILY_ACTION',
        'producer_sha256':sha(__file__),
        'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in source_paths},
        'artifacts':{artifact.name:sha(artifact)},'empirical_inputs':[],
        'single_internal_family_dimension':15,'source_current_count':len(currents),
        'one_declared_Higgs_doublet':{'dimension':2,'Y':'1/2','same_field_and_its_conjugate_used':True,
            'elementary_field_kinetic_term_potential_or_composite_pole_selected':False},
        'invariant_vertices':[{'constituent_q6Y':list(map(str,v['origin'])),
            'uses_conjugate_Higgs':v['conjugate_higgs'],'invariant_tensor_dimension':1,
            'neutral_mass_kernel_rank':v['neutral_rank']} for v in vertices],
        'all12_source_current_vertex_identities_exact':True,
        'Lorentz_rotation_boost_and_CAR_antisymmetry_exact':True,
        'quadratic_action_convention':'L_pair=(1/2) epsilon_AB psi_A^T K(H) psi_B + h.c.; overall sign conventional. H=(0,v) uses neutral-component norm v=v_SM/sqrt(2), so kernel singular values are tree-level mass parameters under canonical kinetic normalization.',
        'three_declared_families_complex_Yukawa_coefficients':27,
        'exact_family_current_component_checks':family_checks,
        'independent_symmetry_allowed_controls':controls,
        'family_scope':'Three copies are declared. If both family factors carry the same unbroken irreducible source A5 action and Higgs is family-neutral, the previous naturality theorem instead forces each Y=yI. That conditional branch gives degeneracy, not the observed hierarchy.',
        'metric_handoff':'These local Lorentz/gauge scalar vertices can be multiplied by sqrt(-g) in the same spin/tetrad action. Their algebra adds no kinetic principal symbol. A source spatial chiral regulator, spin connection, Hamiltonian and vacuum must still be provided.',
        'parameter_scope':'Controls A/B change Yukawa action and its reference. They are not asserted identical complete OPH histories or alternatives at fixed full modular flow.',
        'physical_Standard_Model_masses_derived':False,
        'neutrino_masses_derived':False,
        'physical_source_selected_interacting_action':False}
    (HERE/'results/unified_yukawa_action.json').write_text(json.dumps(public_receipt(result),indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'status':result['status'],'vertices':result['invariant_vertices'],
                      'source_bindings':len(source_paths)},indent=2))


if __name__ == '__main__': main()
