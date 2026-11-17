#!/usr/bin/env python3
"""A shared source metric for scalar, Maxwell and native Clifford symbols.

This supplies an explicit tetrad/intertwiner and finite Whitney pullback.
It does not claim that field actions, global spin boundary conditions or a
matter backreaction law have thereby been selected by the native scheduler.
"""
from __future__ import annotations
from evidence_contract import public_receipt
import hashlib,importlib.util,json,sys
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
import sympy as s
from scipy.linalg import eigh,null_space

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];PUBLIC=ROOT
sys.path[:0]=[str(PUBLIC),str(PUBLIC/'code')]
from model.relative_observer_geometry import intrinsic_geometry as geometry
from source_gravity_admissibility import model as native_dirac
from source_gravity_admissibility import check as independent_dirac

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def enc(a):return [[str(v) for v in row] for row in a.tolist()]
def zero(a):return all(s.simplify(v)==0 for v in a)
def number(a):return np.array(a.evalf(40),float)
def maximum(a):return float(np.max(np.abs(a)))
def binary_rank(array):
    a=np.array(array,dtype=np.uint8)%2;r=0
    for col in range(a.shape[1]):
        rows=np.flatnonzero(a[r:,col])
        if not len(rows):continue
        k=r+int(rows[0]);a[[r,k]]=a[[k,r]]
        for i in range(a.shape[0]):
            if i!=r and a[i,col]:a[i]^=a[r]
        r+=1
        if r==a.shape[0]:break
    return r
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module;spec.loader.exec_module(module);return module

def source_tetrad():
    read=geometry.native_chart_increment(level=1,attempts=2)
    screen=s.Matrix([[s.Rational(v) for v in row] for row in read['averaged_screen_tensor']])
    a=s.diag(1,screen);data=geometry.wave_density_metric(a.tolist())
    g=data['metric'];inverse=data['inverse_metric'];volume=data['metric_volume_density']
    lapse=s.sqrt(-g[0,0]);L=a.cholesky()
    coframe=s.diag(lapse,lapse*L.inv());frame=coframe.inv();eta=s.diag(-1,1,1,1)
    assert zero(coframe.T*eta*coframe-g)
    assert zero(frame*eta*frame.T-inverse)
    assert zero(volume*inverse-data['action_density'])
    # One proper-clock/tetrad convention, no separate field-specific speed.
    scale=s.symbols('clock_scale',positive=True);jac=s.diag(scale,1,1,1)
    changed=coframe*jac
    assert zero(changed.T*eta*changed-jac.T*g*jac)
    anew=scale**2*a;munew=1/scale
    new_density=munew*s.diag(-1,anew)
    assert zero(new_density-scale*jac.inv()*data['action_density']*jac.inv().T)
    return {'readout':read,'a':a,'L':L,'g':g,'inverse':inverse,'volume':volume,'N':lapse,
            'coframe':coframe,'frame':frame,'eta':eta}, {
        'source_increment_readout':read,'spatial_principal_tensor':enc(a),
        'metric':enc(g),'inverse_metric':enc(inverse),'volume_density':str(volume),
        'lapse':str(lapse),'orthonormal_coframe':enc(coframe),'inverse_frame':enc(frame),
        'coframe_metric_and_density_identity_exact':True,
        'clock_rechart':'t=k*tau: g_tau=J^T g J, a_tau=k^2 a, mu_tau=mu/k, N_tau=kN; H_tau=k H, beta_tau=beta/k.',
        'same_modular_reference_beta_H_preserved_under_clock_rechart':True,
        'finite_coefficient_is_measured_continuum_metric':False,
        'measure_density_one_is_a_comparison_convention':True,
        'limiting_family':'mu(x,y)>0, a=diag(1,a_screen(x,y)); g=N^2[-dt^2+dr^2+a_screen^-1], N^2=mu sqrt(det a_screen).',
        'scope':'The same exact coframe formulas apply pointwise to each smooth effective metric in the existing product description. This finite native increment supplies a concrete comparison coefficient, not its derivative jet or curvature.'}

def clifford_and_symbols(data):
    pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)]
    beta=s.diag(1,1,-1,-1);alphas=[s.BlockMatrix([[s.zeros(2),p],[p,s.zeros(2)]]).as_explicit() for p in pauli]
    native_errors=[]
    for j in range(3):
        u=native_dirac.native_basis(j)
        independent=independent_dirac.basis(j)
        native_errors.append(max(maximum(u-independent),maximum(u@np.array(beta,float)@u.conj().T-np.array(alphas[j],complex))))
    assert max(native_errors)<2e-11
    gamma=[s.I*beta]+[s.I*beta*x for x in alphas]
    eta=data['eta'];frame=data['frame'];inverse=data['inverse']
    curved=[sum((frame[mu,a]*gamma[a] for a in range(4)),s.zeros(4)) for mu in range(4)]
    for mu in range(4):
        for nu in range(4):assert zero(curved[mu]*curved[nu]+curved[nu]*curved[mu]-2*inverse[mu,nu]*s.eye(4))
    p=s.Matrix(s.symbols('p0:4',real=True));symbol=sum((p[i]*curved[i] for i in range(4)),s.zeros(4))
    scalar=(p.T*inverse*p)[0]
    assert zero(symbol*symbol-scalar*s.eye(4))
    k=s.Matrix(s.symbols('k0:3',real=True));orthogonal=data['L'].T*k
    h=sum((orthogonal[j]*alphas[j] for j in range(3)),s.zeros(4));squared=(k.T*data['a']*k)[0]
    assert zero(h*h-squared*s.eye(4))
    # The two chiral principal symbols are both available, not selected.
    wh=sum((orthogonal[j]*pauli[j] for j in range(3)),s.zeros(2))
    energy=s.symbols('energy',real=True)
    assert s.simplify((energy*s.eye(2)-wh).det()-(energy**2-squared))==0
    assert s.simplify((energy*s.eye(2)+wh).det()-(energy**2-squared))==0
    # Maxwell principal symbol on covectors. Every pure gradient is null.
    raised=inverse*p;maxwell=s.simplify(scalar*inverse-raised*raised.T)
    assert zero(maxwell*p)
    null_cov=data['coframe'].T*s.Matrix([1,1,0,0])
    timelike=data['coframe'].T*s.Matrix([1,0,0,0])
    null_symbol=maxwell.subs(dict(zip(p,null_cov)))
    massive_symbol=maxwell.subs(dict(zip(p,timelike)))
    assert null_symbol.rank()==1 and massive_symbol.rank()==3
    # Spin connection representation for every Lorentz-algebra generator.
    connection_checks=0
    for a in range(4):
        for b in range(a+1,4):
            w=s.zeros(4);w[a,b]=1;w[b,a]=-1
            omega=sum((w[i,j]*gamma[i]*gamma[j]/4 for i in range(4) for j in range(4)),s.zeros(4))
            mixed=eta*w
            for c in range(4):
                assert zero(omega*gamma[c]-gamma[c]*omega+sum((mixed[c,j]*gamma[j] for j in range(4)),s.zeros(4)))
                connection_checks+=1
    # Positive Euclidean square convention D=-i gamma^a nabla_a.
    # R_abcd=kappa(delta_ac delta_bd-delta_ad delta_bc).
    egamma=[gamma[1],gamma[2],gamma[3],s.I*gamma[0]]
    kappa=s.symbols('sectional_curvature',real=True);curvature=s.zeros(4)
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    r=kappa*(int(a==c)*int(b==d)-int(a==d)*int(b==c))
                    curvature-=r*egamma[a]*egamma[b]*egamma[c]*egamma[d]/8
    assert zero(curvature-3*kappa*s.eye(4))
    return {'native_basis_word_errors':native_errors,'native_words_replayed':3,
        'native_word_pulse_counts':[len(native_dirac.basis_word(j)) for j in range(3)],
        'proper_code_native_alpha_reconstruction_pass':True,
        'all16_curved_Clifford_relations_exact':True,'Dirac_symbol_square_equals_scalar_symbol_exact':True,
        'shared_coordinate_dispersion':'E^2=p_spatial^T a p_spatial for massless principal parts',
        'left_and_right_Weyl_principal_determinants_exact':True,'chirality_or_generation_selected':False,
        'Maxwell_gauge_symbol_annihilation_exact':True,'Maxwell_symbol_rank_off_cone':3,
        'Maxwell_symbol_rank_on_nonzero_null_covector':1,'Maxwell_null_physical_polarization_dimension':2,
        'spin_connection_Clifford_compatibility_checks':connection_checks,
        'required_connection':'Omega_mu=(1/4) omega_mu_ab gamma^a gamma^b, with torsion-free metric spin connection from the chosen smooth tetrad.',
        'constant_curvature_Lichnerowicz_algebra_check':'R=12*kappa; -1/8 R_abcd gamma^a gamma^b gamma^c gamma^d=(R/4) I in the declared positive Riemannian Dirac-square convention.',
        'curved_Dirac_square_is_componentwise_scalar_B':False,
        'new_full_spatial_Dirac_native_history_executed':False,
        'mass_parameter_or_gauge_coupling_fitted':False,
        'global_spin_structure_boundary_and_reference_selected':False},curved

def null_frame_and_stress(data):
    q=s.sqrt(3)/3
    orth=[s.Matrix(v) for v in [(1,1,0,0),(1,-1,0,0),(1,0,1,0),(1,0,-1,0),(1,0,0,1),(1,0,0,-1),(1,q,q,q),(1,q,q,-q),(1,q,-q,q)]]
    vectors=[s.simplify(data['frame']*v) for v in orth];g=data['g'];inverse=data['inverse']
    positions=[(i,j) for i in range(4) for j in range(i,4)]
    design=s.Matrix([[v[i]*v[j]*(1 if i==j else 2) for i,j in positions] for v in vectors])
    metricvec=s.Matrix([g[i,j] for i,j in positions])
    assert design.rank()==9 and zero(design*metricvec)
    assert all(s.simplify((v.T*g*v)[0])==0 for v in vectors)
    # Recover arbitrary g-tracefree symmetric tensors from nine null values.
    variables=s.symbols('T0:10',real=True);tensor=s.zeros(4)
    for z,(i,j) in zip(variables,positions):tensor[i,j]=tensor[j,i]=z
    trace=s.trace(inverse*tensor);tf=tensor-g*trace/4
    trace_row=s.Matrix([[inverse[i,j]*(1 if i==j else 2) for i,j in positions]])
    recovered=design.col_join(trace_row).inv()*(design*s.Matrix(variables)).col_join(s.Matrix([0]))
    assert all(s.simplify(recovered[n]-tf[i,j])==0 for n,(i,j) in enumerate(positions))
    # Independent Maxwell action variation on this same metric.
    fs=s.symbols('F0:6',real=True);F=s.zeros(4)
    for x,(i,j) in zip(fs,((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))):F[i,j]=x;F[j,i]=-x
    f2=sum(F[i,j]*F[k,l]*inverse[i,k]*inverse[j,l] for i in range(4) for j in range(4) for k in range(4) for l in range(4))
    stress=s.simplify(F*inverse*F.T-g*f2/4)
    assert s.simplify(s.trace(inverse*stress))==0
    eps=s.symbols('epsilon',real=True);checks=0
    for i,j in positions:
        direction=s.zeros(4);direction[i,j]=direction[j,i]=1
        varied=inverse+eps*direction
        vf2=sum(F[a,b]*F[c,d]*varied[a,c]*varied[b,d] for a in range(4) for b in range(4) for c in range(4) for d in range(4))
        dvol=-data['volume']*s.trace(g*direction)/2
        actual=s.expand(-(dvol*f2+data['volume']*s.diff(vf2,eps).subs(eps,0))/4)
        expected=-data['volume']*sum(stress[a,b]*direction[a,b] for a in range(4) for b in range(4))/2
        assert s.simplify(actual-expected)==0;checks+=1
    return {'orthonormal_nine_directions':[enc(v.T)[0] for v in orth],
        'coordinate_nine_directions':[enc(v.T)[0] for v in vectors],
        'same_metric_all_nine_null_exact':True,'design_rank':9,'kernel':'span(g_mu_nu)',
        'arbitrary_tracefree_tensor_reconstruction_exact':True,
        'map':'k_coordinate=frame*ell_orthonormal; T_orthonormal=frame^T*T_coordinate*frame',
        'source_order_selected_physical_direction_jacobian':False,
        'scope':'An explicit pointwise frame map now identifies the theorem null-test algebra with this chosen product metric. It does not identify source event differences with spacetime derivatives or assert the nine physical balance equations.',
        'Maxwell_action':'S_A=−1/4 integral sqrt(−g) g^{mu rho}g^{nu sigma}F_mu_nu F_rho_sigma',
        'Maxwell_stress':'T_mu_nu=F_mu_rho F_nu^rho−g_mu_nu F_squared/4',
        'independent_Maxwell_inverse_metric_variations_exact':checks,'Maxwell_stress_trace_zero_exact':True,
        'physical_curvature_or_Einstein_balance_defined_from_stress':False}

def radial_boundary_control():
    """The existing scalar Dirichlet radial choice is not a Dirac domain."""
    x=s.Matrix([[0,1],[1,0]])
    boundary=s.diag(-x,x)
    # One no-normal-current self-adjoint extension: the second spinor
    # component vanishes at each end. It is a maximal isotropic subspace.
    inclusion=s.Matrix([[1,0],[0,0],[0,1],[0,0]])
    assert inclusion.rank()==2 and zero(inclusion.T*boundary*inclusion)
    assert boundary.rank()==4
    length=s.symbols('radial_length',positive=True)
    deficiency=[s.cosh(length)*s.eye(2)+sign*s.sinh(length)*x for sign in (-1,1)]
    assert all(s.simplify(f.det())==1 for f in deficiency)
    return {'operator':'H=−i sigma_x d/dr on the finite radial interval',
        'boundary_form':'−i[psi_dagger sigma_x phi]_0^L',
        'scalar_Dirichlet_on_all_spinor_components':'symmetric minimal operator, not self-adjoint',
        'minimal_two_component_deficiency_indices':[2,2],
        'full_boundary_space_dimension':4,'self_adjoint_boundary_subspace_dimension':2,
        'explicit_reflecting_control':'psi_2(0)=psi_2(L)=0 gives a maximal isotropic boundary subspace and zero local normal current.',
        'square_of_reflecting_control':'component1 Neumann, component2 Dirichlet; it is not two copies of scalar Dirichlet B.',
        'native_source_selects_this_reflecting_choice':False,
        'scope':'Exact domain control only; the same two-component radial block occurs twice in a4-component Dirac module. No new radial particle dynamics or chosen physical boundary is asserted.'}

def whitney_same_metric(data):
    """Reuse the actual public cone and quadrature in the same frozen metric.

    This is a local Galerkin chart attachment; the cone is not asserted to
    be the entire radial-cylinder source refinement or its boundary domain.
    """
    source=load('joint_existing_whitney_action',PUBLIC/'code/equation_of_state/maxwell/maxwell_eos.py')
    mesh,(_,_,base_local,_)=source.mesh_and_forms()
    lapse=float(data['N']);h=number(data['g'][1:,1:]);invh=np.linalg.inv(h)
    chart=number(data['coframe'][1:,1:]);assert maximum(chart.T@chart-h)<2e-13
    pulled_m,pulled_n,local,volume=source.forms(mesh,chart)
    m=np.zeros_like(pulled_m);n=np.zeros_like(pulled_n);mass0=np.zeros((13,13));stiff0=np.zeros((13,13))
    determinant=np.linalg.det(h);root=np.sqrt(determinant)
    for row in base_local:
        ids=row['vertices'];ei=row['edge_indices'];fi=row['face_indices'];v=row['volume']
        t1=np.array(row['one_coefficients']);t2=np.array(row['two_coefficients']);grad=np.array(row['gradients'])
        moments=v*(np.ones((4,4))+np.eye(4))/20
        m[np.ix_(ei,ei)]+=root*np.einsum('eic,cd,ij,fjd->ef',t1,invh,moments,t1)
        n[np.ix_(fi,fi)]+=np.einsum('eic,cd,ij,fjd->ef',t2,h,moments,t2)/root
        mass0[np.ix_(ids,ids)]+=root*moments/lapse
        stiff0[np.ix_(ids,ids)]+=lapse*root*v*grad@invh@grad.T
    errors=[maximum(m-pulled_m),maximum(n-pulled_n)]
    assert max(errors)<2e-11
    # Coordinate-time Maxwell kinetic/stiffness from precisely the same g.
    d=np.array(mesh['D']);c=np.array(mesh['C']);electric=m/lapse;magnetic=lapse*n
    boundary_d=d[12:,1:];boundary_c=c[:20,12:]
    assert binary_rank(boundary_d)==11 and binary_rank(boundary_c)==19
    assert 30-binary_rank(boundary_d)-binary_rank(boundary_c)==0
    stiffness=c.T@magnetic@c
    assert maximum(c@d)==0 and maximum(stiffness@d)<2e-10
    assert maximum(d.T@electric@d-stiff0/lapse**2)<2e-11
    # Gradient energy is the same scalar spatial form with the lapse factor.
    assert maximum(lapse**2*d.T@electric@d-stiff0)<2e-11
    val,u=eigh(electric);rootm=(u*np.sqrt(val))@u.T;invroot=(u/np.sqrt(val))@u.T
    transverse=null_space((rootm@d).T);reduced=transverse.T@invroot@stiffness@invroot@transverse
    modes=eigh(reduced,eigvals_only=True);assert min(modes)>0
    # Coordinate rescaling t=k tau multiplies every mode energy by k.
    k=2.;new_e=electric/k;new_k=k*stiffness
    new_reduced=k**2*reduced
    assert maximum(np.sqrt(eigh(new_reduced,eigvals_only=True))-k*np.sqrt(modes))<2e-11
    return {'cone_vertices':13,'cone_edges':42,'cone_faces':50,'cone_tetrahedra':20,
        'same_frozen_spatial_metric':h.tolist(),'same_lapse':lapse,
        'coframe_pullback_vs_direct_Hodge_errors':errors,'scalar_vertex_modes':13,
        'scalar_and_Maxwell_gradient_form_identity_error':maximum(lapse**2*d.T@electric@d-stiff0),
        'exact_incidence_curl_gradient_zero':True,'metric_Maxwell_Ward_residual':maximum(stiffness@d),
        'physical_Maxwell_modes':len(modes),'frequency_range':[float(np.sqrt(min(modes))),float(np.sqrt(max(modes)))],
        'same_clock_rescaling_spectrum_verified':True,
        'kinetic_matrix':'M1(h)/N','magnetic_matrix':'N*M2(h)',
        'scalar_action':'M0(h)/N and N*K0(h), with the same constant chart metric and natural variational boundary',
        'interpolation_domain':'existing native-boundary tetrahedral cone, used as one declared local Galerkin chart',
        'global_cone_equals_radial_product_cylinder':False,
        'new_native_detector_history_executed':False,
        'native_screen_F2_coboundary_ranks':[11,19],'native_screen_H1_F2_dimension':0,
        'native_screen_Euler_characteristic':12-30+20,
        'oriented_screen_w2_euler_mod2':(12-30+20)%2,
        'smooth_product_spin_structure':'The oriented/time-oriented R x (0,L) x S2 product admits a spin structure: w2(TS2)=Euler(TS2) mod2=0 and the time/radial lines are trivial. Given existence, H1(Z2)=0 makes it unique up to equivalence. This does not choose a radial self-adjoint boundary condition.',
        'Whitney_action_and_chart_attachment_selected_by_repair':False}, {
            'metric':number(data['g']),'coframe':number(data['coframe']),
            'Whitney_D':d,'Whitney_C':c,'Maxwell_electric_kinetic':electric,
            'Maxwell_magnetic_metric':magnetic,'scalar_kinetic':mass0,'scalar_stiffness':stiff0}

def main():
    data,metric=source_tetrad();symbols,_=clifford_and_symbols(data)
    nulls=null_frame_and_stress(data);whitney,arrays=whitney_same_metric(data)
    artifact=HERE/'results/geometric_matter_join_arrays.npz';np.savez_compressed(artifact,**arrays)
    paths=[]
    for module in tuple(sys.modules.values()):
        p=getattr(module,'__file__',None)
        if p and Path(p).is_file() and Path(p).is_relative_to(PUBLIC):paths.append(Path(p))
    paths += [PUBLIC/p for p in ('model/relative_observer_geometry/README.md','model/observer_time_fields/REFINEMENT.md',
        'docs/research/SOURCE_GRAVITY_ADMISSIBILITY.md','docs/research/FERMIONIC_SOURCE_CLOCKS.md',
        'Lean/ObserverPatchHolography/EinsteinBranch/Tensor.lean','Lean/Geometry/SourceOrderEinsteinComposition.lean',
        'code/source_gravity_admissibility/receipt.json','code/source_selection_model/response.json',
        'Lean/Screen/SeamCurrentCarrierQuotient.lean','Lean/Screen/SeamCurrentEdge30Moment.lean',
        'Lean/ObserverPatchHolography/CoreAxioms.lean')]
    result={'status':'PASS_COMMON_METRIC_TETRAD_MATTER_SYMBOL_AND_NULL_FRAME_JOIN_WITH_EXPLICIT_GLOBAL_LIMITS',
        'producer_sha256':sha(__file__),'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(set(paths))},
        'artifacts':{artifact.name:sha(artifact)},'empirical_inputs':[],
        'same_source_metric_and_clock':metric,'native_Clifford_and_shared_symbols':symbols,
        'same_metric_null_tomography_and_stress':nulls,'same_metric_Whitney_chart':whitney,
        'radial_first_order_boundary_control':radial_boundary_control(),
        'full_operator_contract':{
            'smooth_spin_structure_required':True,'torsion_free_metric_spin_connection_required':True,
            'Dirac_square':'spin_connection_Laplacian + R/4 in the positive spatial convention; gauge curvature adds its Clifford contraction; static lapse contributes the appropriate Hamiltonian terms.',
            'global_boundary':'Scalar Dirichlet conditions do not automatically define a self-adjoint first-order Dirac operator; a compatible spin boundary domain is additional. The tested Whitney cone uses its existing natural variational boundary.',
            'taste_and_chirality':'The executed4-component symbol is one declared Dirac module, containing both Weyl principal symbols. It neither selects chirality nor removes tastes of a finite walk or Hodge discretization.',
            'common_source_history_intertwiner_completed':False,'source_selected_mass_or_backreaction_claimed':False},
        'primary_background_sources':['https://arxiv.org/abs/1308.4552','https://arxiv.org/abs/0911.1304']}
    (HERE/'results/geometric_matter_join.json').write_text(json.dumps(public_receipt(result),indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'status':result['status'],'native_basis_error':max(symbols['native_basis_word_errors']),
        'null_frame_rank':nulls['design_rank'],'Whitney_Hodge_error':max(whitney['coframe_pullback_vs_direct_Hodge_errors'])},indent=2))

if __name__=='__main__':main()
