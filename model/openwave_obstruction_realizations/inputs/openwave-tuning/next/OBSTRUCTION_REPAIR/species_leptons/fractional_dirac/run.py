"""Full one-particle covariance, both bands, source current and complete small-CAR controls."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations,product
import hashlib
import json
import numpy as np
from scipy.linalg import expm,block_diag
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
import source as S

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[4]
paths=[HERE/'PLAN.md',HERE/'source.py',HERE.parent/'COMMON_CLOCK.json',
 ROOT/'openwave-tuning/common_dynamics/color_fractional_winding/PROOF.md',
 ROOT/'openwave-tuning/common_dynamics/autonomous_dirac/PROOF.md',
 ROOT/'openwave-push/b/D3_filled_reference_particle/consumed/ROW_23/rc_lib.py',
 HERE.parent/'auxiliary_current/inputs/local_auxiliary_interface.py',
 HERE.parent/'auxiliary_current/MANIFEST.json']
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();pins={str(p.relative_to(ROOT)):sha(p)for p in paths}
clock=json.loads(paths[2].read_text());mass=clock['register_432_midpoints'][1];h=1/8;r=1
I4=np.eye(4)
for i in range(3):
    assert np.max(abs(S.ALPHA[i]@S.BETA+S.BETA@S.ALPHA[i]))==0
    for j in range(3):assert np.max(abs(S.ALPHA[i]@S.ALPHA[j]+S.ALPHA[j]@S.ALPHA[i]-2*(i==j)*I4))==0
assert np.max(abs(S.BETA@S.BETA-I4))==0
# Actual D3 Clifford ordering differs only by the two-qubit swap.
swap=np.zeros((4,4))
for d,sp in product(range(2),repeat=2):swap[2*d+sp,2*sp+d]=1
for a,old in zip(S.ALPHA,[np.kron(q,S.SIGMA[0])for q in S.SIGMA]):assert np.max(abs(a-swap@old@swap.T))==0
assert np.max(abs(S.BETA-swap@np.kron(np.eye(2),S.SIGMA[2])@swap.T))==0

T=[]
for a,b in [(0,1),(0,2),(1,2)]:
    x=np.zeros((3,3),complex);x[a,b]=x[b,a]=.5;T.append(x)
    y=np.zeros((3,3),complex);y[a,b]=-.5j;y[b,a]=.5j;T.append(y)
T.extend([np.diag([.5,-.5,0]).astype(complex),np.diag([1,1,-2]).astype(complex)/(2*np.sqrt(3))])
vec,ax,left=S.current_matrices();charge_error=0.;colour_error=0.
for J in [*vec,*ax,*left]:
    charge_error=max(charge_error,float(np.max(abs(S.Q@J-J@S.Q-J))))
    for t in T:
        C=np.kron(np.eye(8),t);colour_error=max(colour_error,float(np.max(abs(C@J-J@C))))
assert max(charge_error,colour_error)<1e-15
# Both-band frame completeness and actual unitarily transformed currents.
bands=[]
for k in [np.array(q,float)for q in [(0,0,0),(.3,.4,.5),(1.2,-.7,.9),(np.pi,np.pi,0)]]:
    mw=h*mass+r*np.sum(1-np.cos(k));H=mw*S.BETA+sum((np.sin(k[a])*S.ALPHA[a]for a in range(3)),np.zeros((4,4),complex))
    e,U=np.linalg.eigh(H);E=np.sqrt(mw*mw+np.sum(np.sin(k)**2))
    assert np.max(abs(e-np.array([-E,-E,E,E])))<1e-13
    assert np.max(abs(U@U.conj().T-I4))<1e-13
    UU=np.kron(U,np.eye(6));Jall=[]
    for J in [*vec,*ax,*left]:
        jj=UU.conj().T@J@UU
        assert np.max(abs(UU@jj@UU.conj().T-J))<1e-12
        Jall.append(float(np.linalg.norm(jj[:12,12:])))
    bands.append(dict(k=k.tolist(),signed_energies=e.tolist(),positive_negative_current_block_norms=Jall))

# Same physical charge unit on the actual U3 quotient representations.
theta=.37;w=np.exp(-1j*theta/3)*np.eye(3);rw=S.reps(w)
assert np.max(abs(rw-block_diag(np.exp(2j*theta/3)*np.eye(3),np.exp(-1j*theta/3)*np.eye(3))))<1e-15
assert abs(np.linalg.det(w)-np.exp(-1j*theta))<1e-15
centre=np.exp(2j*np.pi/3)*np.eye(3)
for q in [F(2,3),F(-1,3)]:assert np.max(abs(np.exp(2j*np.pi*float(q))*centre-np.eye(3)))<1e-14

# Full192-state covariant background control with all24 independent links.
cells,edges=S.edges(2)
links=[np.exp(-1j*(e+1)/211/3)*expm(1j*(e+1)/97*T[e%8])for e in range(len(edges))]
gauge=[np.exp(-1j*(x+1)/17/3)*expm(1j*(x+1)/31*(T[x%8]+T[(x+3)%8]))for x in range(len(cells))]
changed=[gauge[x]@W@gauge[y].conj().T for (x,y,a),W in zip(edges,links)]
H=S.one_particle(mass,h,links);H2=S.one_particle(mass,h,changed);GG=block_diag(*[S.local_gauge(g)for g in gauge])
cov=float(np.max(abs(H2-GG@H@GG.conj().T)));assert cov<1e-12
assert np.max(abs(H-H.conj().T))<1e-14
fullQ=np.kron(np.eye(8),S.Q);assert np.max(abs(fullQ@H-H@fullQ))<1e-13
weak_cov=0.
for x,g in enumerate(gauge):
    for J in [*vec,*ax,*left]:
        weak_cov=max(weak_cov,float(np.max(abs(S.local_gauge(g)@J@S.local_gauge(g).conj().T-J/np.linalg.det(g)))))
assert weak_cov<1e-12
evals,U=np.linalg.eigh(H);assert np.max(abs(U@U.conj().T-np.eye(192)))<1e-12
# A separate Clifford grading anticommutes with every beta/alpha term, even
# at arbitrary unitary links. Positive Wilson mass gives a uniform gap.
grading=np.kron(np.eye(8),np.kron(np.kron(S.SIGMA[1],np.eye(2)),np.eye(6)))
assert np.max(abs(grading@H+H@grading))<1e-13
assert min(abs(evals))>=h*mass-1e-12

# The free L4 symbol retains all64 momenta, all positive/negative branches.
free=[]
for rr in [0,1]:
    vals=[]
    for j in product(range(4),repeat=3):
        k=2*np.pi*np.array(j)/4;mw=h*mass+rr*np.sum(1-np.cos(k))
        vals.append(float(np.sqrt(mw*mw+np.sum(np.sin(k)**2))))
    minima=sum(abs(e-h*mass)<1e-12 for e in vals)
    assert minima==(8 if rr==0 else 1)
    free.append(dict(Wilson_r=rr,all64_positive_cell_energies=vals,minima_count=minima,
                     negative_branches_retained=True,degeneracy_per_sign_per_momentum=12))

# Auxiliary operator bridge: local upper Dirac sector exactly recovers all
# old colour-spin and V/A operators. The kinetic-off N3 restriction is a
# local operator control, never the vacuum of the full new source.
embed=np.zeros((24,12));embed[:12,:]=np.eye(12)
aux_error=0.
for a,t in enumerate(T):
    for i in range(3):
        new=np.kron(S.SPIN[i],np.kron(np.eye(2),t))
        old=np.kron(np.kron(S.SIGMA[i]/2,np.eye(2)),t)
        aux_error=max(aux_error,float(np.max(abs(embed.T@new@embed-old))))
for J,old in zip([vec[0],*ax[1:]],[np.kron(np.kron(q,S.TAUP),np.eye(3))for q in [np.eye(2),*S.SIGMA]]):
    aux_error=max(aux_error,float(np.max(abs(embed.T@J@embed-old))))
assert aux_error==0

# Complete small CAR control:8 quark modes + neutral/charged lepton references.
small=S.open_small(mass,h);e,Uq=np.linalg.eigh(small)
assert sum(e<0)==4 and sum(e>0)==4
charges8=np.tile([2,-1,2,-1],2) # units of1/3 of the same electron charge
Q8=np.diag(charges8);assert np.max(abs(Q8@small-small@Q8))==0
pos=e>0;neg=~pos
reference_charge=float(np.trace(Uq[:,neg].conj().T@Q8@Uq[:,neg]).real)
assert abs(reference_charge-2)<1e-12
# Resolve the charge degeneracy first, making particle/hole charge accounting
# independent of arbitrary eigenvectors in degenerate energy subspaces.
charge_bands=[]
e8=[];u8=[];q8=[]
for q in [2,-1]:
    ids=np.flatnonzero(charges8==q);ee,uu=np.linalg.eigh(small[np.ix_(ids,ids)])
    assert sum(ee<0)==sum(ee>0)==2
    charge_bands.append(dict(quark_charge_thirds=q,particle_charges=[q]*2,antiparticle_charges=[-q]*2,
                             all_signed_cell_energies=ee.tolist()))
    for j,z in enumerate(ee):
        vv=np.zeros(8,complex);vv[ids]=uu[:,j];e8.append(z);u8.append(vv);q8.append(q)

H10=np.zeros((10,10),complex);H10[:8,:8]=small
q10=np.r_[charges8,0,-3]
Hcar=S.second_quantize(H10);Qcar=S.second_quantize(np.diag(q10))
assert (Qcar@Hcar-Hcar@Qcar).nnz==0
# Current PL tau+ in both cells; colour/spin are the stated invariant slice.
j8=block_diag(np.kron((np.eye(2)-S.SIGMA[0])/2,S.TAUP),np.kron((np.eye(2)-S.SIGMA[0])/2,S.TAUP))
j10=np.zeros((10,10),complex);j10[:8,:8]=j8
l10=np.zeros((10,10),complex);l10[9,8]=1
Jcar=S.second_quantize(j10);Lcar=S.second_quantize(l10)
weak=Jcar@Lcar+Jcar.conj().T@Lcar.conj().T
assert (weak-weak.conj().T).nnz==0
assert (Qcar@weak-weak@Qcar).nnz==0
assert (Qcar@Jcar-Jcar@Qcar-3*Jcar).nnz==0
assert (Qcar@Lcar-Lcar@Qcar+3*Lcar).nnz==0
# Complete particle/hole CAR operator reconstruction on all1024 occupations.
def annihilator(mode):
    rr=[];cc=[];vv=[]
    for w in range(1024):
        if(w>>mode)&1:
            rr.append(w^(1<<mode));cc.append(w);vv.append((-1)**((w&((1<<mode)-1)).bit_count()))
    return sparse.csr_matrix((vv,(rr,cc)),shape=(1024,1024),dtype=complex)
cs=[annihilator(j)for j in range(10)];I1024=sparse.eye(1024,format='csr')
gam=[sum((vv[j].conjugate()*cs[j]for j in range(8)),sparse.csr_matrix((1024,1024),dtype=complex))for vv in u8]
bd=[g if energy>0 else g.conj().T for g,energy in zip(gam,e8)]
max_sparse=lambda M:float(np.max(abs(M.data),initial=0))
carerror=0.
for i,a in enumerate(bd):
    for j,b in enumerate(bd):
        carerror=max(carerror,max_sparse(a@b.conj().T+b.conj().T@a-(I1024 if i==j else 0*I1024)),max_sparse(a@b+b@a))
assert carerror<1e-12
Evac=sum(z for z in e8 if z<0);Qvac=sum(q for q,z in zip(q8,e8)if z<0)
Hre=sum((abs(z)*(a.conj().T@a)for a,z in zip(bd,e8)),Evac*I1024)
Qre=sum(((q if z>0 else -q)*(a.conj().T@a)for a,z,q in zip(bd,e8,q8)),Qvac*I1024)
Qquark=S.second_quantize(np.diag(np.r_[charges8,0,0]))
herr=max_sparse(Hre-Hcar);qerr=max_sparse(Qre-Qquark)
assert max(herr,qerr)<1e-12
filled=[u for u,z in zip(u8,e8)if z<0];filled=np.column_stack(filled)
vac=np.zeros(1024,complex)
for occ in combinations(range(8),4):vac[sum(1<<j for j in occ)]=np.linalg.det(filled[list(occ),:])
assert abs(np.vdot(vac,vac)-1)<1e-12
assert max(np.linalg.norm(a@vac)for a in bd)<1e-12
assert np.linalg.norm(Hcar@vac-Evac*vac)<1e-12
incoming=cs[8].conj().T@vac
evolutions=[]
for coupling in [0,h/64]:
    for celltime in [0.,1.,4.]:
        final=expm_multiply(-1j*celltime*(Hcar+coupling*weak),incoming)
        prob=abs(final)**2;assert abs(prob.sum()-1)<1e-12
        assert np.linalg.norm((Qcar-2*I1024)@final)<1e-12
        evolutions.append(dict(cell_time=celltime,physical_source_time=h*celltime,
                 physical_weak_coupling=0 if coupling==0 else 1/64,
                 all1024_occupation_probabilities=prob.tolist(),
                 complete_probability_sum=float(prob.sum())))
# Full1024 occupation/flux/reference map on this open two-cell U1 control.
records=[];hops=0;undressed_fail=0
for wocc in range(1024):
    leftq=sum(q10[j]for j in range(4)if(wocc>>j)&1)
    rightq=sum(q10[j]for j in range(4,10)if(wocc>>j)&1)
    flux=leftq;reference=leftq+rightq
    assert leftq-flux==0 and rightq+flux-reference==0
    records.append([wocc,int(leftq),int(rightq),int(flux),int(reference)])
coo=Hcar.tocoo()
for a,b,z in zip(coo.row,coo.col,coo.data):
    if a==b or abs(z)==0:continue
    before=records[b];after=records[a];shift=after[1]-before[1]
    if shift:
        hops+=1;assert after[3]==before[3]+shift and after[4]==before[4]
        assert after[1]-before[3]!=0;undressed_fail+=1
assert hops>0 and hops==undressed_fail
# Every weak matrix element preserves total charge. With the lepton reference
# at cell1, a contact at cell0 requires its corresponding flux shift.
wc=weak.tocoo();weak_flux_changes=0
for a,b,z in zip(wc.row,wc.col,wc.data):
    if abs(z)==0:continue
    assert records[a][4]==records[b][4]
    weak_flux_changes+=records[a][3]!=records[b][3]

out=dict(schema='oph.fractional_wilson_dirac_service.v1',input_sha256=pins,
 source=dict(components_per_cell=24,ordering=['Dirac','physical_spin','flavour','colour'],
    common_mass_coefficient_exact=clock['register_432_midpoints_exact'][1],h='1/8',Wilson_r=1,
    physical_quark_charges=['2/3','-1/3'],physical_electron_charge='-1',
    gauge_group='(SU3 x R)/<(zI,2pi)> isomorphic to U3',
    link_representations='up det(W)^-1 W; down W; electron det(W)',
    clock='Hcell=h Hphys; all masses, colour/auxiliary energies and currents have one declared raw source unit',
    field_change='full compact U3 links; not the old193-mode Maxwell carrier',
    central_electric_coefficient_physical='(a/2) Q_link², a=1/2; Q_link generates the same physicaltheta',
    magnetic_coefficient_physical='3 b (3-ReTr W_plaquette), b=3; central small-field term b F²/2',
    colour_electric_coefficient_physical='kappa C_SU3, kappa=1/2,1,2',
    auxiliary_coefficients_physical='g=1/8, Delta=3kappa; multiplied byh inHcell',
    free_filled_reference_background_per_cell='colourless charge-2, representation det(g)^2',
    physical_normal_ordered_charge='dGamma(Q)-2 percell; particleq and hole-q',
    both_bands_and_particle_antiparticle_sectors=True),
 exact_D3_Clifford_swap_intertwiner=True,full_one_particle_dimension=192,
 gauge_covariance_error=cov,local_weak_current_covariance_error=weak_cov,
 charge_raising_commutator_error=charge_error,colour_current_commutator_error=colour_error,
 all192_background_signed_cell_energies=evals.tolist(),both_band_current_controls=bands,
 free_L4_controls=free,auxiliary_upper_component_operator_error=aux_error,
 old_auxiliary_vacuum_transferred=False,
 complete_CAR_control=dict(quark_modes=8,charge_reference_modes=2,full_dimension=1024,
    scope='identity-link open z chain, one physical spin and one colour; fullDirac/flavour/occupation complements retained',
    free_filled_reference_charge_thirds=reference_charge,charge_resolved_bands=charge_bands,
    particle_hole_CAR_error=carerror,particle_hole_H_reconstruction_error=herr,
    particle_hole_charge_reconstruction_error=qerr,complete_contact_evolutions=evolutions,
    full_charge_Ward_identity=True,charged_current_Ward_identity=True,
    every_occupation_flux_reference_record=records,dressed_hop_columns=hops,
    undressed_hop_Gauss_failures=undressed_fail,weak_contact_flux_changing_entries=int(weak_flux_changes),
    interacting_ground_state_claimed=False),
 resource_contract=dict(cells=8,quark_CAR_modes=192,matter_Fock_dimension='2^192',
    auxiliary_dimension='25^8',oriented_links=24,link_Hilbert_space='L2(U3)^tensor24, every Peter-Weyl irrep retained',
    one_particle_logical_code_qubits=8,physical_code_states=192,padding_states=64,padding_energy_cell_units=16,
    full_gauge_or_Fock_cutoff=None,native_full_source_compiler_executed=False,
    vacuum_task='ground full gauge-constrained interacting CAR+auxiliary source in all relevantcharge/baryon sectors; no positive-band three-particle shortcut'),
 physical_axial_ratio_computed=False,physical_hadron_mass_computed=False,physical_clock_seconds_derived=False)
assert pins=={str(p.relative_to(ROOT)):sha(p)for p in paths}
(HERE/'RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k]for k in ['gauge_covariance_error','local_weak_current_covariance_error','auxiliary_upper_component_operator_error','resource_contract']},indent=2))
print('full small CAR dimension1024, dressed hop columns',hops,'weak flux entries',int(weak_flux_changes))
