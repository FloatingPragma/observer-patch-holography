"""Complete dressed-charge records and exact infinite-source error bounds."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt, lcm
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.linalg import eigh, solve
from model.openwave_dirac_maxwell_refinement.operators import (
    pauli_operator, polynomial_member, charged_masses, binomial_half, root_tail,
    COORDINATES, ELECTRIC, MAGNETIC, SPEED_SQUARED)
from model.openwave_transverse_completion.source import flux

CAP = 16
DEGREE = 6
CHARGE = -3
REGULATOR = F(1, 16)
STATE_BUDGET = F(1, 10**10)
RESPONSE_BUDGET = F(1, 10**9)


def sqrt_upper(x):
    x = F(x)
    if x == 0:
        return F(0)
    d = 1 << 64
    return F(isqrt((x.numerator*d*d+x.denominator-1)//x.denominator)+1, d)


def dyadic(z, bits=48):
    d = 1 << bits
    return [(int(round(float(a.real)*d)), int(round(float(a.imag)*d))) for a in z], d


def norm_square(z, denominator=1):
    return F(sum(a*a+b*b for a, b in z), denominator*denominator)


def inner(x, y):
    return (sum(a*c+b*d for (a,b),(c,d) in zip(x,y)),
            sum(a*d-b*c for (a,b),(c,d) in zip(x,y)))


@lru_cache(maxsize=24)
def exact_polynomial(cap, charge, mass, regulator):
    """One complete spin block with Gaussian integer entries / one denominator."""
    mass, regulator = F(mass), F(regulator)
    pad = 2*DEGREE*abs(charge)+2
    op = pauli_operator(cap+pad, charge)
    n = len(op['orbital_basis'])
    a = op['D2'][:n,:n]*512
    # Every entry and every intermediate sum is a Gaussian integer. The row
    # absolute sum is <=4; six products remain strictly below 2**53.
    assert all(z.real == round(z.real) and z.imag == round(z.imag) for z in a.data)
    assert max(np.asarray(abs(a).sum(axis=1)).ravel()) <= 4
    keep = [i for i, (_, q) in enumerate(op['orbital_basis']) if abs(q) <= cap]
    basis = [op['orbital_basis'][i] for i in keep]
    coeff = [binomial_half(k)*mass**(1-2*k)/512**k for k in range(1, DEGREE+1)]
    den = lcm(32, *(z.denominator for z in coeff), (regulator/2).denominator)
    rows = [dict() for _ in basis]
    def add(matrix, scalar):
        factor = int(scalar*den)
        coo = matrix[keep][:,keep].tocoo()
        for i,j,z in zip(coo.row,coo.col,coo.data):
            assert z.real == round(z.real) and z.imag == round(z.imag)
            re, im = rows[i].get(int(j), (0,0))
            rows[i][int(j)] = (re+factor*int(round(z.real)), im+factor*int(round(z.imag)))
    power = sparse.eye(n, dtype=complex, format='csr')
    for k,c in enumerate(coeff,1):
        power = a@power
        assert max(np.asarray(abs(power).sum(axis=1)).ravel()) <= 4**k
        add(power,c)
    add(op['covariant_laplacian_derivatives'][0][:n,:n], regulator/2)
    index = {z:i for i,z in enumerate(basis)}
    for i,(x,q) in enumerate(basis):
        re,im = rows[i].get(i,(0,0))
        rows[i][i] = (re+int(den*(ELECTRIC*sum(v*v for v in flux(x,q,charge))+MAGNETIC)),im)
        for step in (-1,1):
            j=index.get((x,q+step))
            if j is not None:
                re,im=rows[i].get(j,(0,0));rows[i][j]=(re-int(den*MAGNETIC/2),im)
    for i,row in enumerate(rows):
        for j,(re,im) in row.items():
            assert rows[j].get(i,(0,0)) == (re,-im)
    return basis, rows, den


def action(rows, vector):
    return [(sum(a*vector[j][0]-b*vector[j][1] for j,(a,b) in row.items()),
             sum(a*vector[j][1]+b*vector[j][0] for j,(a,b) in row.items())) for row in rows]


def embed(vector, basis, target):
    values=dict(zip(basis, vector))
    return [values.get(z,(0,0)) for z in target]


def exact_residual(rows, den, v, vd, centre, rhs=None, frequency=F(0)):
    """Norm of (H-centre+i frequency)v-rhs, with exact rational arithmetic."""
    centre, frequency = F(centre), F(frequency)
    rhs = [(0,0)]*len(v) if rhs is None else rhs
    scale=lcm(den,centre.denominator,frequency.denominator)
    av=action(rows,v)
    z=[]
    for (ar,ai),(vr,vi),(br,bi) in zip(av,v,rhs):
        z.append((ar*(scale//den)-int(centre*scale)*vr-int(frequency*scale)*vi-scale*br,
                  ai*(scale//den)-int(centre*scale)*vi+int(frequency*scale)*vr-scale*bi))
    return sqrt_upper(norm_square(z,scale*vd))


def finite_spectrum_bound(basis, rows, den, eigen):
    """All eigenvalues, using a complete complex dyadic basis and integer products."""
    es,u=eigen;n=len(basis);hs=1<<28;vs=1<<24
    ar=np.zeros((n,n),dtype=np.int64);ai=np.zeros_like(ar)
    for i,row in enumerate(rows):
        for j,(a,b) in row.items():
            ar[i,j]=a*hs//den;ai[i,j]=b*hs//den
    zr=np.rint(u.real*vs).astype(np.int64);zi=np.rint(u.imag*vs).astype(np.int64)
    en=np.rint(es*hs).astype(np.int64)
    assert n*2*vs**2<2**63
    rowmax=max(sum(abs(int(x))+abs(int(y)) for x,y in zip(a,b)) for a,b in zip(ar,ai))
    assert (rowmax+max(abs(int(x)) for x in en))*vs<2**63
    gr=zr.T@zr+zi.T@zi-np.eye(n,dtype=np.int64)*vs*vs
    gi=zr.T@zi-zi.T@zr
    rr=ar@zr-ai@zi-zr*en;ri=ar@zi+ai@zr-zi*en
    g=sqrt_upper(F(sum(int(x)**2 for x in gr.flat)+sum(int(x)**2 for x in gi.flat),vs**4))
    residual=sqrt_upper(F(sum(int(x)**2 for x in rr.flat)+sum(int(x)**2 for x in ri.flat),(hs*vs)**2))
    # Rounding each real/imag entry down costs <=2n/hs in operator norm.
    rounding=F(2*n,hs);r=residual+rounding*sqrt_upper(1+g)
    k=F(max(abs(int(x)) for x in en),hs)
    assert g<F(1,100)
    error=(2*k*g+(1+g)*r)/(1-g)+rounding
    budget=F(1,100); assert error<=budget
    return F(int(en[1]),hs)-budget, {'dimension':n,'complete_basis_checked':True,
        'spectral_error_upper_exact':str(budget),'second_eigenvalue_lower_exact':str(F(int(en[1]),hs)-budget)}


def ground_member(name,mass,regulator=REGULATOR,cap=CAP):
    mass=F(mass);regulator=F(regulator)
    o=polynomial_member(cap,CHARGE,mass,regulator=regulator)
    n=len(o['basis'])//2;h=o['H_minus_mass'][:n,:n]
    assert np.array_equal(o['H_minus_mass'][:n,n:],np.zeros((n,n)))
    assert np.array_equal(h,o['H_minus_mass'][n:,n:].conj())
    e,u=eigh(h);g=u[:,0]
    basis,rows,den=exact_polynomial(cap,CHARGE,mass,regulator)
    reconstructed=np.zeros_like(h)
    for i,row in enumerate(rows):
        for j,(a,b) in row.items():reconstructed[i,j]=complex(float(F(a,den)),float(F(b,den)))
    assert np.max(abs(h-reconstructed))<3e-14
    second,finite=finite_spectrum_bound(basis,rows,den,(e,u))
    witness=json.loads(Path(__file__).with_name('witnesses.json').read_text())[name+'|'+str(regulator)]
    v=[tuple(map(int,z))for z in witness['ground']];vd=1<<48;norm=norm_square(v,vd)
    assert len(v)==len(basis)
    mu=F(witness['centre'])
    g=np.array([complex(a,b)/vd for a,b in v])/np.sqrt(float(norm))
    padded,prows,pden=exact_polynomial(cap+2*DEGREE*abs(CHARGE),CHARGE,mass,regulator)
    pv=embed(v,basis,padded)
    # A lower square-root bound avoids using an upper bound in a denominator.
    normlow=F(1)-abs(norm-1)
    assert normlow>F(99,100)
    residual=exact_residual(prows,pden,pv,vd,mu)/normlow
    eps=root_tail(mass,DEGREE)
    outside=F((cap+1-abs(CHARGE))**2,8)
    coupling=SPEED_SQUARED/(2*mass)+2*regulator+MAGNETIC/2
    a1=second-eps;assert outside>a1
    infinite_second=a1-coupling*coupling/(outside-a1)
    gap=infinite_second-mu
    assert gap>F(1,10)
    state_error=2*(residual+eps)/gap
    assert state_error<=STATE_BUDGET
    x=np.array([COORDINATES[z[0]][0] for z in basis])
    p_exact=F(sum(a*a+b*b for (a,b),xx in zip(v,x) if xx),vd*vd)/norm
    p_interval=(max(F(0),p_exact-2*STATE_BUDGET),min(F(1),p_exact+2*STATE_BUDGET))
    return {'name':name,'mass':mass,'regulator':regulator,'o':o,'h':h,'e':e,'u':u,'g':g,
        'basis':basis,'rows':rows,'den':den,'padded':padded,'prows':prows,'pden':pden,
        'v':v,'vd':vd,'norm':norm,'mu':mu,'residual':residual,'eps':eps,
        'p_exact':p_exact,'p_interval':p_interval,'x':x,
        'certificate':{'cap':cap,'two_complete_spin_blocks_conjugate':True,
            'source_mass_electron_units_exact':str(mass),'Wilson_energy_coefficient_exact':str(regulator),
            'complete_cap_spectrum':finite,'outside_electric_floor_exact':str(outside),
            'infinite_second_eigenvalue_lower_exact':str(infinite_second),
            'ground_vector_norm_error_upper_exact':str(STATE_BUDGET),
            'exact_padded_polynomial_residual_upper_exact':str(residual),
            'uniform_root_operator_error_exact':str(eps),
            'ground_radius_x_squared_cells_interval_exact':list(map(str,p_interval)),
            'variational_prepared_radius_x_squared_exact':str(p_exact)}}


def response_records(d):
    x=d['x'];h=d['h'];u=d['u'];e=d['e'];g=d['g'];D=-x
    j=1j*(D[None,:]-D[:,None])*h
    contact=-(D[None,:]-D[:,None])**2*h
    weights=abs(u.conj().T@(D*g))**2
    jw=abs(u.conj().T@(j@g))**2;gaps=e-e[0]
    ward=float(np.vdot(g,contact@g).real-2*np.dot(gaps,weights))
    assert abs(ward)<1e-12
    v=d['v'];vd=d['vd'];norm=d['norm'];dv=[(-int(xx)*a,-int(xx)*b)for (a,b),xx in zip(v,x)]
    rhs=embed(dv,d['basis'],d['padded']);rhsnorm=sqrt_upper(norm_square(rhs,vd)/norm)
    result=[]
    for nu in map(F,('1/64','1/16','1/4','1')):
        witness=json.loads(Path(__file__).with_name('witnesses.json').read_text())[d['name']+'|'+str(d['regulator'])]
        xv=[tuple(map(int,z))for z in witness['resolvents'][str(nu)]];xd=vd
        assert len(xv)==len(dv)
        px=embed(xv,d['basis'],d['padded'])
        rr=exact_residual(d['prows'],d['pden'],px,xd,d['mu'],rhs,nu)
        rr+=d['eps']*sqrt_upper(norm_square(px,xd))
        value=2*nu*nu*F(inner(dv,xv)[0],vd*vd)/norm
        # Resolvent identity, state trace distance and ground-centre shift.
        error=2*nu*rhsnorm*rr/(1-abs(norm-1))+2*nu*STATE_BUDGET+2*(d['residual']+d['eps'])
        assert error<=RESPONSE_BUDGET
        diagnostic=float(2*float(nu)**2*np.sum(gaps*weights/(gaps*gaps+float(nu)**2)))
        assert abs(diagnostic-float(value))<1e-11
        result.append({'imaginary_frequency_exact':str(nu),
            'full_infinite_source_response_interval_exact':[str(max(F(0),value-RESPONSE_BUDGET)),str(value+RESPONSE_BUDGET)],
            'response_diagnostic':diagnostic,'complete_resolvent_residual_verified':True})
    if d['regulator']:
        assert all(F(a['full_infinite_source_response_interval_exact'][1])<F(b['full_infinite_source_response_interval_exact'][0]) for a,b in zip(result,result[1:]))
    return {'static_Ward_residual_diagnostic':ward,'current_normalization':'physical electron charge Q=-1, not integer q=-3',
        'all_spectral_gaps_electron_units_diagnostic':gaps.tolist(),
        'all_dipole_weights_diagnostic':weights.tolist(),'all_current_weights_diagnostic':jw.tolist(),
        'second_spin_block_rule':'complex conjugate; every gap and real-dipole weight identical; all cross-spin matrix elements exactly zero',
        'spectral_table_has_individual_outward_error_bounds':False,
        'infinite_response_kernel':'Gamma(nu)=2 nu^2 <Dg,(H-E0)/[(H-E0)^2+nu^2] Dg>',
        'certified_infinite_response':result}


def probe_records(d):
    v=np.array([complex(a,b)/d['vd'] for a,b in d['v']]);v/=np.sqrt(float(d['norm']))
    h=d['h'];basis=d['basis'];x=d['x'];mass=float(d['mass'])
    initial=float(MAGNETIC);prepared=mass+float(np.vdot(v,h@v).real)
    preparation_work=.5*(prepared-initial)
    records=[]
    for k in (0.,np.pi/4,np.pi/2,np.pi):
        # The compensating reference remains at x=0: G=Q(cos(kx)-1).
        generator=1-np.cos(k*x);phase=np.exp(-1j*generator)
        spin_columns=[];spin_phase_energies=[];spin_final_energies=[]
        for state,hamiltonian in ((v,h),(v.conj(),h.conj())):
            columns=[]
            for axis in ('X','Y'):
                factor=1 if axis=='X' else -1j
                for sign in (-1,1):
                    columns.append((state+sign*factor*phase*state)/np.sqrt(8))
            nonselective=.5*np.outer(state,state.conj())+.5*np.outer(phase*state,(phase*state).conj())
            measured=sum((np.outer(w,w.conj())for w in columns),np.zeros_like(h))
            assert np.max(abs(nonselective-measured))<1e-15
            # The literal reported position/flux instrument dephases each
            # pointer branch. Its source-energy change is retained separately.
            dephased=np.diag(np.diag(measured))
            spin_phase_energies.append(mass+float(np.trace(measured@hamiltonian).real))
            spin_final_energies.append(mass+float(np.trace(dephased@hamiltonian).real))
            spin_columns.append(columns)
        probabilities=[sum(float(np.vdot(spin_columns[s][j],spin_columns[s][j]).real)for s in (0,1))/2 for j in range(4)]
        assert abs(sum(probabilities)-1)<1e-14
        phase_after=sum(spin_phase_energies)/2
        after=sum(spin_final_energies)/2
        phase_work=.5*(phase_after-prepared)
        position_work=.5*(after-phase_after)
        work=phase_work+position_work
        original_final=.5*initial+.5*after
        assert abs(original_final-initial-preparation_work-work)<2e-13
        records.append({'wavevector_radians_per_cell':k,
            'physical_probe_generator':'Q cos(k x)+R cos(0)=1-cos(k x); Q=-1,R=+1',
            'conditional_pointer_probabilities_Xminus_Xplus_Yminus_Yplus':probabilities,
            'complete_joint_position_flux_probabilities':[[ (.25*abs(spin_columns[s][j])**2).tolist() for s in (0,1)]for j in range(4)],
            'preparation_failure_probability_exact':'1/2',
            'spin_rule':'unpolarized two-spin mixture; each conjugate source branch is executed with the same probe phase; entries contain spin and preparation one-half factors',
            'original_input_phase_pointer_mean_work_electron_units':phase_work,
            'original_input_terminal_position_flux_mean_work_electron_units':position_work,
            'original_input_probe_mean_work_electron_units':work,
            'final_mean_source_energy_with_failure_electron_units':original_final,
            'actual_terminal_instrument':'phase pointer followed by all spin/position/flux projectors; the final source state is dephased in those reported outcomes'})
    p=d['p_interval'];rlo=sqrt_upper(p[0])-F(1,1<<63);rhi=sqrt_upper(p[1])
    assert rhi<F(1,4)
    # x in {0,1}: exact all-wavevector formula; Taylor inequality is uniform.
    return {'position_flux_labels':[list(z)for z in basis],
        'preparation_input':'matter/reference vacuum and loop electric n=0; mean field energy1/8',
        'prepared_state':'normalized finite dyadic variational trial in each conjugate spin block, unpolarized; true infinite ground within certified norm error',
        'original_preparation_mean_work_electron_units':preparation_work,
        'input_failure_mean_field_energy_exact':'1/8',
        'preparation_external_above_true_field_vacuum':'1/8 minus the full rotor ground energy; separate from charged preparation work',
        'form_factor_exact':'F(k)=1-p_x1 (1-cos k), p_x1=<x^2>; charge response after retaining and separately reporting reference is Q F(k)',
        'reference_inclusive_response_exact':'Q[F(k)-1]; zero at k=0 because the complete state is neutral',
        'radius_x_cells_interval_exact':[str(max(F(0),rlo)),str(rhi)],
        'largest_first_Brillouin_zone_suppression_interval_exact':[str(2*p[0]),str(2*p[1])],
        'independent_radius_measurement':'complete native occupation/flux instrument, not inferred from a current gap',
        'small_wavevector_remainder':'|1-F(k)-k^2<x^2>/2| <= <x^2> k^4/24',
        'order_one_core_onset_resolved':False,
        'reason':'2<x^2><1/16; no wavevector in the first Brillouin zone yields even a one-sixteenth suppression',
        'records':records,'point_probabilities_and_work_are_numerical_diagnostics':True}


def certificate():
    rows=[]
    for name,mass in charged_masses():
        for regulator in (REGULATOR,F(0)):
            d=ground_member(name,mass,regulator)
            rows.append({'family':name,'ground_certificate':d['certificate'],
                         'native_current_response':response_records(d),'complete_core_probe':probe_records(d)})
    return {'schema':'oph.dressed_native_core_response.v1',
        'verdict':'PASS_INFINITE_DRESSED_CURRENT_AND_COMPLETE_CORE_RESOLUTION_CONTROL',
        'source':'same intraband Pauli-square-root plus compact Maxwell and positive Wilson energy member as openwave_dirac_maxwell_refinement',
        'native_charge_unit':'q=3Q=-3; physical dipole and scalar probes use Q=-1',
        'ancestry':'charged masses and electron GeV unit inherited; electric1/32,magnetic1/8,Wilson1/16 declared without core/running targets; zero-Wilson control retained',
        'full_Fock_and_native_complements':'inherited reducing positive one-particle, colourful-vacuum sector; no transfer to pair-producing full Dirac dynamics',
        'source_has_all_integer_fluxes':True,'members':rows,
        'physical_QED_running_or_resolved_electron_core_fit':False,
        'literal_core_onset_inference':'Actual coupling varies with frequency and wavevector, but the independently measured spatial core is subcell and its order-one onset is unresolved; the spectral gap is not assigned as a core radius.'}
