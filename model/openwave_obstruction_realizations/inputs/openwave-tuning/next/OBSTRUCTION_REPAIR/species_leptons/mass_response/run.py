"""Independent exact scalar reduction and compact-shell bound; cached grids are diagnostic."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import math

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
D3 = ROOT / 'openwave-push/b/D3_filled_reference_particle'
FILES = [D3/'lib_d3.py', D3/'field_L12_levels3.npz',
         D3/'consumed/ROW_23/rc_lib.py',
         D3/'consumed/source/coulomb_response_connection_run.py',
         D3/'consumed/D1_dressed_charge/consumed/geometry.json',
         HERE.parent/'COMMON_CLOCK.json',
         ROOT/'openwave-tuning/next/CHAIN3/species_reuse_connection/REGISTER.json']
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
hashes = {str(p.relative_to(ROOT)):sha(p) for p in FILES}

# No trigonometric or matrix-eigensolver assumption enters this exact identity.
m, e, w, b, s2 = sp.symbols('m E omega B s2', positive=True)
plus = -b*(e-m)/(2*e*(e+w-m))
pair = b*(e+m)/(2*e*(e+w+m))
combined = b*m*w/(e*((e+w)**2-m*m))
assert sp.cancel(plus+pair-combined) == 0
assert sp.expand(((e+w)**2-m*m)-(s2+2*e*w+w*w)).subs(e*e,m*m+s2) == 0
# If D(E)=E(s²+2Eω+ω²), D'(E)=s²+4Eω+ω²>0.
den = e*(s2+2*e*w+w*w)
assert sp.expand(sp.diff(den,e)-(s2+4*e*w+w*w)) == 0

clock = json.loads(FILES[5].read_text())
register = json.loads(FILES[6].read_text())
raw = [F(x) for x in clock['register_432_midpoints_exact']]
masses = {'electron':raw[1], 'muon':raw[2], 'tau':raw[0]}
assert masses['electron'] < F(1,12500)
assert masses['muon'] > F(3,250)
assert F(register['effect_intervals'][1][1]) < F(1,12500)
assert F(register['effect_intervals'][2][0]) > F(3,250)
# Shell h/250<=r<=h/125. pi<22/7. The lower bound does not sample a grid.
assert F(1,12500)**2+F(1,125)**2 < F(8001,1000000)**2
assert F(3,250)**2+F(1,125)**2 > F(9,625)**2
shell = F(1,250)*F(49,16*484)*(F(1000000,8001)-F(625,9))
assert shell > F(14,10000)
geom = json.loads(FILES[4].read_text())
assert len(geom['ico_seams']) == 30
assert [len(row) for row in geom['links']] == [2,2,2]
# For u=(1,...,1), u*L(k)u = 4 sum_a(1-cos k_a) <=2|k|².
# Cauchy-Schwarz in L gives (L^-1)_00 >= 1/(u*L*u), a=1/2.

fd = np.load(FILES[1])
k, wt, W, beta, v = (fd[n] for n in ['k','wt','W','beta','v'])
assert k.shape == (5717,3) and W.shape == (5717,24)
assert np.min(W)>0 and np.min(v)>0
ss = np.sum(np.sin(k)**2,axis=1)
B = np.sum(abs(beta)**2*np.cos(k[:,None,:]/2)**2,axis=2)
corners = ss<1e-25
assert int(corners.sum()) == 7

sx = np.array([[0,1],[1,0]],complex)
sy = np.array([[0,-1j],[1j,0]],complex)
sz = np.diag([1.,-1.]).astype(complex)
alpha = np.array([np.kron(pa,sx) for pa in [sx,sy,sz]])
be = np.kron(np.eye(2),sz)
up = np.eye(4)[:,[0,2]]
eye = np.eye(4)
errors=[]
for mu in [1/8, float(masses['muon'])/8, float(masses['electron'])/8]:
    for n in [0,7,83,1001,1726,1727,2199,4386,5716]:
        ee = math.sqrt(mu*mu+ss[n])
        pp = (eye+(mu*be-np.einsum('a,axy->xy',np.sin(k[n]),alpha))/ee)/2
        pm = (eye-(mu*be+np.einsum('a,axy->xy',np.sin(k[n]),alpha))/ee)/2
        mc = 1j*beta[n]*np.exp(.5j*k[n])*np.cos(k[n]/2)
        mats = np.einsum('ja,axy->jxy',mc,alpha)
        denseplus = sum(M.conj().T@pp@M/(mu-ee-ww) for M,ww in zip(mats,W[n]))
        densepair = sum(M@pm@M.conj().T/(mu+ee+ww) for M,ww in zip(mats,W[n]))
        cp = np.trace(up.T@denseplus@up).real/2
        cz = np.trace(up.T@densepair@up).real/2
        scalarplus = np.sum(B[n]*(1-mu/ee)/(2*(mu-ee-W[n])))
        scalarpair = np.sum(B[n]*(1+mu/ee)/(2*(mu+ee+W[n])))
        errors += [abs(cp-scalarplus),abs(cz-scalarpair)]
assert max(errors)<1e-11

records=[]
for h in [1/8,1/16,1/32]:
    for name, mass in [('D3_reference',F(1)),*masses.items()]:
        mu = h*float(mass)
        E = np.sqrt(mu*mu+ss)
        c = v/(2*E)
        p = np.sum(B*W/(E[:,None]*(ss[:,None]+2*E[:,None]*W+W*W)),axis=1)
        omitted = np.sum(B*(1-mu/E[:,None])/(2*(mu-E[:,None]-W)),axis=1)/mu
        for levels in range(4):
            n = 1727+1330*levels
            mask=np.arange(len(k))<n
            cn=mask&corners
            rn=mask&~corners
            total=float(np.sum(wt[mask]*(c[mask]+p[mask])))
            records.append({'h':h,'family':name,'mu':mu,'levels':levels,
                'fractional_rest_shift_diagnostic':total,
                'coulomb':float(np.sum(wt[mask]*c[mask])),
                'plus_and_pair':float(np.sum(wt[mask]*p[mask])),
                'seven_corner_contribution':float(np.sum(wt[cn]*(c[cn]+p[cn]))),
                'excluding_seven_corners':float(np.sum(wt[rn]*(c[rn]+p[rn]))),
                'omitted_pair_control':float(np.sum(wt[mask]*(c[mask]+omitted[mask]))),
                'coupling_off':0.0})

out = {'schema':'oph.filled_reference_mass_response.v1','input_sha256':hashes,
       'source_scope':'D3 free-filled reference, strict order e², rest spin average; conditional family transplant, not interacting poles',
       'symbolic_identity_verified':True,
       'fractional_kernel':'v/(2E)+sum_j B_j omega_j/[E(s²+2E omega_j+omega_j²)]',
       'fractional_kernel_positive_and_strictly_decreasing_in_positive_mass':True,
       'dense_matrix_scalar_max_error':max(errors),
       'shared_register_mass_coefficients':{a:str(q) for a,q in masses.items()},
       'shell':{'r_inner':'h/250','r_outer':'h/125','electron_coefficient_upper':'1/12500',
                'muon_coefficient_lower':'3/250','coulomb_lower':'1/(4 |k|²)',
                'electron_minus_muon_fractional_shift_lower_exact':str(shell),
                'electron_minus_muon_fractional_shift_lower':float(shell),
                'lower_ppm':1e6*float(shell),
                'uniform_for_h':[.125,.0625,.03125],
                'strict_order_relative_ratio_shift_upper_exact':str(-shell),
                'not_an_all_orders_pole_bound':True},
       'quadrature':{'size':len(k),'min_nonzero_component':float(np.min(abs(k[abs(k)>1e-14]))),
                     'unrefined_sine_corners':k[corners].tolist(),
                     'interpretation':'point weights at seven unrefined doubler corners produce artificial 1/m terms; numbers below are diagnostics, not certified integrals'},
       'records':records,
       'physical_lepton_mass_ratio_solved':False,
       'new_native_family_join_compiled':False,
       'perturbative_remainder_bounded':False}
assert hashes == {str(p.relative_to(ROOT)):sha(p) for p in FILES}
(HERE/'RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['dense_matrix_scalar_max_error','shell','quadrature']},indent=2))
for r in records:
    if r['h']==1/8 and r['levels']==3:
        print(r)
