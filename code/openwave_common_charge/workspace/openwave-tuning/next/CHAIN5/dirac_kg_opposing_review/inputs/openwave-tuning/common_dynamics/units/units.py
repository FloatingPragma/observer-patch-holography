"""Shared unit and canonical-charge normalization for the native proposals.

The Maxwell/Coulomb calculation is conditional on the source's canonical,
real-field limit. No observed constants or species masses enter the runner.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np
import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / 'vector_field'))
import linear_response as field


def rational(value, *, positive=False, nonnegative=False):
    if isinstance(value, bool):
        raise ValueError('a physical scalar is required, not a Boolean')
    try:
        answer = F(value)
    except (ValueError, TypeError, OverflowError, ZeroDivisionError) as exc:
        raise ValueError('finite rational physical scalar required') from exc
    if positive and answer <= 0 or nonnegative and answer < 0:
        raise ValueError('scalar outside its stated physical domain')
    return answer


def symbolic(value):
    value = rational(value)
    return s.Rational(value.numerator, value.denominator)


@dataclass(frozen=True)
class Units:
    length: F
    time: F
    action: F
    charge: F

    def __post_init__(self):
        for name in ('length', 'time', 'action', 'charge'):
            object.__setattr__(self, name, rational(getattr(self, name), positive=True))

    def cell_length(self, h):
        return self.length * rational(h, positive=True)

    def elapsed(self, macros, h):
        if type(macros) is not int or macros < 0:
            raise ValueError('nonnegative whole macro count required')
        return self.time * rational(h, positive=True) * macros

    def clock_phase(self, macros, h):
        return self.elapsed(macros, h) / (8 * self.time)

    def momentum(self, continuum_wave_number):
        return self.action * rational(continuum_wave_number) / self.length

    def eigenphase_energy(self, phase_lift, h):
        """Eigenvalue is exp(i phase_lift); the caller supplies its branch."""
        return -self.action * rational(phase_lift) / (self.time * rational(h, positive=True))

    def rest_energy(self, mass_coefficient):
        return self.action * abs(rational(mass_coefficient)) / self.time

    def speed_squared(self, dimensionless_speed_squared):
        return self.length**2 / self.time**2 * rational(dimensionless_speed_squared, positive=True)

    def rest_mass(self, mass_coefficient, dimensionless_speed_squared):
        return self.rest_energy(mass_coefficient) / self.speed_squared(dimensionless_speed_squared)

    def require_same(self, other):
        if not isinstance(other, Units) or self != other:
            raise ValueError('one shared length, time, action and charge map is required')


def field_normalization(alpha, beta, substeps, *, homogeneous):
    """Conditional canonical harmonic coefficients, never a vacuum claim."""
    if homogeneous is not True:
        raise ValueError('the homogeneous field setting is required for this Maxwell normalization')
    alpha = symbolic(rational(alpha, positive=True))
    beta = symbolic(rational(beta, positive=True))
    if type(substeps) is not int or substeps <= 0:
        raise ValueError('positive whole field substep count required')
    rho = s.Rational(2, 3)
    electric = substeps * rho * alpha
    magnetic = substeps * beta
    c2 = electric * magnetic
    g2 = s.sqrt(electric / magnetic)
    return {'electric': electric, 'magnetic': magnetic, 'speed_squared': c2,
            'charge_coupling_squared': g2, 'fine_structure': g2 / (4*s.pi),
            'Coulomb_pair_coefficient_in_hbar_length_over_time_units': electric / (4*s.pi),
            'per_substep_UV_stability_upper': alpha * beta * 158}


def admit_speed(matter_speed_squared, field_coefficients):
    if symbolic(matter_speed_squared) != field_coefficients['speed_squared']:
        raise ValueError('matter and field have different limiting speeds in the shared clock')


def canonical_geometry():
    """Reconstruct unit-flux coordinates, not arbitrary normalized vectors."""
    A0, B0, derivatives, terms, labels = field.bloch_data()
    basis = s.Matrix.hstack(*(B0.col_join(A0.T)).nullspace())
    crossings = s.zeros(3, 36)
    for axis in range(3):
        crossings[axis, 30 + 2*axis] = 1
        crossings[axis, 31 + 2*axis] = 1
    W = basis * (crossings * basis).inv()
    gram = W.T * W
    assert crossings * W == s.eye(3)
    assert crossings * A0 == s.zeros(3, 26)
    assert B0*W == s.zeros(12,3) and A0.T*W == s.zeros(26,3)
    assert gram == s.eye(3)*s.Rational(2,3)
    theta_embedding = W * gram.inv()
    assert theta_embedding.T * W == s.eye(3)
    assert W.T * theta_embedding == s.eye(3)
    P = W*gram.inv()*W.T
    null = s.Matrix.hstack(*A0.nullspace())
    Q = null*(null.T*null).inv()*null.T
    cases = []
    for p in ((1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,2,3)):
        D = sum((p[i]*derivatives[i] for i in range(3)),s.zeros(36,26))
        K = P*D*Q*D.T*P
        v = s.Matrix(p)
        desired = s.eye(3)*v.dot(v)-v*v.T
        assert theta_embedding.T*K*theta_embedding == desired
        cases.append({'momentum':p,'magnetic_form':str(desired)})
    # Independent scalar incidence Schur coefficient. Constant cell charge
    # has norm rho_cell/sqrt(12), so lambda~|k|²/8 implies Green 2/(3|k|²).
    L0 = B0*B0.T
    constant = s.ones(12,1)
    Pconstant = constant*constant.T/12
    inverse = (L0+Pconstant).inv()-Pconstant
    C = []
    for axis, row in enumerate(field.cc.geometry.links()):
        c = s.zeros(12,36)
        for channel, (a,b) in enumerate(row):
            c[b,30+2*axis+channel] = 1
        C.append(c)
    schur = []
    for p in ((1,0,0),(0,1,0),(0,0,1),(1,2,3)):
        c1 = sum((p[i]*C[i] for i in range(3)),s.zeros(12,36))
        c2 = sum((p[i]**2*C[i] for i in range(3)),s.zeros(12,36))
        L1 = s.I*(c1*B0.T-B0*c1.T)
        L2 = c1*c1.T+(c2*B0.T+B0*c2.T)/2
        coefficient = s.simplify((constant.T*(L2-L1*inverse*L1)*constant)[0]/12)
        expected = s.Rational(sum(x*x for x in p),8)
        assert coefficient == expected
        schur.append({'momentum':p,'Laplacian_quadratic_coefficient':str(coefficient)})
    return {'unit_flux_basis': [[str(x) for x in W.row(i)] for i in range(36)],
            'cycle_Gram': str(gram),'symplectic_matrix':str(theta_embedding.T*W),
            'magnetic_checks':cases,'scalar_Schur_checks':schur,
            'real_static_Green_prefactor':'2/3',
            'wrong_orthonormal_charge_crossing':str(s.sqrt(s.Rational(3,2)))}


def reject(call):
    try:
        call()
    except ValueError:
        return True
    raise AssertionError('invalid contract was admitted')


def run():
    inputs = [Path(__file__), Path(field.__file__), Path(field.cc.__file__),
              Path(field.cc.geometry.__file__),
              HERE.parent/'vector_field/MANIFEST.json',
              HERE.parent/'vector_field/LINEAR_RESPONSE.json',
              HERE.parent/'vector_field/CLOCK_JOIN.json',
              HERE.parent/'clock/MANIFEST.json',HERE.parent/'clock/PHASE_CLOCK.json',
              HERE.parent/'continuum/MANIFEST.json',HERE.parent/'energy/MANIFEST.json',
              HERE.parent/'energy/PROOF.md',HERE.parent/'particles/mass_coin.py',
              HERE.parent/'particles/MASS_COIN.json',HERE.parent/'particles/MASS_PROOF.md',
              HERE.parent/'spin/MANIFEST.json',HERE.parent/'spin/PROOF.md']
    pins = {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    geometry = canonical_geometry()
    scalar = field_normalization(F(1,8),F(1,4),1,homogeneous=True)
    dirac = field_normalization(F(1,16),F(3,8),8,homogeneous=True)
    admit_speed(F(1,48),scalar);admit_speed(1,dirac)
    assert scalar['charge_coupling_squared']==s.sqrt(3)/3
    assert dirac['charge_coupling_squared']==s.Rational(1,3)
    assert dirac['per_substep_UV_stability_upper']==s.Rational(237,64)<4
    alternative = field_normalization(F(1,32),F(3,4),8,homogeneous=True)
    assert alternative['speed_squared']==dirac['speed_squared']
    assert alternative['charge_coupling_squared']==dirac['charge_coupling_squared']/2
    # SI-like units are arbitrary positive symbols instantiated with rational
    # controls. The conversion identities hold for every positive choice.
    units = Units(F(7,3),F(5,2),F(11,7),F(13,11))
    assert units.elapsed(3,F(1,8))==3*units.time/8
    assert units.clock_phase(3,F(1,8))==F(3,64)
    assert units.eigenphase_energy(F(-1,8),F(1,8))==units.rest_energy(1)
    for c2 in (F(1,48),F(1)):
        for m in (F(0),F(1),F(3,2)):
            assert units.rest_mass(m,c2)*units.speed_squared(c2)==units.rest_energy(m)
    # A source coefficient changes energy; rescaling units changes the entire
    # source and cannot alter a dimensionless mass ratio or charge coupling.
    assert units.rest_energy(3)/units.rest_energy(1)==3
    qchanged=Units(units.length,units.time,units.action,2*units.charge)
    eps=lambda u,A: symbolic(u.charge**2*u.time/(u.action*u.length))/A
    eps0=eps(units,dirac['electric']);eps1=eps(qchanged,dirac['electric'])
    assert eps1==4*eps0
    c=s.sqrt(symbolic(units.speed_squared(1)))
    fine=symbolic(units.charge**2)/(4*s.pi*eps0*symbolic(units.action)*c)
    assert s.simplify(fine-dirac['fine_structure'])==0
    bad=[
      ('private_field_clock',lambda:units.require_same(Units(units.length,2*units.time,units.action,units.charge))),
      ('one_face_Maxwell_substitution',lambda:field_normalization(F(17,100),F(41,100),1,homogeneous=False)),
      ('wrong_homogeneous_speed',lambda:admit_speed(F(1,48),field_normalization(F(17,100),F(41,100),1,homogeneous=True))),
      ('missing_field_subcycles',lambda:admit_speed(1,field_normalization(F(1,16),F(3,8),1,homogeneous=True))),
      ('zero_time_scale',lambda:Units(1,0,1,1)),
      ('NaN_mass_input',lambda:units.rest_energy(float('nan'))),
      ('negative_refinement',lambda:units.elapsed(1,-1)),
      ('fractional_macro_count',lambda:units.elapsed(F(3,2),1))]
    mutations={name:reject(call) for name,call in bad}
    # Exact free-rest phases exp(±i mu) are source eigenphases. Equal-branch
    # populations beat at twice the single-branch Ramsey phase frequency.
    samples=64;n=np.arange(samples);mu=2*np.pi/samples
    phase=np.exp(-1j*mu*n)
    population=abs((phase+phase.conj())/2)**2
    phase_peak=int(np.argmax(abs(np.fft.fft(phase))))
    population_fft=abs(np.fft.fft(population-population.mean()))
    assert phase_peak==samples-1 and set(np.flatnonzero(population_fft>1e-10))=={2,samples-2}
    # Quantum quasienergy is a phase-lift choice. Two lifts separated by 2pi
    # agree on U and differ in energy by 2pi*hbar/(tau0*h).
    branch_shift=2*s.pi*symbolic(units.action)/(symbolic(units.time)*s.Rational(1,8))
    assert branch_shift>0
    assert pins=={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    encode=lambda r:{k:str(v) for k,v in r.items()}
    result={'schema':'oph.local12.shared_units.v1',
      'checks':{'actual_unit_flux_geometry':True,'canonical_symplectic_pair':True,
                'electric_and_magnetic_forms_same_geometry':True,'independent_static_Schur_prefactor':True,
                'scalar_and_Dirac_shared_speed_settings':True,'speed_does_not_fix_charge_coupling':True,
                'one_mass_energy_momentum_clock_map':True,'charge_unit_rescaling_invariant_coupling':True,
                'eight_invalid_joins_rejected':True,'rest_phase_and_population_beat_distinguished':True,
                'quasienergy_branch_retained':True,'source_bytes_unchanged':True},
      'geometry':geometry,'scalar_homogeneous':encode(scalar),'Dirac_homogeneous':encode(dirac),
      'equal_speed_different_coupling_control':encode(alternative),
      'conversion':{'cell_length':'ell0*h','macro_time':'tau0*h','clock_relative_phase':'n*h/8',
                    'momentum':'hbar*p/ell0','energy':'-hbar*lift(Arg U)/(tau0*h)',
                    'rest_energy':'hbar*abs(m)/tau0','physical_speed_squared':'(ell0/tau0)^2*c_dimensionless^2',
                    'rest_mass':'hbar*tau0*abs(m)/(ell0^2*c_dimensionless^2)',
                    'permittivity_if_Maxwell':'q0^2*tau0/(hbar*ell0*A_e)',
                    'Coulomb_pair_if_real_static_limit':'(hbar*ell0/tau0)*A_e*q1*q2/(4*pi*r)',
                    'dimensionless_charge_coupling':'g_charge^2=A_e/c=sqrt((2/3)*alpha/beta)',
                    'fine_structure_if_Maxwell':'g_charge^2/(4*pi)'},
      'calibration_ancestry':{'architectural':['routing and cycle Gram2/3','mass coefficient m','alpha,beta and field substep count','Clifford/CAR structures and source/readout words'],
                            'conversion_inputs':['ell0','tau0','hbar action unit','q0 charge unit'],
                            'empirical_if_selected':['Setting tau0=hbar*abs(m)/E_target calibrates the target rest energy.','Setting ell0=c_target*tau0/c_dimensionless calibrates a target speed.'],
                            'independent_of_unit_choices':['mass ratios','matter/field speed ratio','g_charge^2','all dimensionless native outcome probabilities'],
                            'absent_claims':['no observed fine-structure constant inserted','no electron species identified','no physical ground or interacting pole inferred']},
      'mutations':mutations,'phase_control':{'phase_Fourier_bin':phase_peak,'population_nonzero_bins':[2,62],
                    'energy_branch_spacing_example':str(branch_shift)},
      'conditional_premises':['Canonical unwrapped field algebra and prepared-state limit on the same source.',
                              'Homogeneous field coefficients and the stated shared macro clock.',
                              'Real static Gauss minimization or an independently proved quantum effective response; integer pointwise flux is not silently relaxed.',
                              'Routed unit-charge phase and coarse unit-flux coordinates retained.',
                              'Floquet, ultraviolet, cutoff and interaction errors controlled in the claimed observable domain.'],
      'scope':'The scalar matched member and Dirac replacement are separately typed. The executed one-face(.17,.41) member is not homogeneous Maxwell. The conditional Coulomb coefficient is not an executed two-charge force or a renormalized quantum coupling. g_charge denotes gauge coupling, not magnetic Lande g.',
      'source_sha256':pins,'official_rows_changed':[]}
    (HERE/'UNITS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'scalar':encode(scalar),'Dirac':encode(dirac),'mutations':mutations}))


if __name__=='__main__':
    run()
