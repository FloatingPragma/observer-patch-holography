"""An assigned SI dictionary on the complete interacting dimer, not its derivation.

Q=e(N-2), S=hbar*S_internal and H_SI=kappa*(K-E0 I) act on the SAME
sixteen-state space and reference. One measured electron energy fixes kappa;
the associated seconds/metres are imposed units, not an independently
recovered physical clock or spatial metric. BIPM fixes h,e,c exactly:
https://www.bipm.org/en/measurement-units/si-defining-constants
"""
from __future__ import annotations

from fractions import Fraction as F
import hashlib
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.empirical_diagnostics import cross_scale
from model.full_fock_field.field import FieldGeometry
from model.interacting_field.lattice import real
from model.particle_response import dimer
from model.particle_response.spectrum import matrix
from model.reference_observation.observation import frozen
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

PLANCK_J_S = F('6.62607015e-34')
ELEMENTARY_CHARGE_C = F('1.602176634e-19')
LIGHT_SPEED_M_S = F(299792458)
MEV_TO_J = 10**6*ELEMENTARY_CHARGE_C
DEFINING_CONSTANTS_SOURCE = 'https://www.bipm.org/en/measurement-units/si-defining-constants'
ROOT = Path(__file__).resolve().parents[2]


class BenchmarkDictionary:
    """Explicit assigned operators and units for the canonical dimer, all axes.

    The positive physical charge is assigned to excess CAR number. A one-charge
    hole is assigned electron charge; the field annihilator itself is not an
    already identified physical electron field. Internal SU(2) is assigned
    angular-momentum units without supplying a Lorentz representation.
    """
    def __init__(self, axis=0, *, root=ROOT):
        if type(axis) is not int or axis not in range(3):
            raise ValueError('coordinate axis 0,1,2 required')
        # Local import permits compare.build() to include this dictionary without
        # duplicating its strict observation contract or creating an import cycle.
        from .compare import load_observations, TARGETS, half_filled_candidates, Q
        self.root = Path(root)
        observations = load_observations(self.root/TARGETS)
        calibration = observations['electron_calibration']
        self.axis = axis
        self.mass_energy_MeV = F(calibration['value'])
        self.mass_uncertainty_MeV = F(calibration['standard_uncertainty'])
        _, gap = half_filled_candidates()
        self.exact_gap = gap
        self.exact_scale_MeV = Q(self.mass_energy_MeV)/gap
        sides = [1,1,1]; sides[axis] = 2
        self.spectrum = dimer.setup(FieldGeometry(tuple(sides),speed=.125,interaction=1.))
        self.dimension = self.spectrum.dimension
        self.gap = math.sqrt(17)/8-1/16
        self.scale_MeV = float(self.mass_energy_MeV)/self.gap
        self.scale_J = self.scale_MeV*float(MEV_TO_J)
        self.hbar_J_s = float(PLANCK_J_S)/(2*math.pi)
        self.seconds_per_model_time = self.hbar_J_s/self.scale_J
        self.metres_per_model_length = (float(LIGHT_SPEED_M_S)*self.seconds_per_model_time/
                                       self.spectrum.geometry.speed)
        self.shifted_generator = frozen(self.spectrum.generator-
                                        self.spectrum.energies[0]*np.eye(16))
        self.hamiltonian_J = frozen(self.scale_J*self.shifted_generator)
        self.reference = self.spectrum.reference
        self.charge_C = frozen(float(ELEMENTARY_CHARGE_C)*(self.spectrum.number-2*np.eye(16)))
        self.internal_spin = dimer.spin_generators(self.spectrum)
        self.spin_J_s = tuple(frozen(self.hbar_J_s*s) for s in self.internal_spin)
        self.projectors = {name:frozen(p) for name,p in dimer.candidate_projectors(self.spectrum).items()}
        self.charge_conjugation = frozen(dimer.charge_conjugation(self.spectrum))
        self._admit()

    def _admit(self):
        """Dimensionless complete-matrix residuals cannot hide behind tiny SI units."""
        h = self.hamiltonian_J/self.scale_J
        q = self.charge_C/float(ELEMENTARY_CHARGE_C)
        spin = [s/self.hbar_J_s for s in self.spin_J_s]
        casimir = sum(s@s for s in spin)
        conjugation = self.charge_conjugation
        residuals = [np.linalg.norm(h-self.shifted_generator,2),
                     np.linalg.norm(h@self.reference,2), np.linalg.norm(h@q-q@h,2),
                     np.linalg.norm(conjugation@q@conjugation.conj().T+q,2),
                     np.linalg.norm(conjugation@h@conjugation.conj().T-h,2),
                     np.linalg.norm(spin[0]@spin[1]-spin[1]@spin[0]-1j*spin[2],2)]
        residuals.extend(np.linalg.norm(h@s-s@h,2) for s in spin)
        for name,p in self.projectors.items():
            charge = (-1 if 'minus' in name else 1)*(2 if 'two' in name else 1)
            gap = self.gap if abs(charge)==1 else (1+math.sqrt(17)/4)/2
            residuals.extend((np.linalg.norm(p@p-p,2),np.linalg.norm(q@p-charge*p,2),
                              np.linalg.norm(h@p-gap*p,2),
                              np.linalg.norm(casimir@p-(.75 if abs(charge)==1 else 0)*p,2)))
        if (self.dimension != 16 or not np.isfinite(residuals).all()
                or max(residuals)>2e-11):
            raise ArithmeticError('assigned SI dictionary does not intertwine the complete model operators')
        return float(max(residuals))

    def model_duration(self, physical_seconds):
        seconds = real(physical_seconds,'assigned SI duration')
        duration = seconds/self.seconds_per_model_time
        if not math.isfinite(duration) or abs(duration)>32:
            raise ValueError('bounded dictionary dynamics requires |model duration|<=32')
        return duration

    def unitary(self, physical_seconds):
        """exp(-i H_SI t_SI/hbar); the ground offset changes only global phase."""
        duration = self.model_duration(physical_seconds)
        return expm(-1j*duration*self.shifted_generator)

    def heisenberg(self, operator, physical_seconds):
        a = matrix(operator,16); u = self.unitary(physical_seconds)
        return u.conj().T@a@u

    def site_charge_C(self, site):
        if type(site) is not int or site not in (0,1):
            raise ValueError('site 0 or 1 required')
        n = sum(self.spectrum.space.number(2*site+s).toarray() for s in (0,1))
        return frozen(float(ELEMENTARY_CHARGE_C)*(n-np.eye(16)))

    def current_A(self, edge=(0,1)):
        """Oriented charge OUT of edge[0], in C/s, on every Fock sector.

        I_xy=-dQ_x/dt=(e/hbar)*dH_SI(A_xy)/dA_xy at A=0.
        The link phase A_xy is dimensionless; this is a supplied Peierls
        correspondence, not a derived dynamical electromagnetic field.
        """
        if type(edge) is not tuple or edge not in ((0,1),(1,0)) or any(type(x) is not int for x in edge):
            raise ValueError('oriented dimer edge (0,1) or (1,0) required')
        q = self.site_charge_C(edge[0])
        return frozen(-1j*(self.hamiltonian_J@q-q@self.hamiltonian_J)/self.hbar_J_s)

    def report(self):
        """Deterministic exact-unit/interval report; floating admission is separate."""
        self._admit()
        root_lo,root_hi = cross_scale.sqrt_bounds(F(17),60)
        gap_lo,gap_hi = root_lo/8-F(1,16),root_hi/8-F(1,16)
        scale = self.mass_energy_MeV/gap_hi,self.mass_energy_MeV/gap_lo
        pi_lo,pi_hi = cross_scale.pi_bounds()
        time = (PLANCK_J_S/(2*pi_hi*scale[1]*MEV_TO_J),
                PLANCK_J_S/(2*pi_lo*scale[0]*MEV_TO_J))
        length = tuple(8*LIGHT_SPEED_M_S*t for t in time)
        compton = (PLANCK_J_S*LIGHT_SPEED_M_S/(2*pi_hi*self.mass_energy_MeV*MEV_TO_J),
                   PLANCK_J_S*LIGHT_SPEED_M_S/(2*pi_lo*self.mass_energy_MeV*MEV_TO_J))
        mass_kg = self.mass_energy_MeV*MEV_TO_J/LIGHT_SPEED_M_S**2
        sigma_kg = self.mass_uncertainty_MeV*MEV_TO_J/LIGHT_SPEED_M_S**2
        covariance = [[str(self.mass_uncertainty_MeV**2)]*2 for _ in range(2)]
        paths = ['model/particle_response/dimer.py','model/particle_response/spectrum.py',
                 'model/full_fock_field/field.py','model/particle_comparison/dictionary.py',
                 'model/particle_comparison/observations.json']
        return {
            'scope':'ASSIGNED_SI_HALF_FILLED_DIMER_BENCHMARK',
            'dimension':16,'all_number_sectors_retained':True,'axis':self.axis,
            'reference':'same unique N=2 pure ground of K=H-N/2; not the faithful Gibbs reference',
            'reference_is_faithful_on_full_algebra':False,
            'operator_dictionary':{'charge':'Q=e*(N-2I)', 'angular_momentum':'S=hbar*S_internal',
                'energy':'H_SI=kappa_J*(K-E0 I)', 'time':'t_SI=(hbar/kappa_J)*s_model',
                'length':'x_SI=(c_SI/c_model)*(hbar/kappa_J)*x_model',
                'local_charge':'Q_x=e*(n_x-I)',
                'oriented_current':'I_xy=-dQ_x/dt_SI=(e/hbar)*dH_SI(A_xy)/dA_xy'},
            'exact_SI_definitions':{'h_J_s':str(PLANCK_J_S),'e_C':str(ELEMENTARY_CHARGE_C),
                'c_m_per_s':str(LIGHT_SPEED_M_S),'MeV_to_J':str(MEV_TO_J),
                'hbar_J_s_formula':'h/(2*pi)','citation':DEFINING_CONSTANTS_SOURCE},
            'calibration':{'observable':'electron_mass_energy_equivalent',
                'value_MeV':str(self.mass_energy_MeV),'standard_uncertainty_MeV':str(self.mass_uncertainty_MeV),
                'citation':'https://physics.nist.gov/cuu/pdf/wall_2022.pdf',
                'one_common_scale_for_all_states':True,'scale_MeV_per_model_energy':self.exact_scale_MeV.as_json(),
                'independently_predicted_absolute_masses':0},
            'assigned_units':{'seconds_per_model_time':cross_scale._interval(time),
                'metres_per_model_length':cross_scale._interval(length),
                'clock_parameter':'assigned Hamiltonian model time; not native incident or recovered modular clock',
                'native_incident_or_modular_clock_intertwiner_supplied':False,
                'relative_standard_uncertainty_first_order':str(self.mass_uncertainty_MeV/self.mass_energy_MeV),
                'uncertainty_rule':'energy scale follows the measured mass; time and length vary inversely, with common calibration fully correlated',
                'uncertainty_includes_physical_dictionary_error':False,
                'intervals_bound_arithmetic_not_physical_model_error':True,
                'units_derive_from_measured_energy_and_imposed_speed_correspondence':True},
            'dimensionless_resolution_check':{
                'reduced_Compton_length_m':cross_scale._interval(compton),
                'reduced_Compton_length_over_assigned_lattice_spacing':
                    (type(self.exact_gap)(F(1,8))/self.exact_gap).as_json(),
                'assigned_lattice_spacing_over_reduced_Compton_length':(8*self.exact_gap).as_json(),
                'meaning':'one assigned lattice spacing exceeds the calibrated reduced Compton length; this dimer supplies no resolved relativistic rest dispersion',
                'universal_no_go_claimed':False},
            'paired_charge_one_assignment':{'names':['electron-like hole','positron-like particle'],
                'projectors':['charge_minus_one_hole','charge_plus_one_particle'],
                'charge_in_e':[-1,1],'internal_spin_squared_in_hbar_squared':['3/4','3/4'],
                'model_gap':self.exact_gap.as_json(),
                'calibrated_mass_energy_MeV':[str(self.mass_energy_MeV)]*2,
                'calibrated_mass_kg':[str(mass_kg)]*2,
                'shared_standard_uncertainty_kg':str(sigma_kg),
                'mass_energy_covariance_MeV_squared':covariance,
                'calibration_correlation_matrix':[[1,1],[1,1]],
                'dimensionless_mass_ratio':'1','normalized_mass_difference':'0'},
            'numerical_operator_admission':{'complete_matrices_checked':True,'dimensionless_tolerance':'2e-11',
                'floating_residuals_are_outward_physical_error_bounds':False},
            'supplied_correspondences':['zero-temperature half-filled reference and chemical potential',
                'physical charge equals excess CAR number times e',
                'internal SU(2) assigned physical angular-momentum units',
                'measured electron gap sets common energy and time units',
                'model coordinate speed assigned to c_SI to set length',
                'dimensionless link phase assigned Peierls electromagnetic coupling'],
            'unestablished':['Lorentz boosts and covariant spin/dispersion',
                'dynamical electromagnetic fields and physical current/detector identification',
                'same-excitation native topological defect and protection',
                'independent physical clock and vacuum/reference selection',
                'interacting chiral spatial continuum and independent measured predictions'],
            'physical_dictionary_supplied':True,'physical_identity_derived':False,
            'physical_clock_independently_calibrated':False,'electromagnetic_dynamics_derived':False,
            'Lorentz_representation_derived':False,'topological_defect_mapping_derived':False,
            'eligible_CP11_physical_comparisons_added':0,'OPH_universe_membership_confirmed':False,
            'source_sha256':bound_pins({p:hashlib.sha256((self.root/p).read_bytes()).hexdigest() for p in paths})}


def build(axis=0, *, root=ROOT):
    return BenchmarkDictionary(axis, root=root).report()


if __name__=='__main__':
    print(cross_scale.canonical(build()),end='')
