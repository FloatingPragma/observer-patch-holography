"""Finite-source Gibbs states and a separately declared continuum ideal gas."""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations_with_replacement
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import kv,zeta

from model.particle_common_contract.provenance import import_closure

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
PARENT=ROOT/'model/conditional_particle_spectrum/results.json'
TEMPERATURES=(1e-12,1e-11,1e-10,1e-4,1e-3,.01,.1,1.,5.,10.,50.,100.,160.,1000.)
QUARKS={'up','down','charm','strange','top','bottom'}
CHARGED={'electron','muon','tau'}


def bindings():
    paths={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md') or p.name=='comparison_inputs.json'}
    paths.update({'model/conditional_particle_spectrum/results.json','model/conditional_particle_spectrum/PROOF.md',
        'model/conditional_particle_spectrum/source.py','model/conditional_particle_spectrum/readout.py',
        'model/common_source/calibrated_matter_results.json','model/native_abundance_epochs/species_identifiability.py'})
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(import_closure(paths))}


@lru_cache(maxsize=1)
def parent():
    value=json.loads(PARENT.read_text())
    for name,digest in value['input_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('calibrated spectrum source changed: '+name)
    return value


@lru_cache(maxsize=1)
def species():
    result=[]
    for row in parent()['elementary_particles']:
        name=row['species'];fermion='positive_modes_per_momentum' in row
        group='quark' if name in QUARKS else 'charged_lepton' if name in CHARGED else 'neutrino' if name.startswith('neutrino_') else 'boson'
        result.append(dict(species=name,group=group,fermion=fermion,
            degrees=row['positive_modes_per_momentum' if fermion else 'physical_modes_per_momentum'],
            mass_GeV=row['rest_coefficient_GeV']['numeric'],mass_exact=row['rest_coefficient_GeV']['exact'],
            mass_convention=row['mass_convention'],independent_mass_prediction=False))
    if len(result)!=17 or sum(r['degrees'] for r in result if r['fermion'])!=90 or sum(r['degrees'] for r in result if not r['fermion'])!=28:
        raise ArithmeticError('complete physical degree inventory differs')
    return result


def occupation(energy,fermion):
    q=math.exp(-energy) if energy<746 else 0.
    return q/(1+q) if fermion else q/(-math.expm1(-energy))


@lru_cache(maxsize=1024)
def integrals(ratio,fermion):
    """Per-degree quantities in powers of T; quadrature errors are estimates."""
    if not math.isfinite(ratio) or ratio<0 or type(fermion) is not bool:
        raise ValueError('finite nonnegative m/T and a Boolean statistics label required')
    factor=1/(2*math.pi**2)
    if ratio==0:
        n=float(zeta(3,1))/math.pi**2*(.75 if fermion else 1)
        rho=math.pi**2/30*(.875 if fermion else 1);pressure=rho/3
        return dict(number=n,energy=rho,pressure=pressure,entropy=rho+pressure,
            heat_capacity=4*rho,pressure_derivative=4*pressure,
            trace_anomaly=0.,quadrature_error_estimates={k:0. for k in ('number','energy','pressure','heat_capacity','pressure_derivative','trace_anomaly')},
            evaluation='analytic massless integrals; binary evaluation of constants')
    if ratio>=200:
        # epsilon>= (x+y)/2; epsilon<=x+y. All bounds are per degree.
        # The final polynomial also dominates the derivative integrals.
        log_upper=math.log(2*factor)+math.log(768+192*ratio+16*ratio**2)-ratio/2
        upper=math.exp(log_upper) if log_upper>-745 else 0.
        return dict(number=0.,energy=0.,pressure=0.,entropy=0.,heat_capacity=0.,
            pressure_derivative=0.,trace_anomaly=0.,quadrature_error_estimates={},
            evaluation='exponentially suppressed contribution omitted within stated analytic bound',
            each_dimensionless_integral_log_error_upper=log_upper,
            each_dimensionless_integral_error_upper=upper,
            analytic_bound_formula='2*(768+192*y+16*y^2)*exp(-y/2)/(2*pi^2), y=m/T>=200',
            underflowed_upper_retained_as_log=upper==0.)
    def value(x,kind):
        energy=math.hypot(x,ratio);f=occupation(energy,fermion)
        if kind=='number':return x*x*f
        if kind=='energy':return x*x*energy*f
        if kind=='pressure':return x**4*f/(3*energy)
        if kind=='trace_anomaly':return x*x*ratio**2*f/energy
        derivative=f*(1-f if fermion else 1+f)
        if kind=='heat_capacity':return x*x*energy*energy*derivative
        return x**4*derivative/3
    result={};errors={}
    for name in ('number','energy','pressure','heat_capacity','pressure_derivative','trace_anomaly'):
        integral,error=quad(value,0,np.inf,args=(name,),epsabs=2e-11,epsrel=2e-11,limit=220)
        result[name]=factor*integral;errors[name]=factor*error
    result['entropy']=result['energy']+result['pressure']
    result.update(quadrature_error_estimates=errors,evaluation='adaptive real quadrature; errors are numerical estimates, not rigorous intervals')
    return result


def gas(temperature):
    if type(temperature) not in (int,float) or not math.isfinite(temperature) or temperature<=0:
        raise ValueError('finite positive temperature in the declared GeV comparison is required')
    rows=[];sums=defaultdict(float)
    for item in species():
        result=integrals(item['mass_GeV']/temperature,item['fermion']);g=item['degrees']
        row=dict(item,ratio_m_over_T=item['mass_GeV']/temperature,
            number_over_T3=g*result['number'],energy_over_T4=g*result['energy'],
            pressure_over_T4=g*result['pressure'],entropy_over_T3=g*result['entropy'],
            heat_capacity_over_T3=g*result['heat_capacity'],
            pressure_derivative_over_T3=g*result['pressure_derivative'],
            trace_anomaly_over_T4=g*result['trace_anomaly'],numerics=result)
        rows.append(row)
        for key in ('number_over_T3','energy_over_T4','pressure_over_T4','entropy_over_T3',
                    'heat_capacity_over_T3','pressure_derivative_over_T3','trace_anomaly_over_T4'):sums[key]+=row[key]
    for row in rows:
        for name,key in [('number','number_over_T3'),('energy','energy_over_T4'),('entropy','entropy_over_T3')]:
            row[name+'_fraction']=row[key]/sums[key]
    aggregates={}
    for group in ('quark','charged_lepton','neutrino','lepton','boson'):
        selected=[r for r in rows if r['group']==group or (group=='lepton' and r['group'] in ('charged_lepton','neutrino'))]
        aggregates[group]={name:sum(r[name+'_fraction'] for r in selected) for name in ('number','energy','entropy')}
    rho=sums['energy_over_T4'];pressure=sums['pressure_over_T4'];entropy=sums['entropy_over_T3']
    sound=sums['pressure_derivative_over_T3']/sums['heat_capacity_over_T3']
    derivative_defect=abs(entropy-sums['pressure_derivative_over_T3'])
    trace_defect=abs(rho-3*pressure-sums['trace_anomaly_over_T4'])
    if derivative_defect>1e-8 or trace_defect>1e-8 or not 0<sound<=1/3+1e-12:
        raise ArithmeticError('ideal-gas thermodynamic identities fail')
    return dict(temperature_GeV=temperature,species=rows,totals=dict(sums),aggregate_fractions=aggregates,
        g_energy=30*rho/math.pi**2,g_entropy=45*entropy/(2*math.pi**2),
        equation_of_state_w=pressure/rho,sound_speed_squared=sound,
        thermodynamic_residuals=dict(pressure_derivative_minus_entropy=derivative_defect,
            trace_anomaly_identity=trace_defect),
        physical_temperature_evolution_predicted=False,
        physical_domain='formal fixed-mass free gas only; no temperature range is certified as physical by this calculation')


def bessel_control(ratio,fermion):
    """Independent massive fugacity series with an analytic truncation bound."""
    if ratio<.1:raise ValueError('control is restricted to m/T>=0.1')
    result=np.zeros(3);tail=np.ones(3)
    for ell in range(1,10001):
        k1,k2=kv(1,ell*ratio),kv(2,ell*ratio)
        terms=ratio**2/(2*math.pi**2)*np.array([k2/ell,3*k2/ell**2+ratio*k1/ell,k2/ell**2])
        result+=(-1 if fermion and ell%2==0 else 1)*terms
        nxt=ell+1;k1,k2=kv(1,nxt*ratio),kv(2,nxt*ratio)
        omitted=ratio**2/(2*math.pi**2)*np.array([k2/nxt,3*k2/nxt**2+ratio*k1/nxt,k2/nxt**2])
        tail=omitted if fermion else omitted/(-math.expm1(-ratio))
        if max(tail)<1e-12:break
    original=integrals(ratio,fermion);direct=np.array([original[k] for k in ('number','energy','pressure')])
    error=float(max(abs(result-direct)))
    if error>2e-10:raise ArithmeticError('independent fugacity series disagrees')
    return dict(mass_over_temperature=ratio,fermion=fermion,terms=ell,
        dimensionless_number_energy_pressure=result.tolist(),analytic_series_tail_upper=tail.tolist(),
        numerical_quadrature_disagreement=error,special_function_roundoff_interval_certified=False)


def finite_source_gibbs(temperature):
    """Literal two-momentum quadratic source, including global boson cap two."""
    if type(temperature) not in (int,float) or not math.isfinite(temperature) or temperature<=0:
        raise ValueError('finite positive source temperature required')
    modes=parent()['mode_inventory']['boson_modes'];fermions=parent()['mode_inventory']['fermion_modes']
    words=[()]+[(i,) for i in range(len(modes))]+list(combinations_with_replacement(range(len(modes)),2))
    energies=np.array([sum(modes[i]['energy_GeV']['numeric'] for i in word) for word in words])
    weights=np.exp(-energies/temperature);partition=float(weights.sum());probabilities=weights/partition
    if len(words)!=1653 or abs(float(probabilities.sum())-1)>2e-14:raise ArithmeticError('complete cap-two thermal law lost mass')
    counts=Counter()
    for word in words:
        key=tuple(sum(modes[i]['species']==name if name!='massless' else modes[i]['species'] in ('photon','gluon') for i in word)
                  for name in ('W','Z','Higgs','massless'))
        counts[key]+=1
    for row in parent()['boson_readout']['complete_energy_projectors']:
        key=tuple(row['occupation_by_energy_class'][name] for name in ('W','Z','Higgs','massless'))
        if counts[key]!=row['projector_rank']:raise ArithmeticError('actual boson energy ranks differ')
    mean_energy=float(probabilities@energies);mean_number=sum(float(p)*len(w) for p,w in zip(probabilities,words))
    fenergies=[row['energy_GeV']['numeric'] for row in fermions]
    foccup=[occupation(e/temperature,True) for e in fenergies]
    f_logz=sum(math.log1p(math.exp(-e/temperature)) for e in fenergies)
    f_energy=sum(e*n for e,n in zip(fenergies,foccup))
    uncapped=sum(occupation(row['energy_GeV']['numeric']/temperature,False) for row in modes)
    return dict(temperature_GeV=temperature,boson_dimension=1653,global_boson_cap=2,
        energy_convention='quasiparticle excitation energy above the selected quadratic vacuum; additive vacuum constant subtracted',
        all_boson_occupation_records=[dict(occupied_modes=list(word),energy_GeV=float(e),probability=float(p),
            log_probability=-float(e)/temperature-math.log(partition),positive_probability_underflowed=bool(p==0))
            for word,e,p in zip(words,energies,probabilities)],
        boson_partition=partition,boson_mean_number=mean_number,boson_mean_energy_GeV=mean_energy,
        boson_entropy=math.log(partition)+mean_energy/temperature,
        uncapped_same_mode_mean_boson_number=uncapped,
        fermion_CAR_modes=len(fermions),fermion_Fock_dimension=2**len(fermions),
        complete_factorized_fermion_occupations=foccup,fermion_log_partition=f_logz,
        complete_factorized_fermion_log_records=[dict(mode=i,
            log_occupied=-e/temperature-math.log1p(math.exp(-e/temperature)),
            log_blank=-math.log1p(math.exp(-e/temperature)),positive_occupied_probability_underflowed=n==0)
            for i,(e,n) in enumerate(zip(fenergies,foccup))],
        fermion_mean_energy_GeV=f_energy,fermion_entropy=f_logz+f_energy/temperature,
        source_temperature_is_supplied=True,continuum_phase_space_or_Universe_abundances=False)


def massless_limit():
    gf=sum(r['degrees'] for r in species() if r['fermion']);gb=sum(r['degrees'] for r in species() if not r['fermion'])
    gn=gb+F(3,4)*gf;gr=gb+F(7,8)*gf
    return dict(fermion_degrees=gf,boson_degrees=gb,number_degrees=str(gn),energy_and_entropy_degrees=str(gr),
        quark_number_fraction=str(F(3,4)*72/gn),lepton_number_fraction=str(F(3,4)*18/gn),
        quark_energy_fraction=str(F(7,8)*72/gr),lepton_energy_fraction=str(F(7,8)*18/gr),
        quark_fraction_among_fermions='4/5',extra_sterile_neutrino_degrees=0,
        exact_global_cap2_boson_high_temperature_mean='3248/1653')


def calculate():
    pins=bindings();rows=[gas(t) for t in TEMPERATURES];reference=gas(10.)
    for row in rows:
        t=row['temperature_GeV'];s=row['totals']['entropy_over_T3']*t**3
        rho=row['totals']['energy_over_T4']*t**4
        row['conditional_adiabatic_expansion']=dict(
            scale_factor_relative_to_T10=(reference['totals']['entropy_over_T3']*10**3/s)**(1/3),
            Hubble_rate_relative_to_T10=math.sqrt(rho/(reference['totals']['energy_over_T4']*10**4)),
            minus_d_log_a_d_log_T=row['totals']['heat_capacity_over_T3']/(3*row['totals']['entropy_over_T3']),
            assumptions=['common local equilibrium and zero chemical potentials','homogeneous adiabatic expansion with a^3*s constant',
                         'flat Friedmann density law for the Hubble ratio'],
            native_expansion_law_derived=False,absolute_time_seconds=None)
    comparisons=[];source=json.loads((HERE/'comparison_inputs.json').read_text())
    for item in source['samples']:
        model=gas(item['temperature_GeV'])
        comparisons.append(dict(**item,ideal_g_energy=model['g_energy'],ideal_g_entropy=model['g_entropy'],
            relative_energy_difference=(model['g_energy']/item['g_rho']-1),
            relative_entropy_difference=(model['g_entropy']/item['g_entropy']-1),
            model_is_interacting_EoS=False,mass_scheme_or_RG_matching_derived=False,
            empirical_validation_or_native_abundance_claim=False))
    controls=[bessel_control(r,f) for r in (.1,1.,10.) for f in (False,True)]
    finite=[finite_source_gibbs(t) for t in (.1,1.,10.,100.)]
    limit=massless_limit()
    checks=dict(all_17_species_counted=True,no_neutrino_or_antiparticle_double_count=True,
        massless_limit_energy_106_75=limit['energy_and_entropy_degrees']=='427/4',
        massless_limit_matches_declared_comparison=limit['quark_number_fraction']=='108/191' and limit['lepton_number_fraction']=='27/191',
        thermodynamic_entropy_pressure_identity=all(r['thermodynamic_residuals']['pressure_derivative_minus_entropy']<1e-8 for r in rows),
        complete_finite_boson_states_retained=all(len(r['all_boson_occupation_records'])==1653 for r in finite),
        finite_cap_distinct_from_uncapped_gas=all(r['boson_mean_number']<r['uncapped_same_mode_mean_boson_number'] for r in finite),
        finite_cap_mean_never_above_two=all(r['boson_mean_number']<=2 for r in finite),
        stable_causal_ideal_gas=all(0<r['sound_speed_squared']<=1/3+1e-12 for r in rows))
    if not all(checks.values()):raise ArithmeticError('thermal history comparison failed')
    if pins!=bindings():raise ArithmeticError('thermal sources changed during calculation')
    return dict(schema='oph.particle_thermal_history.v1',all_checks_pass=True,source_bindings=pins,
        construction=dict(source='calibrated17-species quadratic menu; actual finite cap retained separately',
            continuum_extension='independent free Fock occupations on R3 momentum, thermodynamic volume limit, no zero-mode condensate',
            temperature='supplied common GeV temperature',chemical_potentials='all zero',
            masses='fixed parent action coefficients; MZ-running quark/charged-lepton inputs are not physical thermal poles',
            energy_density='thermal excitation energy with the zero-temperature vacuum contribution subtracted',
            sterile_population='none',pions_or_hadrons_added='none; adding them alongside free quarks would require a separate matching model'),
        high_temperature_limit=limit,temperature_history=rows,finite_source_Gibbs_controls=finite,
        independent_fugacity_series_controls=controls,interacting_reference_comparison=comparisons,
        comparison_reference={k:v for k,v in source.items() if k!='samples'},checks=checks,
        numerical_status='analytic thermodynamic identities; numerical quadrature with error estimates and independent series controls, not certified interval arithmetic',
        physical_admission=dict(native_temperature_clock_energy_dictionary=False,
            native_Gibbs_state_preparation=False,RG_and_thermal_selfenergy_matching=False,
            confinement_and_hadron_replacement=False,electroweak_temperature_dependent_background=False,
            neutrino_decoupling_and_separate_temperatures=False,interaction_corrected_equation_of_state=False,
            cosmological_expansion_attachment=False,physical_abundances_or_freezeout_predicted=False),
        interpretation='The curve is a conditional free comparison with known inadmissible low-temperature colored gas and fixed-background high-temperature extensions. No point is promoted to a native early-Universe population prediction.')


def main():
    parser=argparse.ArgumentParser();m=parser.add_mutually_exclusive_group(required=True)
    m.add_argument('--write',action='store_true');m.add_argument('--check',action='store_true');args=parser.parse_args()
    result=calculate();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')
    elif json.loads(path.read_text())!=result:raise ArithmeticError('thermal history receipt differs')
    print('PASS:17-species conditional thermodynamics and complete finite-cap Gibbs controls')


if __name__=='__main__':main()
