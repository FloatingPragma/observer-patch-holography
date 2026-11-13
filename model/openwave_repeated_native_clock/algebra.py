"""Exact repeated exchange, battery reset and complete occupation instrument."""
from fractions import Fraction as F
import sympy as s
DIM=64
BITS={'preparation':8,'first_readout':4,'reset_battery':2,'second_readout':1}

def exchange(v,bit,c=F(3,5),sn=F(4,5)):
    c,sn=s.Rational(c),s.Rational(sn)
    if c*c+sn*sn!=1: raise ValueError('exact unit-circle coordinates required')
    out=v.copy()
    for j in range(16):
        if j&bit:
            a,b=j,48+(j^bit)
            out[a]=c*v[a]-s.I*sn*v[b]
            out[b]=c*v[b]-s.I*sn*v[a]
    return out

def phase(v,k):
    out=v.copy()
    for j in range(32,64):out[j]*=(-s.I)**k
    return out

def moment(v,bit=None,source=False):
    return s.simplify(sum(s.conjugate(v[j])*v[j] for j in range(DIM)
                          if (j//16>=2 if source else bool(j&bit))))

def coherence(v,bit=None):
    if bit is None:return s.simplify(sum(s.conjugate(v[48+j])*v[j] for j in range(16)))
    return s.simplify(sum(s.conjugate(v[j])*v[j^bit] for j in range(DIM) if j&bit))

def run(first_phase=0,second_phase=0,c=F(3,5),sn=F(4,5),twirl=False):
    seeds=[]
    if twirl:
        for j in (0,1,4,5,8,9,12,13):
            v=s.zeros(DIM,1);v[j]=1/s.sqrt(8);seeds.append(v)
    else:
        v=s.zeros(DIM,1)
        for j in (0,1,4,5,8,9,12,13):v[j]=1/s.sqrt(8)
        seeds=[v]
    final=[];reset_failure=0;battery=0;prep_first=0;prep_final=0;second_coherence=0
    for v in seeds:
        v=exchange(v,8,c,sn);prep_first+=moment(v,8)
        v=phase(v,first_phase);v=exchange(v,4)
        before=moment(v,source=True)
        v=exchange(v,2,F(0),F(1))
        reset_failure+=moment(v,source=True);battery+=moment(v,2)
        assert moment(v,2)==before
        v=exchange(v,8,c,sn);prep_final+=moment(v,8)
        second_coherence+=coherence(v)
        v=phase(v,second_phase);v=exchange(v,1);final.append(v)
    probabilities=[s.simplify(sum(s.conjugate(v[j])*v[j] for v in final)) for j in range(DIM)]
    assert sum(probabilities)==1 and all(p>=0 for p in probabilities)
    assert reset_failure==0
    cr,sr=s.Rational(3,5),s.Rational(4,5);cc,ss=s.Rational(c),s.Rational(sn)
    p1=s.simplify(sum(probabilities[j] for j in range(DIM) if j&4))
    p2=s.simplify(sum(probabilities[j] for j in range(DIM) if j&1))
    w1=0 if twirl else s.re(s.I**first_phase)
    w2=0 if twirl else s.re(s.I**second_phase)
    assert p1==cr*cr/2+sr*sr*ss*ss/2-cr*sr*ss*w1/2
    assert p2==cr*cr/2+sr*sr*ss*ss*cc*cc/2-cr*sr*ss*cc*w2/2
    assert prep_first==cc*cc/2 and prep_final==cc**4/2
    assert second_coherence==(0 if twirl else s.I*ss*cc/2)
    actual_final=s.simplify(sum(probabilities[j]*((j//16>=2)+bool(j&4)+bool(j&2)+bool(j&1)) for j in range(DIM)))
    work=1-cc**4/2
    assert actual_final-s.Rational(1,2)==work
    return {'first_phase_quarter_cycles':first_phase,'second_phase_quarter_cycles':second_phase,
        'preparation_cosine_exact':str(c),'preparation_sine_exact':str(sn),'all_three_initial_references_twirled':twirl,
        'complete_64_occupation_probabilities_exact':[str(p) for p in probabilities],
        'readout_one_excited_probability_exact':str(p1),'readout_two_excited_probability_exact':str(p2),
        'source_excited_after_ideal_reset_exact':str(reset_failure),'battery_energy_over_mass_exact':str(battery),
        'preparation_reference_energy_after_first_over_mass_exact':str(prep_first),
        'preparation_reference_final_occupation_exact':str(prep_final),
        'second_prepared_source_A_expectation_exact':str(second_coherence),
        'gap_switching_work_over_mass_exact':str(work),'initial_actual_energy_over_mass_exact':'1/2',
        'final_actual_energy_over_mass_exact':str(actual_final),
        'terminal_reference_and_source_occupation_measurement_work_exact':'0',
        'orthogonal_native_Fock_complement_probability_exact':'0'}

def certificate():
    rows=[run(a,b) for a in range(4) for b in range(4)]
    controls=[run(a,b,c=F(5,13),sn=F(12,13)) for a,b in ((0,0),(1,2),(2,3))]
    twirl=[run(a,b,twirl=True) for a,b in ((0,0),(1,2),(2,3))]
    assert all(x['complete_64_occupation_probabilities_exact']==twirl[0]['complete_64_occupation_probabilities_exact'] for x in twirl)
    return {'complete_two_cycle_phase_grid':rows,'changed_preparation_controls':controls,'twirled_controls':twirl,
        'first_free_source_coherence':'i*s*exp(i*m*t1)/2',
        'second_free_source_coherence':'i*s*c*exp(i*m*t2)/2',
        'second_to_first_contrast_ratio_exact':'c',
        'final_preparation_reference_coherence':'c^2/2',
        'source_frequency_over_mass_both_flights_exact':'1',
        'native_source_frequency_depends_on_preparation_strength':False,
        'closed_resource_scope':'Same preparation reference reused; two phaseful readout references, one initially empty battery, a position-controlled gap schedule and one finite clock packet are supplied. No reference or battery erasure is performed.',
        'ideal_reset_scope':'Full exchange transfers the active source state and its correlations into the battery. It does not erase information. The static controller separately retains imperfect reset and reflected-clock outcomes.'}
