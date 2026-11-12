"""Induced baryon ordering, complete intrinsic charge records and beta obstruction."""
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
from fractions import Fraction as F
from model.openwave_vectorlike_flavour.flavour import certificate as flavour_certificate
from .source import spectral_certificate,exact_record_witness
from model.lepton_mass_compiler.analyzer import rational_effect_intervals
from model.native_flavor_orbit import intervals as iv

def ratio_intervals():
 e=iv.I(*map(F,rational_effect_intervals()["effect_intervals_k0_k1_k2"][1]))
 def root(a):return iv.I(iv.sqrt_rational(a.lo).lo,iv.sqrt_rational(a.hi).hi)
 out=[]
 for y in (F(2,3),F(1,3)):
  a=e**2+y*y+1
  low=root(2*e**2/(a+root(a**2-4*e**2)))
  out.append(iv.rounded_out(low/e,35))
 return out

def decimal_outward(box):
 result=[]
 for x,rounding in zip((box.lo,box.hi),(ROUND_FLOOR,ROUND_CEILING)):
  with localcontext() as c:
   c.prec=32;c.rounding=rounding;result.append(str(Decimal(x.numerator)/Decimal(x.denominator)))
 return result

def certificate():
 f=flavour_certificate();s=spectral_certificate();states=[exact_record_witness(q) for q in (0,1)]
 with localcontext() as ctx:
  ctx.prec=40;mu=Decimal(f['quark_branches'][0]['light_gap_in_electron_units']);md=Decimal(f['quark_branches'][1]['light_gap_in_electron_units']);unit=Decimal(f['calibration']['rows'][0]['electron_anchored_GeV'])
  energies=[];mu_box,md_box=ratio_intervals();anchor=F(str(unit));baseline=iv.I(*map(F,s['lowest_interval_rational']))
  for charge,state in zip((0,1),states):
   nu=charge+1;nd=3-nu;shift=nu*(mu-1)+nd*(md-1)
   interval=baseline+nu*(mu_box-1)+nd*(md_box-1)
   energies.append({'charge':charge,'light_up_count':nu,'light_down_count':nd,'source_gap_interval_rational':interval.strings(),'source_gap_interval':decimal_outward(interval),'conditional_GeV_interval_rational':(interval*anchor).strings(),'conditional_GeV_interval':decimal_outward(interval*anchor)})
  delta=md-mu;assert 0<delta<1
  return {'verdict':'PASS_NATIVE_BARYON_ORDERING_AND_SIGNED_INTRINSIC_PROFILE_WITH_BETA_OBSTRUCTION',
   'source':'Full15-native-coordinate/two-spin CAR source on two graph cells; positive graph hopping, native colour, colourful-spin and colourful-weak Casimirs; m_up,m_down are declared number energies taken from exact free vectorlike gaps',
   'parameters':{'colour_U':1,'spin_J':'.25','weak_J':'.25','hopping_t':'.125'},
   'architecture_ancestry':'The equal positive spin/weak-current weights are an exploratory extension of the spin-only source, motivated by testing baryon and two-baryon structure. They are not a blind physical prediction or fitted hadron masses.',
   'free_to_interacting_boundary':'Using free branch gaps as composite number-energy coefficients is a declared matching rule, not a proof of an interacting-band projection of the vectorlike Weyl/Higgs action.',
   'exact_free_mass_ratio_intervals':{'up':mu_box.strings(),'down':md_box.strings()},'calibration_interval_scope':'Rational analyzer and spectral intervals propagated exactly; decimal endpoints rounded outward. The empirical electron anchor is held at its declared central value; experimental uncertainty and physical model error are excluded.',
   'exact_equal_mass_spectrum_certificate':s,'complete_intrinsic_records':states,'baryon_energies':energies,
   'neutron_like_minus_proton_like_source_exact_definition':'m_down_light-m_up_light','neutron_like_minus_proton_like_electron_units':str(delta),'neutron_like_minus_proton_like_conditional_GeV':str(delta*unit),'neutron_like_minus_proton_like_conditional_GeV_interval_rational':((md_box-mu_box)*anchor).strings(),
   'exact_beta_threshold_obstruction':{'identity':'Delta E=m_down_light-m_up_light<bare_light_analyzer_mass=m_e for every G,M>0 in this common-family construction','electron_threshold_electron_units':'1','n_to_p_e_antinu_ground_channel_energetically_allowed':False,'massless_neutrino_is_already_the_most_permissive_case':True},
   'intrinsic_profile_dictionary':'Each original occupation outcome uses number centre X=(x1+x2+x3)/3 and actual native charges. This is a specified translation-invariant graph observable, not a measured relativistic nucleon form factor.',
   'physical_limitations':['absolute masses fail nucleon matching','energy ordering and inner/outer signs do not establish physical neutron/proton identification','no observed beta decay or physical charge radius follows','source graph, number-centre dictionary and equal current weights are supplied']}
