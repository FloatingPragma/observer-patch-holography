"""Portable SVG layout for the authenticated eighteen-channel candidate map."""
from html import escape
import json
import textwrap
from .builder import number, decimal, box

WIDTH=4000
MARGIN=120
INK='#152c40'
MUTED='#506477'
TEAL='#087f80'
BLUE='#285ca5'
AMBER='#a56816'
RED='#a43e39'
LIGHT='#f3f6f8'

# Paths remain relative to the SVG, so a checkout and a GitHub raw download work.
def local(path): return '../../'+path

LABELS={
 'electron':('e⁻','electron','e⁺; Q = +1','Q = −1; target spin 1/2','Charged action mode 0; left and eᶜ CAR components.','Charged infrared threshold and measured observable.'),
 'muon':('μ⁻','muon','μ⁺; Q = +1','Q = −1; target spin 1/2','Charged action mode 1 in the same specified family.','Unstable resonance, weak decay and width.'),
 'tau':('τ⁻','tau','τ⁺; Q = +1','Q = −1; target spin 1/2','Charged action mode 2 in the same specified family.','Unstable resonance, weak decay and width.'),
 'up':('u','up quark','anti-u; Q = −2/3; 3 ↔ 3̄','Q = +2/3; target spin 1/2','Up-type exact current weights; family action unspecified.','Interacting color sector and mass-scheme matching.'),
 'charm':('c','charm quark','anti-c; Q = −2/3; 3 ↔ 3̄','Q = +2/3; target spin 1/2','Same up-type CAR interface; family coordinate unspecified.','Family/mixing action and an interacting observable.'),
 'top':('t','top quark','anti-t; Q = −2/3; 3 ↔ 3̄','Q = +2/3; target spin 1/2','Same up-type CAR interface; no chosen top mass action.','Unstable-quark observable, decay and matching.'),
 'down':('d','down quark','anti-d; Q = +1/3; 3 ↔ 3̄','Q = −1/3; target spin 1/2','Down-type exact current weights; family action unspecified.','Gauge-invariant mass convention and color response.'),
 'strange':('s','strange quark','anti-s; Q = +1/3; 3 ↔ 3̄','Q = −1/3; target spin 1/2','Same down-type CAR interface; family unspecified.','Family/mixing action and an interacting observable.'),
 'bottom':('b','bottom quark','anti-b; Q = +1/3; 3 ↔ 3̄','Q = −1/3; target spin 1/2','Same down-type CAR interface; family unspecified.','Family/mixing action and an interacting observable.'),
 'electron_neutrino':('νₑ','electron flavor','antineutrino channel; Dirac/Majorana unselected','Q = 0; target spin 1/2','Neutral left-doublet CAR component; no νᶜ in M15.','Mass mechanism and flavor-to-mass mixing.'),
 'muon_neutrino':('νμ','muon flavor','antineutrino channel; Dirac/Majorana unselected','Q = 0; target spin 1/2','Same neutral source component; family unspecified.','Oscillation response; flavor is not a mass eigenstate.'),
 'tau_neutrino':('ντ','tau flavor','antineutrino channel; Dirac/Majorana unselected','Q = 0; target spin 1/2','Same neutral source component; family unspecified.','Oscillation response; flavor is not a mass eigenstate.'),
 'photon':('γ','photon','self-conjugate target','Q = 0; target spin 1','Q = T₃ + Y stabilizes a supplied neutral Higgs direction.','Selected vacuum and massless transverse response.'),
 'gluon':('g','gluon multiplet','real adjoint; eight color directions','Q = 0; target spin 1','Eight exact color-current directions, one target channel.','Gauge-invariant color dynamics; no free colored pole.'),
 'W_plus':('W⁺','charged weak boson','W⁻ (separate target card)','Q = +1; target spin 1','Weak raising current T⁺; [Q,T⁺] = T⁺.','Broken-phase gauge action and charged resonance.'),
 'W_minus':('W⁻','charged weak boson','W⁺ (separate target card)','Q = −1; target spin 1','Weak lowering current T⁻; [Q,T⁻] = −T⁻.','Broken-phase gauge action and conjugate resonance.'),
 'Z':('Z','neutral weak boson','self-conjugate target','Q = 0; target spin 1','Neutral span(T₃,Y); mixing coefficient supplied.','Selected mixing, vacuum and neutral resonance.'),
 'Higgs':('h','Higgs radial mode','self-conjugate radial target','Q = 0; target spin 0','Radial response along a supplied neutral doublet.','Selected vacuum, scalar kinetic response and resonance.'),
}

class Canvas:
    def __init__(self):
        self.items=[]
        self.bottom=0
    def rect(self,x,y,w,h,fill='white',stroke=None,rx=20):
        self.items.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}"'+(f' stroke="{stroke}" stroke-width="2"' if stroke else '')+'/>')
        self.bottom=max(self.bottom,y+h)
    def text(self,x,y,text,size=30,color=INK,weight=400,mono=False,anchor='start'):
        if size<24: raise ValueError('atlas text below readable minimum')
        family='monospace' if mono else 'Arial, Helvetica, sans-serif'
        self.items.append(f'<text x="{x}" y="{y}" fill="{color}" font-family="{family}" font-size="{size}" font-weight="{weight}" text-anchor="{anchor}">{escape(str(text))}</text>')
        self.bottom=max(self.bottom,y+size*.25)
    def para(self,x,y,w,text,size=29,color=INK,weight=400,leading=None):
        lead=leading or round(size*1.32)
        # Conservative character budget leaves space for wide glyphs.
        lines=textwrap.wrap(text,width=max(12,int(w/(size*.58))),break_long_words=False,break_on_hyphens=False)
        for line in lines:
            self.text(x,y,line,size,color,weight);y+=lead
        return y
    def line(self,x1,y1,x2,y2,color=TEAL,dashed=False,arrow=False):
        self.items.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="4"'+(' stroke-dasharray="12 9"' if dashed else '')+(' marker-end="url(#arrow)"' if arrow else '')+'/>')
    def link(self,x,y,text,path,size=27,color=BLUE):
        self.items.append('<a href="'+escape(path,quote=True)+'">')
        self.text(x,y,text,size,color)
        self.items.append('</a>')
    def section(self,y,index,title,subtitle):
        self.text(MARGIN,y,index,30,TEAL,700)
        self.text(MARGIN+90,y,title,47,INK,700)
        self.para(MARGIN+90,y+49,WIDTH-2*MARGIN-90,subtitle,28,MUTED)
        return y+112


def mass_text(ref):
    if ref['target'].endswith('neutrino'):return 'SM flavor channel: no individual mass eigenvalue'
    if ref['target']=='photon':return 'SM mass parameter: 0 GeV; experimental bound < 10⁻²⁷ GeV'
    if ref['target']=='gluon':return 'SM benchmark: massless gauge parameter (theoretical)'
    suffix='; 90% CL' if ref.get('confidence_level')=='90%' else ''
    return f"SM input mc²: {ref['value']} ± {ref['uncertainty']} {ref['unit']}{suffix}"


def convention(ref):
    t=ref['target']
    if t in ('up','down','strange'):return 'MSbar running mass at μ = 2 GeV'
    if t in ('charm','bottom'):return 'MSbar running mass at μ = m(μ)'
    if t=='top':return 'Direct event-kinematic extraction; not a pole definition'
    if t=='Z':return 'Breit-Wigner parameter; differs from the pole convention'
    if t in ('W_plus','W_minus'):return 'PDG2026 average excludes the CDF2022 result'
    if t=='Higgs':return 'PDG2026 mass average; quoted scale factor 1.5'
    if t=='tau':return 'PDG2026; CODATA2022 comparison is shown separately'
    if t in ('electron','muon'):return 'CODATA2022; ratios use the published ratio entries'
    if t.endswith('neutrino'):return 'Mass mixing and absolute eigenvalues require separate data'
    return 'Reference conventions do not assign a source particle'


def joined_fields(c,y,data):
    """Finite constructions, with the hypotheses visible beside each arrow."""
    if 'finite_joins' not in data:return y
    y=c.section(y,'03','Three constructive joins on specified sources',
        'These are finite mathematical attachments. The conditional Standard Model benchmark supplies additional physical action, state and matching inputs.')
    w=(WIDTH-2*MARGIN-60)/3
    cards=[
      ('Complete controller → full CAR',
       'A supplied attenuator r retains all sixteen flags, including loss. The complete C-duty inverse also realizes the prescribed positive family matrix Y. For any Y > 0, the full-occupation pairing is real-linear and has norm 2 tr(Y).',
       'Specified cubic kinetic action + Gaussian reference gives E² = c²p² + mⱼ² in its controlled free limit. All eight valleys and the 39 massless internal spectators are retained.',
       'CHARGED_FIELD_ADMISSION.md'),
      ('Quantum Higgs → same Gauss law',
       'The finite quantum-Higgs rank-one vertex has norm 2 and commutes with all twelve Gauss generators. Complete memory feedback admits the same finite gauge/Higgs H₀, without dividing by the active duty probability.',
       'This supplies a finite interacting Hamiltonian join. It supplies no selected Higgs vacuum, interacting chiral continuum, charged pole, lifetime or physical clock.',
       'INTERACTING_CHIRAL_JOIN.md'),
      ('Compact texture → localized fields',
       'The declared opposite-winding texture has a four-state spatial band. Its exact invariant CAR subalgebra admits the same current and Higgs pairing; energy readout uses the paired quadratic ground.',
       'At phase 0: four positive charged frequency groups m_light ± 1/256 and m_heavy ± 1/256, each enclosed with error < 1/25000. Positive Nambu row counts 8,8,4,4; static source-clock gap > 1/40000. These are additive frequencies, not Dirac rest masses or GeV.',
       'DEFECT_FIELD_JOIN.md')]
    height=510
    for i,(title,a,b,proof) in enumerate(cards):
        x=MARGIN+i*(w+30);c.rect(x,y,w,height,'#e9f3f5','#9cc6ca')
        yy=c.para(x+26,y+48,w-52,title,34,TEAL,700)
        yy=c.para(x+26,yy+16,w-52,a,28,MUTED)
        yy=c.para(x+26,yy+16,w-52,b,28,MUTED)
        if yy>y+height-40:raise ValueError('finite-join panel overflow: '+title)
        c.link(x+26,y+height-22,'Owning finite theorem and exact controls',local('model/particle_candidate_atlas/'+proof),25)
    return y+height+70


def matching_table(c,y,data):
    if 'finite_joins' not in data:return y
    rows=data['finite_joins']['measured_mass_join']['branches']
    c.text(MARGIN,y,'Source coefficients and necessary correction budgets',35,TEAL,700)
    c.para(MARGIN,y+47,WIDTH-2*MARGIN,'Displayed decimals are rounded; exact intervals are linked. The Higgs amplitude v is declared. Each correction floor is an infimum needed to reach both quoted CODATA ratio windows with one shared electron scale. It is not a computed interaction, a fit-selected source change or a proof of a physical match.',28,MUTED)
    y+=135
    xs=[MARGIN+24,MARGIN+650,MARGIN+1230,MARGIN+1820,MARGIN+2410,MARGIN+3060]
    heads=['P-law branch','light coefficient / v','middle coefficient / v','heavy coefficient / v','relative matching floor','additive ‖E‖ / mₑ floor']
    c.rect(MARGIN,y,WIDTH-2*MARGIN,64,'#18394b')
    for x,t in zip(xs,heads):c.text(x,y+42,t,25,'white',700)
    names={'forward':'Supplied forward','measured_endpoint':'Measured-alpha endpoint',
           'mssm_beta4_proxy':'MSSM proxy β=4','sm_color_beta7':'SM color β=7'}
    for i,row in enumerate(rows):
        yy=y+64+i*64;c.rect(MARGIN,yy,WIDTH-2*MARGIN,64,'white' if i%2==0 else '#edf3f6',rx=0)
        bounds=row['matching_budgets']
        # All displayed decimals derive from exact outward intervals in audit.json.
        vals=[names[row['branch']]]+[number(v,12) for v in row['masses_over_declared_Higgs_amplitude']]
        floor=bounds['necessary_relative_matching_error_lower']
        additive=bounds['minimal_additive_E_over_source_light_mass_interval'][0]
        vals += [number([floor,floor],11),number([additive,additive],11)]
        for x,t in zip(xs,vals):c.text(x,yy+42,t,27,INK)
    y+=64+4*64+40
    c.para(MARGIN,y,WIDTH-2*MARGIN,'No source-derived correction upper bound is available. A zero floor means the quoted windows admit the source ratios; it does not identify particles. Family-dependent kinetic or self-energy terms can respect gauge symmetry, while their source origin and values require an independent calculation.',28,MUTED)
    return y+130


def calibrated_controller(c,y,data):
    if 'finite_joins' not in data:return y
    row=data['finite_joins']['all_three_input_controller']
    def fmt(key,n):return number([row[key],row[key]],n)
    c.rect(MARGIN,y,WIDTH-2*MARGIN,252,'#e9eef8','#b5c6e3')
    c.text(MARGIN+28,y+49,'Separate three-input calibration: the complete attenuator can encode the charged benchmark',35,BLUE,700)
    c.para(MARGIN+28,y+97,WIDTH-2*MARGIN-56,'Using the electron, muon and PDG tau central masses as three supplied inputs, inverse Fourier algebra selects a pass amplitude, phase and common energy scale. Every one of the sixteen flags is retained. This construction has zero held-out mass predictions.',29,MUTED)
    c.text(MARGIN+28,y+170,
        f"r ≈ {fmt('pass_amplitude',14)};  θ ≈ {row['phase_radians_approx']:.14f};  common scale ≈ {fmt('common_mass_scale_GeV',10)} GeV. Rational-controller errors < 10⁻³⁶ GeV.",29,BLUE)
    c.text(MARGIN+28,y+216,'The scale is not a source-selected Higgs vacuum value. No fitted-point native tape, reference-covariance fit or physical mass-matching theorem is claimed.',27,AMBER)
    return y+310


def obstruction_panels(c,y,data):
    if 'finite_joins' not in data:return y
    y=c.section(y,'07','Exact obstructions and unproved physical attachments',
        'Red statements exclude the written identification in its stated source. They do not exclude every action, quantum field theory or OPH realization.')
    items=[
      ('Exact Koide on the same deformed masses',
       'For χ > 0 and positive roots, Q = (1 + exp(−2χ))/3 < 2/3. The same three coefficients cannot also obey exact balanced Koide. A changed response operator requires its own map.'),
      ('Literal defaults as observed charged ratios',
       'Current-gated C gives 1 : 119.900581 : 896.516443. Coherent phase 0 gives 1 : 1 : (34 + 24√2). Neither directly matches the charged-lepton ratios under one common unit. Unknown interactions are outside this exclusion.'),
      ('One valley by deleting the other seven',
       'The actual sine action has eight zeros of alternating winding. A valley projector is not reducing under the admitted local gauge action. The interacting chiral regulator and mirror treatment require additional control.'),
      ('Fixed quantum link as an identity background',
       'The source link is nilpotent. Its state-vector distance from I is at least (√5−1)/2; its expectation differs from 1 by at least 1/2. A nonzero mean does not remove order-one bidirectional relative fluctuations.'),
      ('Static defect frequency as a particle rest mass',
       'The proved compact branches add Eⱼ ± m. They do not obey the homogeneous Dirac square. Creating their CAR excitation does not create or move the prescribed background winding; no autonomous defect dispersion is derived.'),
      ('Pixel self-consistency selects χ',
       'The two admitted scale laws retain strict self-maps after χ → χ + 10⁻¹⁴. Each perturbed law has a nearby self-consistent root and a disjoint Q interval. Uniqueness for one fixed law does not select that law. No explicit-shift uniqueness or fitted root is claimed.')]
    w=(WIDTH-2*MARGIN-60)/3;h=285
    for i,(title,body) in enumerate(items):
        x=MARGIN+(i%3)*(w+30);yy=y+(i//3)*(h+24)
        c.rect(x,yy,w,h,'#fff1ee','#d9b1ac')
        pos=c.para(x+25,yy+45,w-50,title,31,RED,700)
        pos=c.para(x+25,pos+12,w-50,body,27,MUTED)
        if pos>yy+h-12:raise ValueError('obstruction panel overflow: '+title)
    y+=2*(h+24)+32
    c.rect(MARGIN,y,WIDTH-2*MARGIN,160,'#fff8ed','#e4d5bd')
    c.text(MARGIN+28,y+48,'Amber: physical conditions not established by these finite attachments',33,AMBER,700)
    c.para(MARGIN+28,y+93,WIDTH-2*MARGIN-56,'Selected interacting state and physical clock; chiral gauge continuum; charged thresholds or unstable resonances; decay widths and mixing; neutral masses; measured-observable corrections; an actual source-defect creation map and species identity. The conditional Standard Model benchmark assumes its action, matched renormalized parameters and measured references independently.',28,MUTED)
    return y+220


def render(data):
    c=Canvas();refs={r['target']:r for r in data['reference']['masses']}
    c.rect(0,0,WIDTH,255,'#102b3d',rx=0)
    c.text(120,100,'From complete source records to particle candidates',72,'white',700)
    c.text(120,159,'Exact finite attachments • configured action spectra • measured-reference comparisons',34,'#c3e1e6')
    c.line(122,215,245,215,TEAL,arrow=True);c.text(270,225,'proved attachment under stated premises',27,'white')
    c.line(1050,215,1170,215,AMBER,dashed=True,arrow=True);c.text(1195,225,'candidate identity or physical matching',27,'white')
    c.text(2550,225,'18 target channels  |  15 internal weights  |  0 source identifications',27,'#ffe0a3')
    top=[
      ('01','Complete native records','Twelve repair ports; all failure, disabled and quantum-complement outcomes retained.','Specified source',TEAL),
      ('02','Selected current grammar','The admitted band/current construction supplies exact color, weak and hypercharge currents.','Exact finite algebra',TEAL),
      ('03','Fifteen internal weights','An explicit full-Fock intertwiner attaches odd CAR creators. Three families are configured.','Exact representation',TEAL),
      ('04','Separate charged actions','Coherent Y=C(z)²/3 and capacity-dependent gZ_C are distinct configured operators.','Conditional actions',BLUE),
      ('05','Coherent branch: CAR + Dirac','Coherent coefficient + specified cubic action + Gaussian reference: 3+1 free fields, all eight valleys.','Controlled free limit',BLUE),
      ('06','Measured particle identity','Interacting vacuum, chiral regulator, physical clock and observable matching require independent control.','Physical match unproved',AMBER),
    ]
    w=(WIDTH-2*MARGIN-5*24)/6;y=295;h=410
    for i,(num,title,body,badge,color) in enumerate(top):
        x=MARGIN+i*(w+24)
        c.rect(x,y,w,h,'white','#d7e1e8');c.text(x+24,y+44,num,28,color,700)
        yy=c.para(x+24,y+96,w-48,title,34,INK,700)
        yy=c.para(x+24,yy+16,w-48,body,27,MUTED)
        c.text(x+24,y+h-27,badge,25,color,700)
        if i<5:c.line(x+w+2,y+190,x+w+22,y+190,AMBER if i==4 else TEAL,i==4,True)
    y=c.section(783,'01','The exact internal dictionary','All-left Weyl convention: uᶜ,dᶜ,eᶜ are conjugate components. Internal weights alone do not identify spacetime spin or a natural species.')
    irreps=[('Q_L','(3, 2, Y=1/6)','1',6),('uᶜ','(3̄, 1, Y=−2/3)','-4',3),('dᶜ','(3̄, 1, Y=1/3)','2',3),('L_L','(1, 2, Y=−1/2)','-3',2),('eᶜ','(1, 1, Y=1)','6',1)]
    iw=(WIDTH-2*MARGIN-4*25)/5
    for i,(name,label,q,count) in enumerate(irreps):
        x=MARGIN+i*(iw+25);c.rect(x,y,iw,154,'#e9f3f5')
        c.text(x+25,y+47,name+'   '+label,33,INK,700)
        c.text(x+25,y+94,f'{count} weights; q₆Y = {q}',28,TEAL)
        c.text(x+25,y+132,'source masks: '+', '.join(map(str,sorted(data['source_weight_groups'][q]))),25,MUTED)
    y+=212
    c.rect(MARGIN,y,WIDTH-2*MARGIN,220,'#e9eef8','#b5c6e3')
    c.text(MARGIN+28,y+50,'Complete conditional Standard Model benchmark',38,BLUE,700)
    c.para(MARGIN+28,y+98,WIDTH-2*MARGIN-56,'Supplying the renormalized chiral gauge action, its state, matched renormalized parameters, measured references and a neutrino-mass extension gives the conventional particle benchmark. Its masses, spins, charges and conjugates below are conditional Standard Model inputs or assignments. They are separate from the OPH source calculations.',30,MUTED)
    c.text(MARGIN+28,y+190,'Independent Yukawa parameters are permitted. Exact Koide is not imposed on measured poles; source dynamics do not select these benchmark inputs.',27,BLUE)
    y+=275
    y=c.section(y,'02','Eighteen target channels and their antiparticles','Spin and electric charge below are target/reference labels. The charged free-field theorem supplies model spin 1/2; physical source masses and species identities are unproved.')
    grid_w=(WIDTH-2*MARGIN-2*28)/3;card_h=395
    for idx,row in enumerate(data['candidate_rows']):
        col=idx%3
        if col==0 and idx:y+=card_h+24
        x=MARGIN+col*(grid_w+28);t=row['target'];symbol,name,anti,quantum,source,missing=LABELS[t]
        accent=TEAL if idx<3 else BLUE if idx<12 else '#7a5a95'
        c.rect(x,y,grid_w,card_h,'white','#d6e0e7');c.rect(x,y,9,card_h,accent,rx=4)
        c.text(x+28,y+56,symbol,50,accent,700);c.text(x+160,y+51,name,37,INK,700)
        c.text(x+160,y+89,quantum,26,MUTED)
        yy=c.para(x+28,y+130,grid_w-56,'Conjugate: '+anti,26,MUTED)
        yy=c.para(x+28,yy+8,grid_w-56,mass_text(refs[t]),29,INK,700)
        yy=c.para(x+28,yy+3,grid_w-56,convention(refs[t]),25,MUTED)
        yy=c.para(x+28,yy+10,grid_w-56,'Source: '+source,28,TEAL)
        yy=c.para(x+28,yy+4,grid_w-56,'Required: '+missing,27,AMBER)
        if yy>y+card_h-42:raise ValueError('candidate card overflow: '+t)
        c.text(x+28,y+card_h-20,'Physical source mass unproved • no defect assignment',25,AMBER)
    y+=card_h+76
    c.rect(MARGIN,y,WIDTH-2*MARGIN,185,'#e9eef8','#b5c6e3')
    c.text(MARGIN+28,y+48,'Three neutrino mass eigenfields, distinct from the three flavor cards',35,BLUE,700)
    c.para(MARGIN+28,y+96,WIDTH-2*MARGIN-56,'νₑ, νμ, ντ = U × (ν₁, ν₂, ν₃). The extension specifies the mixing matrix U and symbolic masses m₁, m₂, m₃; no numerical individual masses are assigned here. Oscillations constrain mass-squared differences, while the absolute scale and Dirac/Majorana choice require separate inputs.',28,MUTED)
    c.text(MARGIN+28,y+157,'Mass-eigenstate antiparticles are distinct in a Dirac extension; a Majorana extension has self-conjugate neutral fields. Neither choice is derived by the source dictionary.',27,AMBER)
    y+=240
    y=joined_fields(c,y,data)
    y=c.section(y,'04','Five charged-action comparisons','The electron unit is optional calibration. The coherent phase fits mu/e only; tau is held out. P variants have distinct source-coordinate and architecture ancestry.')
    c.rect(MARGIN,y,WIDTH-2*MARGIN,72,'#18394b')
    xs=[MARGIN+24,MARGIN+650,MARGIN+1090,MARGIN+1550,MARGIN+1970,MARGIN+2400,MARGIN+2830,MARGIN+3280]
    heads=['Configured branch','mu/e action','tau/e action','Δ(mu/e)/u_ref','Δ(tau/e)/u_ref','mu GeV*','tau GeV*','Δtau/u_PDG']
    for x,text in zip(xs,heads):c.text(x,y+47,text,27,'white',700)
    names={'calibrated_coherent':'Coherent: mu/e fitted','forward':'P: supplied forward','measured_endpoint':'P: measured-alpha endpoint','mssm_beta4_proxy':'P: MSSM proxy β=4','sm_color_beta7':'P: SM color β=7'}
    for i,row in enumerate(data['lepton_branches']):
        yy=y+72+i*74;c.rect(MARGIN,yy,WIDTH-2*MARGIN,74,'white' if i%2==0 else '#edf3f6',rx=0)
        comp=row['comparisons']
        values=[names[row['branch']],number(row['ratios'][1],9),number(row['ratios'][2],9),
            number(comp['muon_electron_ratio']['residual_over_quoted_uncertainty'],5),
            number(comp['tau_electron_ratio']['residual_over_quoted_uncertainty'],5),
            number(row['optional_charged_GeV'][1],9),number(row['optional_tau_GeV'],9),
            number(comp['tau_GeV_PDG2026']['residual_over_quoted_uncertainty'],6)]
        for x,text in zip(xs,values):c.text(x,yy+47,text,27,INK,500)
    y+=72+5*74+42
    c.text(MARGIN,y,'CODATA2022: mu/e = 206.7682827 ± 0.0000046; tau/e = 3477.23 ± 0.23.',29,INK)
    c.text(MARGIN,y+42,'*Supplied electron unit = 0.00051099895069 GeV. Tau: CODATA2022 1.77686 ± 0.00012 GeV; PDG2026 1.77693 ± 0.00009 GeV.',29,INK)
    c.text(MARGIN,y+84,'*The supplied electron-unit uncertainty is not propagated in these conditional GeV columns.',27,MUTED)
    c.text(MARGIN,y+126,'Residual / quoted uncertainty is descriptive, not statistical significance. Comparisons share data; no physical theory-error budget is supplied.',28,AMBER)
    y+=202
    y=calibrated_controller(c,y,data)
    y=matching_table(c,y,data)
    y=c.section(y,'05','Why the two Koide constructions differ','For positive roots rⱼ = 1 + ρ cos(δ + 2πj/3), the exact quotient is Q = Σrⱼ²/(Σrⱼ)² = (1 + ρ²/2)/3.')
    half=(WIDTH-2*MARGIN-30)/2
    panels=[
      ('Balanced coherent source',[
        'Y = C(z)²/3; ρ = √2 and Q = 2/3 in the positive-root chamber.',
        'Default z = 1: action ratios 1 : 1 : (34 + 24√2) ≈ 1 : 1 : 67.94. The native default does not give three unequal charged masses.',
        f"The comparison chooses θ ≈ {data['phase_radians_approx']:.12f} from mu/e. The finite theorem does not derive this phase.",
        'All original flags are retained. A generator frequency becomes physical energy only with a clock; repair cost is not automatically a mass.'
      ],BLUE),
      ('Conditional capacity response (MCPR)',[
        'χ = (u²/512)/(1 − u/77); ρ = √2 exp(−χ); δ = 2/9 + ζ(P,u). Thus Q < 2/3 for u > 0 and positive roots.',
        'The complete C-duty inverse realizes gZ_C. Its positive-Y full-CAR lift is conditional. P₂,P₃,P₁₀ overlap; they are not mass eigenprojectors.',
        '512 = 2^(12−3) and 77 = (5+2)(12−1) are finite counts. Selecting these registers as the amplitude law and choosing the 2/9 baseline are stipulated, target-informed inputs.',
        'A distinct 16-flag attenuator realizes supplied r = exp(−χ), with every loss flag retained. The native unread diffusion does not select this controller.'
      ],TEAL)]
    for i,(title,paras,color) in enumerate(panels):
        x=MARGIN+i*(half+30);c.rect(x,y,half,440,'white','#d6e0e7');c.text(x+28,y+48,title,36,color,700)
        yy=y+94
        for p in paras:yy=c.para(x+28,yy,half-56,p,29,MUTED)+15
        if yy>y+438:raise ValueError('Koide panel overflow')
    y+=485
    if 'finite_joins' in data:
        c.rect(MARGIN,y,WIDTH-2*MARGIN,205,'#e9f3f5','#9cc6ca')
        c.text(MARGIN+28,y+49,'Exact native-memory diffusion is not the scalar damping ansatz',35,TEAL,700)
        c.para(MARGIN+28,y+95,WIDTH-2*MARGIN-56,'For H = R + R⁻¹ and J = i(R − R⁻¹), both [D(H),R] and [D(J),R] are nonzero in all four duties: always, C, E, CE. The complete 12×128 raw records and exact all-lag group inverse give this result at p=1/12, q=1/6, gain=1.',28,MUTED)
        c.text(MARGIN+28,y+172,'Scope: the unread h·G₁ coefficient at the fixed default rates. The complete dilation and controlled weak limit are retained; general self-energy maps are unproved.',27,MUTED)
        c.link(WIDTH-840,y+49,'Diffusion proof and exact witnesses',local('model/particle_candidate_atlas/COMPATIBILITY.md'),26)
        y+=265
    else:
        y+=20
    y=c.section(y,'06','What the discrepancy audit establishes','Exact arithmetic constrains the configured inputs. A valid finite theorem with an unmet physical premise remains valid; it does not establish the physical identification.')
    thirds=(WIDTH-2*MARGIN-2*30)/3
    audit_cards=[
      ('Rounding is not the mismatch','Independent 100-digit spectra lie inside the rational boxes. At fixed inputs, P mu/e widths range from about 5×10⁻²³ to 3.3×10⁻¹⁴. Forward and proxy residuals are about 1.8×10⁻⁵ and 1.5×10⁻⁵. Display rounding cannot account for those differences.'),
      ('Input ancestry is part of the claim','Forward and measured-endpoint coordinates are declared decimals; their physical endpoint uncertainty is unspecified. Closure boxes certify roots within chosen census and scale laws. Measured-alpha ancestry and target-informed architecture prevent a prospective source-only claim.'),
      ('Central values are estimates','Exact central agreement is not required. The balanced held-out tau lies inside both quoted windows, so its offset alone requires no correction. Interacting kinetic, self-energy, infrared/resonance and clock errors are uncomputed. No known physical cause or joint significance is claimed.'),
    ]
    for i,(title,body) in enumerate(audit_cards):
        x=MARGIN+i*(thirds+30);c.rect(x,y,thirds,325,'#fff8ed','#e4d5bd')
        yy=c.para(x+26,y+48,thirds-52,title,33,AMBER,700)
        c.para(x+26,yy+15,thirds-52,body,28,MUTED)
    y+=380
    c.para(MARGIN,y,WIDTH-2*MARGIN,'The lightest root is a near-cancellation: small amplitude or phase changes can move heavy/light ratios substantially. A common rescaling cannot repair a ratio. A protected charge labels a whole sector, not a unique excitation or lifetime; unstable μ, τ, W and Z require resonance observables.',28,MUTED)
    y+=120
    y=obstruction_panels(c,y,data)
    y=c.section(y,'08','Protected source objects, without forced species links','The eleven objects belong to specified source sectors. Topological winding, an internal charge and a mass eigenmode are distinct data. No bijection to the eighteen targets is derived.')
    cols=3;dw=(WIDTH-2*MARGIN-2*28)/3;dh=160
    aliases={
      'bare_native_source':'Bare source: scalar commutant only',
      'ordered_exchange_charge':'Ordered exchange: sharp Z₂ charge',
      'selected_spin_charge':'Selected spin cube: internal U(1)',
      'single_occupancy_spin_memory':'Single-occupancy spin code: M₂₅₆',
      'joint_positive_field_canonical_sector':'Positive-field canonical pairs',
      'common_hidden_spin_factor':'Common hidden spin qutrit factor',
      'cube_eta_range':'Cube η symmetry: rank-1820 range',
      'compact_opposite_texture_cores':'Prescribed textures: winding ±1',
      'one_hole_spin_readout':'One-hole repair sector',
      'centered_full_CAR_charge_response':'Centered full-CAR charge response',
      'kinetic_Weyl_nodes_and_free_fields':'Eight kinetic nodes: winding ±1',
    }
    for i,r in enumerate(data['protected_objects']):
        x=MARGIN+(i%cols)*(dw+28);yy=y+(i//cols)*(dh+20)
        c.rect(x,yy,dw,dh,'white','#d6e0e7')
        c.para(x+22,yy+43,dw-44,aliases[r['id']],29,TEAL,700)
        c.para(x+22,yy+92,dw-44,r['kind'].replace('_',' ')+'; particle link unproved.',26,MUTED)
    x=MARGIN+2*(dw+28);yy=y+3*(dh+20)
    c.rect(x,yy,dw,dh,'#f6f0e8','#e4d5bd')
    c.para(x+22,yy+43,dw-44,'Beyond the eighteen targets',29,AMBER,700)
    c.para(x+22,yy+92,dw-44,'No hadron or quantum-graviton mass is assigned by this map.',26,MUTED)
    y+=4*(dh+20)+50
    c.line(MARGIN,y,WIDTH-MARGIN,y,'#bbcad4');y+=59
    c.text(MARGIN,y,'Proofs, source inputs and primary measurements',37,INK,700)
    link_rows=[
      [('Full-Fock current attachment','model/standard_model_source_map/FERMION_ATTACHMENT.md'),('Full-CAR coherent action','model/charged_lepton_response/CAR_KOIDE.md'),('Free relativistic branches','model/charged_lepton_response/RELATIVISTIC_MASS.md'),('Phase and amplitude provenance','model/charged_lepton_response/SOURCE_PHASE_SELECTION.md')],
      [('Capacity law and branches','model/charged_lepton_response/P_DEPENDENCE.md'),('Mass-response error theorem','model/charged_lepton_response/MASS_RESPONSE_BRIDGE.md'),('Finite conditional-premise theorem','code/particles/leptons/derive_charged_face_incidence_conditional_theorem.py'),('Bare/effective compatibility','model/particle_candidate_atlas/COMPATIBILITY.md')],
      [('Complete attenuator','model/particle_candidate_atlas/BARE_EFFECTIVE_CONTROLLER.md'),('Generic full-CAR field join','model/particle_candidate_atlas/CHARGED_FIELD_ADMISSION.md'),('Necessary mass-matching budgets','model/particle_candidate_atlas/MEASURED_MASS_JOIN.md'),('Pixel closure nonselection','model/particle_candidate_atlas/CAPACITY_NON_SELECTION.md')],
      [('Finite quantum-Higgs join','model/particle_candidate_atlas/INTERACTING_CHIRAL_JOIN.md'),('Compact defect-field join','model/particle_candidate_atlas/DEFECT_FIELD_JOIN.md'),('Complete native response','model/particle_candidate_atlas/NATIVE_RESPONSE.md'),('Exact intervals and input hashes','model/particle_candidate_atlas/audit.json')],
      [('Conditional SM benchmark','model/particle_candidate_atlas/CONDITIONAL_STANDARD_MODEL.md'),('PDG electroweak theory','https://pdg.lbl.gov/2026/reviews/rpp2026-rev-standard-model.pdf'),('PDG neutrino mixing','https://pdg.lbl.gov/2026/reviews/rpp2026-rev-neutrino-mixing.pdf'),('PDG quantum chromodynamics','https://pdg.lbl.gov/2026/reviews/rpp2026-rev-qcd.pdf')],
      [('Three-input controller calibration','model/particle_candidate_atlas/CALIBRATED_TRIPLET.md'),('Physical mass-response conditions','model/charged_lepton_response/MASS_RESPONSE_BRIDGE.md'),('Exact original-probability readout','model/charged_lepton_response/PHASE_AND_PREMISES.md'),('Source-to-candidate interfaces','model/charged_lepton_response/CANDIDATE_MAP.md')],
      [('CODATA2022 lepton references',data['reference']['sources']['CODATA2022']),('PDG2026 leptons',data['reference']['sources']['PDG2026_leptons']),('PDG2026 quark conventions',data['reference']['sources']['PDG2026_quarks']),('PDG2026 gauge/Higgs data',data['reference']['sources']['PDG2026_gauge_higgs'])],
    ]
    for row in link_rows:
        y+=48
        for i,(title,path) in enumerate(row):c.link(MARGIN+i*940,y,title,path if path.startswith('https:') else local(path),27)
    y+=55
    prefix='Authenticated canonical input manifest' if data['canonical_child_export_consumed'] else 'LAYOUT PREVIEW: canonical child authentication pending; input manifest'
    c.text(MARGIN,y,f"{prefix}: {data['source_binding_count']} paths. Exact intervals and hashes: atlas audit.json.",26,MUTED)
    c.text(MARGIN,y+42,'Source experiments and large operator proofs are retained evidence. The atlas executes bounded interval, closure and spectral checks; physical assignments are unproved.',26,MUTED)
    height=int(max(y+105,c.bottom+45))
    metadata=escape(json.dumps({'schema':data['schema'],'input_sha256':data['input_sha256'],
        'physical_particle_map':[],'physical_mass_predictions':[]},sort_keys=True))
    return ('<?xml version="1.0" encoding="UTF-8"?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{height}" viewBox="0 0 {WIDTH} {height}" role="img" aria-labelledby="title desc">'
        '<title id="title">Complete source records and eighteen particle candidate channels</title>'
        '<desc id="desc">Exact finite attachments and conditional charged-action comparisons. Physical identities and defect-to-particle assignments are unproved.</desc>'
        '<metadata>'+metadata+'</metadata>'
        '<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto"><path d="M0,0 L0,6 L6,3 z" fill="#a56816"/></marker></defs>'
        f'<rect width="{WIDTH}" height="{height}" fill="{LIGHT}"/>'+'\n'.join(c.items)+'</svg>\n')
