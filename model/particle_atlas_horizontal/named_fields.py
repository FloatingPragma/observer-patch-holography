"""Named-particle mechanisms, drawn from the checked field-action certificate."""
import math


def mechanism(c, family, x, y, w):
    from .drawing import TEAL, BLUE, ORANGE, INK, GOLD, MUTED
    if family in ('photon', 'W_plus', 'Z'):
        c.circle(x+600,y,155,'#fff8e6',GOLD)
        c.text(x+600,y+20,'H₀',78,ORANGE,700,anchor='middle')
        c.text(x+600,y+240,'reference field',40,INK,anchor='middle')
        c.line(x+880,y,x+1800,y,TEAL,7,arrow=True)
        if family=='photon':
            c.text(x+1900,y-130,'Q H₀ = 0',66,TEAL,700)
            c.wave(x+1900,y+40,w-2300,105,3,TEAL)
            c.text(x+1900,y+245,'unchanged background → no Higgs rest term',40,TEAL)
        else:
            c.text(x+1900,y-130,'generator changes H₀',56,ORANGE,700)
            c.wave(x+1900,y+40,w-2300,105,3,ORANGE,True)
            c.text(x+1900,y+245,'quadratic cost → rest frequency',43,ORANGE)
    elif family=='Higgs':
        cx=x+1300
        pts=[(cx+900*t,y-400*t*t+190) for t in [i/50 for i in range(-50,51)]]
        c.path('M'+' L'.join(f'{a},{b}' for a,b in pts),ORANGE,8)
        c.circle(cx+320,y+139,60,ORANGE)
        c.line(cx-230,y-40,cx+500,y-40,ORANGE,5,arrow=True)
        c.text(cx,y+300,'radial oscillation h',48,ORANGE,700,anchor='middle')
        c.line(x+2550,y,x+3350,y-170,BLUE,6,arrow=True)
        c.line(x+2550,y,x+3350,y+170,TEAL,6,arrow=True)
        c.text(x+3450,y-150,'|D H₀|² → W, Z',53,BLUE,700)
        c.text(x+3450,y+190,'Y H₀ → left ↔ right',47,TEAL,700)
    elif family in ('charged_leptons','quarks'):
        for i,(label,col) in enumerate((('L',BLUE),('R',ORANGE))):
            xx=x+1500+i*2700;c.circle(xx,y,180,'white',col)
            c.text(xx,y+30,label,98,col,700,anchor='middle')
        c.line(x+1800,y-60,x+3800,y-60,TEAL,8,arrow=True)
        c.line(x+3800,y+60,x+1800,y+60,TEAL,8,arrow=True)
        c.text(x+2850,y-200,'Higgs background × Yukawa coefficient',43,INK,anchor='middle')
        c.text(x+2850,y+290,'left–right coupling at p = 0 gives a rest frequency',41,TEAL,anchor='middle')
        if family=='quarks':
            for i,col in enumerate(('#b94c45','#4a8752','#4265a3')):
                c.circle(x+350+i*120,y+150,38,col)
    elif family=='gluon':
        cx=x+1500
        for i in range(8):
            a=math.pi*i/4;xx=cx+260*math.cos(a);yy=y+220*math.sin(a)
            c.line(cx,y,xx,yy,TEAL,5);c.circle(xx,yy,48,TEAL)
        c.text(x+2300,y-110,'[Tₐ, Tᵦ] ≠ 0',74,TEAL,700)
        c.text(x+2300,y+50,'color waves interact with color waves',43,INK)
        c.text(x+2300,y+220,'Tₐ H₀ = 0: no Higgs color mass',43,ORANGE)
    else:
        for i,label in enumerate(('νe','νμ','ντ')):
            yy=y-190+i*190;c.text(x+450,yy+20,label,65,BLUE,700)
            for j in range(3):c.line(x+850,yy,x+2500,y-190+j*190,'#a8bdc9',3)
        for j in range(3):
            yy=y-190+j*190;c.text(x+2600,yy+18,f'ν{j+1}',61,TEAL,700)
            c.wave(x+3100,yy,w-3500,37,3+j,TEAL)
        c.text(x+w/2,y+330,'weak flavor → mass phases → weak flavor',43,INK,anchor='middle')


def render(c,data,y):
    from .drawing import W,M,INK,MUTED,TEAL,ORANGE,BLUE,GOLD,href
    c.items.append('<g id="named-field-nature">')
    c.text(M,y,'THE NAMED PARTICLES: WHAT THEIR FIELDS DO',88,INK,700,serif=True)
    c.para(M,y+110,W-2*M,'Gold names identify the comparison particles. The equations describe the declared Standard Model field completion; the finite source-to-field bridges have their own stated domains. No mass reference is silently promoted to a source prediction.',47,INK)
    cw=(W-2*M)/4
    for i,panel in enumerate(data['field_nature']['panels']):
        x=M+i*cw;c.rect(x,y+250,cw-65,650,'#f0ece2',radius=22)
        c.para(x+55,y+340,cw-180,panel['title'],58,INK,700,max_y=y+440)
        yy=c.para(x+55,y+490,cw-180,panel['text'],40,INK,max_y=y+700)
        c.para(x+55,yy+30,cw-180,panel['boundary'],37,ORANGE,max_y=y+865)
    order=('photon','Higgs','W_plus','Z','gluon','charged_leptons','quarks','neutrinos')
    names={'photon':'PHOTON γ','Higgs':'HIGGS BOSON h','W_plus':'W+ AND W− BOSONS','Z':'Z BOSON','gluon':'GLUONS g','charged_leptons':'ELECTRON · MUON · TAU','quarks':'THE SIX QUARKS','neutrinos':'THREE NEUTRINO FLAVORS'}
    members={'photon':'Neutral · spin 1 · two transverse helicities','Higgs':'Neutral · spin 0 · radial scalar excitation','W_plus':'Charges +1 and −1 · spin 1 · three polarizations each','Z':'Neutral · spin 1 · three polarizations','gluon':'Eight color directions · spin 1 · color adjoint','charged_leptons':'electron e− / positron e+\nmuon μ− / antimuon μ+\ntau τ− / antitau τ+','quarks':'up u · charm c · top t: Q = +2/3\ndown d · strange s · bottom b: Q = −1/3\nAll: spin 1/2 · color triplet','neutrinos':'electron neutrino νe\nmuon neutrino νμ\ntau neutrino ντ\nNeutral · spin 1/2 · mixed mass fields'}
    short={'photon':'m² = 0 in the quadratic action','Higgs':'mₕ² = 2 λ v²','W_plus':'mW² = g²v² / 4','Z':'mZ² = (g² + g′²)v² / 4','gluon':'No Higgs-generated gluon mass','charged_leptons':'Mℓ = vYℓ / √2','quarks':'Mu,d = vYu,d / √2','neutrinos':'Dirac Yν or Majorana κ must be specified'}
    for i,family in enumerate(order):
        yy=y+1060+i*1130
        row=next(r for r in data['field_nature']['families'] if r['family']==family)
        c.rect(M,yy,W-2*M,1080,'#fffdf7',radius=22)
        c.rect(M,yy,4800,1080,'#f8e4b7',radius=22)
        c.para(M+75,yy+130,4630,names[family],74,INK,700,max_y=yy+270)
        cy=yy+340
        for line in members[family].split('\n'):
            cy=c.para(M+75,cy,4580,line,44,INK,max_y=yy+690)+15
        c.para(M+75,yy+740,4570,short[family],51,ORANGE,700,max_y=yy+860)
        c.para(M+75,yy+930,4580,row['antiparticle_meaning'],35,INK,max_y=yy+1045)
        x=M+5000
        c.label(x,yy+75,'PROGRAM OPERATION → FIELD MECHANISM',TEAL)
        c.para(x,yy+160,6320,row['source_operation'],42,INK,max_y=yy+380)
        mechanism(c,family,x,yy+670,6300)
        x=M+11600
        c.label(x,yy+75,'WHY THIS KIND OF EXCITATION',TEAL)
        cy=c.para(x,yy+165,7380,row['field_nature'],46,INK,max_y=yy+470)
        c.para(x,cy+80,7380,row['mass_reason'],49,TEAL,700,max_y=yy+850)
        c.link(x,yy+995,'Exact field algebra and hypotheses','model/particle_atlas_horizontal/FIELD_NATURE.md',36)
        x=M+19400
        c.label(x,yy+75,'APPARATUS → PHYSICS RECORD',BLUE)
        cy=c.para(x,yy+165,7650,row['readout'],45,INK,max_y=yy+465)
        c.para(x,cy+90,7650,row['caveat'],42,ORANGE,max_y=yy+840)
        selected=row['links'][-2:]
        for j,link in enumerate(selected):
            c.link(x,yy+925+j*65,link['label'],href(link),34)
    end=y+1060+len(order)*1130+80
    finite_outputs(c,data,end)
    c.items.append('</g>')
    return end+1420


def finite_outputs(c,data,y):
    from .drawing import W,M,INK,TEAL,ORANGE
    report=data['observer_outputs']
    gap=report['neutral_scalar']['occupation_spectral_response']['gap']
    c.text(M,y,'COMPUTED SOURCE OUTPUTS: ONE FINITE PROGRAM AT A TIME',75,TEAL,700,serif=True)
    c.para(M,y+100,W-2*M,'These calculations use declared coefficients without particle-mass data. They establish finite interacting records and controlled errors. Their energies are in program units; identifying any line with a named particle requires the spatial, spectral and detector bridge.',44,INK)
    cards=[
        ('SCALAR GAP FROM COMPLETE RECORDS',f'The 307-state scalar source has two bright levels and 305 retained dark states. For the declared test Y = diag(1,2,3), d = 1, the levels are −7, 8 and 1. Its bright gap is Δ = {gap}.','J₀₁(2t) / [4J₀₁(t)] = cos²(Δt/2)','The same preparation weight cancels. Ground cooling is unnecessary for this protocol; both complete experiments and every outcome are retained.','HIGGS_RECORDS.md'),
        ('CHARGED ADDITION AND REMOVAL','The same source gives 1650 addition states and 36 removal states. An 18 × 18 Gram matrix resolves the addition response while retaining 1614 dark states. A charged probe compensates the inserted or removed charge.','H = V + dN_H + ε_F Q*','Neutral data cannot determine ε_F: this term vanishes in the neutral sector but shifts charged lines. It is an explicit program coefficient.','CHARGED_HIGGS_RECORDS.md'),
        ('FAITHFUL STATE, SAME CLOCK','For the same test Hamiltonian, choose β = ln 2. The Gibbs state retains every one of the 307 levels. A conjugate charged reference makes its purification gauge invariant.','−log ρ = βH + log Z · I','Modular parameter s corresponds to Heisenberg time −βs. This establishes a finite state/clock identity, without assigning physical seconds or GeV.','THERMAL_CLOCK.md'),
        ('CONTROLLED GAUGE-FIELD RECORDS','The declared gauge-coordinate program retains the native currents and all admitted fermion occupations. A spinor/quantum-flux experiment computes every pair of detector records.','cutoff error + certified time-step error','At flux cutoff 30, the exact rational numerical half-L1 error is below 1.155 × 10⁻³⁵. Representation-cutoff error is separate; no interacting chiral continuum is claimed.','FIELD_ROUTE.md'),
        ('READ THE FINITE COUPLINGS','Complete native endpoint records recover ε/|J| = 1. The common scalar source compiles γ = 1/4. Complete scalar counts also determine its pair-vertex norm and number-cost magnitude.','specified program → observable response','These operational coefficients have proved meanings. Matching them to running strong/electroweak couplings or a continuum Higgs phase requires the same-source response and scale map.','GAUGE_INPUT.md'),
        ('INVERT SPECTRA TO COEFFICIENTS','Use the neutral gap Δ and the separation δ of two visible charged lines. Their difference of squares reconstructs a source family coefficient; the additive charged offset cancels.','σⱼ(Y)² = (Δ² − δⱼ²) / (4κ)','κ is 1 for a doublet channel, 2 for a singlet. The actual source coefficient is recovered with original nonzero residues. This is not a physical lepton-mass identification.','BARE_EFFECTIVE_BOUNDARY.md')]
    cw=(W-2*M)/len(cards)
    for i,(title,text,formula,boundary,proof) in enumerate(cards):
        x=M+i*cw;c.rect(x,y+240,cw-65,1080,'#e7f0ec',radius=20)
        c.para(x+60,y+340,cw-190,title,49,TEAL,700,max_y=y+470)
        c.para(x+60,y+530,cw-190,text,45,INK,max_y=y+775)
        c.para(x+60,y+865,cw-190,formula,48,TEAL,700,max_y=y+950)
        c.para(x+60,y+1020,cw-190,boundary,39,ORANGE,max_y=y+1230)
        c.link(x+60,y+1280,'Source operators, complete records and bounds','model/particle_observer_outputs/'+proof,34)
