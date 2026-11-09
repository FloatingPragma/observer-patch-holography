"""Large outer panels: finite defect operations and observer interpretation."""
import math


class ReadableCanvas:
    """Give the outer explanatory panels larger type than the dense atlas."""
    def __init__(self, canvas):self.canvas=canvas
    def __getattr__(self,name):return getattr(self.canvas,name)
    def text(self,x,y,text,size=48,*args,**kwargs):
        return self.canvas.text(x,y,text,size*(1.6 if size<=80 else 1.25),*args,**kwargs)
    def para(self,x,y,width,text,size=48,*args,**kwargs):
        return self.canvas.para(x,y,width,text,size*(1.6 if size<=80 else 1.25),*args,**kwargs)
    def link(self,x,y,label,path,size=36):
        return self.canvas.link(x,y,label,path,size*1.4)
    def bars(self,*args,**kwargs):
        from .drawing import Canvas
        return Canvas.bars(self,*args,**kwargs)


def carrier(c, x, y, scale, data, stage):
    from .drawing import TEAL,ORANGE,LINE,INK
    vertices=data['carrier']['vertices'];faces=data['carrier']['faces']
    projected=[]
    for vertex in vertices:
        a,b,z=vertex['coordinates_display']
        projected.append((x+scale*(.82*a+.57*z),y-scale*(.3*a+.85*b-.43*z),.48*a-.52*b+.7*z))
    for i in sorted(range(len(faces)),key=lambda j:sum(projected[k][2] for k in faces[j])):
        points=[projected[k] for k in faces[i]]
        negative=stage['all_face_holonomy_signs'][i]<0
        c.path('M'+' L'.join(f'{a:.2f},{b:.2f}' for a,b,_ in points)+' Z',
               ORANGE if negative else LINE,4,fill='#eab45d' if negative else '#f1f3ed',opacity=.72)
    for i in stage['defect_faces']:
        points=[projected[k] for k in faces[i]]
        xx=sum(p[0] for p in points)/3;yy=sum(p[1] for p in points)/3
        c.text(xx,yy+18,f'{i}: −1',45,INK,700,anchor='middle')
    for a,b,_ in projected:c.circle(a,b,12,TEAL)


def box(c,x,y,title,subtitle,height):
    from .drawing import INK,MUTED
    c.rect(x,y,9120,height,'#fffdf7',radius=35)
    c.para(x+100,y+150,8880,title,92,INK,700,max_y=y+400)
    c.para(x+100,y+470,8880,subtitle,56,MUTED,max_y=y+790)


def render(c,data,height):
    from .drawing import INK,MUTED,TEAL,BLUE,ORANGE,RED,GOLD,loop_diagram,physics_inset
    c=ReadableCanvas(c)
    h=data['observer_outputs']['holonomy_transport'];q=data['observer_outputs']['consensus_records']
    left=360;right=38360
    c.rect(0,0,10000,height,'#e8e6dd',radius=0)
    c.rect(38000,0,10000,height,'#e2ece8',radius=0)
    for x,title,line in ((left,'WHERE RECORDS BEGIN','FINITE SOURCE EVENTS AND PROTECTED HOLONOMY'),
                         (right,'WHAT OBSERVERS RECONSTRUCT','RECORDS → MEANINGS → EFFECTIVE PHYSICS')):
        c.para(x,240,9150,title,142,INK,700,max_y=630)
        c.para(x,720,9150,line,58,TEAL,700,max_y=940)
    # Left: actual finite source and the distinction between kinds of holonomy.
    box(c,left,1200,'An observer is a finite process','Read data. Apply a specified update. Store a record. Let that record influence the next operation. No surrounding physical space is needed to define these steps.',3100)
    for i,(title,detail) in enumerate((('READ','a local load or memory'),('REPAIR','an allowed finite operation'),('RECORD','the actual outcome'))):
        xx=left+230+i*2930;c.rect(xx,2240,2550,820,'#e6f0ee',TEAL,25)
        c.text(xx+1275,2460,title,84,TEAL,700,anchor='middle')
        c.para(xx+100,2630,2330,detail,59,INK,max_y=2930)
        if i<2:c.line(xx+2580,2640,xx+2860,2640,TEAL,7,arrow=True)
    c.para(left+100,3330,8880,'Example: a repair of loads (1,0) records either a transfer or a hold. The controller uses that record to select its next seam. An address in this program is not automatically a measured spatial position.',59,INK,max_y=3800)
    c.link(left+100,4090,'Literal primitive repair and record policy','model/observer/native.py',48)

    box(c,left,4600,'Go around a loop; compare what comes back','A loop can return the occupied ports to their starting pattern while transporting the internal memory differently. The native example exchanges two memory names.',3700)
    loop_diagram(c,left+150,5900,8800)
    c.text(left+180,6600,'Holonomy = the net transformation around the loop',64,TEAL,700)
    c.para(left+180,6800,8730,'The exchange X has two character sectors, +1 and −1. Allowed native repair maps commute with the protected Q label. They can change other data while preserving that label. This proves a memory invariant; a spatial defect needs its own realization.',60,INK,max_y=7370)
    c.para(left+180,7610,8730,'The coin-summed native load-current contrasts are a gradient. Their face-loop sum is zero. A nontrivial memory holonomy and a nonzero spatial curvature are distinct observables.',57,ORANGE,max_y=7990)
    c.link(left+180,8170,'Native memory protection and complete current check','model/particle_observer_outputs/HOLONOMY_TRANSPORT.md',46)

    box(c,left,8600,'A concrete finite curvature defect','Declare a sign on each of the actual carrier’s 30 links. Multiply the signs around each triangular face. Product +1 is flat; product −1 is a curvature defect. All twenty face records are retained.',3500)
    stage=h['graphical_path']['stages'][1]
    carrier(c,left+2500,10500,650,h,stage)
    c.text(left+4850,10040,'Bf = ∏ Zedge',92,ORANGE,700)
    c.para(left+4850,10290,4050,'A negative product is a mismatch after transport around that loop. Gauge relabeling at a vertex cannot change it.',58,INK,max_y=10860)
    c.para(left+4850,11100,4050,'Carrier-address diagram; this mesh is not a declaration of observer spacetime.',50,MUTED,max_y=11550)
    c.link(left+120,11860,'Full finite link space, gauge orbits and loop observables','model/particle_observer_outputs/HOLONOMY_TRANSPORT.md',46)

    box(c,left,12400,'Create a pair, move it, return the work','A reversible local gate trades one work quantum 2J for two J-cost defects. A move transfers one negative face sign to its neighbor. Bringing the endpoints together can annihilate them and restore the work.',4600)
    stages=h['graphical_path']['stages'];chosen=(stages[0],stages[1],stages[5])
    for i,(stage,label) in enumerate(zip(chosen,('FLAT + WORK','PAIR CREATED','ENDPOINT MOVED'))):
        xx=left+1580+i*2990
        c.text(xx,13450,label,49,TEAL,700,anchor='middle')
        carrier(c,xx,14310,440,h,stage)
        c.text(xx,15300,f"defects: {stage['defect_energy_in_J']}J; work: {stage['work_energy_in_J']}J",47,INK,anchor='middle')
        if i<2:c.line(xx+1210,14300,xx+1720,14300,ORANGE,8,arrow=True)
    c.para(left+150,15700,8800,'The uncoupled defect/work ledger reads 2J throughout creation, motion and annihilation. J is a declared source coefficient. The coupled field needs its separate controller-work account; relativistic mass needs the energy–momentum and detector map.',56,INK,max_y=16160)
    c.text(left+150,16400,'Complete birth readout: 9/25 no pair; 16/25 pair',66,TEAL,700)
    c.para(left+150,16600,8800,'These probabilities follow from the declared (3/5, 4/5) coherent pulse. Both histories, every zero label and outside-Gauss outcomes are retained.',48,MUTED,max_y=16880)

    box(c,left,17300,'Why a local repair cannot erase an odd defect','Fix a region and restrict updates to those preserving its boundary holonomy. Internal edge signs cancel in pairs: the product of face holonomies inside equals the loop product on the boundary.',4600)
    c.rect(left+480,18520,3900,1500,'#e8f1ec',TEAL,60)
    c.text(left+2430,18740,'REGION R',65,TEAL,700,anchor='middle')
    c.circle(left+1800,19350,210,ORANGE);c.text(left+1800,19380,'−1',100,'white',700,anchor='middle')
    c.path(f'M{left+2100},19350 Q{left+3200},18850 {left+3550},19500',TEAL,10,arrow=True)
    c.text(left+2420,19830,'moving inside preserves parity',52,TEAL,anchor='middle')
    c.line(left+4600,19360,left+5200,19360,RED,8,dashed=True,arrow=True)
    c.circle(left+5780,19350,210,ORANGE);c.text(left+5780,19380,'−1',100,'white',700,anchor='middle')
    c.para(left+6400,19030,2400,'Crossing the boundary changes the protected sector.',55,RED,max_y=19700)
    c.text(left+170,20320,'QR = ∏f∈R Bf = ∏edge∈∂R Zedge',85,TEAL,700)
    c.para(left+170,20600,8760,'If QR = −1, a permitted update cannot make every face inside flat. It can move defects or remove pairs. Pair annihilation across the boundary is outside that protecting rule. Gauge symmetry alone does not make every individual defect permanent.',58,INK,max_y=21200)
    c.link(left+170,21650,'Exact protection condition and complete gate truth tables','model/particle_observer_outputs/HOLONOMY_TRANSPORT.md',46)

    box(c,left,22200,'Let the protected defect change a detector law','A declared interaction makes the actual regional holonomy control the scalar vertex: H_Q = (3I − Q)V / 2. The strengths 1 and 2 are program choices. Every field input and birth record is retained.',3600)
    for i,(label,value) in enumerate((('NO PAIR / Q = +1','occupied: 1/4'),('PAIR / Q = −1','occupied: 3/4'))):
        xx=left+200+i*4510;c.rect(xx,23400,4210,650,'#e4efea',TEAL,25)
        c.text(xx+2105,23630,label,61,TEAL,700,anchor='middle')
        c.text(xx+2105,23900,value,78,INK,700,anchor='middle')
    c.para(left+170,24400,8760,'These are conditional detector laws at one shared tick; the complete joint experiment keeps birth weights 9/25 and 16/25. The same controller works with native protected memory in a separate construction. Neither attachment assigns a physical particle name.',60,INK,max_y=25100)
    c.link(left+170,25450,'Exact defect → interaction → complete record attachment','model/particle_defect_attachment/PROOF.md',48)

    box(c,left,26200,'The apparatus and controller also exchange energy','Protected moves preserve the full coupled action. Birth changes Q and hence the scalar coupling. Its work account must include the field, even when the field’s mean energy happens to stay zero.',3000)
    c.text(left+180,27400,'controller work:  0       +√56       −√56',72,TEAL,700)
    c.text(left+180,27700,'original weights:  9/25     8/25       8/25',66,INK)
    c.para(left+180,28070,8740,'These weights belong to a separate two-time energy measurement. That measurement changes the scalar preparation; it is not the occupation experiment above. The uncoupled 2J defect ledger alone does not describe this coupled apparatus.',60,ORANGE,max_y=28700)
    c.link(left+180,29000,'Full preparation-work accounting and measurement disturbance','model/particle_defect_attachment/PROOF.md#preparation-work-is-not-hidden-in-the-defect-cost',45)

    # Right: actual finite instrument, meanings and effective descriptions.
    box(c,right,1200,'An apparatus converts a question into a record','Preparation supplies a state. The selected apparatus couples a property of that state to an orthogonal pointer label. A detector stores one label for this run; the complete law predicts repeated-run frequencies.',3500)
    c.text(right+170,2290,'p(r | preparation, apparatus) = Tr(ρ Er)',89,TEAL,700)
    c.para(right+170,2530,8720,'This is the Born rule for that instrument. Er is the effect implemented by its interaction, and all effects sum to the identity. Changing the apparatus changes the question and can change the probabilities.',62,INK,max_y=3170)
    probs=q['complete_ideal_experiments']['momentum']['original_event_probabilities']
    c.bars(right+600,3960,7900,1100,probs,['0','1','2','3'],TEAL)
    c.text(right+170,4400,'One click is an event; a histogram tests a probability law.',55,INK,700)
    c.link(right+170,4570,'Exact four-bin instrument and all observer records','model/particle_observer_outputs/CONSENSUS_RECORDS.md',46)

    box(c,right,5000,'Two observers can use different local codes','Copy the orthogonal detector record into each observer’s memory. These controlled copies preserve the complete instrument and any untouched quantum reference. They do not clone an unknown quantum state.',3600)
    for i,(name,codes) in enumerate((('OBSERVER A',q['alice_encoding']),('OBSERVER B',q['bob_encoding']))):
        xx=right+220+i*4470;c.rect(xx,6090,4200,1300,'#e7eef5',BLUE,30)
        c.text(xx+150,6280,name,72,BLUE,700)
        c.text(xx+150,6520,'event       0    1    2    3',59,INK)
        c.text(xx+150,6720,'local code  '+'    '.join(map(str,codes)),59,INK)
        c.para(xx+150,6990,3900,'A dictionary maps each local code back to its event meaning.',52,INK,max_y=7280)
    c.para(right+180,7720,8750,'The event key includes the preparation and apparatus setting. “Bin 0” in a position experiment and “bin 0” in a momentum experiment need not mean the same observable.',61,INK,max_y=8200)
    c.link(right+180,8430,'Complete copy unitary and interpretation maps','model/particle_observer_outputs/CONSENSUS_RECORDS.md',46)

    box(c,right,8900,'Compare meanings; retain disagreements','Observers translate their codes, compare assignments and store the comparison. Agreement is an output of those maps. A wrong dictionary can produce a mismatch without changing the detector physics.',3800)
    for i,(label,col) in enumerate((('AGREE',TEAL),('MISMATCH',ORANGE),('FAILURE',RED))):
        xx=right+200+i*2960;c.rect(xx,10050,2670,680,'#f3f0e8',col,24)
        c.text(xx+1335,10430,label,65,col,700,anchor='middle')
    wrong=q['wrong_dictionary_comparison']
    c.para(right+180,11100,8740,'The deliberately swapped dictionary disagrees with probability '+wrong['disagreement']+'. Correcting the dictionary repairs the interpretation; it does not alter the original Born weights or select a preferred click.',60,INK,max_y=11660)
    c.para(right+180,11900,8740,'Compatible local assignments on a finite event cover glue to one shared assignment. Every mismatch and failed transmission is kept. Agreement among a few records does not establish a universal field theory.',57,ORANGE,max_y=12410)

    box(c,right,13000,'An effective theory describes the agreed events','The simulation produces finite records. Observers introduce positions, momenta, particles and spacetime coordinates as a mathematical description of regularities in compatible records.',4900)
    steps=('FINITE EVENT RECORDS','AGREED INTERPRETATION','EFFECTIVE EQUATIONS','PREDICT / COMPARE MORE RECORDS')
    for i,title in enumerate(steps):
        yy=14150+i*680;c.rect(right+300,yy,8520,500,'#e4efea',TEAL,25)
        c.text(right+4560,yy+300,title,64,TEAL,700,anchor='middle')
        if i<3:c.line(right+4560,yy+530,right+4560,yy+635,TEAL,8,arrow=True)
    c.para(right+170,17160,8760,'A 3+1-dimensional relativistic field theory is a proposed effective description, valid where its source, state, observable map and error estimates justify it. The consensus step alone does not establish that limit.',59,INK,max_y=17680)

    box(c,right,18200,'Continuum equations are mathematical idealizations','Each executable carrier, memory and event record is finite. A smooth field or four-dimensional spacetime can organize those records without an additional continuous universe being generated.',3300)
    c.text(right+180,19300,'finite records → controlled approximation → smooth model',60,TEAL,700)
    c.text(right+180,19620,'prediction error ≤ cutoff + volume + readout + matching errors',56,TEAL,700)
    c.para(right+180,20030,8740,'The approximation must control the same complete observables. Increasing resolution, taking a mathematical limit or adopting a coordinate convention does not silently supply chirality, a vacuum, a mass or an experimental identity.',59,INK,max_y=20740)
    c.link(right+180,21220,'Finite event gluing and observable error contraction','model/particle_observer_outputs/CONSENSUS_RECORDS.md',46)

    physics_inset(c,right,21800,9120)
    c.para(right+130,23800,8840,'E_rest = mc² is the rest-energy relation in a justified relativistic sector. E = c|p| describes a massless Lorentz-vacuum mode. These equations organize energy and momentum readouts; they do not identify every repair counter with mass.',57,INK,max_y=24430)
    c.para(right+130,24700,8840,'Uncertainty belongs to the actual pair of instruments and generators: ΔX ΔP ≥ ½ |⟨[X,P]⟩|. The familiar ℏ/2 bound uses the canonical relation on its valid domain. Measuring one does not make the other distribution meaningless.',56,INK,max_y=25340)
    c.text(right+130,25750,'SIMULATION-GENERATED OBSERVABLES',84,TEAL,700)

    box(c,right,26200,'A local clock depends on what is accessible','Specify both the global state and the local observable algebra. Restrict the actual state before calculating its modular generator. A global thermal clock does not automatically restrict to the same local dynamics.',3000)
    c.text(right+180,27400,'global ρ → actual restriction → ρregion → Kregion = −log ρregion',60,TEAL,700)
    c.para(right+180,27800,8740,'For one family of the finite scalar source, Kregion contains an induced pair term with coefficient log(13440/41873). A gauge-invariant probe detects it even when every local number probability is fixed. The term comes from the reduced interacting state.',60,INK,max_y=28450)
    c.para(right+180,28610,8740,'This is an internal family region, not a spatial ball. Geometric clocks, stress and area need their own common-field attachment.',55,ORANGE,max_y=28900)
    c.link(right+180,29000,'Exact regional modular operator and source-admitted entropy response','model/particle_regional_geometry/PROOF.md',44)
