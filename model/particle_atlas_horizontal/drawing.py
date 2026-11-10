"""Deterministic native SVG primitives and scientific explanatory diagrams."""
from html import escape
import math
import textwrap

W = 28000
M = 360
INK = '#142c3c'
MUTED = '#526a78'
PAPER = '#fbf8f0'
TEAL = '#087e85'
ORANGE = '#bb661c'
BLUE = '#3861a5'
RED = '#a13e39'
LINE = '#ccd5d5'
GOLD = '#e7a149'


class Canvas:
    def __init__(self):
        self.items = []
        self.text_boxes = []

    def rect(self, x, y, w, h, fill='white', stroke=None, radius=20):
        self.items.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{radius}" fill="{fill}"'+(f' stroke="{stroke}" stroke-width="3"' if stroke else '')+'/>')

    def text(self, x, y, text, size=48, color=INK, weight=400, anchor='start', serif=False):
        font = 'Georgia, Times New Roman, serif' if serif else 'Arial, Helvetica, sans-serif'
        self.items.append(f'<text x="{x}" y="{y}" font-family="{font}" font-size="{size}" fill="{color}" font-weight="{weight}" text-anchor="{anchor}">{escape(str(text))}</text>')

    def para(self, x, y, width, text, size=48, color=MUTED, weight=400, max_y=None):
        lines = textwrap.wrap(str(text), width=max(10,int(width/(size*.57))), break_long_words=False, break_on_hyphens=False)
        leading = size*1.38
        bottom = y+(len(lines)-1)*leading+size*.3
        if max_y is not None and bottom > max_y:
            raise ValueError(f'text overflow {bottom}>{max_y}: {text[:80]}')
        for line in lines:
            self.text(x,y,line,size,color,weight); y += leading
        self.text_boxes.append((x,bottom,width,str(text)))
        return y

    def path(self, d, color=INK, width=5, fill='none', dashed=False, arrow=False, opacity=1):
        self.items.append(f'<path d="{d}" stroke="{color}" stroke-width="{width}" fill="{fill}" opacity="{opacity}" stroke-linecap="round" stroke-linejoin="round"'+(' stroke-dasharray="22 17"' if dashed else '')+(' marker-end="url(#arrow)"' if arrow else '')+'/>')

    def line(self,x1,y1,x2,y2,color=LINE,width=4,dashed=False,arrow=False):
        self.path(f'M{x1},{y1} L{x2},{y2}',color,width,dashed=dashed,arrow=arrow)

    def circle(self,x,y,r,fill='white',stroke=None,width=5):
        self.items.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}"'+(f' stroke="{stroke}" stroke-width="{width}"' if stroke else '')+'/>')

    def link(self,x,y,label,path,size=36):
        href = path if path.startswith(('https://','#')) else '../../'+path
        self.items.append(f'<a href="{escape(href,quote=True)}" target="_blank">')
        self.text(x,y,label+' ↗',size,BLUE)
        self.items.append('</a>')

    def label(self,x,y,text,color=TEAL):
        self.text(x,y,str(text).upper(),34,color,700)

    def wave(self,x,y,width,height,cycles=4,color=TEAL,packet=False):
        pts=[]
        for i in range(181):
            t=i/180
            amp=math.exp(-((t-.5)/.26)**2) if packet else 1
            pts.append((x+width*t,y+height*amp*math.sin(2*math.pi*cycles*t)))
        self.path('M'+' L'.join(f'{a:.2f},{b:.2f}' for a,b in pts),color,8)

    def bars(self,x,y,width,height,probs,labels,color=TEAL,show_values=True):
        from fractions import Fraction
        n=len(probs);step=width/n
        self.line(x,y,x+width,y,LINE,4)
        for i,(p,label) in enumerate(zip(probs,labels)):
            q=Fraction(p);h=height*float(q)
            self.rect(x+i*step+step*.18,y-h,step*.64,max(h,2),color,radius=8)
            self.text(x+(i+.5)*step,y+56,label,37,MUTED,anchor='middle')
            if show_values:self.text(x+(i+.5)*step,y-h-22,str(q),38,color,700,anchor='middle')


def graph(c,x,y,scale=1):
    """A repair graph is drawn as a graph, not an assumed physical triangle."""
    pts=[(x,y-170*scale),(x-220*scale,y+155*scale),(x+220*scale,y+155*scale)]
    for i in range(3):
        a,b=pts[i],pts[(i+1)%3]
        c.line(*a,*b,TEAL,7,arrow=True)
    for i,(a,b) in enumerate(pts):
        c.circle(a,b,64*scale,'white',TEAL,6)
        c.text(a,b+18*scale,('0','5','11')[i],46*scale,INK,700,anchor='middle')
    return pts


def loop_diagram(c,x,y,w):
    named=((0,2),(1,2),(1,0),(2,0))
    step=w/4
    for i,(a,b) in enumerate(named):
        cx=x+(i+.5)*step;pts=graph(c,cx,y,.6)
        for k,col,txt in ((a,ORANGE,'A'),(b,BLUE,'B')):
            px,py=pts[k];c.circle(px+28,py-30,30,col)
            c.text(px+28,py-17,txt,33,'white',700,anchor='middle')
        if i<3:c.line(cx+step*.34,y,cx+step*.61,y,TEAL,5,arrow=True)
    c.text(x+w/2,y+235,'same occupied ports; the two names swap',42,INK,anchor='middle')
    c.text(x+w/2,y+305,'one history: 1/8; all other outcomes retained',38,TEAL,anchor='middle')
    c.text(x+w/2,y+365,'X exchanges 0/1; its ± characters persist',39,ORANGE,anchor='middle')


def charge_diagram(c,x,y,w):
    for i,(name,col) in enumerate((('+',TEAL),('−',BLUE))):
        xx=x+w*.27+i*w*.46
        c.rect(xx-220,y-210,440,360,'#e5f0ef' if i==0 else '#e9edf6',col,32)
        c.text(xx,y-15,'Q = '+name+'1',72,col,700,anchor='middle')
        c.path(f'M{xx-130},{y+75} C{xx-200},{y-90} {xx+200},{y-90} {xx+130},{y+75}',col,6,arrow=True)
    c.text(x+w/2,y+260,'repair stays within each sector',45,INK,anchor='middle')
    c.text(x+w/2,y+328,'[K(record), Q] = 0',49,TEAL,anchor='middle')


def bridge_diagram(c,x,y,w):
    c.circle(x+w*.23,y,90,'#e5f0ef',TEAL)
    c.text(x+w*.23,y+18,'Q',66,TEAL,700,anchor='middle')
    c.path(f'M{x+w*.32},{y} L{x+w*.65},{y}',ORANGE,8,dashed=True,arrow=True)
    c.circle(x+w*.76,y,125,'#f8e6cf',ORANGE)
    c.text(x+w*.76,y+22,'n',72,ORANGE,700,anchor='middle')
    c.text(x+w/2,y-185,'local realization',62,ORANGE,700,anchor='middle')
    c.text(x+w/2,y+250,'two different invariants',45,INK,anchor='middle')
    c.text(x+w/2,y+320,'an embedding must connect them',43,ORANGE,anchor='middle')


def texture_diagram(c,x,y,w):
    centers=[(x+w*.32,y),(x+w*.77,y)]
    for j,(cx,cy) in enumerate(centers):
        for a in range(12):
            t=2*math.pi*a/12
            r0,r1=(110,230) if j==0 else (230,110)
            c.line(cx+r0*math.cos(t),cy+r0*math.sin(t),cx+r1*math.cos(t),cy+r1*math.sin(t),ORANGE if j==0 else BLUE,6,arrow=True)
        c.circle(cx,cy,65,ORANGE if j==0 else BLUE)
        c.text(cx,cy+20,'+1' if j==0 else '−1',48,'white',700,anchor='middle')
    c.path(f'M{centers[0][0]+80},{y+90} C{x+w*.5},{y+420} {x+w*.6},{y+420} {centers[1][0]-60},{y+90}',LINE,4,dashed=True)
    c.text(x+w/2,y+330,'compact pair; the far field closes',43,INK,anchor='middle')
    c.text(x+w/2,y+390,'opposite winding does not identify an antiparticle',35,ORANGE,anchor='middle')


def mode_diagram(c,x,y,w):
    for i in range(9):
        yy=y-200+i*55;c.line(x+200,yy,x+w-180,yy,LINE,2)
    for i in range(16):
        xx=x+200+i*(w-380)/15;c.line(xx,y-200,xx,y+240,LINE,2)
    c.circle(x+w*.5,y,150,'#fae5c7',ORANGE,3)
    c.wave(x+230,y+20,w-460,205,4,TEAL,True)
    c.text(x+w/2,y+365,'a localized field mode',49,TEAL,anchor='middle')


def field_diagram(c,x,y,w):
    left=x+190;bottom=y+235
    c.line(left,bottom,left+w-360,bottom,MUTED,5,arrow=True)
    c.line(left,bottom,left,y-280,MUTED,5,arrow=True)
    pts=[]
    for i in range(101):
        p=i/100;pts.append((left+(w-450)*p,bottom-125*math.sqrt(1+14*p*p)))
    c.path('M'+' L'.join(f'{a:.1f},{b:.1f}' for a,b in pts),ORANGE,9)
    c.line(left,bottom,left+w-450,bottom-125*math.sqrt(14),TEAL,6)
    c.text(left-65,y-275,'E',50,INK)
    c.text(x+w-180,bottom+50,'|p|',50,INK)
    c.text(left+30,bottom-160,'mc²',49,ORANGE,700)
    c.text(x+w/2,y+375,'E² = c²p² + m²c⁴',56,INK,700,anchor='middle')


def apparatus_diagram(c,x,y,w):
    c.wave(x+100,y,w*.36,110,3,TEAL,True)
    c.path(f'M{x+w*.48},{y-210} L{x+w*.48},{y+210} L{x+w*.68},{y} Z',BLUE,6,'#e5ecf8')
    c.line(x+w*.7,y,x+w*.88,y-130,TEAL,7,arrow=True)
    c.line(x+w*.7,y,x+w*.88,y+130,ORANGE,7,arrow=True)
    c.rect(x+w*.9,y-200,70,400,INK,radius=10)
    c.circle(x+w*.915,y-130,35,GOLD)
    c.text(x+w*.5,y+345,'couple → route → register',47,INK,anchor='middle')


def record_diagram(c,x,y,w):
    for i in range(2):
        xx=x+160+i*w*.48
        c.rect(xx,y-230,w*.36,410,'white',LINE,22)
        c.label(xx+35,y-165,'Observer '+('A' if i==0 else 'B'))
        for j,t in enumerate(('run 24','setting X','outcome 0')):
            c.text(xx+35,y-65+j*75,t,45,INK)
        c.circle(xx+w*.3,y+123,26,TEAL)
    c.line(x+w*.44,y,x+w*.52,y,TEAL,6,arrow=True)
    c.text(x+w/2,y+325,'shared record, conditional state',46,INK,anchor='middle')


DIAGRAMS = [loop_diagram,charge_diagram,bridge_diagram,texture_diagram,mode_diagram,field_diagram,apparatus_diagram,record_diagram]


def render(data):
    c=Canvas();story=data['story'];families=data['families'];measure=data['measurement']
    c.items.append('<g transform="translate(10000,0)">')
    stages=story['stages']
    if len(stages)!=8:raise ValueError('eight hero stages required for the illustrated sequence')
    c.rect(0,0,W,420,INK,radius=0)
    c.text(M,190,'FINITE EVENTS, SHARED RECORDS, EFFECTIVE PHYSICS',135,'#fff9ed',700,serif=True)
    c.text(M,305,'OPH PARTICLE ATLAS  /  THE HORIZONTAL EXPLAINER',48,'#a9d6d4',700)
    c.text(W-M,175,'02',175,GOLD,700,anchor='end',serif=True)
    c.text(W-M,300,'Follow the arrows. Read the assumptions beside them.',46,'white',anchor='end')
    c.text(M,520,'SOLID: a proved step under stated premises',41,TEAL,700)
    c.text(M+4200,520,'DASHED: a physical connection to derive',41,ORANGE,700)
    c.text(M+8700,520,'BLUE: a supplied field model or calibration',41,BLUE,700)
    c.text(M+14300,520,'RED: an excluded shortcut',41,RED,700)
    c.para(M,620,W-2*M,'SIMULATION-GENERATED OBSERVABLES: the executable produces finite events and records. Observers compare their meanings and use field equations to describe the agreed events. The middle panels show candidate descriptions and their dependencies, not an objective continuum world forming.',47,INK)
    col=(W-2*M)/8
    for i,s in enumerate(stages):
        x=M+i*col
        if i%2:c.rect(x-25,770,col,1420,'#f3efe4',radius=0)
        c.text(x+50,930,f'{i+1:02}',125,'#d7cdb7',700,serif=True)
        c.para(x+270,840,col-340,s['title'],65,INK,700,max_y=1005)
        DIAGRAMS[i](c,x+45,1270,col-115)
        status=s.get('status','conditional')
        color=ORANGE if any(k in status.lower() for k in ('unproved','missing','bridge','open')) else TEAL if 'proved' in status.lower() or 'exact' in status.lower() else BLUE
        c.label(x+50,1730,s.get('status_label',status),color)
        explanation=s.get('explanation',s.get('plain_explanation',s.get('text','')))
        c.para(x+50,1820,col-160,explanation,44,INK,max_y=2100)
        links=s.get('proof_links',[])
        selected=links[:1]+[p for p in links if '.pdf' in p['path']][:1]
        if len(selected)<2:selected=links[:2]
        for j,link in enumerate(selected):
            path=link['path']+('#'+link['anchor'] if link.get('anchor') and '#' not in link['path'] else '')
            c.link(x+50,2130+j*52,link['label'][:65],path,32)
        if i<7:
            dashed=story['edges'][i]['status']=='open'
            c.line(x+col-150,1230,x+col+40,1230,ORANGE if dashed else TEAL,8,dashed,True)

    experiment(c,data,2370)
    theorem_rails(c,data,4530)
    family_lanes(c,data,6220)
    end=calibration_map(c,data,6220)
    endpoint_closeup(c,data,end+100)
    from . import named_fields, source_connections, operational_readouts, refinement_readouts, ground_readouts, common_contract
    end=named_fields.render(c,data,end+2460)
    end=source_connections.render(c,data,end+100)
    end=operational_readouts.render(c,data,end+100)
    end=refinement_readouts.render(c,data,end+100)
    end=ground_readouts.render(c,data,end+100)
    end=common_contract.render(c,data,end+100)
    footer(c,data,end+100)
    height=end+1030
    c.items.append('</g>')
    from . import observer_wings
    observer_wings.render(c,data,height)
    source_connections.render_candidate_wings(c,data)
    defs='<defs><marker id="arrow" markerWidth="12" markerHeight="10" refX="10" refY="5" orient="auto" markerUnits="strokeWidth"><path d="M0,0 L10,5 L0,10 Z" fill="context-stroke"/></marker></defs>'
    # Explicit titles and descriptive text support screen readers and text search.
    total_width=W+20000
    head=f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="{total_width}" height="{height}" viewBox="0 0 {total_width} {height}" role="img" aria-labelledby="title desc">'
    desc='<title id="title">OPH particle atlas: from a loop to a detector</title><desc id="desc">A horizontal illustrated dependency map with a protected-memory example, conditional spatial defect and field connections, preparation and quantum measurement diagrams, all Standard Model families, finite theorem premises and calibration-removal strategies. Dashed arrows mark unproved physical connections.</desc>'
    return head+desc+defs+f'<rect width="{total_width}" height="{height}" fill="{PAPER}"/>'+''.join(c.items)+'</svg>\n'


def experiment(c,data,y):
    c.text(M,y,'THE EXPERIMENT: WHAT AN OBSERVER ACTUALLY DOES',82,INK,700,serif=True)
    c.para(M,y+95,W-2*M,'A preparation makes a reproducible state. An apparatus couples a chosen observable to a pointer. A detector leaves a record. Comparing records makes a shared observation. These are physical processes; neither an observer’s wish nor agreement chooses the outcome.',48,MUTED)
    y+=260
    x=M
    c.label(x,y,'PREPARE')
    c.text(x,y+90,'Make the same input',65,INK,700)
    c.rect(x,y+145,2600,1070,'#e7f1ef',radius=30)
    c.rect(x+100,y+250,500,400,INK,radius=50)
    c.text(x+350,y+420,'state',66,'white',700,anchor='middle')
    c.text(x+350,y+520,'recipe',66,'white',700,anchor='middle')
    c.wave(x+660,y+450,1700,160,5,TEAL,True)
    c.text(x+110,y+775,'ψ = (1, 1, 0, 0) / √2',66,INK,700)
    c.para(x+110,y+895,2360,'Supplied four-bin projective instruments. Coherent amplitudes occupy two bins. The two apparatus settings describe alternative experiments on the same preparation.',47,INK)
    c.text(x+110,y+1160,'All displayed probabilities are exact.',40,TEAL,700)
    fork=x+3000
    c.line(x+2650,y+580,fork,y+580,TEAL,8,arrow=True)
    c.path(f'M{fork},{y+580} C{fork+350},{y+580} {fork+120},{y+270} {fork+650},{y+270}',TEAL,7,arrow=True)
    c.path(f'M{fork},{y+580} C{fork+350},{y+580} {fork+120},{y+930} {fork+650},{y+930}',BLUE,7,arrow=True)
    base=x+3750
    m=data['measurement']
    for k,(title,col,probs,lab) in enumerate((('POSITION SETTING',TEAL,m['position_probabilities'],m['position_labels']),('MOMENTUM SETTING',BLUE,m['momentum_probabilities'],m['momentum_labels']))):
        yy=y+k*680
        c.label(base,yy,title,col)
        c.text(base,yy+90,'Resolve places' if k==0 else 'Mix phases, then resolve modes',60,INK,700)
        if k==0:
            for j in range(4):
                c.path(f'M{base+40},{yy+205+j*55} L{base+820},{yy+205+j*55}',col,6)
        else:
            for j in range(4):
                c.path(f'M{base+40},{yy+205+j*55} C{base+280},{yy+205+j*55} {base+500},{yy+370-j*55} {base+820},{yy+370-j*55}',col,6,opacity=.75)
            c.rect(base+280,yy+165,300,270,'#e5eaf6',col,10)
            c.text(base+430,yy+325,'F₄',78,col,700,anchor='middle')
        for j in range(4):
            c.rect(base+870,yy+180+j*70,105,53,col if j!=2 else '#d7dfe3',radius=8)
        c.text(base+1000,yy+220,'detector bins',42,MUTED)
        c.bars(base+1480,yy+435,1750,600,probs,lab,col)
        c.para(base+3400,yy+210,2400,'A place-sensitive pointer records one position bin.' if k==0 else 'A Fourier analyzer routes coherent spatial phases to different pointer bins.',47,INK)
        c.text(base+3400,yy+465,'p = (1/2, 1/2, 0, 0)' if k==0 else 'p = (1/2, 1/4, 0, 1/4)',49,col,700)
        c.line(base+5900,yy+300,base+6250,yy+300,col,7,arrow=True)
        c.text(base+6380,yy+140,'IF '+('x0' if k==0 else 'k0')+' CLICKS',43,col,700)
        conditional=m['after_position_x0']['momentum_probabilities'] if k==0 else m['after_momentum_k0']['position_probabilities']
        c.bars(base+6380,yy+435,1800,650,conditional,m['momentum_labels'] if k==0 else m['position_labels'],ORANGE)
        c.text(base+6380,yy+565,'other setting: uniform outcomes',41,ORANGE)
    endx=M+16400
    c.label(endx,y,'COMPARE RECORDS')
    c.text(endx,y+90,'One event; a shared account',62,INK,700)
    record_diagram(c,endx-40,y+440,4200)
    c.para(endx,y+900,4200,'The state assigned after a click is conditional on that recorded outcome. Keeping every branch preserves the full probability law. Consensus checks compatible records and their history; it does not select a preferred click.',48,INK)
    c.text(M,y+1480,'BORN RULE',45,TEAL,700)
    c.text(M+3600,y+1480,'CONDITION ON A RECORD',45,BLUE,700)
    c.text(M+9600,y+1480,'UNCERTAINTY IS A DISTRIBUTION',45,ORANGE,700)
    c.text(M,y+1570,'p(r) = Tr(ρ Eᵣ),   Σ Eᵣ = I',64,INK)
    c.text(M+3600,y+1570,'ρᵣ = Kᵣ ρ Kᵣ† / p(r),  p(r) > 0',59,INK)
    c.para(M,y+1655,9300,'Finite Born uniqueness uses effect-dependent weights and compatible refinements. Native closed histories supply a positive functional; access to arbitrary complex detector effects is an additional bridge.',35,MUTED)
    c.para(M+9600,y+1550,11100,'A sharply localized state has a spread of momentum outcomes. “Undefined” here means no single sharp value, not absence of a momentum distribution. The continuum bound ΔX ΔP ≥ ħ/2 requires the same canonical X and P and their domain; a finite four-bin illustration does not prove that physical identification.',47,INK)
    links=data.get('measurement_content',{}).get('proof_links',[])
    selected=[links[i] for i in (0,1,4,6)]
    for j,link in enumerate(selected):
        c.link(M+j*5250,y+1780,link['label'][:80],href(link),35)
    physics_inset(c,M+21800,y-170,5400)


def href(link):
    return link.get('path',link.get('url',''))+('#'+link['anchor'] if link.get('anchor') else '')


def physics_inset(c,x,y,w):
    c.rect(x,y,w,1730,'#f0e6d6',radius=30)
    c.text(x+75,y+100,'Why waves? Why mass?',80,INK,700,serif=True)
    c.wave(x+130,y+350,w-260,145,7,TEAL,True)
    c.para(x+75,y+560,w-150,'A state carries amplitudes and relative phases. Alternatives add before detection. Constructive and destructive interference change the probability of a click.',46,INK)
    # Two energy curves share one limiting propagation speed.
    xx=x+230;yy=y+1280;ww=w*.52
    c.line(xx,yy,xx+ww,yy,MUTED,4,arrow=True)
    c.line(xx,yy,xx,yy-430,MUTED,4,arrow=True)
    points=[(xx+ww*i/100,yy-90*math.sqrt(1+16*(i/100)**2)) for i in range(101)]
    c.path('M'+' L'.join(f'{a:.1f},{b:.1f}' for a,b in points),ORANGE,8)
    c.line(xx,yy,xx+ww,yy-360,TEAL,6)
    c.text(xx-70,yy-440,'E',47,INK)
    c.text(xx+ww+15,yy+40,'|p|',47,INK)
    c.text(xx+40,yy-135,'mc²',47,ORANGE,700)
    c.text(xx+ww+180,yy-280,'massive',45,ORANGE,700)
    c.text(xx+ww+180,yy-195,'massless',45,TEAL,700)
    c.para(x+75,y+1435,w-150,'Rest mass sets the zero-momentum energy mc² in a justified relativistic sector and common clock. A repair bill or trapped mode frequency needs an energy dictionary before it determines that mass.',45,INK)
    c.link(x+75,y+1660,'Free-field mass theorem','model/charged_lepton_response/RELATIVISTIC_MASS.md',35)


def theorem_rails(c,data,y):
    c.text(M,y,'WHAT THE FINITE THEOREMS DO',78,INK,700,serif=True)
    c.para(M,y+95,W-2*M,'Each theorem is a bridge with an entry condition. Read across: what the construction supplies → what follows → what must connect to the same physical state and clock.',47,MUTED)
    roles={r['id']:r for r in data['theorems']['roles']}
    g=data['families']['gauge_theorems'][0]
    rails=[
      ('SYMMETRY & INTERACTIONS',g['premises'],g['result'],g['boundary'],g['links'],'u(1) ⊕ su(2) ⊕ su(3)'),
      ('POSITION & MOMENTUM',roles['carrier_rank_three']['input'],roles['carrier_rank_three']['output']+' A separate ordered-spin field construction supplies X and its translation P.',roles['carrier_rank_three']['not_implied'],roles['carrier_rank_three']['links'][:1]+roles['field_position_momentum']['links'][:2],'rank 3 carrier → frame; same-field [X, P] = iħ'),
      ('TIME & ENERGY',roles['source_clocks']['input'],roles['source_clocks']['output'],roles['source_clocks']['not_implied'],roles['source_clocks']['links'],'H = −log P;  ρ_L ∝ P^L;  τ = −Ls')]
    for i,(title,inp,out,bound,refs,formula) in enumerate(rails):
        yy=y+190+i*480
        c.rect(M,yy,W-2*M,450,'#eaf0ee' if i%2==0 else '#eef0f6',radius=15)
        c.label(M+45,yy+60,title)
        c.para(M+45,yy+125,3500,inp,40,INK,max_y=yy+365)
        c.line(M+3800,yy+118,M+4200,yy+118,TEAL,7,arrow=True)
        c.para(M+4380,yy+70,7000,out,43,INK,max_y=yy+245)
        c.text(M+4380,yy+325,formula,49,TEAL,700)
        c.line(M+11600,yy+118,M+12000,yy+118,ORANGE,7,True,True)
        c.para(M+12150,yy+70,7200,bound,43,INK,max_y=yy+290)
        if i==0:
            c.text(M+4380,yy+252,'Finite quantum links + Higgs: [G_a, H] = 0; full fermion currents preserved.',38,TEAL,700)
            c.link(M+12150,yy+330,'Exact Gauss interaction and its chiral-field boundary',
                   'model/particle_candidate_atlas/INTERACTING_CHIRAL_JOIN.md',35)
        if i==2:
            c.text(M+12150,yy+330,'The clock and ordered-spin field here are different declared sources.',38,ORANGE)
        for j,r in enumerate(refs[:3]):c.link(M+4380+j*6500,yy+417,r['label'][:75],href(r),32)
        theorem_icon(c,i,M+22000,yy+200,4500)


def theorem_icon(c,index,x,y,w):
    if index==0:
        for j,(n,name,color) in enumerate(((1,'U(1)',ORANGE),(3,'SU(2)',BLUE),(8,'SU(3)',TEAL))):
            cx=x+450+j*1500
            for k in range(n):
                ang=2*math.pi*k/n;c.circle(cx+100*math.cos(ang),y-30+100*math.sin(ang),27,color)
            c.text(cx,y+180,name,55,color,700,anchor='middle')
        c.text(x+w/2,y-145,'brackets + inner symmetry, not counts alone',37,INK,anchor='middle')
    elif index==1:
        cx=x+700
        for i in range(12):
            a=2*math.pi*i/12;c.circle(cx+160*math.cos(a),y+160*math.sin(a),18,TEAL)
            c.line(cx+160*math.cos(a),y+160*math.sin(a),cx+160*math.cos(a+2*math.pi*5/12),y+160*math.sin(a+2*math.pi*5/12),LINE,2)
        c.line(cx+260,y,x+1700,y,TEAL,6,arrow=True)
        for dx,dy,col in ((550,0,ORANGE),(0,-170,BLUE),(-250,150,TEAL)):
            c.line(x+2200,y,x+2200+dx,y+dy,col,6,arrow=True);c.text(x+2240+dx,y+dy,'x' if dx==550 else 'y' if dx==0 else 'z',43,col)
        c.text(x+3350,y+30,'G = 4 P_slow',44,INK,700,anchor='middle')
    else:
        for i,(label,col) in enumerate((('repair attempts n',TEAL),('modular parameter s',BLUE),('physical duration',ORANGE))):
            yy=y-120+i*120;c.text(x,yy+15,label,40,col,700)
            c.line(x+1650,yy,x+w-80,yy,col,5,i==2,True)
            for j in range(6):c.line(x+1770+j*360,yy-18,x+1770+j*360,yy+18,col,4)


def family_lanes(c,data,y):
    c.text(M,y,'EVERY FAMILY: SOURCE → FIELD → DYNAMICS → OBSERVATION',80,INK,700,serif=True)
    c.para(M,y+100,20700,'These are full dependency chains, including their missing links. The right-hand masses and spins belong to the supplied Standard Model benchmark. No row is a proved identification of an individual OPH defect with a measured species.',47,MUTED)
    starts=[M,M+2600,M+7000,M+11500,M+16600]
    widths=[2300,4000,4100,4700,4200]
    heads=['FAMILY / CONJUGATE','SOURCE FEATURE','FIELD & THEOREM','PROPAGATION / MASS','NAMED PARTICLES: MASS REFERENCE']
    for x,title in zip(starts,heads):c.label(x,y+250,title)
    refs={r['target']:r for r in data['v1']['reference']['masses']}
    for i,lane in enumerate(data['families']['lanes']):
        yy=y+310+i*510
        c.rect(M-35,yy,21100,480,'white' if i%2==0 else '#f1eee6',radius=0)
        c.rect(M-35,yy,14,480,TEAL if i<4 else BLUE,radius=0)
        c.rect(starts[4]-70,yy+5,widths[4]+140,470,'#fff0cf',GOLD,radius=15)
        c.para(starts[0]+20,yy+75,widths[0]-80,lane['title'],53,INK,700,max_y=yy+175)
        targets=lane.get('targets',[])
        syms='  ·  '.join(t.get('symbol',t.get('id','')) for t in targets)
        c.para(starts[0]+20,yy+230,widths[0]-80,syms,58,TEAL)
        spin='/'.join(dict.fromkeys(str(t['spin']) for t in targets))
        charge='/'.join(dict.fromkeys(str(t['charge']) for t in targets))
        color='/'.join(dict.fromkeys(str(t['color']) for t in targets))
        c.para(starts[0]+20,yy+315,widths[0]-80,f'Spin {spin} · Q {charge} · {color}',37,MUTED,max_y=yy+385)
        anti='; '.join(dict.fromkeys(t['antiparticle'] for t in targets))
        if any(t['id'].endswith('neutrino') for t in targets):anti='antineutrino channels; Dirac/Majorana choice'
        c.para(starts[0]+20,yy+425,widths[0]-80,'Anti: '+anti,34,BLUE,max_y=yy+480)
        for j,s in enumerate(lane['stages'][:3]):
            x=starts[j+1]
            c.para(x,yy+64,widths[j+1],s['title'],47,INK,700)
            c.para(x,yy+158,widths[j+1],s.get('text',''),41,MUTED,max_y=yy+330)
            if s.get('formula'):c.para(x,yy+368,widths[j+1],s['formula'],37,BLUE,max_y=yy+399)
            links=s.get('links',[])
            selected=links[:1]+[p for p in links if '.pdf' in p['path']][:1]
            for k,link in enumerate(selected):c.link(x,yy+421+k*49,link.get('label','Proof')[:76],href(link),30)
            c.line(x-260,yy+240,x-70,yy+240,ORANGE if j else TEAL,6,j>0,True)
        x=starts[4];last=lane['stages'][-1]
        ids=lane.get('target_ids',[t['id'] for t in targets])
        if len(ids)<=3:
            for j,t in enumerate(ids):
                r=refs[t]
                names={'electron':'electron (e−)','muon':'muon (μ−)','tau':'tau (τ−)',
                       'W_plus':'W+ boson','W_minus':'W− boson','Z':'Z boson','Higgs':'Higgs boson (H)'}
                name=names.get(t,t+' quark' if t in ('up','down','charm','strange','top','bottom') else t)
                text=(f"{name}: {r['value']} ± {r.get('uncertainty','?')} GeV" if r.get('value') is not None else name+': mass not specified')
                if t.endswith('neutrino'):text=t.replace('_',' ')+': flavor, not a mass eigenstate'
                if t=='photon':text='photon (γ): gauge mass 0; bound < 10⁻²⁷ GeV'
                if t=='gluon':text='gluons (g): 8 color fields; zero gauge mass, no free pole'
                if t in ('up','down','strange'):text+='; MSbar(2 GeV), 90% CL'
                if t in ('charm','bottom'):text+='; MSbar(mq), 90% CL'
                if t in ('electron','muon'):text+='; CODATA2022'
                if t in ('tau','top','W_plus','W_minus','Z','Higgs'):text+='; PDG2026'
                c.text(x,yy+70+j*66,text,41,INK,700)
        c.para(x,yy+265,widths[4],last.get('text',''),37,MUTED,max_y=yy+416)
        links=last.get('links',[])
        if links:
            link=next((p for p in links if '.pdf' in p['path']),links[0]);c.link(x,yy+461,link.get('label','Measurement definition')[:82],href(link),31)


def calibration_map(c,data,y):
    start=22100;width=W-start-M
    c.para(start,y,width,'CALIBRATIONS & DECLARED DESIGN CHOICES',67,INK,700,max_y=y+170)
    c.para(start,y+210,width,'Choose the hardware and program; compute every observer record. Architecture need not be unique or forced. A setting chosen from a measured target calibrates that comparison; a held-out output tests the construction.',43,MUTED,max_y=y+340)
    items=data['families'].get('calibration_paths',[])
    cols=2;cw=width/cols
    for i,item in enumerate(items):
        row,col=divmod(i,cols);xx=start+col*cw;yy=y+390+row*815
        c.rect(xx+15,yy,cw-50,780,'#e9edf6',radius=25)
        title=item.get('title',item.get('name',item.get('id','Input')))
        c.para(xx+65,yy+70,cw-160,title,45,BLUE,700,max_y=yy+150)
        if item.get('tested_result'):
            c.text(xx+65,yy+167,item['tested_result'],28,ORANGE,700)
        c.label(xx+65,yy+200,item.get('input_classification','SUPPLIED'))
        c.para(xx+65,yy+260,cw-160,item['input'],36,INK,max_y=yy+365)
        c.label(xx+65,yy+410,'DERIVE / TEST',ORANGE)
        c.para(xx+65,yy+465,cw-160,item['removal_strategy'],36,INK,max_y=yy+615)
        c.para(xx+65,yy+655,cw-160,'Pass: '+item['success_condition'],32,MUTED,max_y=yy+735)
        if item.get('links'):c.link(xx+65,yy+760,item['links'][0]['label'][:70],href(item['links'][0]),29)
    table_y=y+390+math.ceil(len(items)/cols)*815+50
    c.text(start,table_y,'Charged-ratio cross-checks: the same frozen inputs',48,INK,700)
    c.text(start,table_y+65,'Branch / ancestry',34,MUTED,700)
    c.text(start+2350,table_y+65,'μ / e',34,MUTED,700)
    c.text(start+3850,table_y+65,'τ / e',34,MUTED,700)
    names={'calibrated_coherent':'Koide: μ/e fitted, τ held out','forward':'P: forward (μ/e outside window)',
           'measured_endpoint':'P: measured-α endpoint','mssm_beta4_proxy':'P: MSSM proxy β=4','sm_color_beta7':'P: color β=7'}
    for i,row in enumerate(data['v1']['lepton_branches']):
        yy=table_y+130+i*60
        c.text(start,yy,names.get(row['branch'],row['branch']),35,INK)
        c.text(start+2350,yy,f"{float(row['comparisons']['muon_electron_ratio']['center_decimal']):.9f}",37,BLUE)
        c.text(start+3850,yy,f"{float(row['comparisons']['tau_electron_ratio']['center_decimal']):.9f}",37,BLUE)
    c.para(start,table_y+470,width,'CODATA2022: μ/e = 206.7682827 ± 0.0000046; τ/e = 3477.23 ± 0.23. No physical correction error budget is supplied. Three-input attenuation fits all three masses; it has no held-out prediction.',33,MUTED)
    c.link(start,table_y+615,'Exact intervals, fitted phase and GeV comparison','docs/research/OPH_PARTICLE_CANDIDATE_ATLAS.svg',34)
    return max(table_y+680,11400)


def endpoint_closeup(c,data,y):
    c.items.append('<g id="endpoints">')
    c.text(M,y,'ZOOM IN: WHERE THE STORY STARTS AND ENDS',80,INK,700,serif=True)
    c.para(M,y+100,W-2*M,'A primitive observer is a process that reads a record and changes its next operation. An apparatus is a controlled interaction that turns a property of a prepared state into a stable, comparable record.',48,MUTED)
    left=M;right=14300
    c.rect(left,y+230,13300,2010,'#f1eee5',radius=22)
    c.rect(right,y+230,W-right-M,2010,'#eaf0ee',radius=22)
    c.text(left+60,y+330,'THE START: REPAIR, MEMORY, AND A CREATION BOUNDARY',56,INK,700)
    c.para(left+60,y+430,13100,'The example begins with two occupied ports and a declared memory. The local repair conserves total load. The controller sees the current record and chooses its next seam; it does not need a spatial picture.',42,INK)
    # Literal first native repair: seam (0,1), spectator at port 11.
    xx=left+80;yy=y+760
    c.label(xx,yy-140,'1. READ / REPAIR')
    c.circle(xx+140,yy,85,ORANGE);c.circle(xx+700,yy,85,'white',TEAL)
    c.text(xx+140,yy+15,'1',50,'white',700,anchor='middle')
    c.text(xx+700,yy+15,'0',50,INK,700,anchor='middle')
    c.line(xx+260,yy,xx+570,yy,TEAL,6,arrow=True)
    c.text(xx+90,yy+140,'port 0',37,INK);c.text(xx+650,yy+140,'port 1',37,INK)
    c.para(xx,yy+250,2100,'Balance the chosen endpoints. An odd total has two placements; each complete branch is kept.',40,INK)
    c.text(xx,yy+475,'j = n₀ − n′₀',51,TEAL,700)
    c.text(xx,yy+565,'record = (coin, current j)',38,MUTED)
    xx=left+2500
    c.label(xx,yy-140,'2. KEEP BOTH RECORDS')
    for i,row in enumerate(data['endpoints']['primitive_branches']):
        cy=yy-20+i*260
        c.rect(xx,cy-70,2670,205,'white',radius=15)
        c.text(xx+50,cy,f"coin {row['record'][0]}: (1,0) → ({','.join(map(str,row['output_seam_loads']))})",47,INK,700)
        c.text(xx+50,cy+80,f"record ({','.join(map(str,row['record']))}); probability {row['weight']}",38,TEAL)
    c.para(xx,yy+540,2670,'Positive current selects seam (0,5); a hold selects (0,1). This is the literal two-event policy.',39,INK)
    xx=left+5600
    c.label(xx,yy-140,'3. TRANSPORT A MEMORY')
    c.rect(xx,yy-70,2750,470,'#fffaf1',radius=16)
    c.text(xx+70,yy+40,'moving loop history → X',52,TEAL,700)
    c.text(xx+70,yy+150,'|+⟩ → |+⟩',55,TEAL)
    c.text(xx+70,yy+255,'|−⟩ → −|−⟩',55,BLUE)
    c.text(xx+70,yy+355,'weight 1/8; other outcomes retained',36,ORANGE)
    c.para(xx,yy+490,2770,'The illustrated three-seam loop uses a separate supplied schedule. It exposes a sign of the coherent memory; every repair preserves its Q sector.',40,INK)
    xx=left+8900
    c.label(xx,yy-140,'4. LOCALIZE / CREATE?',ORANGE)
    c.line(xx,yy+130,xx+520,yy+130,ORANGE,8,True,True)
    c.circle(xx+1150,yy+130,165,'#f6ddba',ORANGE)
    c.text(xx+1150,yy+148,'+1',65,ORANGE,700,anchor='middle')
    c.circle(xx+1950,yy+130,165,'#e1e8f6',BLUE)
    c.text(xx+1950,yy+148,'−1',65,BLUE,700,anchor='middle')
    c.para(xx,yy+420,3500,'The localized winding pair is a supplied field background. A source process that creates this pair from the prepared native state has no proof here.',43,ORANGE)
    c.para(left+80,y+1600,12950,'Protection says what cannot be erased by the allowed updates. It does not show how a spatial core was born. Creating a fermion in a bound mode fills that mode; it does not create the background winding. A drawing cannot identify these two operations.',46,INK,max_y=y+1820)
    c.text(left+80,y+1905,'[K_record,Q] = 0  ⇒  no transition between sharp Q sectors',57,TEAL,700)
    c.link(left+80,y+2070,'Actual primitive observer and record controller','model/observer/native.py',37)
    c.link(left+6600,y+2070,'Memory loop, protection and spatial boundary','model/defect_protection/README.md',37)
    c.link(left+80,y+2160,'Why occupying a field mode does not create the winding','model/particle_atlas_horizontal/STORY.md',37)
    # A supplied two-path / four-screen-port analogue of a double slit.
    c.text(right+70,y+330,'THE END: TWO PATHS, ONE SCREEN RECORD',57,INK,700)
    c.para(right+70,y+430,W-right-M-150,'A double-slit-style analogue uses two coherent paths and four screen bins. Both settings prepare (|0⟩+|1⟩)/√2. The apparatus can leave the paths indistinguishable or attach a distinct memory tag to each.',42,INK)
    x=right+100;cy=y+1120
    c.circle(x+140,cy,80,TEAL);c.text(x,cy+180,'prepare',38,INK)
    c.path(f'M{x+250},{cy} Q{x+650},{cy-300} {x+1200},{cy-300}',TEAL,7,arrow=True)
    c.path(f'M{x+250},{cy} Q{x+650},{cy+300} {x+1200},{cy+300}',TEAL,7,arrow=True)
    c.rect(x+1270,cy-530,80,180,INK,radius=0)
    c.rect(x+1270,cy-240,80,480,INK,radius=0)
    c.rect(x+1270,cy+350,80,180,INK,radius=0)
    c.text(x+1480,cy-285,'path 0',38,INK);c.text(x+1480,cy+315,'path 1',38,INK)
    c.wave(x+2100,cy-300,1600,65,4,TEAL)
    c.wave(x+2100,cy+300,1600,65,4,TEAL)
    c.rect(x+4050,cy-400,1900,800,'white',BLUE,22)
    c.text(x+5000,cy-155,'PATH MARKER',46,BLUE,700,anchor='middle')
    c.text(x+5000,cy-30,'same tag / distinct tags',40,INK,anchor='middle')
    c.text(x+5000,cy+135,'then Fourier mixer F₄',43,BLUE,anchor='middle')
    c.line(x+6100,cy,x+6630,cy,TEAL,7,arrow=True)
    h=data['double_path'];labels=h['screen_labels']
    c.text(x+6950,cy-330,'No path record: amplitudes interfere',42,TEAL,700)
    c.bars(x+6950,cy+20,4400,520,h['unmarked']['screen_probabilities'],labels,TEAL)
    c.text(x+6950,cy+320,'Distinct path tags: paths add as probabilities',40,ORANGE,700)
    c.bars(x+6950,cy+670,4400,520,h['orthogonal_path_tags']['screen_probabilities'],labels,ORANGE)
    c.para(x,y+1770,6250,'An unread tag is enough. The particle and marker become correlated. Ignoring that marker removes the cross term, so the screen pattern changes without a person looking.',44,INK,max_y=y+1900)
    c.para(right+70,y+1940,W-right-M-150,'The detector couples a screen bin to a pointer and writes a record: run, preparation, apparatus setting, path tag if present, and screen bin. With tags, all eight (path, screen) records have probability 1/8. A particular run stores one outcome; repeated runs test the Born distribution.',43,INK,max_y=y+2110)
    c.text(right+70,y+2140,'Supplied four-bin apparatus; no laboratory particle identity is inferred.',36,ORANGE,700)
    c.link(right+7200,y+2140,'Exact states, record laws and partial which-path information','model/particle_atlas_horizontal/MEASUREMENT.md',32)
    c.items.append('</g>')


def footer(c,data,y):
    c.rect(0,y-30,W,930,INK,radius=0)
    c.text(M,y+100,'READ THE MAP WITHOUT LOSING THE SCIENCE',75,'white',700,serif=True)
    panels=[('A wave is a coherent state','Alternative amplitudes can interfere before a record distinguishes them. The detector samples the probabilities of its specified instrument; agreement does not supply a missing amplitude law.'),('A mass is an energy relation','In a relativistic massive sector, energy remains at zero momentum. Repair effort becomes a mass only after deriving its Hamiltonian, common clock, normalization and measured observable.'),('A theorem has a domain','Finite internal symmetry, a protected memory bit, a winding background and a continuum particle are different objects. The bridges are named so they can be proved or rejected on the same source.')]
    cw=(W-2*M)/3
    for i,(title,text) in enumerate(panels):
        x=M+i*cw;c.text(x,y+240,title,61,GOLD,700)
        c.para(x,y+330,cw-150,text,46,'#d4e0e4')
    c.text(M,y+650,f"Authenticated v1 scientific inputs: {len(data['v1']['input_sha256'])}  •  No source mass or species claim is promoted by this visual format.",40,'#c4d4db')
    c.link(M,y+760,'Numerical atlas and all mass conventions','docs/research/OPH_PARTICLE_CANDIDATE_ATLAS.svg',40)
    c.link(M+7000,y+760,'Exact illustration and provenance','model/particle_atlas_horizontal/README.md',40)
    c.link(M+14000,y+760,'Full dependency and assumption map','model/particle_atlas_horizontal/STORY.md',40)
