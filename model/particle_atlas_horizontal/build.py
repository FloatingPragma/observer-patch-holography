"""Reproducible horizontal SVG, pan/zoom reader and bounded illustration audit."""
import argparse
import hashlib
from html import escape
import json
from pathlib import Path
from urllib.parse import unquote
import xml.etree.ElementTree as ET

from . import drawing, measurement, endpoints

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SVG = ROOT/'docs/research/OPH_PARTICLE_ATLAS_V2.svg'
HTML = ROOT/'docs/research/OPH_PARTICLE_ATLAS_V2.html'
AUDIT = HERE/'audit.json'
V1 = ROOT/'model/particle_candidate_atlas/audit.json'
OUTPUTS = ROOT/'model/particle_observer_outputs/results.json'
CONTENT = ('story.json','families.json','measurement_content.json','theorems.json','field_nature.json')


def canonical(value):
    return json.dumps(value,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+'\n'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def links(value):
    if isinstance(value,dict):
        if 'path' in value and 'label' in value:yield value
        for v in value.values():yield from links(v)
    elif isinstance(value,list):
        for v in value:yield from links(v)


def load():
    data={'v1':json.loads(V1.read_text()),'measurement':measurement.certificate(),
          'endpoints':endpoints.certificate(), 'double_path':measurement.double_path_certificate()}
    for name in CONTENT:data[Path(name).stem]=json.loads((HERE/name).read_text())
    data['observer_outputs']=json.loads(OUTPUTS.read_text())
    for rel,expected in data['observer_outputs']['input_sha256'].items():
        if sha(ROOT/rel)!=expected:raise ValueError('stale observer-output input: '+rel)
    pins=data['v1']['input_sha256']
    for rel,expected in pins.items():
        path=(ROOT/rel).resolve()
        if not path.is_relative_to(ROOT) or sha(path)!=expected:
            raise ValueError('stale v1 scientific input: '+rel)
    if data['v1']['physical_mass_predictions'] or data['v1']['physical_particle_map']:
        raise ValueError('v1 physical scope changed; editorial review required')
    validate_content(data)
    paths=set(pins)|set(data['observer_outputs']['input_sha256'])|{str(V1.relative_to(ROOT)),str(OUTPUTS.relative_to(ROOT))}
    paths.update(str(p.relative_to(ROOT)) for p in HERE.iterdir()
                 if p.is_file() and (p.suffix in ('.py','.md','.js') or p.name in CONTENT))
    for link in links([data[n.removesuffix('.json')] for n in CONTENT]):
        name=unquote(link['path'].split('#')[0])
        if name.startswith('https://'):continue
        path=(ROOT/name).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file():
            raise ValueError('missing or escaped reference: '+name)
        paths.add(name)
    data['input_sha256']={p:sha(ROOT/p) for p in sorted(paths)}
    return data


def validate_content(data):
    ids=[t for lane in data['families']['lanes'] for t in lane['target_ids']]
    expected={r['target'] for r in data['v1']['candidate_rows']}
    if len(ids)!=18 or len(set(ids))!=18 or set(ids)!=expected:
        raise ValueError('complete elementary target coverage failed')
    named=data['field_nature']['families']
    if len(named)!=18 or {r['target'] for r in named}!=expected:
        raise ValueError('named-field coverage failed')
    if data['field_nature']['physical_particle_assignments'] or any(
            r['source_predicted_physical_mass'] is not None for r in named):
        raise ValueError('named-field explanation promoted a physical prediction')
    if data['observer_outputs']['physical_particle_assignments'] or data['observer_outputs']['source_predicted_GeV_masses']:
        raise ValueError('finite output explanation promoted a physical prediction')
    expected_stages=('loop','charge','local_bridge','texture','mode','field','apparatus','record')
    if tuple(s['id'] for s in data['story']['stages'])!=expected_stages:
        raise ValueError('hero identity chain changed')
    edges=data['story']['edges']
    if [(r['from'],r['to']) for r in edges]!=list(zip(expected_stages,expected_stages[1:])):
        raise ValueError('hero edges must cover the same sequence')
    if any(edges[i]['status']!='open' for i in (1,2,4)):
        raise ValueError('unproved physical bridge promoted')
    for lane in data['families']['lanes']:
        if len(lane['stages'])!=4 or any(t['source_physical_identification'] for t in lane['targets']):
            raise ValueError('family physical scope changed')
    for row in data['families']['calibration_paths']:
        if not all(row.get(k) for k in ('input','removal_strategy','success_condition','links')):
            raise ValueError('calibration path lost its required input or test')


def reader(svg,data):
    full_links=[];seen=set()
    extra=[{'label':label,'path':'model/particle_observer_outputs/'+name} for name,label in (
        ('HOLONOMY_TRANSPORT.md','Finite holonomy birth, work exchange and regional protection'),
        ('CONSENSUS_RECORDS.md','Record copying, observer meanings and finite event gluing'),
        ('BARE_EFFECTIVE_BOUNDARY.md','Same-clock spectral coefficient inverse'),
        ('CONTROLLER_READOUT.md','Controller attenuation and phase-class readout'),
        ('GAUGE_INPUT.md','Operational source couplings and conditional gauge ray'))]
    extra += [{'label':label, 'path':'model/'+owner+'/PROOF.md'} for owner,label in (
        ('particle_defect_attachment','Protected memory and holonomy coupled to complete field records'),
        ('particle_field_causality','Same-field causal influence and complete local-history bounds'),
        ('particle_gravity_response','Same-action energy transfer, work and geometric premises'))]
    extra += [{'label':label, 'path':'model/'+owner+'/PROOF.md'} for owner,label in (
        ('particle_scalar_ground','Full-carrier scalar ground, singlet preparation and complete occupation records'),
        ('particle_regional_geometry','Induced regional modular pair interaction and entropy response'),
        ('particle_field_preparation','Finite gauge-word preparation and same-clock Wilson-gap bounds'),
        ('particle_ground_charged','Complete charged spectra about the full scalar ground'),
        ('particle_prepared_field_readout','One prepared spatial field with all 61 local detector outcomes'),
        ('particle_local_gauge_state','Mixed gauge-state family with fixed-parameter local roughness bounds'),
        ('particle_local_wilson_density','Local Wilson low-mode weight and one-body sign approximation'),
        ('particle_local_history_family','Complete polynomial field records at fixed nonzero time and growing volume'),
        ('particle_sign_response_boundary','Native inverse-gap electric response and the spectral shortcut boundary'),
        ('particle_smooth_sign_field','Smooth interacting action with complete local quantum records'),
        ('particle_joint_translation_readout','Joint field/reference translations, position backaction and finite uncertainty'),
        ('particle_physical_causality','Refinement-uniform physical causal bound for complete local records'),
        ('particle_field_state_response','Complete interacting thermal preparation and flux records'),
        ('particle_geometric_response','Same-action spatial strain, work and pressure with declared geometry'),
        ('particle_spatial_ground_state','Gap-free full spatial low-energy preparation'),
        ('particle_ground_spectral_records','Charge-conjugate thresholds and complete source-energy records'),
        ('particle_gauge_characteristics','Native trace-dual gauge and matter characteristics'),
        ('particle_common_contract','One selected source, state, reference and clock'),
        ('particle_gauss_ground_state','Physical Gauss preparation on the full spatial carrier'),
        ('particle_charged_residue','Original charged-threshold weights'),
        ('particle_spatial_modular_response','Physical regional support and modular response'),
        ('particle_link_cutoff_limit','Fixed-volume electric cutoff and complete-record convergence'),
        ('particle_physical_boundary_state','Physical preparation with boundary charges and complete records'),
        ('particle_yang_mills_gap_contract','Same-clock local correlations and a continuum spectral-gap criterion'))]
    extra.append({'label':'All-time gap control from conserved energy and its refinement boundary',
                  'path':'model/particle_field_preparation/REFINEMENT.md'})
    for item in links([data[k] for k in ('story','families','measurement_content','theorems','field_nature')]+[extra]):
        href=item['path']+('#'+item['anchor'] if item.get('anchor') else '')
        if href in seen:continue
        seen.add(href)
        url=href if href.startswith('https://') else '../../'+href
        full_links.append(f'<li><a href="{escape(url,quote=True)}" target="_blank">{escape(item["label"])}</a></li>')
    return '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>OPH particle atlas v2: illustrated reader</title>
<style>*{box-sizing:border-box}body{margin:0;height:100vh;display:flex;flex-direction:column;overflow:hidden;background:#fbf8f0;color:#142c3c;font:15px Arial,sans-serif}header{min-height:76px;flex-shrink:0;background:#142c3c;color:white;display:flex;gap:10px;align-items:center;padding:12px 18px;flex-wrap:wrap;position:relative;z-index:3}header strong{font:22px Georgia,serif;margin-right:12px}button,a.download{border:1px solid #789ba7;background:#244454;color:white;border-radius:6px;padding:9px 13px;cursor:pointer;text-decoration:none;font-size:14px}button:hover,a.download:hover{background:#397078}#hint{font-size:12px;color:#b9d2db;margin-left:auto}#viewport{flex:1;min-height:0;overflow:hidden;touch-action:none;cursor:grab;position:relative}#art{transform-origin:0 0;position:absolute;left:0;top:0}#art svg{display:block;max-width:none}#zoom{min-width:48px;text-align:center}#refs{display:none;position:fixed;right:0;top:0;width:min(560px,90vw);height:100vh;overflow:auto;background:#fffaf0;border-left:1px solid #ccd5d5;z-index:4;padding:22px}#refs a{color:#285b9f}#refs li{margin:11px 0}#refs h2{font:26px Georgia,serif}#nav{position:fixed;bottom:20px;right:22px;z-index:3;display:flex;gap:8px}#nav button{box-shadow:0 3px 12px #0002}@media(max-width:850px){header{min-height:118px}#hint{display:none}}@media print{header,#nav,#refs{display:none}#viewport{overflow:visible;height:auto}#art{position:static;transform:none!important}#art svg{width:100%;height:auto}}</style>
<header><strong>Simulation-generated observables</strong><button data-row="origin">Defect origins</button><button data-row="story">Derivation chain</button><button data-row="observers">Observers & meaning</button><button data-row="experiment">The experiment</button><button data-row="theorems">Theorem bridges</button><button data-row="endpoints">Endpoints & slits</button><button data-row="families">Every family</button><button data-row="particles">Mass references</button><button data-row="fieldnature">Photon, Higgs & matter</button><button data-row="requirements">Six requirements</button><button data-row="readouts">Actual field readouts</button><button data-row="contract">One program</button><button data-row="ground">Ground states &amp; cones</button><button data-row="refinement">Refinement & pressure</button><button data-row="calibrations">Remove calibrations</button><button id="overview">Overview</button><button id="references">Proofs & papers</button><a class="download" href="OPH_PARTICLE_ATLAS_V2.svg" download>SVG ↓</a><span style="flex-basis:100%;font-size:13px;color:#f2d79a">Named fields: photon · Higgs · W and Z · gluons · electron, muon, tau · quarks · neutrinos. Gold boxes show comparison particles; arrows show the construction and its proof boundaries.</span><span id="hint">Drag to pan · wheel to zoom · click blue proof links</span></header>
<main id="viewport" aria-label="Draggable and zoomable particle atlas"><div id="art">'''+svg+'''</div></main><aside id="refs"><button id="closeRefs">Close</button><h2>Exact proofs and paper sections</h2><p>Full references behind the illustrated steps. Supplied hypotheses and unproved bridges are labeled in the spread.</p><ul>'''+''.join(full_links)+'''</ul></aside>
<nav id="nav"><button id="left" aria-label="Pan left">←</button><button id="right" aria-label="Pan right">→</button><button id="out" aria-label="Zoom out">−</button><button id="in" aria-label="Zoom in">+</button><button id="zoom">100%</button></nav>
<script>
const viewport=document.querySelector('#viewport'),art=document.querySelector('#art');let scale=.25,x=0,y=0,drag=null;
function draw(){art.style.transform=`translate(${x}px,${y}px) scale(${scale})`;document.querySelector('#zoom').textContent=Math.round(scale*100)+'%'}
function zoom(f,cx=viewport.clientWidth/2,cy=viewport.clientHeight/2){const b=art.querySelector('svg').viewBox.baseVal;const floor=Math.min(.035,.97*viewport.clientWidth/b.width,.97*viewport.clientHeight/b.height);const next=Math.min(1.8,Math.max(floor,scale*f));x=cx-(cx-x)*next/scale;y=cy-(cy-y)*next/scale;scale=next;draw()}
function row(name){const places={origin:[0,0],observers:[38000,0],story:[10320,680],experiment:[10320,2320],theorems:[10320,4500],families:[10320,6150],particles:[26800,6150],calibrations:[32000,6150],endpoints:[10320,11500],fieldnature:[10320,13900],requirements:[10320,25400],readouts:[10320,29200],refinement:[10320,33800],ground:[10320,38950],contract:[10320,43850]};scale=Math.min(.43,viewport.clientWidth/(['origin','observers'].includes(name)?9900:5600));x=35-places[name][0]*scale;y=35-places[name][1]*scale;draw()}
viewport.addEventListener('wheel',e=>{e.preventDefault();const r=viewport.getBoundingClientRect();zoom(Math.exp(-e.deltaY*.0015),e.clientX-r.left,e.clientY-r.top)},{passive:false});
viewport.addEventListener('pointerdown',e=>{if(e.target.closest('a'))return;drag={px:e.clientX,py:e.clientY,x,y};viewport.setPointerCapture(e.pointerId);viewport.style.cursor='grabbing'});
viewport.addEventListener('pointermove',e=>{if(!drag)return;x=drag.x+e.clientX-drag.px;y=drag.y+e.clientY-drag.py;draw()});
function stop(){drag=null;viewport.style.cursor='grab'}viewport.addEventListener('pointerup',stop);viewport.addEventListener('pointercancel',stop);
document.querySelectorAll('[data-row]').forEach(b=>b.onclick=()=>row(b.dataset.row));document.querySelector('#in').onclick=()=>zoom(1.3);document.querySelector('#out').onclick=()=>zoom(1/1.3);document.querySelector('#left').onclick=()=>{x+=viewport.clientWidth*.72;draw()};document.querySelector('#right').onclick=()=>{x-=viewport.clientWidth*.72;draw()};
document.querySelector('#overview').onclick=()=>{const b=art.querySelector('svg').viewBox.baseVal;scale=Math.min(viewport.clientWidth/b.width,viewport.clientHeight/b.height)*.97;x=(viewport.clientWidth-b.width*scale)/2;y=8;draw()};
document.querySelector('#references').onclick=()=>document.querySelector('#refs').style.display='block';document.querySelector('#closeRefs').onclick=()=>document.querySelector('#refs').style.display='none';
document.addEventListener('keydown',e=>{if(e.key==='ArrowRight'){x-=viewport.clientWidth*.5;draw()}if(e.key==='ArrowLeft'){x+=viewport.clientWidth*.5;draw()}if(e.key==='Escape')document.querySelector('#refs').style.display='none'});row('origin');
</script></html>\n'''


def build():
    data=load();svg=drawing.render(data);html=reader(svg,data)
    parsed=ET.fromstring(svg)
    for node in parsed.iter('{http://www.w3.org/2000/svg}a'):
        href=node.get('href','')
        if not href:raise ValueError('empty graphic link')
        if href.startswith(('https://','#')):continue
        path=(SVG.parent/unquote(href.split('#')[0])).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file():raise ValueError('broken graphic link: '+href)
    audit={
        'schema':'oph.horizontal_particle_atlas.v2',
        'input_sha256':data['input_sha256'],
        'source_scientific_experiments_replayed':False,
        'source_physical_identifications':[],
        'reference_data_changed':False,
        'hero_stage_ids':[s['id'] for s in data['story']['stages']],
        'family_targets':{r['id']:r['target_ids'] for r in data['families']['lanes']},
        'measurement_illustration':data['measurement'],
        'primitive_endpoint_illustration':data['endpoints'],
        'double_path_illustration':data['double_path'],
        'finite_theorem_roles':data['theorems'],
        'named_field_meanings':data['field_nature'],
        'finite_observer_output_connections':data['observer_outputs']['finite_output_connections'],
        'common_construction_specification_resolved':True,
        'finite_physical_preparation_including_boundary_charges_resolved':True,
        'six_physical_obligations_resolved':False,
        'physical_requirement_status':data['observer_outputs']['physical_requirement_status'],
        'layout':{'width':parsed.get('width'),'height':parsed.get('height')},
        'svg_sha256':hashlib.sha256(svg.encode()).hexdigest(),
        'html_sha256':hashlib.sha256(html.encode()).hexdigest(),
    }
    if any(sha(ROOT/p)!=h for p,h in data['input_sha256'].items()):
        raise ValueError('inputs changed during illustration build')
    return svg,html,canonical(audit)


def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args()
    outputs=build()
    for path,text in zip((SVG,HTML,AUDIT),outputs):
        if a.write:path.write_text(text)
        elif not path.is_file() or path.read_text()!=text:raise ValueError('stale horizontal atlas: '+str(path))
    print('PASS: horizontal SVG, standalone reader and bounded audit')


if __name__=='__main__':main()
