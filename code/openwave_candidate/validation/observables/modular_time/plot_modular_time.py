"""Plot source-bound modular correlations; no native evolution rerun."""
from pathlib import Path
import csv
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

here=Path(__file__).resolve().parent
receipt=json.loads((here/'MODULAR_TIME.json').read_text())
rows=list(csv.DictReader((here/'correlations.csv').open()))
s=[float(r['s']) for r in rows]
fig,axes=plt.subplots(2,1,figsize=(8,6),sharex=True,layout='constrained')
axes[0].plot(s,[float(r['C_real']) for r in rows],label='Re C(s)',lw=1.8)
axes[0].plot(s,[float(r['C_imag']) for r in rows],label='Im C(s)',lw=1.8)
axes[0].set(ylabel='Observer correlation',title='Modular flow of a conditioned native observer state')
axes[0].legend(loc='upper right',ncol=2);axes[0].grid(alpha=.2)
axes[1].plot(s,[float(r['observable_01_real']) for r in rows],label='Re A(s)[0,1]',lw=1.5)
axes[1].plot(s,[float(r['observable_01_imag']) for r in rows],label='Im A(s)[0,1]',lw=1.5)
axes[1].set(xlabel='Dimensionless modular parameter s',ylabel='Observable element')
axes[1].legend(loc='upper right',ncol=2);axes[1].grid(alpha=.2)
fig.supxlabel('Native preparation duration = 0.1; modular parameter has no assigned physical-time conversion.',fontsize=9)
fig.savefig(here/'modular_correlations.png',dpi=180)
fig.savefig(here/'modular_correlations.svg')
files=[Path(__file__).resolve(),here/'MODULAR_TIME.json',here/'correlations.csv']
(here/'PLOT_SOURCES.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2)+'\n')
