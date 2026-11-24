"""Compare photometric enclosure in AGC 114905 with the eligible outer-orbit sample."""
import json

import matplotlib.pyplot as plt
import numpy as np

import disk_domain

HERE = disk_domain.fit.HERE


def draw():
    bank = json.loads((HERE / 'disk_domain.json').read_text())
    rings = json.loads((HERE.parent / 'udg_transfer/rings.json').read_text())['panels']['figure8']['rings']
    ink, teal, copper = '#27354a', '#176c72', '#a04c24'
    plt.rcParams.update({'font.size': 9, 'pdf.fonttype': 42})
    fig, (gas, orbits) = plt.subplots(1, 2, figsize=(8.7, 3.2), layout='constrained')
    radius = np.linspace(0, 19, 150)
    gas.plot(radius, [100*disk_domain.gas_fraction_inside(r) for r in radius], color=teal)
    measured = [r['radius_kpc'] for r in rings]
    gas.scatter(measured, [100*disk_domain.gas_fraction_inside(r) for r in measured],
                s=24, color=ink, zorder=3, label='Radii of measured speeds')
    gas.axhline(90, color='0.55', linestyle=':', linewidth=1)
    gas.annotate('90% of gas at 14.30 kpc', (14.30,90), xytext=(5.0,98), fontsize=8,
                 arrowprops={'arrowstyle':'-', 'color':'0.5'})
    gas.annotate('66.9%', (measured[-1],66.9), xytext=(11.0,59), color=ink,
                 arrowprops={'arrowstyle':'-', 'color':'0.5'})
    gas.set(xlabel='Cylindrical radius (kpc)', ylabel='Enclosed gas mass (%)',
            title='AGC 114905 gas profile', xlim=(0,19), ylim=(0,104))
    gas.legend(loc='lower right', frameon=False, fontsize=8)
    names = ['NGC3838', 'NGC2685', 'NGC6798']
    for i,name in enumerate(names):
        obj = bank['exterior_orbits'][name]
        orbits.errorbar(obj['observed_kms'], i+0.10, xerr=obj['observed_error_kms'],
                        fmt='o', capsize=3, color=ink, label='Measured ± quoted error' if i==0 else None)
        lo,hi = obj['phase_conditioned_speed_range_kms']
        orbits.plot([lo,hi],[i-0.10,i-0.10],color=copper,lw=3,solid_capstyle='butt',
                    marker='|',markersize=6,label='Predicted cut-phase range' if i==0 else None)
    limits=[bank['exterior_orbits'][k]['enclosure_lower_bound'] for k in names]
    orbits.set(yticks=range(3),yticklabels=[f'{k[:3]} {k[3:]}\n≥{int(v*100)}% inside' for k,v in zip(names,limits)],
               ylim=(2.85,-0.55), xlim=(124,205), xlabel='Outer circular speed (km/s)',
               title='Deep orbits with ≥80% baryon enclosure')
    orbits.grid(axis='x', alpha=0.2)
    orbits.legend(loc='lower right', frameon=False, fontsize=8)
    for axes in (gas, orbits):
        axes.spines[['right','top']].set_visible(False)
    fig.savefig(HERE / 'comparison.pdf', metadata={'CreationDate':None, 'ModDate':None})
    plt.close(fig)


if __name__ == '__main__':
    draw()
