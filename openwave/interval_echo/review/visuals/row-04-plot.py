"""Render the row 4 pair-persistence figure from its archived detector counts."""
import argparse
import json
from pathlib import Path
from zipfile import ZipFile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
TEAL = '#087666'
RUST = '#a34323'
INK = '#182b3b'
MUTED = '#4d6372'


def receipt(location, name):
    """Read a named row 4 receipt from the archive or its extracted directory."""
    if location.is_dir():
        return json.loads((location / name).read_bytes())
    with ZipFile(location) as archive:
        return json.loads(archive.read(name))


def canvas(title, subtitle, size=(10.5, 5.15)):
    fig = plt.figure(figsize=size, facecolor='white')
    fig.text(.065, .925, title, fontsize=20, weight='bold', color=INK)
    fig.text(.065, .855, subtitle, fontsize=11, color=MUTED)
    return fig


def style(ax):
    ax.set_axisbelow(True)
    ax.grid(axis='y', color='#dbe3e7', linewidth=.7)
    for side in ['top', 'right']:
        ax.spines[side].set_visible(False)
    for side in ['left', 'bottom']:
        ax.spines[side].set_color('#aabac2')
    ax.tick_params(labelsize=10, colors=INK, length=3)


def save(fig, out, filename):
    fig.savefig(out / filename, dpi=165, facecolor='white', metadata={'Software': 'Matplotlib'})
    plt.close(fig)


def persistence(src, out):
    data = receipt(src, 'row_04_receipt.json')
    bound = receipt(src, 'row_04_window.json')
    assert data['addresses'][-1] == 'quiet'
    fig = canvas('A local receiver still finds the prepared pair',
                 'Row 4  ·  Cell 13 of a 20-cell sphere  ·  4,096 freshly prepared releases per point')
    ax = fig.add_axes([.085, .24, .85, .52])
    style(ax)
    exposures = [r['exposure'] for r in data['arms'][0]['readings']]
    pos = np.arange(len(exposures))
    for arm in data['arms']:
        native = arm['strength'] == data['native_kappa']
        offset, color, marker = (-.105, TEAL, 'o') if native else (.105, RUST, 's')
        name = 'Native binding' if native else 'Binding removed'
        reads = arm['readings']
        assert [r['exposure'] for r in reads] == exposures
        hits = np.array([sum(r['counts'][:-1]) for r in reads])
        n = data['shots_per_setting']
        assert all(sum(r['counts']) == n for r in reads)
        ax.scatter(pos + offset, 100 * hits/n, s=68, color=color, marker=marker,
                   label=f'{name} (κ = {arm["strength"]:g})', edgecolors='white', zorder=3)
        for x, h in zip(pos + offset, hits):
            ax.annotate(f'{h:,}', (x, 100*h/n), xytext=(0, -17 if native else 10),
                        textcoords='offset points', ha='center', color=color, weight='bold')
    ax.set(xlim=(-.47, 2.47), ylim=(0, 104), yticks=[0, 25, 50, 75, 100],
           ylabel='Releases with a local pair detection (%)', xlabel='Free exposure at readout')
    ax.set_xticks(pos, labels=exposures)
    ax.legend(loc='center right', bbox_to_anchor=(1, .54), frameon=True, facecolor='white', framealpha=.96,
              edgecolor='#e0e6e8', fontsize=9)
    fig.text(.065, .052, f'The finite model bounds the native local pair probability at or above '
             f'{bound["all_exposures_lower_bound_percent"]}% throughout exposures 0–{bound["last_exposure"]}.',
             fontsize=10, color=MUTED)
    save(fig, out, 'row-04-persistence.png')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--row04', type=Path, default=HERE.parent/'row-04/row-04-evidence.zip')
    parser.add_argument('--output', type=Path, default=HERE)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=True, parents=True)
    plt.rcParams.update({'font.family':'DejaVu Sans', 'text.color': INK, 'axes.labelcolor':INK, 'axes.labelsize':10})
    persistence(args.row04, args.output)


if __name__ == '__main__':
    main()
