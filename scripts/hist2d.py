"""
Author: Shujia Huang
Date: 2017-06-13 09:15:04
"""
import sys
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("Agg")

from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt


def scale_format(value):

    if value < 0.0001:
        value = '%.5f' % value

    elif value < 0.001:
        value = '%.4f' % value

    elif value < 0.01:
        value = '%.3f' % value

    elif value < 0.1:
        value = '%.2f' % value

    elif value < 0.5:
        value = '%.1f' % value

    else:
        value = '%.0f' % value

    return value


def draw_hist2d(argv):

    MIN_MAF = 0.003
    ax = plt.gca()

    data = pd.read_table(argv[0])
    labels = argv[1].strip().split(':')
    out_fig_file = argv[2]

    data = pd.DataFrame(data.values, columns=labels, dtype=float)

    data = data[data[labels[0]] + data[labels[1]] > 0]
    data = data[data[labels[1]] < 1]

    x, y = data[labels[0]], data[labels[1]]

    min_value = max([min(x), min(y)])
    if min_value == 0:
        min_value = min([min(x[x > 0]), min(y[y > 0])])

    x[x == 0] = min_value
    y[y == 0] = min_value
    im = ax.hist2d(x=np.log10(x), y=np.log10(y), bins=100,
                   norm=LogNorm(), cmap=plt.cm.hsv)

    min_a = min([min(np.log10(x)), min(np.log10(y))])
    max_a = max([max(np.log10(x)), max(np.log10(y))])
    ax.plot([min_a, max_a], [min_a, max_a], 'k--', linewidth=1)

    ax.set_xlim([np.log10(MIN_MAF), max_a])
    ax.set_ylim([np.log10(MIN_MAF), max_a])

    locs, _ = plt.xticks()
    # ax.set_xticklabels(['${10}^{%.1f}$'%i for i in locs])
    ax.set_xticklabels([scale_format(10 ** i) for i in locs])

    locs, _ = plt.yticks()
    # ax.set_yticklabels(['${10}^{%.1f}$'%i for i in locs])
    ax.set_yticklabels([scale_format(10 ** i) for i in locs])

    ax.set_xlabel(labels[0], fontsize=14)
    ax.set_ylabel(labels[1], fontsize=14)

    plt.colorbar(im[-1])

    plt.tight_layout()
    plt.savefig(out_fig_file)

    #plt.show()


if __name__ == '__main__':

    # usage: python hist2d.py basvarc_plot.txt.gz "Real:Basevar" test.pdf
    draw_hist2d(sys.argv[1:])
