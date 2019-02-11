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

    if value < 0.01:
        value = '%.3f' % value

    elif value < 0.1:
        value = '%.2f' % value

    elif value < 0.5:
        value = '%.1f' % value

    else:
        value = '%.0f' % value

    return value


def draw_hist2d(argv):

    ax = plt.gca()

    data = pd.read_table(argv[0])
    labels = argv[1].strip().split(':')
    out_fig_file = argv[2]

    data = pd.DataFrame(data.values, columns=labels, dtype=float)

    data = data[data[labels[0]] + data[labels[1]] > 0]
    data = data[data[labels[1]] < 1]

    x, y = data[labels[0]], data[labels[1]]

    # im = ax.hist2d(x=x, y=y, bins=100, norm=LogNorm(), cmap=plt.cm.hsv)
    # ax.plot([min(x), max(x)], [min(y), max(y)], 'k--', linewidth=1)
    x[x == 0] = 0.0001
    y[y == 0] = 0.0001
    im = ax.hist2d(x=np.log10(x), y=np.log10(y), bins=100,
                   norm=LogNorm(), cmap=plt.cm.hsv)

    min_a = min([min(np.log10(x)), min(np.log10(y))])
    max_a = max([max(np.log10(x)), max(np.log10(y))])
    ax.plot([min_a, max_a], [min_a, max_a], 'k--', linewidth=1)

    locs, labs = plt.xticks()
    # ax.set_xticklabels(['${10}^{%.1f}$'%i for i in locs])
    ax.set_xticklabels([scale_format(10 ** i) for i in locs])

    locs, labs = plt.yticks()
    # ax.set_yticklabels(['${10}^{%.1f}$'%i for i in locs])
    ax.set_yticklabels([scale_format(10 ** i) for i in locs])

    ax.set_xlabel(labels[0], fontsize=14)
    ax.set_ylabel(labels[1], fontsize=14)

    plt.colorbar(im[-1])

    plt.tight_layout()
    plt.savefig(out_fig_file)

    # plt.show()


if __name__ == '__main__':

    draw_hist2d(sys.argv[1:])


