"""
Author : Shujia Huang
Date : 2016-04-06
"""
import sys
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.optimize import curve_fit
from scipy.misc import factorial

from scipy import stats

import geneview as gv
##############################################
def jointplot(argv):

    data = pd.read_table(argv[0], header=None)
    data = data * 100
    labels = argv[1].strip().split(':')
    out_fig_file = argv[2]

    data = pd.DataFrame(data.values, columns=labels)
    data = data[data[labels[0]] + data[labels[1]]>1]

    data[labels[1]][data[labels[1]]>50] = 100.0 - data[labels[1]][data[labels[1]]>50]
    data[labels[0]][data[labels[0]]>50] = 100.0 - data[labels[0]][data[labels[0]]>50]
    g = gv.jointplot(x=data[labels[1]], y=data[labels[0]], kind="hex",
                     joint_kws=dict(gridsize=30), color='#FF0000',
                     ratio=3, dropna=True,
                     marginal_kws=dict(bins=20))
    plt.savefig(out_fig_file)
    plt.show()


def jointplot_old(argv):

    data = pd.read_table(argv[0], header=None)
    data = pd.DataFrame(data.values, columns=['1KG', 'NIFTY'])
    data = data[data['1KG'] + data['NIFTY'] > 1]
    data['NIFTY'][data['NIFTY']>50] = 100.0 - data['NIFTY'][data['NIFTY']>50]
    data['1KG'][data['1KG']>50] = 100.0 - data['1KG'][data['1KG']>50]
    g = gv.jointplot(x=data['NIFTY'], y=data['1KG'], kind="hex",
                                          joint_kws=dict(gridsize=30))#, bins='log'))
    plt.savefig('jointplot.pdf')
    plt.show()


def hist2d(argv):
    data = pd.read_table(argv[0], header=None)
    data = pd.DataFrame(data.values, columns=['1KG', 'NIFTY'])
    data = data[data['1KG'] + data['NIFTY'] > 0]
    data['NIFTY'][data['NIFTY']>50] = 100.0 - data['NIFTY'][data['NIFTY']>50]
    data['1KG'][data['1KG']>50] = 100.0 - data['1KG'][data['1KG']>50]
    gv.hist2d(x=data['NIFTY'], y=data['1KG'], bins=100)
    plt.show()


def depth_bar(argv):
    # load data
    bin_size = 100
    data = {}
    with open(argv[0]) as I:
        # chr12 60000   6
        for line in I:
            col = line.strip().split()
            depth = int(col[2])
            if depth < 100 or depth > 40000:
                continue
            depth = depth // bin_size if depth % bin_size == 0 else depth // bin_size + 1
            depth *= bin_size
            data[depth] = data.get(depth, 0) + 1

    # covert to array data
    data = np.array([[k, d] for k, d in sorted(data.items(), key = lambda k:k[0])])

    # draw the main figure
    ax = plt.gca()

    ind = data[:,0]
    width = 100
    ax.bar(ind, data[:,1], width, edgecolor=None)
    plt.plot(ind, data[:,1], 'r--', linewidth=0.8)
    print y
    plt.savefig('bar.pdf')
    plt.show()


# poisson function, parameter lamb is the fit parameter
def poisson(k, lamb):
    return (lamb**(k)/factorial(k)) * np.exp(-lamb)


def depth_hist(argv):
    data = []
    num = 0
    with open(argv[0]) as I:
        for line in I:
            num += 1
            if num % 1000000 == 0:
                print '[INFO]Loading lines', num
            col = line.strip().split()
            depth = int(col[2])
            if depth > 25000:
                continue
            data.append(depth * 10)

    num_bins = 50
    data = np.array(data)
    ax = plt.gca()
    n, bins, patches = ax.hist(data, num_bins, normed=1, alpha=0.5)

    """
    # calculate binmiddles
    bin_middles = 0.5*(bins[1:] + bins[:-1])
    # fit with curve_fit
    parameters, cov_matrix = curve_fit(poisson, bin_middles, n)
    x_plot = np.linspace(0, 20, data.max())
    plt.plot(x_plot, poisson(x_plot, *parameters), 'r--')
    #mu = data.mean()
    #plt.plot(bins, np.array([poisson(x, mu) for x in bins]), 'r--')
    """
    mu = data.mean()
    sigma = data.std()
    y = mlab.normpdf(bins, mu - 300 * 10, sigma)
    ax.plot(bins, y, 'r--')
    ax.set_xlim(0, 25000 * 10)
    plt.savefig('hist.pdf')
    #plt.show()


if __name__ == '__main__':
    jointplot(sys.argv[1:])
