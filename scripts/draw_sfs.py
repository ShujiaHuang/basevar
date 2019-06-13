"""
Draw the site frequence spectrum.

Copyright (c) Shujia Huang
Date: 2015-09-20
"""
import sys
import os

import numpy as np
import matplotlib.pyplot as plt

STEP_BIN = round(1.0 * 5/1008, 4)

def main(argv):

    # Load data first
    h, data = loadData(argv[0])
    summary(data)

    known = data[data[:,0]==1.0][:,1:2] # [AF]
    novel = data[data[:,0]==0.0][:,1:2] # [AF]

    # Distribution
    known_af = calcuDist(known)
    novel_af = calcuDist(novel)

    af_bin_order, af = fixData(known_af, novel_af)

    #draw(argv[1], [['MLEAF', mleaf_bin_order, mleaf],
    #               ['AF', af_bin_order, af]])
    draw2(argv[1], ['AF', af_bin_order, af])

def summary(data):
    """
    Output the data summary
    """
    known = data[data[:,0] == 1.0][:,1:2] # [AF]
    novel = data[data[:,0] == 0.0][:,1:2] # [AF]

    known_sum, known_summary = _summary(known[:,0]) # Just AF
    novel_sum, novel_summary = _summary(novel[:,0]) # Just AF

    key_list_inorder = sorted(list(set(known_summary.keys() +
                                   novel_summary.keys())))
    print '# Just for AF:\n\n#Frequence\tKnow(%)\tNovel(%)'
    for k in key_list_inorder:
        if k not in known_summary: known_summary[k] = 0
        if k not in novel_summary: novel_summary[k] = 0
        known_r = round(known_summary[k] * 100/ float(known_sum + novel_sum), 2)
        novel_r = round(novel_summary[k] * 100/ float(known_sum + novel_sum), 2)
        print (k + ': ' + str(known_summary[k]) + '(' + str(known_r) + '%)' +
               '\t' + str(novel_summary[k]) + '(' + str(novel_r) + '%)')

def _summary(data):
    """
    """
    stat = {}
    sumup  = 0
    for d in data:
        if d > 0.05:
            stat['>0.05'] = stat.get('>0.05', 0) + 1
        #elif d > 0.01:
        #    stat['0.01-0.05'] = stat.get('0.01-0.05', 0) + 1
        elif d > 0.005:
            stat['0.005-0.05'] = stat.get('0.005-0.05', 0) + 1
        else:
            stat['0.005'] = stat.get('0.005', 0) + 1
        sumup += 1

    return sumup, stat

def draw(figfile, data):
    """
    """
    plt.style.use('ggplot')

    a_sz    = len(data)
    col_num = 2
    row_num = a_sz / col_num if a_sz % col_num == 0 else (a_sz / col_num) + 1
    # Init the figure
    fig, axes = plt.subplots(ncols = col_num,
                             nrows = row_num,
                             sharey  = False,
                             figsize = (12 * col_num, 8 * row_num),
                             tight_layout = True)
    ax = axes.ravel() # Get all the sub figures and assign them to an array
    # The Y axis' name
    for i, d in enumerate(data):
        sub_plot(ax[i], d[0], d[1], d[2])

    # Svae figure
    fig.savefig(figfile)

def draw2(figfile, data):
    """
    ``data``: just one row
    """
    plt.style.use('ggplot')

    a_sz    = len(data)
    col_num = 1
    row_num = 1
    # Init the figure
    fig, axes = plt.subplots(ncols = col_num,
                             nrows = row_num,
                             sharey = False,
                             figsize = (12 * col_num, 8 * row_num),
                             tight_layout = True)
    ax = axes
    sub_plot(ax, data[0], data[1], data[2])

    # Svae figure
    fig.savefig(figfile)

def sub_plot(ax, title_name, x_axis, data):
    """
    Argv:
        `data`: [known, novel]
    """
    color1 = plt.rcParams['axes.color_cycle'][1] # color for known
    color2 = plt.rcParams['axes.color_cycle'][0] # color for novel

    # Set figure
    ax.tick_params(top = False, right = False, labelsize = 16)
    ind = np.arange(len(data)) # the x locations for the groups
    width = 1.0

    # Drawing
    p1 = ax.bar(ind, data[:,0], width, color = 'g')
    #p2 = ax.bar(ind, data[:,1], width, color = color2, bottom = data[:,0])

    n = len(ind)
    ax.set_xticks(ind[0:n:2] + width / 2.)
    ax.set_xticklabels(x_axis[0:n:2],
                       rotation = 'vertical',
                       fontsize = 16)
#ax.legend((p1[0], p2[0]), ('Known', 'Novel'), fontsize = 16, loc = 'upper right')
    ax.set_xlabel('MAF', fontsize = 18, color = 'k')
    ax.set_ylabel('Proportion of variants', fontsize = 18, color = 'k')
    ax.grid(b = True)

def fixData(known, novel):
    key_list_inorder = sorted(list(set(known.keys() + novel.keys())))
    data = []
    sum = 0
    for k in key_list_inorder:
        if k not in known: known[k] = 0
        if k not in novel: novel[k] = 0
        data.append([known[k], novel[k]])
        sum += (known[k] + novel[k])

    return key_list_inorder, np.array(data) / float(sum)

def calcuDist(data):
    """
    """
    # Set bin
    data[:,0] = [calcubin(i) for i in data[:,0]] # AF

    af_dist_dict = {}
    for i in data[:,0]:
        af_dist_dict[i] = af_dist_dict.get(i, 0) + 1

    return af_dist_dict

def calcubin(datum, amplifier = 10000):
    """
    """
    # Amplifier the data make them easy to calculate.
    step = int(STEP_BIN * amplifier)
    datum = int(datum * amplifier)
    bin_size = datum // step + 1 if datum % step > 0 else datum // step

    return bin_size * step / float(amplifier) #

def loadData(file):
    """
    Load the data from input file
    """
    data   = []
    header = None
    I = open(file)
    while 1:

        lines = I.readlines(100000)
        if not lines: break

        for line in lines:
            #CHROM  POS Known(1)/Novel(0)   AF
            col = line.strip('\n').split()
            if col[0][0] == '#':
                header = col[2:]
            else:
                af = float(col[3])
                if af > 0.0 and af < 1.0: # Ingnore homo reference
                    data.append([float(col[2]), float(col[3])
                                 if float(col[3]) < 0.5 else 1.0 - float(col[3])])
    I.close()

    return header, np.array(data)

if __name__ == '__main__':

    main(sys.argv[1:])