"""
This is the main program of BaseVar. It's the toppest layer of
BaseVar's tool sets.

Autor: Shujia Huang
Date: 2016-10-06 16:38:00
"""
import sys
import time


def basetype():
    from caller.executor import BaseTypeBamRunner

    bt = BaseTypeBamRunner()
    bt.run()

    return


def basetype_mpileup():
    from caller.executor import BaseTypeRunner

    bt = BaseTypeRunner()
    bt.run()

    return


def vqsr():
    from caller.executor import VQSRRuner
    vq = VQSRRuner()
    vq.run()

    return


def nearby_indel():

    from caller.executor import NearbyIndelRunner
    nbi = NearbyIndelRunner()
    nbi.run()

    return


def merge():
    from caller.executor import MergeRunner

    mg = MergeRunner()
    mg.run()

    return

def coverage():
    from caller.executor import CoverageRunner
    cvg = CoverageRunner()
    cvg.run()


if __name__ == '__main__':

    runner = {'basetype': basetype,
              'merge': merge,
              'coverage': coverage,
              'nbi': nearby_indel,
              'VQSR': vqsr
              }

    if len(sys.argv) == 1 or (sys.argv[1] not in runner):
        sys.stderr.write('[Usage] python %s [option]\n\n' % sys.argv[0])
        sys.stderr.write('\n\t'.join(['Option:'] + runner.keys()) + '\n\n')
        sys.exit(1)

    command = sys.argv[1]
    runner[command]()

    sys.stderr.write('** %s ALL DONE %s **\n' % (command, time.asctime()))
    sys.stderr.write('>> For the flowers bloom in the desert <<\n')