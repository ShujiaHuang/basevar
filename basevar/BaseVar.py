"""
This is the main program of BaseVar. It's the toppest layer of
BaseVar's tool sets.

Autor: Shujia Huang
Date: 2016-10-06 16:38:00
"""
import sys
import time


def basetype():
    from caller.executor import BaseTypeRunner

    bt = BaseTypeRunner()
    bt.run()


def vqsr():
    from caller.executor import VQSRRuner
    vq = VQSRRuner()
    vq.run()


def merge():
    from caller.executor import MergeRunner

    mg = MergeRunner()
    mg.run()


if __name__ == '__main__':

    runner = {'basetype': basetype,
              'merge': merge,
              'VQSR': vqsr}

    if len(sys.argv) == 1 or (sys.argv[1] not in runner):
        sys.stderr.write('[Usage] python %s [option]\n\n' % sys.argv[0])
        sys.stderr.write('\n\t'.join(['Option:'] + runner.keys()))
        sys.exit(1)

    command = sys.argv[1]
    runner[command]()

    sys.stderr.write('** %s ALL DONE %s **' % (command, time.asctime()))
    sys.stderr.write('>> For the flowers bloom in the desert <<')