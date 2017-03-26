"""
This is the main program of BaseVar. It's the toppest layer of
BaseVar's tool sets.

Autor: Shujia Huang
Date: 2016-10-06 16:38:00
"""
import sys
import time


def basetype():
    from caller.executor import Runner

    bt = Runner()
    bt.basetype()


def coverage():
    from caller.executor import Runner
    cvg = Runner()
    cvg.coverage()


if __name__ == '__main__':

    runner = {'basetype': basetype,
              'coverage': coverage}

    if len(sys.argv) == 1 or (sys.argv[1] not in runner):
        print >> sys.stderr, '[Usage] python [option] %s' % sys.argv[0]
        print >> sys.stderr, '\n\t'.join(['Option:'] + runner.keys())
        sys.exit(1)

    command = sys.argv[1]
    runner[command]()

    print >> sys.stderr, '** %s ALL DONE %s **' % (command, time.asctime())
    print >> sys.stderr, '>> For the flowers bloom in the desert <<'
