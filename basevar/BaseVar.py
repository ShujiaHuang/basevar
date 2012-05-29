"""
This is the main program of BaseVar. It's the toppest of all
the BaseVar's tool sets.
Autor: Shujia Huang
Date: 2016-10-06 16:38:00
"""
import sys
import time


def basetype():
    from executor import RunBaseType

    bt = RunBaseType()
    bt.run()


def bqprob():
    from genotype import BQprob
    BQprob()


if __name__ == '__main__':

    runner = {'bqprob': bqprob, 'basetype': basetype}

    if len(sys.argv) == 1 or (sys.argv[1] not in runner):
        print >> sys.stderr, '[Usage] python [option] %s' % sys.argv[0]
        print >> sys.stderr, '\n\t'.join(['Option:'] + runner.keys())
        sys.exit(1)

    command = sys.argv[1]
    runner[command]()

    print >> sys.stderr, '** %s ALL DONE %s **' % (command, time.asctime())
    print >> sys.stderr, '>> For the flowers bloom in the desert <<'
