"""
This is a caller process
"""
import sys
import time
import multiprocessing


class CallerProcess(multiprocessing.Process):
    """
    simple class to represent a single process, which is run as part of
    a multi-process job.

    It's a shield for ``func``
    """

    def __init__(self, func, *args, **kwargs):
        """Multiprocess Constructor.
        """
        multiprocessing.Process.__init__(self)
        self.single_process = func(*args, **kwargs)

    def run(self):
        """ Run the BaseVar process"""
        self.single_process.run()


def process_runner(processes):
    """run and monitor the process"""

    for p in processes:
        p.start()

    # listen for signal while any process is alive
    while True in [p.is_alive() for p in processes]:
        try:
            time.sleep(1)

        except KeyboardInterrupt:
            sys.stderr.write('KeyboardInterrupt detected, terminating '
                             'all processes...\n')
            for p in processes:
                p.terminate()

            sys.exit(1)

    # Make sure all process are finished
    for p in processes:
        p.join()

    return
