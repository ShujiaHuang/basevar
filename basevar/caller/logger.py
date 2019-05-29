"""decorator for logger
"""

import time
from functools import wraps
import inspect


def logger(func):
    @wraps(func)
    def wrapper(*args, **kwargs):

        command = func.__name__

        start_time = time.time()

        print("\n** %s Start at %s **\n" % (command, time.asctime()))
        func(*args, **kwargs)
        elapsed_time = time.time() - start_time

        print ("[INFO] function = %s\n" % command)
        print ("[INFO] called_from_line: %d\n" % inspect.currentframe().f_back.f_back.f_lineno)
        print ('** %s done at %s, %d seconds elapsed **\n' % (command, time.asctime(), elapsed_time))

    return wrapper

