import sys

from logbook import Logger, StreamHandler

LOG_NAME = "BaseVar"

log_handler = StreamHandler(sys.stderr)
log_handler.push_application()
logger = Logger(LOG_NAME)
