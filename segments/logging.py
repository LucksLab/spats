'''
 This module contains helper classes and functions related to logging and log message formatting.
'''

import logging
import os
import traceback

from logging import DEBUG, INFO, WARNING, ERROR
VERBOSE = 5

class _GlobalConfig(object):
    def __init__(self):
        self.level = INFO
_logConfig = _GlobalConfig()


def addLogger(objself, loggerContext = None, level = None):
    '''
     @param objself:  Some object to which to add the logger
     @param loggerContext:  Optionally specify the name/context of this logger
     @param level:  Optionally specify the level at which to filter log messages
    '''
    objself.loggerContext = str(loggerContext) if loggerContext is not None else objself.__class__.__name__
    objself.loggerLevel = level or _logConfig.level
    objself.logger = logging.getLogger(objself.loggerContext)
    objself.logger.setLevel(objself.loggerLevel)



class LoggingClass(object):
    '''
     Base class to use for clases where you want a logger.
     Use the following in the __init__ method of derived classes:
         LoggingClass.__init__(self, loggerContext, level)
    '''

    # TAI:  add loggerContext into the defaultFormat and use __name__ as the defaultContext for most classes
    #defaultFormat = "%(asctime)s [%(levelname)s] [%(name)s::%(funcName)s()@%(lineno)d]: %(message)s"
    # updated defaultFormat to handle `self.info`, etc patterns below
    #defaultFormat = "%(asctime)s [%(levelname)s] %(message)s"
    defaultFormat = "%(asctime)s [%(levelname)s %(name)s] %(message)s"
    threadFormat = "%(asctime)s [%(thread)|%(levelname)s %(name)s] %(message)s"
    noStampFormat = "[%(levelname)s %(name)s] %(message)s"

    def __init__(self, loggerContext = None, level = None, initLogging = False, loggingFormat = defaultFormat):
        '''
         @param loggerContext:  Optionally specify the name/context of this logger
         @param level:  Optionally specify the level at which to filter log messages
        '''
        if initLogging:
            initializeLogging(format = loggingFormat)
        if loggerContext == self:
            # for nicer display
            loggerContext = "{}@{}".format(self.__class__.__name__, id(self))
        addLogger(self, loggerContext, level)

    def __getstate__(self):
        ## Include this because loggers in older versions of python cannot be pickled.
        state = self.__dict__.copy()
        del state['logger']
        return state

    def _show(self, level, args):
        self.logger.log(level, ' '.join([str(x) for x in args]))

    def setLogLevel(self, level):
        # note that loggers are by default on a class basis, not an instance basis
        # you can control this with loggerContext
        self.logger.setLevel(level)
        self.loggerLevel = level

    def resetLogLevel(self):
        self.logger.setLogLevel(self.loggerLevel)

    def setLevelVerbose(self):
        self.setLogLevel(VERBOSE)

    def isLevelVerbose(self):
        return (self.loggerLevel == VERBOSE)

    def verbose(self, *args):
        self._show(VERBOSE, args)

    def verbosef(self, fmtStr, *args):
        self._show(VERBOSE, [ fmtStr.format(*args) ])

    def setLevelDebug(self):
        self.setLogLevel(DEBUG)

    def isLevelDebug(self):
        return (self.loggerLevel == DEBUG)

    def debug(self, *args):
        self._show(DEBUG, args)

    def debugf(self, fmtStr, *args):
        self._show(DEBUG, [ fmtStr.format(*args) ])

    def setLevelInfo(self):
        self.setLogLevel(INFO)

    def info(self, *args):
        self._show(INFO, args)

    def infof(self, fmtStr, *args):
        self._show(INFO, [ fmtStr.format(*args) ])

    def isLevelInfoOrHigher(self):
        return (self.loggerLevel <= INFO)

    def setLevelWarn(self):
        self.setLogLevel(WARNING)

    def warn(self, *args):
        self._show(WARNING, args)

    def warnf(self, fmtStr, *args):
        self._show(WARNING, [ fmtStr.format(*args) ])

    def setLevelError(self):
        self.setLogLevel(ERROR)

    def error(self, *args):
        self._show(ERROR, args)

    def errorf(self, fmtStr, *args):
        self._show(ERROR, [ fmtStr.format(*args) ])

    def exc(self, e, *args):
        self._show(ERROR, ["Exception: {}".format(e)] + list(args) + [ "\n", traceback.format_exc() ])

    def dumpStack(self, limit = 128, *args):
        entries = traceback.extract_stack(limit = limit)[:-1]
        lines = [ "  [{}:{}] @ {}:    {}".format('/'.join(e[0].split('/')[-2:]), e[1], e[2], e[3][:80]) for e in entries ]
        self._show(INFO, ["\n".join((list(args) or [ "Stack trace:" ]) + lines)])
    printStackTrace = dumpStack

    def fatal(self, *args):
        self._show(ERROR, args)
        self._show(ERROR, ["Fatal error, aborting w/o a stack trace."])
        os._exit(-1)

    def logToFile(self, logfile):
        '''
         Convenience method to add and set up logging to the specified logFile.
         Note this can also be achieved for the default logger with the initializeLoggging() function below.
         @param logfile:  Path to the logfile into which to write all log messages.
        '''
        loghandler = logging.FileHandler(logfile)
        loghandler.setFormatter(logging.Formatter(LoggingClass.defaultFormat))
        self.logger.addHandler(loghandler)


def initializeLogging(level = INFO, format = LoggingClass.defaultFormat, filename = None):
    '''
     Call this once for a program to initialize logging with the defaultFormat defined above.
     @param level:  Optionally specify the level at which to filter log messages; defaults to INFO
     @param format:  Optionally specify a format string for all log messages; defaults to LoggingClass.defaultFormat
     @param filename:  Optionally specify the path of a file into which to dump the log messages
    '''
    _logConfig.level = level
    if (filename is not None):
        logging.basicConfig(level=level, format=format, filename=filename)
    else:
        logging.basicConfig(level=level, format=format)
