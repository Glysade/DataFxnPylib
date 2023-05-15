"""
=======
log.py
=======

Copyright (C) 2017-2022 Glysade, LLC

Ruse logging module
"""

import sys
import datetime

class Level:
    """
    Global logging levels
    """
    All         = 0
    Debug       = 1
    Info        = 2
    Warning     = 3
    Error       = 4
    Critical    = 5

logLevel = Level.All


def debug(*args):
    """
    Log a debug message

    :param args: messages
    """

    if Level.Debug < logLevel:
        return
    _to_log("%-9s" % "Debug", args)


def info(*args):
    """
    Log an informational message

    :param args: messages
    """

    if Level.Info < logLevel:
        return
    _to_log("%-9s" % "Info", args)


def warning(*args):
    """
    Log an warning message

    :param args: messages
    """

    if Level.Warning < logLevel:
        return
    _to_log("%-9s" % "Warning", args)


def error(*args):
    """
    Log an error message

    :param args: messages
    """

    if Level.Error < logLevel:
        return
    _to_log("%-9s" % "Error", args)


def critical(*args):
    """
    Log an critical message

    :param args: messages
    """

    _to_log("%-9s" % "Critical", args)


def set_level(level):
    """
    Set the global log level

    :param level:
    """
    global logLevel
    logLevel= level
    print('level = ',logLevel)


def _to_log(level, args):
    """
    Writes log message

    :param level: log level as str
    :param args: messages
    """
    msg = _to_string(args)
    with open('task.log','a') as fh:
        logline = level + ": " + datetime.datetime.now().isoformat().split('.')[0] + " " + msg + "\n"
        fh.write(logline)
    print(msg,file=sys.stderr)
    sys.stderr.flush()


def _to_string(args):
    """
    Convert a list of messages to a single string

    :param args: messages
    :return: joined string
    """
    msg = ""
    try:
        msg = " ".join(str(x) for x in args)
    except:
        pass
    return msg
