"""
=======
util.py
=======

Copyright (C) 2017-2022 Glysade, LLC

Helper functions for Ruse
"""

import io
import itertools
import json
import multiprocessing
import os
import psutil
from http import HTTPStatus

import requests
from typing import List, TypeVar, Dict

T = TypeVar('T')


class RuseServerException(Exception):
    """
    An exception class for ruse service errors
    """
    pass


def num_processors() -> int:
    """
    Find the number of processors on the local machine

    :return: number of processors
    """
    return multiprocessing.cpu_count()


def is_exe(fpath):
    """
    :param fpath:
    :return: true if thid file is executable
    """
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program: str) -> str:
    """
    Finds the full path of an executable on the PATH environment (checking that it is, in fact, executable)

    :param program: the executable to find
    :return: full path, or None if the executable is not on the PATH
    """

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def find_on_path(alternatives: List[str]) -> str:
    """
    Return the full path of the first of the input executables found on the PATH. Throws ValueError if none of the
    alternatives can be found

    :param alternatives:
    :return:
    """
    """ find one of alternative programs on the path"""
    for program in alternatives:
        exe = which(program)
        if exe:
            return exe

    raise ValueError("Unable to find any of {} on path!".format(alternatives))


def ruse_server() -> str:
    """
    :return: A (string) url to the ruse server
    """

    ruse_url = os.environ['RUSE_SERVER'] if 'RUSE_SERVER' in os.environ else 'http://localhost:9000'
    return ruse_url


def ruse_public_server() -> str:
    """
    :return: A (string) url to the ruse server
    """

    ruse_url = os.environ['RUSE_PUBLIC_SERVER'] if 'RUSE_PUBLIC_SERVER' in os.environ else 'http://localhost:9000'
    return ruse_url


def ruse_work_dir() -> str:
    """
    :return: full path to the svcrunner work directory
    """

    return os.environ['WorkDirPath'] if 'WorkDirPath' in os.environ \
        else os.path.join(os.path.dirname(__file__),
                          os.path.normpath(
                              '../../../ruse_docker/docker_files/ruse/master/files'))


def ruse_server_present() -> bool:
    """
    :return: True if the ruse server is present
    """
    try:
        response = requests.get(ruse_server() + '/tasks')
        return response.status_code in [requests.codes.ok, requests.codes.see_other]
    except requests.ConnectionError:
        return False


def run_external_process(exes: List[str], args: List[str]) -> int:
    """
    A helper method to run an external process.  See :class:`ruse.util.run_external_process.RunExternalProcess`

    :param exes: A list of possible alternative executable names
    :param args: program arguments
    :return: program exit code
    """
    from ruse.util.run_external_process import RunExternalProcess
    process = RunExternalProcess()
    return process.run_external_process(exes, args)


# don't import type for getter- creates circular references
#  full signature causes error in Sphnix
# def getter(json_obj: 'ruse.util.util.data_table.JsonType', key: str):
def getter(json_obj, key: str):
    """
    Given a key extract value from a input Json python dictionary

    :param json_obj: dictionary
    :param key: the key to extract
    :return: value of key, or None if key is not present
    """

    return json_obj[key] if key in json_obj else None


def dict_to_list(in_dict: Dict[T, T]) -> List[T]:
    """
    Converts a dictionary to a list

    :param in_dict: input dictionary
    :return: output list
    """

    return list(itertools.chain.from_iterable([[k, in_dict[k]] for k in in_dict]))


# full signature causes error in Sphinx
# def input_field_value(json_in: 'ruse.util.util.data_table.JsonType', field: str, optional:bool = False) -> 'ruse.util.util.data_table.DataTypes':
def input_field_value(json_in, field: str,
                      optional: bool = False):
    """
    Retrieve the data attribute for an input field from submitted json

    :param json_in: python data structure from submitted Json
    :param field: input field name
    :param optional: If true return None if the input field is not present
    :return: input field data
    """

    if optional and field not in json_in['inputFields']:
        return None
    return json_in['inputFields'][field]['data']


# full signature causes error in Sphinx
# def input_field_type(json_in: 'ruse.util.util.data_table.JsonType', field: str) -> str:
def input_field_type(json_in, field: str) -> str:
    """
    Retrieve the contentType attribute for an input field from submitted json

    :param json_in: python data structure from submitted Json
    :param field: input field name
    :return: input field content type
    """

    return json_in['inputFields'][field]['contentType']


# full signature causes error in Sphnix
# def has_input_field(json_in: 'ruse.util.util.data_table.JsonType', field: str) -> bool:
def has_input_field(json_in, field: str) -> bool:
    """
    :param json_in:  python data structure from submitted Json
    :param field: input field name
    :return: True if the input field is present in the Json
    """
    return field in json_in['inputFields']


def get_task_id() -> int:
    """
    Gets the id of the current task from the current working directory
    :return: Current task id
    """

    return int(os.path.split(os.getcwd())[-1])


# full signature causes error in Sphnix
# def get_task_result(task_id: int) -> 'ruse.util.util.data_table.JsonType':
def get_task_result(task_id: int):
    """
    Retrieves result.json from the ruse server for a given task id

    :param task_id: task id
    :return: output json for task
    """
    server = ruse_server()
    url = server + "/task/%d/result.json" % task_id
    response = requests.get(url)
    if response.status_code == HTTPStatus.OK:
        return json.loads(response.text)

    return None


def get_memory_in_gb():
    mem = psutil.virtual_memory()
    return mem.total / (1024.0 * 1024.0 * 1024.0)
