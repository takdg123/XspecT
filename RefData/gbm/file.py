import os.path
import copy
import re

import datetime

from gbm.time import Met
from .detectors import Detector


class GbmFile:
    REGEX_PATTERN = r'^glg_(?P<data_type>.+)_(?P<detector>[bnl][0-9abl]|all)_(?P<trigger>(?:bn)?)(?P<uid>(?:\d{9}|\d{6}' \
                    r'_\d\dz|\d{6}))(?P<meta>(?:_.+)?)_v(?P<version>\d\d)\.(?P<extension>.+)$'

    def __init__(self):
        self.directory = ''
        self.trigger = False
        self.data_type = None
        self._detector = None
        self.uid = None
        self.meta = None
        self.version = 0
        self.extension = 'fit'

    def _init_by_dict(self, values):
        for key, val in values.items():
            # handle properties differently
            try:
                p = getattr(self, key)
                if isinstance(p, property):
                    p.__set__(self, val)
                else:
                    self.__setattr__(key, val)
            except AttributeError:
                raise ValueError("{} is not a valid attribute.".format(key))

    @property
    def detector(self):
        if not self._detector:
            return 'all'
        return self._detector

    @detector.setter
    def detector(self, value):
        if value == 'all':
            self._detector = None
        elif isinstance(value, Detector):
            self._detector = value
        else:
            if isinstance(value, str):
                d = Detector.from_str(value)
                self._detector = d if d else value
            elif isinstance(value, int):
                d = Detector.from_num(value)
                if d:
                    self._detector = d
                else:
                    raise ValueError("Invalid detector value")

    def version_str(self):
        if isinstance(self.version, int):
            v = "{:02d}".format(self.version)
        else:
            v = self.version
        return v

    def basename(self):
        if self.trigger:
            u = 'bn' + self.uid
        else:
            u = self.uid

        if self.meta:
            return str.format("glg_{}_{}_{}{}_v{}.{}",
                              self.data_type, self.detector, u, self.meta, self.version_str(), self.extension)
        return str.format("glg_{}_{}_{}_v{}.{}",
                          self.data_type, self.detector, u, self.version_str(), self.extension)

    def path(self):
        return os.path.join(self.directory, self.basename())

    def __str__(self):
        return self.path()

    def __repr__(self):
        return self.basename()

    @classmethod
    def create(cls, **kwargs):
        obj = cls()
        obj._init_by_dict(kwargs)
        return obj

    @classmethod
    def from_path(cls, path):
        m = re.match(cls.REGEX_PATTERN, os.path.basename(path), re.I | re.S)

        result = None
        if m:
            result = cls.create(**m.groupdict())
            result.directory = os.path.dirname(path)

        return result

    def detector_list(self):
        result = []
        for d in Detector:
            x = copy.copy(self)
            x.detector = d
            result.append(x)
        return result

    @classmethod
    def list_from_paths(cls, path_list, unknown=None):
        result = []
        for p in path_list:
            f = GbmFile.from_path(p)
            if f:
                result.append(f)
            else:
                if unknown is not None:
                    unknown.append(p)
                else:
                    raise ValueError('Unrecognized file name')
        return result


def scan_dir(path, hidden=False, recursive=False, absolute=False, regex=None):
    """
    Scans the given directory for files.

    Parameters
    ==========
    :param path: The root directory to scan.
    :param hidden: Set true if you want to include hidden files.
    :param recursive: Set true if you want to scan subdirectories within the given path.
    :param absolute: Set true if you want the absolute path of each file returned.
    :param regex: Set if you want to only return files matching the given regular expression.
    :return: List of files within the directory each file contains their absolute full path.
    """
    result = []
    for f in os.listdir(path):
        if not hidden:
            if f.startswith('.'):
                continue
        file_path = os.path.join(path, f)
        if absolute:
            file_path = os.path.abspath(file_path)
        if os.path.isfile(file_path):
            if regex and re.search(regex, f) is None:
                continue
            result.append(file_path)
        elif recursive:
            files = scan_dir(file_path, hidden, recursive, absolute, regex)
            result.extend(files)
    return result


def all_exists(file_list, parent_dir=None):
    """
    Does all the files in the list exist in the filesystem?

    Parameters
    ==========
    :param file_list: List of file names to check
    :param parent_dir: parent directory (optional)
    :return: True is all files exist
    """
    if not file_list:
        return False
    for f in file_list:
        if parent_dir is not None:
            path = os.path.join(parent_dir, f.basename())
        else:
            path = str(f)
        if not os.path.exists(path):
            return False
    return True


def has_detector(file_list, detector):
    """
    Does the file list contain a file for the given detector?

    Parameters
    ==========
    :param file_list: List of file names
    :param detector: Detector being searched
    :return: True if the list of file names include the given detector
    """
    for f in file_list:
        if f.detector == detector:
            return True
    return False


def is_complete(file_list):
    """
    Does the file list contain a file for every detector?

    Parameters
    ==========
    :param file_list: List of files that represent a detector set
    :return: True is the file list contains a file for every detector
    """
    for d in Detector:
        if not has_detector(file_list, d):
            return False
    return True


def max_version(file_list):
    """
    Returns the maximum _version of file name in the given list

    Parameters
    ==========
    :param file_list: list of file names
    :return: Largest _version number in the list
    """
    result = None
    for f in file_list:
        try:
            v = int(f.version)
            if result is None or v > result:
                result = v
        except ValueError:
            pass

    return result


def min_version(file_list):
    """
    Returns the minimum _version of file name in the given list

    Parameters
    ==========
    :param file_list: list of file names
    :return: Smallest _version number in the list
    """
    result = None
    for f in file_list:
        try:
            v = int(f.version)
            if result is None or v < result:
                result = v
        except ValueError:
            pass

    return result


def ymd_path(base, name):
    if isinstance(name, str) or isinstance(name, GbmFile):
        v = name if isinstance(name, str) else name.basename()
        m = re.match(r'.*_(?:(?:bn)?)(\d{6})(?:(\d{3})|(_\d\d)?)_.*', v, re.I | re.S)
        if m:
            d = datetime.datetime.strptime(m.group(1), "%y%m%d")
            return os.path.join(base, d.strftime('%Y-%m-%d'), os.path.basename(name))
    elif isinstance(name, Met):
        return os.path.join(base, name.datetime.strftime('%Y-%m-%d'))
    elif isinstance(name, datetime.datetime) or isinstance(name, datetime.date):
        return os.path.join(base, name.strftime('%Y-%m-%d'))
    raise ValueError("Can't parse a YMD value")
