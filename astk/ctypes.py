# -*- coding: utf-8 -*-

"""
astk.ctypes
~~~~~~~~~~~~~~~~~
This module provide custom  parameter types.
"""
from pathlib import Path
from typing import TypeVar, Dict, Iterable, Tuple, Sequence


FilePath = TypeVar('FilePath', str, Path)
