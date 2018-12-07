#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# utils.py
#
# This file is part of pyplanes, a software distributed under the MIT license.
# For any question, please contact mathieu@matael.org.
#
# Copyright (c) 2018 The pyplanes authors
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#

import os
import yaml

# from .eqf import EqFluidJCA
# from .elastic import Elastic
# from .pem import PEM
# from .screen import Screen
from .fluid import Fluid

__MEDIUMCLASSES_MAP = {
    _.MEDIUM_TYPE: _ for _ in [Fluid]
}


def from_yaml(filename, force=None):
    """Reads medium definition from YAML file filename. Raises an IOError if the file
    is not found, ValueError if medium type is not known and a LookupError if the
    parameter definition is incomplete.

    One may set the optional argument 'force' to a medium type to force loading this one.

    Parameters
    ----------
    filename: str
        path to the file defining the medium

    force: type
        forced type (class) of the output

    Returns
    -------
        A FEMMedium's subclass instance
    """

    if not os.path.exists(filename):
        raise IOError('Unable to locate file {}'.format(filename))

    with open(filename) as fh:
        yaml_data = yaml.load(fh)

    if yaml_data.get('medium_type') is None and force is None:
        raise LookupError('Unspecified medium type')

    if force is None:
        medium_class = __MEDIUMCLASSES_MAP.get(yaml_data['medium_type'])
    else:
        medium_class = force
    if medium_class is None:
        raise ValueError('Medium type is not known')
    else:
        medium = medium_class()
        medium.from_dict(yaml_data)
        return medium
