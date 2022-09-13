# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 11:29:03 2022

@author: bener807
"""

import json
import numpy
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, 'C:/python/useful_definitions/')
import useful_defs as udfs
udfs.set_nes_plot_style()


def json_read_dictionary(file_name):
    """Read .json file."""
    with open(file_name, 'r') as handle:
        j = json.load(handle)
    return j


if __name__ == '__main__':
    # Read data
    file_name = 'output_files/tofu_drf_kin.json'
    j = json_read_dictionary(file_name)

    udfs.plot_matrix(j['matrix'], (j['x'], j['y']), log=True,
                     xlabel=f'E ({j["x_unit"]})',
                     ylabel=f'$t_{{TOF}}$ ({j["y_unit"]})')
