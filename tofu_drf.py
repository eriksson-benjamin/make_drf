# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 14:04:31 2022

@author: bener807
"""


import uproot
import os
import numpy as np
import json


def json_write_dictionary(file_name, to_save, check=True):
    """Write dictionary to .json file."""
    inp = ' '
    if check:
        if os.path.exists(file_name):
            inp = input((f'{file_name} already exists.\nDo you want to '
                         f'overwrite it? [y/n] '))
            if inp != 'y':
                print('File was not overwritten.')
                return 0
    with open(file_name, 'w') as handle:
        json.dump(to_save, handle)


def read_data(energy, S1, light_yield):
    """Return flight times, deposited erg, S2 mask, and kinematic cuts mask."""
    file = f'{energy}keV_S1%3A{S1+1}_ToFuMatrix.root'
    with uproot.open(f'data/{file}') as handle:
        # Grab tree
        tree = handle['tree3D;1']

        # Grab data from tree
        t_tof = tree['timeWres'].array(library='np')
        if light_yield:
            E_S1 = tree['S1EdMeVee'].array(library='np')
            E_S2 = tree['S2EdMeVee'].array(library='np')
        else:
            E_S1 = tree['S1Ed'].array(library='np')
            E_S2 = tree['S2Ed'].array(library='np')
        i_S2 = tree['S2c'].array(library='np')
        kin_S1 = tree['S1kin'].array(library='np')
        kin_S2 = tree['S2kin'].array(library='np')

    return 0.4 * t_tof, E_S1, E_S2, i_S2, kin_S1, kin_S2


def mask_S1(E_S1, threshold):
    """Return mask for S1 energy threshold."""
    mask = (E_S1 >= threshold)

    return mask


def mask_S2(E_S2, t_tof, i_S2, detector, threshold):
    """
    Return mask for S1 energy threshold.

    Parameters
    ----------
    E_S2 : ndarray,
         1D array of deposited energies in S2s.
    t_tof : ndarray,
          1D array of S1-S2 flight times.
    i_S2 : ndarray,
         1D array of numbers to identify which S2 interaction ocurred in.
    detector : int,
             Number between 0-31 to identify which detector to analyze.
    threshold : float,
              Threshold for given S2.
    """
    # Mask for detector number
    det_mask = i_S2 == detector

    # Mask for energy threshold
    erg_mask = E_S2 >= threshold

    mask = det_mask & erg_mask
    return E_S2[mask], t_tof[mask], i_S2[mask]


def return_thresholds(light_yield):
    """Return detector energy thresholds (MeV)."""
    if light_yield:
        f_name = 'input_files/thresholds_MeVee.txt'
    else:
        f_name = 'input_files/thresholds_MeV.txt'
    thresholds = np.loadtxt(f_name, usecols=[1])

    # Convert to keV(ee)
    thresholds *= 1000
    return thresholds[:5], thresholds[5:]


def mask_kincut(kin_S1, kin_S2):
    """Return mask for kinematic cuts."""
    mask = ((kin_S1 == 1) & (kin_S2 == 1))

    return mask


def masked(mask, t_tof, E_S1, E_S2, i_S2):
    """Return masked events."""
    return t_tof[mask], E_S1[mask], E_S2[mask], i_S2[mask]


def save_json(drf, t_bin_centres, E_bin_centres, info, name, file_name,
              light_yield):
    """Save DRF as json file."""
    # Convert from numpy array to list
    drf_list = [list(col) for col in drf.T]
    E_list = list(np.array(E_bin_centres, dtype='float'))
    t_list = list(t_bin_centres)

    # Create dictionary to save
    to_save = {'matrix': drf_list, 'x_unit': 'keV', 'y_unit': 'ns',
               'x': E_list, 'y': t_list, 'name': name, 'info': info}

    # Write to .json file
    json_write_dictionary(file_name, to_save)


def main(kinematic_cuts, light_yield):
    """Calculate detector response function."""
    # Bins
    E_bin_centres = np.arange(1000, 18000 + 50, 50)
    t_bin_centres = np.arange(0, 200, 0.4)
    t_bin_edges = np.arange(-0.2, 200.2, 0.4)

    # Initialize DRF matrix
    drf_matrix = np.zeros([len(t_bin_centres), len(E_bin_centres)])

    # S1/S2 thresholds
    S1_thr, S2_thr = return_thresholds(light_yield)

    for drf_col, energy in enumerate(E_bin_centres):
        print(f'Processing: {energy} keV')
        for S1 in range(5):
            # Get data
            t_tof, E_S1, E_S2, i_S2, kin_S1, kin_S2 = read_data(energy, S1,
                                                                light_yield)

            # Create mask for S1 energy thresholds
            mask_1 = mask_S1(E_S1, S1_thr[S1])

            # Create mask for kinematic cuts
            if kinematic_cuts:
                mask_2 = mask_kincut(kin_S1, kin_S2)
            else:
                mask_2 = np.array(len(mask_1) * [True])

            # Sort out data outside mask
            mask = mask_1 & mask_2
            t_tof, E_S1, E_S2, i_S2 = masked(mask, t_tof, E_S1, E_S2, i_S2)

            # Sort through S2's, take thresholds into account.
            for S2 in range(32):
                E_masked, t_masked, i_masked = mask_S2(E_S2, t_tof, i_S2, S2,
                                                       S2_thr[S2])

                # Histogram and add to DRF matrix
                t_counts, _ = np.histogram(t_masked, bins=t_bin_edges)
                drf_matrix[:, drf_col] += t_counts

    return drf_matrix, t_bin_centres, E_bin_centres


if __name__ == '__main__':
    kinematic_cuts = True
    light_yield = True
    drf_matrix, t_bin_centres, E_bin_centres = main(kinematic_cuts,
                                                    light_yield)
    info = ('DRF for TOFu, individual thresholds and energy dependent time '
            'resolution applied to each S1 and S2. Kinematic cuts are '
            'applied with scaling factors a=1, b=1, c=1. Units of light '
            'yield used when generating DRF.')
    # info = ('DRF for TOFu, individual thresholds and energy dependent time '
    #         'resolution applied to each S1 and S2. No kinematic cuts are '
    #         'applied. Units of light yield used when generating DRF.')
    name = 'TOFu DRF'
    file_names = ['tofu_drf.json', 'tofu_drf_kin.json',
                  'tofu_drf_scaled_kin.json']

    save_json(drf_matrix, t_bin_centres, E_bin_centres, info, name,
              file_names[1], light_yield)
