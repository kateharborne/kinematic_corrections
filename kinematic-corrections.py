#!/usr/bin/env python

import numpy as np
import h5py
import os
from argparse import ArgumentParser
import csv
import pandas

# Define arguments ###
parser = ArgumentParser()
parser.add_argument("-i", "--infile", action="store", dest="inputFile", type=str, default=None,
                    help="The path to the CSV file containing the kinematic properties of observaed galaxies.")
parser.add_argument("-o", "--outfile", action="store", dest="outputFile", type=str, default=None,
                    help="The path to the output CSV file with added corrected kinematic values")

# Parse arguments ###
args = parser.parse_args()
kin_uc = pandas.read_csv(args.inputFile)

if "psf_over_re" not in kin_uc: # if psf_over_re is not computed in the file, calculate from Re and fwhm
    if {"Re", "fwhm"}.issubset(kin_uc):  
        kin_uc = kin_uc.assign(psf_over_re = (kin_uc.fwhm / 2.355) / kin_uc.Re)
    else: 
        print("\nERROR: No blurring information available! \n "
              "Please provide measurement radius (\"Re\") and fwhm (\"fwhm\")\n"
              "or the fraction sigma_PSF / Re (\"psf_over_re\") in your --infile.\n"
              " --- ABORTING ---\n")
        exit()  # Check that blurring information has been supplied.

if all([item not in kin_uc for item in {"obs_lr", "obs_elr", "obs_vsig"}]):
    print("\nERROR: No kinematic information available! \n "
          "Please provide observed lambdaR (\"obs_lr\")\n"
          "or observed elliptical lambdaR (\"obs_elr\")\n"
          "or observed V/sigma (\"obs_vsig\") in your --infile.\n"
          " --- ABORTING ---\n")
    exit()  # Check that kinematic information has been supplied.

if "e" not in kin_uc:
    print("\nERROR: No ellipticity information available! \n "
          "Please provide measurement ellipticity (\"e\") in your --infile.\n"
          " --- ABORTING ---\n")
    exit()  # Check that kinematic information has been supplied.

if "n" not in kin_uc:
    print("\nERROR: No sersic index information available! \n "
          "Please provide sersic index (\"n\") in your --infile.\n"
          " --- ABORTING ---\n")
    exit()  # Check that kinematic information has been supplied.

if "Reff_fac" not in kin_uc:
    print("\nWARNING: No measurement radius factor information available! \n "
          "Corrected kinematics will be calculated assuming observations made at 1 Reff.\n"
          "If this is not the case, please provide measurement radius factor (\"Reff_fac\")\n"
          "in your --infile for other appropriate correction.\n")
    kin_uc = kin_uc.assign(Reff_fac = 1)

# Correction constants
lr_f1_a = 7.48; lr_f1_b = 4.08; lr_f1_c = 1.60; lr_f1_d = 2.89
lr_f1_e = lr_f1_a / (1 + np.exp(lr_f1_d))
lr_f2_obse = 0.10033595; lr_f2_logn = -0.22313183; lr_f2_reff = -0.12270405; lr_f2_C = 0.2188

elr_f1_a = 7.49; elr_f1_b = 4.01; elr_f1_c = 1.57; elr_f1_d = 2.84
elr_f1_e = elr_f1_a / (1 + np.exp(elr_f1_d))
elr_f2_loge = 0.01638699; elr_f2_logn = -0.19020806; elr_f2_reff = -0.13429812; elr_f2_C = 0.27955548

vsig_f1_a = 7.55; vsig_f1_b = 4.42; vsig_f1_c = 1.55; vsig_f1_d = 2.73
vsig_f1_e = vsig_f1_a / (1 + np.exp(vsig_f1_d))
vsig_f2_obse = -0.10142709; vsig_f2_logn = 0.02407725; vsig_f2_reff = -0.05636119; vsig_f2_C = 0.11757231

if "obs_lr" in kin_uc:
    corr_lr_f1 = (lr_f1_a / (1 + np.exp((lr_f1_b * (kin_uc.psf_over_re)**lr_f1_c) + lr_f1_d))) - lr_f1_e
    corr_lr_f2 = (lr_f2_obse * kin_uc.e) + (lr_f2_logn * np.log10(kin_uc.n)) + (lr_f2_reff * kin_uc.Reff_fac) + lr_f2_C
    dlr = corr_lr_f1 + (kin_uc.psf_over_re * corr_lr_f2)
    kin_uc = kin_uc.assign(corr_lr = 10**(np.log10(kin_uc.obs_lr) - dlr))

if "obs_elr" in kin_uc:
    corr_elr_f1 = (elr_f1_a / (1 + np.exp((elr_f1_b * (kin_uc.psf_over_re)**elr_f1_c) + elr_f1_d))) - elr_f1_e
    corr_elr_f2 = (elr_f2_loge * np.log10(kin_uc.e)) + (elr_f2_logn * np.log10(kin_uc.n)) + (elr_f2_reff * kin_uc.Reff_fac) + elr_f2_C
    delr = corr_elr_f1 + (kin_uc.psf_over_re * corr_elr_f2)
    kin_uc = kin_uc.assign(corr_elr = 10**(np.log10(kin_uc.obs_elr) - delr))

if "obs_vsig" in kin_uc:
    corr_vsig_f1 = (vsig_f1_a / (1 + np.exp((vsig_f1_b * (kin_uc.psf_over_re)**vsig_f1_c) + vsig_f1_d))) - vsig_f1_e
    corr_vsig_f2 = (vsig_f2_obse * kin_uc.e) + (vsig_f2_logn * np.log10(kin_uc.n)) + (vsig_f2_reff * kin_uc.Reff_fac) + vsig_f2_C
    dvsig = corr_vsig_f1 + (3 * kin_uc.psf_over_re * corr_vsig_f2)
    kin_uc = kin_uc.assign(corr_vsig = 10**(np.log10(kin_uc.obs_vsig) - dvsig))

kin_uc.to_csv(args.outputFile, index=False)
print("New data file with added corrected kinematics written at: ", args.outputFile)