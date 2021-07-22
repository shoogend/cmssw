
#######################################################################
# Make _cfi.py files, each corresponding to one of the datasets       #
# in the CMSSW_11_2_0_pre5 RelVal samples used by L1TrkAlgo.          #
#                                                                     #
# The list of these datasets is hard-coded below, as is the name of   #
# the directory dirName to which they will be written.                #
# The directory will be created if not yet existing.                  #
#                                                                     #
# Run from linux with:                                                #
#   python getMC_11_2_0.py                                            #
#                                                                     #
# N.B. You must have a valid GRID proxy to use this script            #
#      and must have done cmsenv !                                    #
#######################################################################

from makeDataCfiFiles import *
import os

if (not os.path.exists("getCMSdata_cfi.py")):
    print "ERROR: You must run this script from the directory containing it!"
    exit(1)

# Name of output directory (will be created if not existing)
dirName = "../../RelVal_1130/python/"

# Names of datasets to be processed.
dataList = [
"/RelValTTbar_14TeV/CMSSW_11_3_0_pre6-113X_mcRun4_realistic_v6_2026D49noPU-v1/GEN-SIM-DIGI-RAW",
"/RelValTTbar_14TeV/CMSSW_11_3_0_pre3-PU_113X_mcRun4_realistic_v3_2026D49PU200_rsb-v1/GEN-SIM-DIGI-RAW",
"/RelValSingleMuFlatPt2To100/CMSSW_11_3_0_pre6-113X_mcRun4_realistic_v6_2026D49noPU-v1/GEN-SIM-DIGI-RAW"
]

makeDataCfiFiles(dirName, dataList);
