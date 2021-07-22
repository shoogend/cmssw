
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
dirName = "../../RelVal_1120/python/"

# Names of datasets to be processed.
dataList = [
'/RelValTTbar_14TeV/CMSSW_11_2_0_pre5-110X_mcRun4_realistic_v3_2026D49noPU-v1/GEN-SIM-DIGI-RAW',
'/RelValTTbar_14TeV/CMSSW_11_2_0_pre5-PU25ns_110X_mcRun4_realistic_v3_2026D49PU200-v1/GEN-SIM-DIGI-RAW',

'/RelValSingleElectronFlatPt1p5To8/CMSSW_11_2_0_pre6-112X_mcRun4_realistic_v2_2026D49noPU_L1T-v1/GEN-SIM-DIGI-RAW',
'/RelValSingleElectronFlatPt1p5To8/CMSSW_11_2_0_pre6-PU25ns_112X_mcRun4_realistic_v2_2026D49PU200_L1T-v2/GEN-SIM-DIGI-RAW',
'/RelValSingleEFlatPt2To100/CMSSW_11_2_0_pre6-112X_mcRun4_realistic_v2_2026D49noPU_L1T-v1/GEN-SIM-DIGI-RAW',
'/RelValSingleEFlatPt2To100/CMSSW_11_2_0_pre6-PU25ns_112X_mcRun4_realistic_v2_2026D49PU200_L1T-v2/GEN-SIM-DIGI-RAW',

'/RelValSingleMuFlatPt1p5To8/CMSSW_11_2_0_pre6-112X_mcRun4_realistic_v2_2026D49noPU_L1T-v1/GEN-SIM-DIGI-RAW',
'/RelValSingleMuFlatPt1p5To8/CMSSW_11_2_0_pre6-PU25ns_112X_mcRun4_realistic_v2_2026D49PU200_L1T-v2/GEN-SIM-DIGI-RAW',
'/RelValSingleMuFlatPt2To100/CMSSW_11_2_0_pre6-112X_mcRun4_realistic_v2_2026D49noPU_L1T-v1/GEN-SIM-DIGI-RAW',
'/RelValSingleMuFlatPt2To100/CMSSW_11_2_0_pre6-PU25ns_112X_mcRun4_realistic_v2_2026D49PU200_L1T-v2/GEN-SIM-DIGI-RAW',

'/RelValDisplacedMuPt1p5To8/CMSSW_11_2_0_pre6-112X_mcRun4_realistic_v2_2026D49noPU_L1T-v1/GEN-SIM-DIGI-RAW',
'/RelValDisplacedMuPt1p5To8/CMSSW_11_2_0_pre6-PU25ns_112X_mcRun4_realistic_v2_2026D49PU200_L1T-v2/GEN-SIM-DIGI-RAW',
'/RelValDisplacedMuPt2To100/CMSSW_11_2_0_pre6-112X_mcRun4_realistic_v2_2026D49noPU_L1T-v1/GEN-SIM-DIGI-RAW',
'/RelValDisplacedMuPt2To100/CMSSW_11_2_0_pre6-PU25ns_112X_mcRun4_realistic_v2_2026D49PU200_L1T-v2/GEN-SIM-DIGI-RAW'
]

makeDataCfiFiles(dirName, dataList);
