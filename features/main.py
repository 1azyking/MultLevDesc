# The main function for data processing
import elementInfo, featureInfo, crystalInfo, groupInfo

from objects import BASISDIR
from collections import Counter
import numpy as np


if __name__ == '__main__':
    # Retrieve all element information
    listElementInfo = elementInfo.readAllElementInfo()

    # Retrieve Group
    listNegaGroups = []
    listPosiGroups = []
    # listNegaGroups = groupInfo.readGroupFromOutFiles(listElementInfo, True)  # Retrieve all Group information from the out file
    # listPosiGroups = groupInfo.readGroupFromOutFiles(listElementInfo, False)  # Retrieve all Group information from the out file
    listNegaGroups, listPosiGroups = groupInfo.readAllGroupInfo()
    # groupInfo.saveAllGroupInfo(listNegaGroups, True)
    # groupInfo.saveAllGroupInfo(listPosiGroups, False)

    # Retrieve all NLO crystal information
    listCrystalInfo = crystalInfo.readAllCrystalInfo(listElementInfo)
    #crystalInfo.saveAllCrystalInfo(listCrystalInfo)

    # Construct descriptors corresponding to crystals
    listFeatureInfo = featureInfo.createAllFeatureInfo(listPosiGroups, listNegaGroups, listElementInfo, listCrystalInfo)
    featureInfo.saveAllFeatureInfo(listFeatureInfo)
