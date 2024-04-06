# Read HOMO, LUMO, dipole moments, etc. for all metallic oxides (MOs) and acid radicals (ARs)
import math
import os
import objects
import method
import numpy as np


# Get all the bonds in the out file
def getBondsFromOutFile(strFile, listElementInfo):

    listBonds = []

    outFile = open(strFile,"r")
    outFile.seek(0,0)
    nTitleLines = 0
    bReadCoord = False
    listAtoms = []
    while 1:
        strLine = outFile.readline()
        if not strLine:
            break

        strLine = ' '.join(strLine.split())
        if 'Input orientation:' in strLine:
            bReadCoord = True
            nTitleLines = 4
            continue

        if bReadCoord and nTitleLines > 0:
            nTitleLines = nTitleLines -1
            continue

        if bReadCoord:
            listTmpInfos = strLine.split(" ")

            #The first digit must be a number or the read ends
            if not listTmpInfos[0].isnumeric():
                break

            curAtom = objects.CAtom(listTmpInfos[1])
            curAtom.m_strName = listTmpInfos[1] 
            curAtom.m_strX = listTmpInfos[3]
            curAtom.m_strY = listTmpInfos[4]
            curAtom.m_strZ = listTmpInfos[5]
            listAtoms.append(curAtom)

    #Calculate all bonds, filter by atom type and distance
    listBonds = method.getBondsFromAtoms(listAtoms, listElementInfo)
    return listBonds


#Calculate bondValencePara
def calFlexibility(listGroupInfo,listElementInfo, isNega):
    strPolarPath = "PATH/TO/posi-out-polar/"
    if isNega:
        strPolarPath = "PATH/TO/nega-out-polar/"

    listfile = os.listdir(strPolarPath)
    for curFile in listfile:
        strFileName = curFile.strip(".out")
        curFile = strPolarPath + curFile
        if not curFile.endswith("out"):
            continue

        nGroupIndex = -1
        for item in listGroupInfo:
            nGroupIndex = nGroupIndex + 1
            if item.m_strID == strFileName:
                break
        if nGroupIndex < 0:
            break

        listBonds = getBondsFromOutFile(curFile, listElementInfo)

        listFlexibilities = []
        listBondLengths = []
        for curBond in listBonds:
            listFlexibilities.append(float(curBond.m_strFlexibility))
            listBondLengths.append(float(curBond.m_strLength))

        dMeanBondLength = 152 #default: 0 (pm)
        if len(listBondLengths) > 0:
            dMeanBondLength = np.mean(listBondLengths) * 100
        dVolume = (4 * math.pi * dMeanBondLength * dMeanBondLength * dMeanBondLength) / 3

        nFlexNum = len(listFlexibilities)
        if nFlexNum == 0:
            listGroupInfo[nGroupIndex].m_strAverFlexibility = "0.0"
            listGroupInfo[nGroupIndex].m_strVolume = str(dVolume)
        else:
            dTotalFlexibility = 0.0
            for flexItem in listFlexibilities:
                dTotalFlexibility = dTotalFlexibility + flexItem
            listGroupInfo[nGroupIndex].m_strAverFlexibility = str(dTotalFlexibility/nFlexNum)
            listGroupInfo[nGroupIndex].m_strVolume = str(dVolume)

    return  listGroupInfo


#Get the GroupInfo directly from the full text
def readPostiveGroupInfo(strFile):
    listGroupInfo = []
    file = open(strFile,"r")
    while 1:
        strLine = file.readline()
        if not strLine:
            break

        strLine = strLine.strip()
        listInfos = strLine.split(",")
        if len(listInfos) < 25:
            continue

        curGroup = objects.CPostiveGroupInfo(listInfos[0])
        curGroup.m_strName = listInfos[1]
        curGroup.m_strAtomicNumber = listInfos[2]
        curGroup.m_strCharge = listInfos[3]
        curGroup.m_strMultiplicity = listInfos[4]
        curGroup.m_strHomo = listInfos[5]
        curGroup.m_strLumo = listInfos[6]
        curGroup.m_strLumoHomo = listInfos[7]
        curGroup.m_strDipoleTotal = listInfos[8]
        curGroup.m_strQuadrupoleXX = listInfos[9]
        curGroup.m_strQuadrupoleYY = listInfos[10]
        curGroup.m_strQuadrupoleZZ = listInfos[11]
        curGroup.m_strQuadrupoleXY = listInfos[12]
        curGroup.m_strQuadrupoleXZ = listInfos[13]
        curGroup.m_strQuadrupoleYZ = listInfos[14]
        curGroup.m_strAnisoQuadrupole = listInfos[15]
        curGroup.m_listPolarFreq = listInfos[16].split(";")
        curGroup.m_listIsoPolar = listInfos[17].split(";")
        curGroup.m_listAnisoPolar = listInfos[18].split(";")
        curGroup.m_listHyperPolarX = listInfos[19].split(";")
        curGroup.m_listHyperPolarY = listInfos[20].split(";")
        curGroup.m_listHyperPolarZ = listInfos[21].split(";")
        curGroup.m_listTotalHyperPolar = listInfos[22].split(";")
        curGroup.m_listVectorHyperPolar = listInfos[23].split(";")
        curGroup.m_strAverFlexibility = listInfos[24]
        curGroup.m_strVolume = listInfos[25]
        listGroupInfo.append(curGroup)

    file.close()

    return listGroupInfo


#Read all the information about the Group from the txt file.
#One can use readGroupFromOutFiles to read from the out file first, 
#and then use this function to read the information after it is complete.
def readAllGroupInfo():
    strPositiveFile = objects.BASISDIR + "public/positive-group-in.csv"
    listPosiGroups = readPostiveGroupInfo(strPositiveFile)
    strPositiveFile = objects.BASISDIR + "public/negative-group-in.csv"
    listNegaGroups = readPostiveGroupInfo(strPositiveFile)

    return listNegaGroups,listPosiGroups

#Save all information
def saveAllGroupInfo(listGroupInfo, isNega):
    strOutFile = objects.BASISDIR + "public/positive-group-out.csv"
    if isNega:
        strOutFile = objects.BASISDIR + "public/negative-group-out.csv"
    file = open(strOutFile,"a+")
    file.write("\n--------Groups---------------\n")
    file.write(objects.CGroupInfo.joinTitle(","))
    file.write("\n")
    for curGroup in listGroupInfo:
        file.write(curGroup.joinToString(","))
        file.write("\n")
    file.close()

# Export results directly to a csv file when called individually
if __name__ == '__main__':
    print('getGroupInfo')




