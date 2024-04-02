# Read HOMO, LUMO, dipole moments, etc. for all cathode and anode groups
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


# Get GroupInfo directly from Gaussian output file, # Get HOMO, LUMO
def readEnergyFromOut(strPath):

    listAlphaOcc = []
    listAlphaVirt = []
    listBetaOcc = []
    listBetaVirt = []

    listGroupInfo = []
    listfile = os.listdir(strPath)
    for curFile in listfile:

        curFileName = curFile.strip(".out")
        curGroup = objects.CPostiveGroupInfo(curFileName)

        curFile = strPath + curFile
        file = open(curFile,"r")
        listAlphaOcc.clear()
        listAlphaVirt.clear()
        listBetaOcc.clear()
        listBetaVirt.clear()
        while 1:
            strLine = file.readline()
            if not strLine:
                break

            strLine = ' '.join(strLine.split())
            if 'Charge = ' in strLine:
                listTmpInfo = strLine.split(" ")
                curGroup.m_strCharge = listTmpInfo[2]
                curGroup.m_strMultiplicity = listTmpInfo[5]
            if 'Alpha occ. eigenvalues' in strLine:
                listAlphaOcc = listAlphaOcc + strLine.split(' ')
            if 'Alpha virt. eigenvalues' in strLine:
                listAlphaVirt = listAlphaVirt + strLine.split(' ')
            elif 'Beta occ. eigenvalues' in strLine:
                listBetaOcc = listBetaOcc + strLine.split(' ')
            elif 'Beta virt. eigenvalues' in strLine:
                listBetaVirt = listBetaVirt + strLine.split(' ')

        strMaxOcc = "-999999"
        listAlphaOcc.reverse()
        for strEnergy in listAlphaOcc:
            strEnergy = strEnergy.strip()
            strEnergy = strEnergy.strip('\n')
            if method.is_number(strEnergy):
                strMaxOcc = strEnergy
                break

        strMinVirt = "999999"
        for strEnergy in listAlphaVirt:
            strEnergy = strEnergy.strip()
            strEnergy = strEnergy.strip('\n')
            if method.is_number(strEnergy):
                strMinVirt = strEnergy
                break

        curGroup.m_strHomo = strMaxOcc
        curGroup.m_strLumo = strMinVirt
        curGroup.m_strLumoHomo = str(float(strMinVirt) - float(strMaxOcc))

        strMaxOcc = "-999999"
        listBetaOcc.reverse()
        for strEnergy in listBetaOcc:
            strEnergy = strEnergy.strip()
            strEnergy = strEnergy.strip('\n')
            if method.is_number(strEnergy):
                strMaxOcc = strEnergy
                break

        strMinVirt = "999999"
        for strEnergy in listBetaVirt:
            strEnergy = strEnergy.strip()
            strEnergy = strEnergy.strip('\n')
            if method.is_number(strEnergy):
                strMinVirt = strEnergy
                break

        if float(curGroup.m_strLumoHomo) > float(strMinVirt)-float(strMaxOcc):
            curGroup.m_strHomo = strMaxOcc
            curGroup.m_strLumo = strMinVirt
            curGroup.m_strLumoHomo = str(float(strMinVirt) - float(strMaxOcc))

        #Hartree-->eV
        curGroup.m_strHomo = str(float(curGroup.m_strHomo) * 27.2116)
        curGroup.m_strLumo = str(float(curGroup.m_strLumo) * 27.2116)
        curGroup.m_strLumoHomo = str(float(curGroup.m_strLumoHomo) * 27.2116)

        listGroupInfo.append(curGroup)

        file.close()

    return listGroupInfo


#Get GroupInfo directly from Gaussian output file
#dipole moment, polarizability, hyperpolarizability, etc
def readDiPoleAndPolarFromOut(strPath, listGroupInfo):

    if len(listGroupInfo) < 1:
        return []

    listfile = os.listdir(strPath)
    for curFile in listfile:

        nIndex = -1
        strFileName = curFile.strip(".out")
        for item in listGroupInfo:
            nIndex = nIndex + 1
            if item.m_strID == strFileName:
                break
        if nIndex < 0:
            break

        nIsDiPole = 0
        nIsQuadrupole = 0
        nIsAlpha = 0
        nIsStaticBeta = 0
        nIsDynamicBeta = 0
        nIsStartRead = 0
        nIsReadInputOrient = 1
        curFile = strPath + curFile
        file = open(curFile,"r")
        while 1:
            strLine = file.readline()
            if not strLine:
                break

            strLine = ' '.join(strLine.split())

            #Read after the line: Population analysis using the SCF density
            if nIsStartRead < 1 and 'Population analysis using the SCF density' in strLine:
                nIsStartRead = 1
                continue

            if nIsStartRead < 1:
                continue

            if nIsReadInputOrient > 0:

                # Dipole moment(Debye)
                if 'Dipole moment (field-independent basis, Debye)' in strLine:
                    nIsDiPole = 1
                    continue

                # Quadrupole moment(Debye/Ang)
                if 'Quadrupole moment (field-independent basis, Debye-Ang)' in strLine:
                    nIsQuadrupole = 1
                    continue

                #Frequencies(Hartree-->nm)
                if 'Frequencies=' in strLine:
                    listTmpInfo = strLine.split(' ')
                    del listTmpInfo[0]
                    del listTmpInfo[0]
                    listGroupInfo[nIndex].m_listPolarFreq.append("0.00000")
                    for item in listTmpInfo:
                        listGroupInfo[nIndex].m_listPolarFreq.append(str(1239.8424122 / (float(item) * 27.2116)))

                if 'Dipole polarizability, Alpha (input orientation)' in strLine:
                    nIsAlpha = 1
                    continue

                if 'Beta(0;0,0):' in strLine:
                    nIsStaticBeta = 1
                    continue

                if 'Beta(-2w;w,w)' in strLine:
                    nIsDynamicBeta = 1
                    continue

                if nIsDiPole == 1:
                    if 'Tot=' in strLine:
                        strLine = strLine.replace("D","E")#The exponential representation in Gaussian is converted to E
                        listTmpInfo = strLine.split(" ")
                        listGroupInfo[nIndex].m_strDipoleTotal = listTmpInfo[7]
                        nIsDiPole = -1

                if nIsQuadrupole == 1:
                    if 'XX=' in strLine:
                        strLine = strLine.replace("D","E")
                        listTmpInfo = strLine.split(" ")
                        listGroupInfo[nIndex].m_strQuadrupoleXX = listTmpInfo[1]
                        listGroupInfo[nIndex].m_strQuadrupoleYY = listTmpInfo[3]
                        listGroupInfo[nIndex].m_strQuadrupoleZZ = listTmpInfo[5]
                    elif 'XY' in strLine:
                        strLine = strLine.replace("D","E")
                        listTmpInfo = strLine.split(" ")
                        listGroupInfo[nIndex].m_strQuadrupoleXY = listTmpInfo[1]
                        listGroupInfo[nIndex].m_strQuadrupoleXZ = listTmpInfo[3]
                        listGroupInfo[nIndex].m_strQuadrupoleYZ = listTmpInfo[5]
                        nIsQuadrupole = 0

                if nIsAlpha == 1:
                    strLine = strLine.replace("D","E")
                    listTmpInfo = strLine.split(" ")
                    if 'aniso' in strLine:
                        listGroupInfo[nIndex].m_listAnisoPolar.append(listTmpInfo[1])
                    elif 'iso' in strLine:
                        listGroupInfo[nIndex].m_listIsoPolar.append(listTmpInfo[1])
                    elif len(strLine) < 2:
                        nIsAlpha = 0

                if nIsStaticBeta == 1:
                    strLine = strLine.replace("D","E")
                    listTmpInfo = strLine.split(" ")
                    if 'x' in strLine and len(listTmpInfo[0]) < 3:
                        listGroupInfo[nIndex].m_listHyperPolarX.append(listTmpInfo[1])
                    elif 'y' in strLine and len(listTmpInfo[0]) < 3:
                        listGroupInfo[nIndex].m_listHyperPolarY.append(listTmpInfo[1])
                    elif 'z' in strLine and '||' not in strLine and len(listTmpInfo[0]) < 3:
                        listGroupInfo[nIndex].m_listHyperPolarZ.append(listTmpInfo[1])
                    elif 'zzz' in strLine:
                        nIsStaticBeta = 0

                if nIsDynamicBeta == 1:
                    strLine = strLine.replace("D","E")
                    listTmpInfo = strLine.split(" ")
                    if 'x' in strLine and len(listTmpInfo[0]) < 3:
                        listGroupInfo[nIndex].m_listHyperPolarX.append(listTmpInfo[1])
                    elif 'y' in strLine and len(listTmpInfo[0]) < 3:
                        listGroupInfo[nIndex].m_listHyperPolarY.append(listTmpInfo[1])
                    elif 'z' in strLine and '||' not in strLine and len(listTmpInfo[0]) < 3:
                        listGroupInfo[nIndex].m_listHyperPolarZ.append(listTmpInfo[1])
                    elif len(strLine) < 2:
                        nIsDynamicBeta = 0
                        nIsReadInputOrient = 0 #End here
            else:
                if 'Beta(0;0,0):' in strLine:
                    nIsStaticBeta = 1
                    continue

                if 'Beta(-2w;w,w)' in strLine:
                    nIsDynamicBeta = 1
                    continue

                if nIsStaticBeta == 1:
                    strLine = strLine.replace("D","E")
                    listTmpInfo = strLine.split(" ")
                    if '||' in strLine  and 'z' not in strLine and len(listTmpInfo[0]) < 3:
                        listGroupInfo[nIndex].m_listVectorHyperPolar.append(listTmpInfo[1])
                    elif 'zzz' in strLine:
                        nIsStaticBeta = 0

                if nIsDynamicBeta == 1:
                    strLine = strLine.replace("D","E")
                    listTmpInfo = strLine.split(" ")
                    if '||' in strLine and 'z' not in strLine and len(listTmpInfo[0]) < 3:
                        listGroupInfo[nIndex].m_listVectorHyperPolar.append(listTmpInfo[1])
                    elif len(strLine) < 2:
                        nIsDynamicBeta = 0
                        break

        file.close()

    return listGroupInfo


#Directly get GroupInfo from Gaussian output files
def readGroupFromOutFiles(listElementInfo, isNega):
    strEnergyPath = "PATH/TO/posi-out-energy/"
    strPolarPath = "PATH/TO/posi-out-polar/"
    if isNega:
        strEnergyPath = "PATH/TO/nega-out-energy/"
        strPolarPath = "PATH/TO/nega-out-polar/"

    listGroupInfo = readEnergyFromOut(strEnergyPath)
    listGroupInfo = readDiPoleAndPolarFromOut(strPolarPath, listGroupInfo)

    #Calculate the bondValencePara and Flexibility for each Group
    listGroupInfo = calFlexibility(listGroupInfo,listElementInfo, isNega)

    nIndex = -1
    for item in listGroupInfo:
        nIndex = nIndex + 1
        listGroupInfo[nIndex].check()
        listGroupInfo[nIndex].calDerivatives()

    return listGroupInfo

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
    strPositiveFile = objects.BASISDIR + "positive-group-in.csv"
    listPosiGroups = readPostiveGroupInfo(strPositiveFile)
    strPositiveFile = objects.BASISDIR + "negative-group-in.csv"
    listNegaGroups = readPostiveGroupInfo(strPositiveFile)

    return listNegaGroups,listPosiGroups

#Save all information
def saveAllGroupInfo(listGroupInfo, isNega):
    strOutFile = objects.BASISDIR + "positive-group-out.csv"
    if isNega:
        strOutFile = objects.BASISDIR + "negative-group-out.csv"
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




