# Obtain all crystal information from an external file
import math
import os
import sys
import objects
import method
import numpy as np

# Read all SpaceGroup-related information from a txt file
# Note: Some SpaceGroup symbols may be expressed differently, which are not currently handled
def readAllSpaceGroupInfo():

    dictSpaceGroupInfo = {}

    # Select the directory where the txt file is located
    strSpaceGroupFile = objects.BASISDIR + "public/spacegroup-in.csv"
    file = open(strSpaceGroupFile,"r")
    while 1:
        strLine = file.readline()
        if not strLine:
            break

        strLine = strLine.strip()
        listInfos = strLine.split(",")

        for curInfo in listInfos:
            listdetails = curInfo.split(";")
            if len(listdetails) < 2:
                continue
            # Note: name serves as an index
            if listdetails[1] not in dictSpaceGroupInfo.keys():
                dictSpaceGroupInfo[listdetails[1]] = listdetails[0]

    file.close()

    return dictSpaceGroupInfo


# Read the unit cell parameters of all crystals from cif files and save to a file
# listCifFiles is used to filter crystals; when it's empty, all crystals in the folder are read
def readCrystalParaFromCif(listCifFiles):

    strCifPath = objects.BASISDIR + "cif/"
    listCifParas = []

    # Traverse all files in the folder
    if len(listCifFiles) < 1:
        listCifFiles = os.listdir(strCifPath)

    # Obtain encoding based on the system
    strCodeType = "utf-8"
    sysType = sys.platform
    if sysType == "Windows":
        strCodeType = "utf-8"
    elif sysType == "Linux":
        strCodeType = "gbk"

    for curFile in listCifFiles:
        listLengthPara = []
        listAngelPara = []
        strVolume = ""
        strSpaceID = ""

        strFileName = curFile.strip(".cif")
        curFile = strCifPath + curFile
        # Use utf-8 to avoid 'gbk' codec can't decode byte 0xa9 error
        # Alternatively, convert the file to ANSI beforehand; Linux defaults to utf-8, while Windows defaults to gbk
        file = open(curFile,"r",encoding=strCodeType)
        while 1:
            strLine = file.readline()
            if not strLine:
                break

            # Handle cases where multiple spaces are used as separators
            strLine = ' '.join(strLine.split())
            if '_cell_length_' in strLine:
                strLine = strLine.strip()
                strLine = strLine.strip(")")
                strLine = strLine.replace("(", "")
                listTmpInfo = strLine.split(" ")
                listLengthPara.append(listTmpInfo[1])
            elif '_cell_angle_' in strLine:
                strLine = strLine.strip()
                strLine = strLine.strip(")")
                strLine = strLine.replace("(", "")
                listTmpInfo = strLine.split(" ")
                listAngelPara.append(listTmpInfo[1])
            elif '_cell_volume' in strLine:
                strLine = strLine.strip()
                strLine = strLine.strip(")")
                strLine = strLine.replace("(", "")
                listTmpInfo = strLine.split(" ")
                strVolume = listTmpInfo[1]
            elif '_space_group_IT_number' in strLine:
                strLine = strLine.strip()
                strLine = strLine.strip(")")
                strLine = strLine.replace("(", "")
                listTmpInfo = strLine.split(" ")
                strSpaceID = listTmpInfo[1]

        file.close()

        listLengthDiffs = []
        listLengthDiffs.append(math.fabs((float(listLengthPara[0].strip()) - float(listLengthPara[1].strip()))))
        listLengthDiffs.append(math.fabs((float(listLengthPara[0].strip()) - float(listLengthPara[2].strip()))))
        listLengthDiffs.append(math.fabs((float(listLengthPara[1].strip()) - float(listLengthPara[2].strip()))))
        # listLengthDiffs.append(math.fabs(float(listCurPara[0].strip())))
        # listLengthDiffs.append(math.fabs(float(listCurPara[1].strip())))
        # listLengthDiffs.append(math.fabs(float(listCurPara[2].strip())))
        listLengthDiffs.sort()

        # Normalize
        # if ((listLengthDiffs[2] - listLengthDiffs[0]) < 1e-3):
        #     listLengthDiffs[0] = 0.0
        #     listLengthDiffs[1] = 0.0
        #     listLengthDiffs[2] = 0.0
        # else:
        #     listLengthDiffs[1] = (listLengthDiffs[1] - listLengthDiffs[0])/(listLengthDiffs[2] - listLengthDiffs[0])
        #     listLengthDiffs[2] = 1.0
        #     listLengthDiffs[0] = 0.0

        listAngelDiffs = []
        listAngelDiffs.append(math.fabs((float(listAngelPara[0].strip()) - float(listAngelPara[1].strip()))))
        listAngelDiffs.append(math.fabs((float(listAngelPara[0].strip()) - float(listAngelPara[2].strip()))))
        listAngelDiffs.append(math.fabs((float(listAngelPara[1].strip()) - float(listAngelPara[2].strip()))))
        # listAngelDiffs.append(math.fabs(float(listCurPara[3].strip())))
        # listAngelDiffs.append(math.fabs(float(listCurPara[4].strip())))
        # listAngelDiffs.append(math.fabs(float(listCurPara[5].strip())))
        listAngelDiffs.sort()

        # Normalize
        # if ((listAngelDiffs[2] - listAngelDiffs[0]) < 1e-3):
        #     listAngelDiffs[0] = 0.0
        #     listAngelDiffs[1] = 0.0
        #     listAngelDiffs[2] = 0.0
        # else:
        #     listAngelDiffs[1] = (listAngelDiffs[1] - listAngelDiffs[0])/(listAngelDiffs[2] - listAngelDiffs[0])
        #     listAngelDiffs[2] = 1.0
        #     listAngelDiffs[0] = 0.0

        listCurPara = []
        listCurPara.append(strFileName)
        listCurPara.append(strSpaceID)
        listCurPara.append(str(listLengthDiffs[0]))
        listCurPara.append(str(listLengthDiffs[2]))
        listCurPara.append(str(listAngelDiffs[0]))
        listCurPara.append(str(listAngelDiffs[2]))
        listCurPara.append(strVolume)

        # dMeanDiff = (listLengthDiffs[0]+listLengthDiffs[1]+listLengthDiffs[2]) / 3
        # listCurPara.append(str(dMeanDiff))
        # dMeanDiff = (listAngelDiffs[0]+listAngelDiffs[1]+listAngelDiffs[2]) / 3
        # listCurPara.append(str(dMeanDiff))

        # dSisso1 = listLengthDiffs[0]*listLengthDiffs[0]*listLengthDiffs[0] - listAngelDiffs[2]
        # listCurPara.append(str(dSisso1))

        # Add variance
        # listCurPara.append(str(np.var(listLengthDiffs)))
        # listCurPara.append(str(np.var(listAngelDiffs)))

        listCifParas.append(listCurPara)

    # Write to output file
    strCrystalFile = objects.BASISDIR + "cif-out.csv"
    file = open(strCrystalFile,"a+")
    file.write("\n--------Crystal Parameters---------------\n")
    for item in listCifParas:
        file.write(",".join(item))
        file.write("\n")

    file.close()


# Read the Wyckoff matrices of all crystals from cif files and save to a file
# listCifFiles is used to filter crystals; when it's empty, all crystals in the folder are read
def readCrystalWyckoffFromCif(listCifFiles,listElementInfo):

    # Obtain encoding based on the system
    strCodeType = "utf-8"
    sysType = sys.platform
    if sysType == "Windows":
        strCodeType = "utf-8"
    elif sysType == "Linux":
        strCodeType = "gbk"

    listElements = ["H","Li","Be","B","C","N","O","F","Na","Mg","Al","Si","P","S","Cl","K","Ca","Sc",
                    "Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Rb","Sr","Y",
                    "Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Cs","Ba","La","Lu",
                    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi"]

    dictEleIndex = {}
    for eleItem in listElementInfo:
        dictEleIndex[eleItem.m_strName] = int(eleItem.m_strID)

    # Read all possible Wyckoff positions from the file
    dictWyckoff = {}
    listWyckoffVector = []
    strWychkoffFile = objects.BASISDIR + "public/WyckoffPosition.csv"
    WychkoffFile = open(strWychkoffFile,"r",encoding=strCodeType)
    while 1:
        strLine = WychkoffFile.readline()
        if not strLine:
            break
        strLine = strLine.strip()
        arrTmpPositions = strLine.split(",")
        dictWyckoff[arrTmpPositions[0]] = 0
        if arrTmpPositions[0] not in listWyckoffVector:
            listWyckoffVector.append(arrTmpPositions[0])
    WychkoffFile.close()

    # Traverse all files in the folder
    strCifPath = objects.BASISDIR + "cif/"
    if len(listCifFiles) < 1:
        listCifFiles = os.listdir(strCifPath)

    listResult = []
    for curFile in listCifFiles:
        nMultiplicityIndex = -1
        nLineNum = -1
        dictCurWychoff = dictWyckoff.copy()
        dictCurCount = {}
        curFile = strCifPath + curFile
        # Use utf-8 to avoid 'gbk' codec can't decode byte 0xa9 error
        # Alternatively, convert the file to ANSI beforehand; Linux defaults to utf-8, while Windows defaults to gbk
        file = open(curFile,"r",encoding=strCodeType)
        while 1:
            strLine = file.readline()
            if not strLine:
                break

            # Handle cases where multiple spaces are used as separators
            strLine = ' '.join(strLine.split())
            if '_atom_site_label' in strLine:
                nLineNum = 0
                continue
            if nLineNum > -1 and '_atom_' in strLine and '_atom_site_symmetry_multiplicity' not in strLine:
                nLineNum += 1
                continue
            if nLineNum > -1 and '_atom_site_symmetry_multiplicity' in strLine:
                nLineNum += 1
                nMultiplicityIndex = nLineNum
                continue

            # Reading started
            if nMultiplicityIndex > -1 and '_atom_' not in strLine:
                arrAtomInfos = strLine.split()
                # Reading ended upon encountering an empty line or 'loop'
                if (len(arrAtomInfos) < nMultiplicityIndex+2):
                    break
                strAtomType = arrAtomInfos[0].strip('0123456789')
                # strResultIndex = arrAtomInfos[0].strip('0123456789')
                # strResultIndex += "-"
                # strResultIndex += arrAtomInfos[nMultiplicityIndex+1]
                strResultIndex = arrAtomInfos[nMultiplicityIndex]
                dictCurWychoff[strResultIndex] += dictEleIndex[strAtomType]#arrAtomInfos[nMultiplicityIndex]
                if strResultIndex in dictCurCount.keys():
                    dictCurCount[strResultIndex] += 1
                else:
                    dictCurCount[strResultIndex] = 1
        file.close()

        # Retrieve proportions
        # for key in dictCurWychoff.keys():
        #     dictCurWychoff[key] = dictCurWychoff[key] / nCurTotalAtomNum
        for key in dictCurWychoff.keys():
            if key in dictCurCount.keys():
                dictCurWychoff[key] = dictCurWychoff[key] / dictCurCount[key]

        listResult.append(dictCurWychoff)

    strOutFile = objects.BASISDIR + "Wyckoff-out.csv"
    outfile = open(strOutFile,"a+")
    outfile.write("ICSD")
    for item in listWyckoffVector:
            outfile.write(",")
            outfile.write(item)

    outfile.write("\n")
    outfile.close()

    # Write to output file
    outfile = open(strOutFile,"a+")
    for nIndex in range(0,len(listCifFiles)):
        outfile.write(listCifFiles[nIndex])
        # Check if this column is zero; write if not zero
        for item in listWyckoffVector:
                outfile.write(",")
                outfile.write(str(listResult[nIndex][item]))
        outfile.write("\n")

    outfile.close()

# Extract bond length distribution, flexibility of negative groups in the crystal, and magnitude of bond length variation from CIF files.
def calBondAndFlexibility(strCifFile, listElementInfo, crystalObject):

    # Retrieve all atoms in the negative groups within the crystal, excluding positive groups.
    listNegaAtoms = []
    for item in crystalObject.m_listNegativeGroups:
        strFormula = item
        dictElements = {}
        method.getElementFromFormula(strFormula,dictElements)
        for key in dictElements.keys():
            if key not in listNegaAtoms:
                listNegaAtoms.append(key)

    # Use utf-8 to avoid 'gbk' codec can't decode byte 0xa9 error
    # Alternatively, convert the file to ANSI beforehand
    listParas = [999.0, 999.0, 999.0, 999.0, 999.0, 999.0]
    bIsReadCoord = False
    nAtomNameCol = -1
    nAtomXCol = -1
    nAtomYCol = -1
    nAtomZCol = -1
    listAtoms = []
    nLineIndex = -1
    file = open(strCifFile,"r",encoding='utf-8')
    while 1:
        strLine = file.readline()
        if not strLine:
            break

        # Handle cases where multiple spaces are used as separators
        strLine = ' '.join(strLine.split())

        # Read unit cell parameters
        if '_cell_length_' in strLine or '_cell_angle_' in strLine:
            strLine = strLine.strip()
            strLine = strLine.strip(")")
            strLine = strLine.replace("(", "")
            listTmpInfo = strLine.split(" ")
            listParas.append(listTmpInfo[1])

        # Read all negative groups within the unit cell
        strLine.strip()
        if "_atom_site" in strLine:

            if not bIsReadCoord:
                bIsReadCoord = True

            nLineIndex =  nLineIndex + 1

            if strLine == "_atom_site_label":
                nAtomNameCol = nLineIndex
            elif strLine == "_atom_site_fract_x":
                nAtomXCol = nLineIndex
            elif strLine == "_atom_site_fract_y":
                nAtomYCol = nLineIndex
            elif strLine == "_atom_site_fract_z":
                nAtomZCol = nLineIndex

            continue

        if bIsReadCoord:
            strLine = strLine.replace(")", "")
            strLine = strLine.replace("(", "")
            listTmpInfo = strLine.split(" ")
            if len(listTmpInfo) < 4 or "#" in strLine:
                break

            strCurAtomName = listTmpInfo[nAtomNameCol].strip("0123456789")
            if  strCurAtomName in listNegaAtoms:
                curAtom = objects.CAtom(listTmpInfo[nAtomNameCol]) # Atomic name
                curAtom.m_strName = strCurAtomName
                curAtom.m_strX = listTmpInfo[nAtomXCol] # Atomic type
                curAtom.m_strY = listTmpInfo[nAtomYCol] # Atomic type
                curAtom.m_strZ = listTmpInfo[nAtomZCol] # Atomic type
                if curAtom.m_strX == ".":
                    curAtom.m_strX = "0.0"
                if curAtom.m_strY == ".":
                    curAtom.m_strY = "0.0"
                if curAtom.m_strZ == ".":
                    curAtom.m_strZ = "0.0"
                listAtoms.append(curAtom)

    file.close()

    # Crystal constants
    listLengthParas = []
    listLengthParas.append(float(listParas[0]))
    listLengthParas.append(float(listParas[1]))
    listLengthParas.append(float(listParas[2]))
    listAngelParas = []
    listAngelParas.append(float(listParas[3]))
    listAngelParas.append(float(listParas[4]))
    listAngelParas.append(float(listParas[5]))

    # Calculate SymmFunc for each atom
    dCutOffRadiu = 6.0
    listNearAtoms = []
    listNearDist = []
    nAtomNum = len(listAtoms)
    for nIndex in range(0,nAtomNum):
        listNearAtoms.clear()
        for nSubIndex in range(0,nAtomNum):
            if nIndex == nSubIndex:
                continue
            dDistance = method.calDistance(listAtoms[nIndex], listAtoms[nSubIndex],listLengthParas, listAngelParas)
            if dDistance < dCutOffRadiu:
                listNearAtoms.append(listAtoms[nSubIndex])
                listNearDist.append(dDistance)

        listAtoms[nIndex].calSymmetryFunction(listNearAtoms,listNearDist, dCutOffRadiu)

    crystalObject.m_listNegaAtoms = listAtoms

    # Get Flexibility
    listBonds = method.getBondsFromAtoms(listAtoms, listElementInfo, listLengthParas, listAngelParas)
    crystalObject.m_listNegaBonds = listBonds

    return  crystalObject


# Retrieve CrystalInfo directly from the complete text
def readCrystalInfoFromTxt(strFile):
    listCrystalInfo = []
    file = open(strFile,"r")
    while 1:
        strLine = file.readline()
        if not strLine:
            break

        strLine = strLine.strip()
        listInfos = strLine.split(",")
        if len(listInfos) < 19:
            continue

        curCrystal = objects.CCrystalObject(listInfos[0].strip())
        curCrystal.m_strName = listInfos[0].strip()
        curCrystal.m_strSpaceGroupName = listInfos[1].strip()
        curCrystal.m_strSpaceGroupID = listInfos[2].strip()
        curCrystal.m_strLengthMinDiff = listInfos[3].strip()
        curCrystal.m_strLengthMaxDiff = listInfos[4].strip()
        curCrystal.m_strAngleMinDiff = listInfos[5].strip()
        curCrystal.m_strAngleMaxDiff = listInfos[6].strip()
        curCrystal.m_strVolume = listInfos[7].strip()
        curCrystal.m_listPostiveGroups = listInfos[8].strip().split(";")
        curCrystal.m_listNegativeGroups = listInfos[9].strip().split(";")
        curCrystal.setFormula(listInfos[10].strip())

        curCrystal.m_strGapLevel = listInfos[11].strip()
        curCrystal.m_strBandGap = listInfos[12].strip()

        curCrystal.m_strBiRefLevel = listInfos[13].strip()
        curCrystal.m_listBeReflength = listInfos[14].strip().split(";")
        curCrystal.m_listBeRef = listInfos[15].strip().split(";")

        curCrystal.m_strMaxDijLevel = listInfos[16].strip()
        curCrystal.m_lisMaxDijlength = listInfos[17].strip().split(";")
        curCrystal.m_listMaxDij = listInfos[18].strip().split(";")

        listCrystalInfo.append(curCrystal)

    file.close()

    return listCrystalInfo


# Save the distribution of elements in the crystal
def saveElementDistribution(listCrystalInfo,listElementInfo):


    listTargetEle = ["ICSD","Formula",'H', 'Li', 'Be', 'B', 'C', 'N', 'O','F',
                     'Na','Mg','Al','Si','P','S','Cl',
                     'K','Ca','Sc','Ti','V','Cr', 'Mn','Co','Ni','Cu', 'Zn','Ga','Ge','As', 'Se','Br',
                     'Rb', 'Sr','Y','Zr','Nb', 'Mo','Ru','Pd','Ag','Cd','In', 'Sn','Sb','Te','I',
                     'Cs','Ba','La', 'Hf','Ta','W','Os','Pt','Au','Hg','Tl','Pb','Bi']

    # Write to output file
    strCrystalFile = objects.BASISDIR + "elementDistribution-out.csv"
    file = open(strCrystalFile,"a+")
    file.write("\n--------elementDistribution---------------\n")
    file.write("\n")
    file.write(",".join(listTargetEle))

    for item in listCrystalInfo:
        listTmpNum = []
        listTmpNum.append(item.m_strID)
        listTmpNum.append(item.getFormula())
        for nIndex in range(2,len(listTargetEle)):
            if listTargetEle[nIndex] in item.m_dictElements.keys():
                #nNum = item.m_dictElements[ listTargetEle[nIndex] ]
                listTmpNum.append(str(1))
            else:
                listTmpNum.append("0")

        file.write("\n")
        file.write(",".join(listTmpNum))
        print(listTmpNum)

    file.close()


# Read all Crystal-related information from a txt file
def readAllCrystalInfo(listElementInfo):
    # Select the directory where the txt file is located
    strCrystalFile = objects.BASISDIR + "in_and_out/crystal-in.csv"
    listCrystalInfo = readCrystalInfoFromTxt(strCrystalFile)

    listCenterSpaceID = [2,10,11,12,13,14,15,47,48,49,50,51,52,53,54,55,56,57,58,59,60,
                         61,62,63,64,65,66,67,68,69,70,71,72,73,74,83,84,85,86,87,88,123,
                         124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,
                         147,148,162,163,164,165,166,167,175,176,191,192,193,194,200,201,202,203,204,
                         205,206,221,222,223,224,225,226,227,228,229,230]

    # Read all SpaceGroup information. If SpaceGroupID is not present in the original file, read it from the respective file
    # nIndex = -1
    # dictSpaceGroupInfo = readAllSpaceGroupInfo()
    # for curCrystal in listCrystalInfo:
    #     nIndex = nIndex + 1
    #     strCifFile = objects.BASISDIR + "cif/" +curCrystal.m_strID+".cif"
    #     if not os.path.exists(strCifFile):
    #         print("The %s.cif is not exist\n"%(curCrystal.m_strID))
    #         continue
        # Calculate Bond and Flexibility ****
        #listCrystalInfo[nIndex] = calBondAndFlexibility(strCifFile, listElementInfo, curCrystal)

        # Fill in the SpaceGroup ID
        # strCurSpaceGroupName = curCrystal.m_strSpaceGroupName.strip("Z")
        # strCurSpaceGroupName = strCurSpaceGroupName.strip("S")
        # if len(curCrystal.m_strSpaceGroupID) < 1 and strCurSpaceGroupName in dictSpaceGroupInfo.keys():
        #     listCrystalInfo[nIndex].m_strSpaceGroupID = dictSpaceGroupInfo[strCurSpaceGroupName]

        # Mark centrosymmetric space groups
        # if listCrystalInfo[nIndex].m_strSpaceGroupID == "":
        #     listCrystalInfo[nIndex].m_strSpaceGroupID = "0"
        # elif int(listCrystalInfo[nIndex].m_strSpaceGroupID) in listCenterSpaceID:
        #     listCrystalInfo[nIndex].m_strSpaceGroupID = "0"

    # Retrieve crystal parameters for all crystal CIFs
    # listCrystalCifs = []
    # for curCrystal in listCrystalInfo:
    #     listCrystalCifs.append(curCrystal.m_strID+".cif")
    # readCrystalParaFromCif(listCrystalCifs)
    #readCrystalWyckoffFromCif(listCrystalCifs,listElementInfo)

    # Save the distribution of elements in the crystal
    # saveElementDistribution(listCrystalInfo,listElementInfo)

    return listCrystalInfo

# Save all information
def saveAllCrystalInfo(listCrystalInfo):
    strCrystalFile = objects.BASISDIR + "in_and_out/crystal-out.csv"
    file = open(strCrystalFile,"a+")
    file.write("\n--------Crystals---------------\n")
    for curCrystal in listCrystalInfo:
        file.write(curCrystal.joinToString(","))
        file.write("\n")
    file.close()


if __name__ == '__main__':

    curCrystal = objects.CCrystalObject("111")
    curCrystal.setFormula("Ba2NO3(OH)3")
    print (curCrystal.m_dictElements)

