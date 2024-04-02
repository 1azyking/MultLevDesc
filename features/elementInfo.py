# Retrieve all element information from an external file

import os
import sys
sys.path.append("..") # Import from parent directory
import objects
import method

# Directly extract ElementInfo from the complete text
def readElementInfoFromTxt(strFile):
    listElementInfo = []
    file = open(strFile,"r")
    while 1:
        strLine = file.readline()
        if not strLine:
            break

        strLine = strLine.strip();
        listInfos = strLine.split(",")
        if len(listInfos) < 20:
            continue

        curElement = objects.CElementObject(listInfos[0]); # By atomic number
        curElement.m_strName = listInfos[1] # By element symbol
        curElement.m_strMass = listInfos[2]
        curElement.m_strPeriodNum = listInfos[3]
        curElement.m_strGroupNum = listInfos[4]
        curElement.m_strRadui = listInfos[5]
        curElement.m_strEleNegativity = listInfos[6]
        curElement.m_strSValenceEleNum = listInfos[7]
        curElement.m_strPValenceEleNum = listInfos[8]
        curElement.m_strDValenceEleNum = listInfos[9]
        curElement.m_strFValenceEleNum = listInfos[10]
        curElement.m_strSUnfilledState = listInfos[11]
        curElement.m_strPUnfilledState = listInfos[12]
        curElement.m_strDUnfilledState = listInfos[13]
        curElement.m_strFUnfilledState = listInfos[14]
        curElement.m_strionizationEner = listInfos[15]
        curElement.m_strEleAffinity = listInfos[16]
        curElement.m_strMeltingPoint = listInfos[17]
        curElement.m_strBoilingPoint = listInfos[18]
        curElement.m_strDensity = listInfos[19]
        listElementInfo.append(curElement)

    file.close()

    return  listElementInfo


# Read all Element-related information from a txt file
def readAllElementInfo():
    # Select the directory where the txt file is located
    strElementFile = objects.BASISDIR + "public/element-in.csv"
    listElementInfo = readElementInfoFromTxt(strElementFile)
    return  listElementInfo

# Save all information
def saveAllElementInfo(listElementInfo):
    strElementFile = objects.BASISDIR + "element-out.csv"
    file = open(strElementFile,"a+")
    file.write("\n--------Elements---------------\n")
    for curElement in listElementInfo:
        file.write(curElement.joinToString(","))
        file.write("\n")
    file.close()


