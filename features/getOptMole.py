import os
import sys
import method
import elementInfo
import shutil
import numpy as np

# Obtain optimized structures directly from Gaussian output files
def getOptMoleFromOutFiles():

    listElementInfo = elementInfo.readAllElementInfo()

    strTZVPD = "@/public4/home/sc55809/mySoft/basisset/def2-TZVPD/X.gbs"

    #Select the directory where the .out files are located
    strOutPath = "PATH/TO/YOUR/OUT/FILES"
    listfile = os.listdir(strOutPath)
    for curFile in listfile:
        curFile = strOutPath + curFile
        if curFile.endswith("gjf"):
            continue

        outFile = open(curFile,"r")
        strgjfFile = curFile.replace("out","gjf")
        gjfFile = open(strgjfFile,"r")
        gjfFile.seek(0,0)
        strnewFile = curFile.replace("out","new")
        newFile = open(strnewFile,"a+")

        #Write the title of the gjf file to new, while storing all atom names
        ntitleNum = 7
        listAtoms = []
        nIndex = 0
        while 1:
            strLine = gjfFile.readline()
            if not strLine:
                break

            strLine = strLine.strip()
            nIndex = nIndex + 1
            if 'chk' in strLine:
                ntitleNum = ntitleNum + 1

            if nIndex <= ntitleNum:
                newFile.write(strLine)
                newFile.write("\n")
                continue

            if nIndex > ntitleNum:
                if len(strLine) > 0:
                    listCoord = strLine.split()
                    listAtoms.append(listCoord[0])
                else:
                    break

        nIndex = 0
        nStartLine = 0
        outFile.seek(0,0)
        while 1:
            strLine = outFile.readline()
            if not strLine:
                break

            nIndex = nIndex + 1
            if 'Input orientation:' in strLine:
                nStartLine = nIndex

        nIndex = 0
        nStartLine = nStartLine + 5
        nEndLine = nStartLine+len(listAtoms)
        listAtomCoords = []
        outFile.seek(0,0)
        while 1:
            strLine = outFile.readline()
            if not strLine:
                break

            nIndex = nIndex + 1
            if nIndex >= nStartLine and nIndex < nEndLine:
                listInfos = strLine.split()
                strCoord = listInfos[3] + "    " + listInfos[4] + "    " + listInfos[5]
                listAtomCoords.append(strCoord)

        nIndex = 0
        for curAtom in listAtoms:
            strLine = curAtom + "               "
            strLine = strLine + listAtomCoords[nIndex]
            nIndex = nIndex + 1
            newFile.write(strLine)
            newFile.write("\n")

        nIndex = 0
        gjfFile.seek(0,0)
        nStartLine = ntitleNum+len(listAtoms)
        while 1:
            strLine = gjfFile.readline()
            if not strLine:
                break

            if nIndex > nStartLine:
                newFile.write(strLine)

        newFile.write("\n")

        newFile.write("532nm 1064nm 10600nm\n")
        newFile.write("\n")
        dictElements = {}
        strFormula = curFile.strip(".out").strip("0123456789")
        method.getElementFromFormula(strFormula,dictElements)
        for curElement in listElementInfo:
            if curElement.m_strName in dictElements.keys():
                curTZVPD = strTZVPD.replace("X", curElement.m_strName)
                newFile.write(curTZVPD)
                newFile.write("\n")
        newFile.write("\n")

        outFile.close()
        gjfFile.close()
        newFile.close()

#Reading atomic coordinates from QM input file of qm4d
def getAtomPositions():

    strOutPath = "PATH/TO/YOUR/OUT/FILES"
    strOutFile = strOutPath + "result.csv"
    outFile = open(strOutFile,"a+")

    for curStep in range(0,501,5):
        curFile = strOutPath + "testeth.com_bak" + str(curStep)
        gjfFile = open(curFile,"r")
        gjfFile.seek(0,0)

        nIndex = 0
        while 1:
            strLine = gjfFile.readline()
            if not strLine:
                break

            nIndex = nIndex + 1
            if nIndex == 8:
                strLine = strLine.strip()
                strLine = strLine.replace(" ", ",")
                outFile.write(strLine)
                outFile.write("\n")
                break

        gjfFile.close()

    outFile.close()

#Rename
def renameTogjf():
    strOutPath = "PATH/TO/YOUR/OUT/FILES"
    listfile = os.listdir(strOutPath)
    for curFile in listfile:
        curFile = strOutPath + curFile
        strgjfFile = curFile.replace(".new",".gjf")
        os.rename(curFile,strgjfFile)

def renameTo():
    strOutPath = "PATH/TO/YOUR/OUT/FILES"
    listfile = os.listdir(strOutPath)
    for curFile in listfile:
        curFile = strOutPath + curFile
        strgjfFile = curFile.replace(".new",".gjf")
        os.rename(curFile,strgjfFile)

def copyFiles():
    strInPath = "PATH/TO/YOUR/CIF/FILES"
    strOutPath = "PATH/TO/NEW/FILES"

    listFileNames = []
    strFilePath = "PATH/TO/needCopy.txt"
    file = open(strFilePath,"r")
    while 1:
        strLine = file.readline()
        if not strLine:
            break

        listFileNames.append(strLine.strip())

    for curFile in listFileNames:
        strInFile = strInPath + curFile + ".cif"
        strOutFile = strOutPath + curFile + ".cif"
        shutil.copyfile(strInFile, strOutFile)


def getInformation():

    strPath = "PATH/TO/YOUR/OUT/FILES"
    listCifFiles = os.listdir(strPath)

    #Obtain related coding method based on the system
    strCodeType = "utf-8"
    sysType = sys.platform
    if sysType == "Windows":
        strCodeType = "utf-8"
    elif sysType == "Linux":
        strCodeType = "gbk"

    listInfos = []
    for curFile in listCifFiles:
        listCurInfo = []
        strFileName = curFile.strip(".cif")
        listCurInfo.append(strFileName)
        curFile = strPath + curFile
        # Use utf-8, otherwise it will report error 'gbk' codec can't decode byte 0xa9 in
        # You can also convert the file to ANSI beforehand, linux defaults to utf-8, windows defaults to gbk.
        file = open(curFile,"r",encoding=strCodeType)
        while 1:
            strLine = file.readline()
            if not strLine:
                break

            strLine = ' '.join(strLine.split())
            if '_space_group_IT_number' in strLine:
                arrTmpInfo = strLine.split()
                listCurInfo.append(arrTmpInfo[1])
                listInfos.append(listCurInfo)
                break

        file.close()
        if len(listCurInfo) < 2:
            listCurInfo.append("9999")
            listInfos.append(listCurInfo)


    #Write the output file(new.csv)
    strOutFile = "PATH/TO/SAVE/new.csv"
    file = open(strOutFile,"a+")
    for item in listInfos:
        file.write(",".join(item))
        file.write("\n")

    file.close()


def FormatOutput():

    strInputPsfFile = "PATH/TO/11.csv"
    strOutPsfFile = "PATH/TO/11.out"

    listAllItems = []
    inFile = open(strInputPsfFile,"r")
    while 1:
        strLine = inFile.readline()
        if not strLine:
            break

        strLine = strLine.strip()
        strLine = strLine.strip("\n")
        listItems = strLine.split(",")
        for nIndex in range(0, len(listItems)):
            listItems[nIndex] = listItems[nIndex].strip()
            if (listItems[nIndex] != ""):
                listAllItems.append(listItems[nIndex])

    #output
    nItemNumInLine = 9
    outFile = open(strOutPsfFile, "w+")
    for nIndex in range(0, len(listAllItems)):
        if (nIndex % nItemNumInLine) == 0:
            outFile.write("\n")

        outFile.write("%9s"%(listAllItems[nIndex]))

    inFile.close()
    outFile.close()


if __name__ == '__main__':
    #getOptMoleFromOutFiles()
    renameTo()
    #getAtomPositions()
    #copyFiles()
    #getInformation()












