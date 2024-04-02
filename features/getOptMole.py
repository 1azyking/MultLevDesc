import os
import sys
import method
import elementInfo
import numpy as np

#从Gaussian输出文件中直接获取优化后的结构
def getOptMoleFromOutFiles():

    #获取所有元素信息
    listElementInfo = elementInfo.readAllElementInfo()

    strTZVPD = "@/public4/home/sc55809/mySoft/basisset/def2-TZVPD/X.gbs"

    #选择out文件所在目录
    strOutPath = "D:/Data/"
    #遍历文件夹中所有文件
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

        #将gjf文件的title写入new，同时存储所有原子名称
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

        #首先从out中获取最后一个Input orientation:
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

        #向new文件中写入原子和坐标
        nIndex = 0
        for curAtom in listAtoms:
            strLine = curAtom + "               "
            strLine = strLine + listAtomCoords[nIndex]
            nIndex = nIndex + 1
            newFile.write(strLine)
            newFile.write("\n")

        #写入gjf文件中剩余部分
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

        #添加Polar计算的部分

        #获取所有的元素
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

#从qm4d的QM输入文件中读取原子坐标
def getAtomPositions():

    #选择out文件所在目录
    strOutPath = "D:/New-cation/"
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

#改名
def renameTogjf():
    #选择out文件所在目录
    strOutPath = "D:/New-cation/"
    #遍历文件夹中所有文件
    listfile = os.listdir(strOutPath)
    for curFile in listfile:
        curFile = strOutPath + curFile
        strgjfFile = curFile.replace(".new",".gjf")
        os.rename(curFile,strgjfFile)

#改名
def renameTo():
    #选择out文件所在目录
    strOutPath = "D:/Data/"
    #遍历文件夹中所有文件
    listfile = os.listdir(strOutPath)
    for curFile in listfile:
        curFile = strOutPath + curFile
        strgjfFile = curFile.replace(".new",".gjf")
        os.rename(curFile,strgjfFile)

#根据txt复制文件
import shutil
def copyFiles():
    #选择out文件所在目录
    strInPath = "D:/Data/cif/"
    strOutPath = "D:/Data/new/"

    listFileNames = []
    strFilePath = "D:/Data/needCopy.txt"
    file = open(strFilePath,"r")
    while 1:
        strLine = file.readline()
        if not strLine:
            break

        listFileNames.append(strLine.strip())

    #遍历文件夹中所有文件
    for curFile in listFileNames:
        strInFile = strInPath + curFile + ".cif"
        strOutFile = strOutPath + curFile + ".cif"
        shutil.copyfile(strInFile, strOutFile)


#从文本中读取一些信息
def getInformation():

    #选择out文件所在目录
    strPath = "D:/Data/new/"
    listCifFiles = os.listdir(strPath)

    #根据系统获得编码方式
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
        # 使用utf-8，否则会报错'gbk' codec can't decode byte 0xa9 in
        # 也可以事先将文件转为ANSI,linux默认是utf-8,windows默认gbk
        file = open(curFile,"r",encoding=strCodeType)
        while 1:
            strLine = file.readline()
            if not strLine:
                break

            #处理多个空格分隔的情况
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


    #写入输出文件
    strOutFile = "D:/Data/new.csv"
    file = open(strOutFile,"a+")
    for item in listInfos:
        file.write(",".join(item))
        file.write("\n")

    file.close()


#格式化输，从csv读入二维数组。按行优先形式每隔多少个输出一行
def FormatOutput():

    strInputPsfFile = "D:/11.csv"
    strOutPsfFile = "D:/11.out"

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












