#从外部文件中获得所有晶体信息
import math
import os
import sys
import objects
import method
import numpy as np

#从txt文件中读取所有SpaceGroup相关信息
#注意有些空间群符号表达方式不同，这些暂时没有处理
def readAllSpaceGroupInfo():

    dictSpaceGroupInfo = {}

    #选择txt文件所在目录
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
            #注意,name作为索引
            if listdetails[1] not in dictSpaceGroupInfo.keys():
                dictSpaceGroupInfo[listdetails[1]] = listdetails[0]

    file.close()

    return dictSpaceGroupInfo


#从cif文件中读取所有晶体的晶胞参数，并保存到文件
#listCifFiles用于过滤晶体，如果为空时，读取文件夹下所有晶体
def readCrystalParaFromCif(listCifFiles):

    strCifPath = objects.BASISDIR + "cif/"
    listCifParas = []

    #遍历文件夹中所有文件
    if len(listCifFiles) < 1:
        listCifFiles = os.listdir(strCifPath)

    #根据系统获得编码方式
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
        # 使用utf-8，否则会报错'gbk' codec can't decode byte 0xa9 in
        # 也可以事先将文件转为ANSI,linux默认是utf-8,windows默认gbk
        file = open(curFile,"r",encoding=strCodeType)
        while 1:
            strLine = file.readline()
            if not strLine:
                break

            #处理多个空格分隔的情况
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

        #归一化
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

        #归一化
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

        #添加方差
        #listCurPara.append(str(np.var(listLengthDiffs)))
        #listCurPara.append(str(np.var(listAngelDiffs)))

        listCifParas.append(listCurPara)

    #写入输出文件
    strCrystalFile = objects.BASISDIR + "cif-out.csv"
    file = open(strCrystalFile,"a+")
    file.write("\n--------Crystal Parameters---------------\n")
    for item in listCifParas:
        file.write(",".join(item))
        file.write("\n")

    file.close()


#从cif文件中读取所有晶体的Wyckoff矩阵，并保存到文件
#listCifFiles用于过滤晶体，如果为空时，读取文件夹下所有晶体
def readCrystalWyckoffFromCif(listCifFiles,listElementInfo):

    #根据系统获得编码方式
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

    #从文件中读取所有可能的Wychkoff position
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

    #遍历文件夹中所有文件
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
        # 使用utf-8，否则会报错'gbk' codec can't decode byte 0xa9 in
        # 也可以事先将文件转为ANSI,linux默认是utf-8,windows默认gbk
        file = open(curFile,"r",encoding=strCodeType)
        while 1:
            strLine = file.readline()
            if not strLine:
                break

            #处理多个空格分隔的情况
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

            #开始读取了
            if nMultiplicityIndex > -1 and '_atom_' not in strLine:
                arrAtomInfos = strLine.split()
                #读取结束时遇到空行或loop
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

        #取比例
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

    #写入输出文件
    outfile = open(strOutFile,"a+")
    for nIndex in range(0,len(listCifFiles)):
        outfile.write(listCifFiles[nIndex])
        #检查此列是否为0,不为0时写入
        for item in listWyckoffVector:
                outfile.write(",")
                outfile.write(str(listResult[nIndex][item]))
        outfile.write("\n")

    outfile.close()

#从cif文件中获取键长分布，晶体中阴性基团的Flexibility，键长的变化幅度
def calBondAndFlexibility(strCifFile, listElementInfo, crystalObject):

    #获得晶体中阴性基团中的所有原子，也就是此时不考虑阳性基团
    listNegaAtoms = []
    for item in crystalObject.m_listNegativeGroups:
        strFormula = item
        dictElements = {}
        method.getElementFromFormula(strFormula,dictElements)
        for key in dictElements.keys():
            if key not in listNegaAtoms:
                listNegaAtoms.append(key)

    # 使用utf-8，否则会报错'gbk' codec can't decode byte 0xa9 in
    # 也可以事先将文件转为ANSI
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

        #处理多个空格分隔的情况
        strLine = ' '.join(strLine.split())

        #读取晶胞参数
        if '_cell_length_' in strLine or '_cell_angle_' in strLine:
            strLine = strLine.strip()
            strLine = strLine.strip(")")
            strLine = strLine.replace("(", "")
            listTmpInfo = strLine.split(" ")
            listParas.append(listTmpInfo[1])

        #读取晶胞内所有阴性基团的
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
                curAtom = objects.CAtom(listTmpInfo[nAtomNameCol]) #原子名称
                curAtom.m_strName = strCurAtomName
                curAtom.m_strX = listTmpInfo[nAtomXCol] #原子类型
                curAtom.m_strY = listTmpInfo[nAtomYCol] #原子类型
                curAtom.m_strZ = listTmpInfo[nAtomZCol] #原子类型
                if curAtom.m_strX == ".":
                    curAtom.m_strX = "0.0"
                if curAtom.m_strY == ".":
                    curAtom.m_strY = "0.0"
                if curAtom.m_strZ == ".":
                    curAtom.m_strZ = "0.0"
                listAtoms.append(curAtom)

    file.close()

    #晶体常数
    listLengthParas = []
    listLengthParas.append(float(listParas[0]))
    listLengthParas.append(float(listParas[1]))
    listLengthParas.append(float(listParas[2]))
    listAngelParas = []
    listAngelParas.append(float(listParas[3]))
    listAngelParas.append(float(listParas[4]))
    listAngelParas.append(float(listParas[5]))

    #计算每个原子的SymmFunc
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

    #获取Flexibility
    listBonds = method.getBondsFromAtoms(listAtoms, listElementInfo, listLengthParas, listAngelParas)
    crystalObject.m_listNegaBonds = listBonds

    return  crystalObject


#从完整的文本中直接获取CrystalInfo
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


#保存晶体中元素的分布情况
def saveElementDistribution(listCrystalInfo,listElementInfo):


    listTargetEle = ["ICSD","Formula",'H', 'Li', 'Be', 'B', 'C', 'N', 'O','F',
                     'Na','Mg','Al','Si','P','S','Cl',
                     'K','Ca','Sc','Ti','V','Cr', 'Mn','Co','Ni','Cu', 'Zn','Ga','Ge','As', 'Se','Br',
                     'Rb', 'Sr','Y','Zr','Nb', 'Mo','Ru','Pd','Ag','Cd','In', 'Sn','Sb','Te','I',
                     'Cs','Ba','La', 'Hf','Ta','W','Os','Pt','Au','Hg','Tl','Pb','Bi']

    #写入输出文件
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


#从txt文件中读取所有Crystal相关信息
def readAllCrystalInfo(listElementInfo):
    #选择txt文件所在目录
    strCrystalFile = objects.BASISDIR + "crystal-in.csv"
    listCrystalInfo = readCrystalInfoFromTxt(strCrystalFile)

    listCenterSpaceID = [2,10,11,12,13,14,15,47,48,49,50,51,52,53,54,55,56,57,58,59,60,
                         61,62,63,64,65,66,67,68,69,70,71,72,73,74,83,84,85,86,87,88,123,
                         124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,
                         147,148,162,163,164,165,166,167,175,176,191,192,193,194,200,201,202,203,204,
                         205,206,221,222,223,224,225,226,227,228,229,230]

    #读入所有SpaceGroup，如果原文件中没有SpaceGroupID，则从相应文件中读取
    # nIndex = -1
    # dictSpaceGroupInfo = readAllSpaceGroupInfo()
    # for curCrystal in listCrystalInfo:
    #     nIndex = nIndex + 1
    #     strCifFile = objects.BASISDIR + "cif/" +curCrystal.m_strID+".cif"
    #     if not os.path.exists(strCifFile):
    #         print("The %s.cif is not exist\n"%(curCrystal.m_strID))
    #         continue
        #计算Bond和Flexibility ****
        #listCrystalInfo[nIndex] = calBondAndFlexibility(strCifFile, listElementInfo, curCrystal)

        # 填写空间群ID
        # strCurSpaceGroupName = curCrystal.m_strSpaceGroupName.strip("Z")
        # strCurSpaceGroupName = strCurSpaceGroupName.strip("S")
        # if len(curCrystal.m_strSpaceGroupID) < 1 and strCurSpaceGroupName in dictSpaceGroupInfo.keys():
        #     listCrystalInfo[nIndex].m_strSpaceGroupID = dictSpaceGroupInfo[strCurSpaceGroupName]

        #标记中心对称空间群
        # if listCrystalInfo[nIndex].m_strSpaceGroupID == "":
        #     listCrystalInfo[nIndex].m_strSpaceGroupID = "0"
        # elif int(listCrystalInfo[nIndex].m_strSpaceGroupID) in listCenterSpaceID:
        #     listCrystalInfo[nIndex].m_strSpaceGroupID = "0"

    #获取所有晶体cif的晶体参数
    # listCrystalCifs = []
    # for curCrystal in listCrystalInfo:
    #     listCrystalCifs.append(curCrystal.m_strID+".cif")
    # readCrystalParaFromCif(listCrystalCifs)
    #readCrystalWyckoffFromCif(listCrystalCifs,listElementInfo)

    #保存晶体的元素分布
    # saveElementDistribution(listCrystalInfo,listElementInfo)

    return listCrystalInfo

#保存所有的信息
def saveAllCrystalInfo(listCrystalInfo):
    strCrystalFile = objects.BASISDIR + "crystal-out.csv"
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

