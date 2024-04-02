# 读取所有阴阳基团的HOMO，LUMO，偶极矩等
import math
import os
import objects
import method
import numpy as np


# 获取out文件中的所有键
def getBondsFromOutFile(strFile, listElementInfo):

    listBonds = []

    #获取所有原子
    outFile = open(strFile,"r")
    outFile.seek(0,0)
    nTitleLines = 0
    bReadCoord = False
    listAtoms = []
    while 1:
        strLine = outFile.readline()
        if not strLine:
            break

        strLine = ' '.join(strLine.split()) #处理多个空格间隔的情况
        if 'Input orientation:' in strLine:
            bReadCoord = True
            nTitleLines = 4
            continue

        if bReadCoord and nTitleLines > 0:
            nTitleLines = nTitleLines -1
            continue

        if bReadCoord:
            listTmpInfos = strLine.split(" ")

            #第一位必须为数字，否则读取结束
            if not listTmpInfos[0].isnumeric():
                break

            curAtom = objects.CAtom(listTmpInfos[1])
            curAtom.m_strName = listTmpInfos[1] #在getBondsFromAtoms函数中会被修改
            curAtom.m_strX = listTmpInfos[3]
            curAtom.m_strY = listTmpInfos[4]
            curAtom.m_strZ = listTmpInfos[5]
            listAtoms.append(curAtom)

    #计算所有键,通过原子类型和距离进行过滤
    listBonds = method.getBondsFromAtoms(listAtoms, listElementInfo)
    return listBonds


#计算bondValencePara
def calFlexibility(listGroupInfo,listElementInfo, isNega):
    #选择out文件所在目录
    strPolarPath = "D:\Data/posi-out-polar/"
    if isNega:
        strPolarPath = "D:\Data/nega-out-polar/"

    #获取out文件中所有原子
    #遍历文件夹中所有文件
    listfile = os.listdir(strPolarPath)
    for curFile in listfile:
        strFileName = curFile.strip(".out")
        curFile = strPolarPath + curFile
        if not curFile.endswith("out"):
            continue

        #获得当前文件对应的已经创建的Group类
        nGroupIndex = -1
        for item in listGroupInfo:
            nGroupIndex = nGroupIndex + 1
            if item.m_strID == strFileName:
                break
        if nGroupIndex < 0:
            break

        #获取所有键
        listBonds = getBondsFromOutFile(curFile, listElementInfo)

        #计算所有键相关的Flexibility，
        listFlexibilities = []
        listBondLengths = []
        for curBond in listBonds:
            listFlexibilities.append(float(curBond.m_strFlexibility))
            listBondLengths.append(float(curBond.m_strLength))

        #平均键长
        dMeanBondLength = 152 #默认是O的，单位pm。1埃=100pm
        if len(listBondLengths) > 0:
            dMeanBondLength = np.mean(listBondLengths) * 100
        dVolume = (4 * math.pi * dMeanBondLength * dMeanBondLength * dMeanBondLength) / 3

        #更新Group中的相应变量
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


#从Gaussian输出文件中直接获取GroupInfo,#获得HOMO，LUMO
def readEnergyFromOut(strPath):

    listAlphaOcc = []
    listAlphaVirt = []
    listBetaOcc = []
    listBetaVirt = []

    #遍历文件夹中所有文件
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

        #检查Alpha电子的
        strMaxOcc = "-999999"
        listAlphaOcc.reverse()
        for strEnergy in listAlphaOcc:
            strEnergy = strEnergy.strip()
            strEnergy = strEnergy.strip('\n')
            if method.is_number(strEnergy):
                strMaxOcc = strEnergy
                break

        #注意，虚轨道不需要反转
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

        #检查Beta电子的
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

        #由Hartree转为eV
        curGroup.m_strHomo = str(float(curGroup.m_strHomo) * 27.2116)
        curGroup.m_strLumo = str(float(curGroup.m_strLumo) * 27.2116)
        curGroup.m_strLumoHomo = str(float(curGroup.m_strLumoHomo) * 27.2116)

        listGroupInfo.append(curGroup)

        file.close()

    return listGroupInfo


#从Gaussian输出文件中直接获取GroupInfo,#获得偶极矩，极化率，超极化率等
def readDiPoleAndPolarFromOut(strPath, listGroupInfo):

    #Group应该在readEnergyAndDiPoleFromOut中被创建
    if len(listGroupInfo) < 1:
        return []

    #遍历文件夹中所有文件
    listfile = os.listdir(strPath)
    for curFile in listfile:

        #获得当前文件对应的已经创建的类
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

            #处理多个空格分隔的情况
            strLine = ' '.join(strLine.split())

            #从Population analysis using the SCF density行后开始读取
            if nIsStartRead < 1 and 'Population analysis using the SCF density' in strLine:
                nIsStartRead = 1
                continue

            if nIsStartRead < 1:
                continue

            if nIsReadInputOrient > 0:

                # 读取偶极矩，单位Debye
                if 'Dipole moment (field-independent basis, Debye)' in strLine:
                    nIsDiPole = 1
                    continue

                # 读取四极矩，单位Debye/Ang
                if 'Quadrupole moment (field-independent basis, Debye-Ang)' in strLine:
                    nIsQuadrupole = 1
                    continue

                #获得计算中采用的频率,由Hartree转为nm
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
                        strLine = strLine.replace("D","E")#Gaussian中的指数表示转为E
                        listTmpInfo = strLine.split(" ")
                        listGroupInfo[nIndex].m_strDipoleTotal = listTmpInfo[7]
                        nIsDiPole = -1

                if nIsQuadrupole == 1:
                    if 'XX=' in strLine:
                        strLine = strLine.replace("D","E")#Gaussian中的指数表示转为E
                        listTmpInfo = strLine.split(" ")
                        listGroupInfo[nIndex].m_strQuadrupoleXX = listTmpInfo[1]
                        listGroupInfo[nIndex].m_strQuadrupoleYY = listTmpInfo[3]
                        listGroupInfo[nIndex].m_strQuadrupoleZZ = listTmpInfo[5]
                    elif 'XY' in strLine:
                        strLine = strLine.replace("D","E")#Gaussian中的指数表示转为E
                        listTmpInfo = strLine.split(" ")
                        listGroupInfo[nIndex].m_strQuadrupoleXY = listTmpInfo[1]
                        listGroupInfo[nIndex].m_strQuadrupoleXZ = listTmpInfo[3]
                        listGroupInfo[nIndex].m_strQuadrupoleYZ = listTmpInfo[5]
                        nIsQuadrupole = 0

                if nIsAlpha == 1:
                    strLine = strLine.replace("D","E")#Gaussian中的指数表示转为E
                    listTmpInfo = strLine.split(" ")
                    if 'aniso' in strLine:
                        listGroupInfo[nIndex].m_listAnisoPolar.append(listTmpInfo[1])
                    elif 'iso' in strLine:
                        listGroupInfo[nIndex].m_listIsoPolar.append(listTmpInfo[1])
                    elif len(strLine) < 2:
                        nIsAlpha = 0

                if nIsStaticBeta == 1:
                    strLine = strLine.replace("D","E")#Gaussian中的指数表示转为E
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
                    strLine = strLine.replace("D","E")#Gaussian中的指数表示转为E
                    listTmpInfo = strLine.split(" ")
                    if 'x' in strLine and len(listTmpInfo[0]) < 3:
                        listGroupInfo[nIndex].m_listHyperPolarX.append(listTmpInfo[1])
                    elif 'y' in strLine and len(listTmpInfo[0]) < 3:
                        listGroupInfo[nIndex].m_listHyperPolarY.append(listTmpInfo[1])
                    elif 'z' in strLine and '||' not in strLine and len(listTmpInfo[0]) < 3:
                        listGroupInfo[nIndex].m_listHyperPolarZ.append(listTmpInfo[1])
                    elif len(strLine) < 2:
                        nIsDynamicBeta = 0
                        nIsReadInputOrient = 0 #到这就结束了，文本下面是偶极矩方向的值
            else:
                if 'Beta(0;0,0):' in strLine:
                    nIsStaticBeta = 1
                    continue

                if 'Beta(-2w;w,w)' in strLine:
                    nIsDynamicBeta = 1
                    continue

                if nIsStaticBeta == 1:
                    strLine = strLine.replace("D","E")#Gaussian中的指数表示转为E
                    listTmpInfo = strLine.split(" ")
                    if '||' in strLine  and 'z' not in strLine and len(listTmpInfo[0]) < 3:
                        listGroupInfo[nIndex].m_listVectorHyperPolar.append(listTmpInfo[1])
                    elif 'zzz' in strLine:
                        nIsStaticBeta = 0

                if nIsDynamicBeta == 1:
                    strLine = strLine.replace("D","E")#Gaussian中的指数表示转为E
                    listTmpInfo = strLine.split(" ")
                    if '||' in strLine and 'z' not in strLine and len(listTmpInfo[0]) < 3:
                        listGroupInfo[nIndex].m_listVectorHyperPolar.append(listTmpInfo[1])
                    elif len(strLine) < 2:
                        nIsDynamicBeta = 0
                        break #读取结束

        file.close()

    return listGroupInfo


#从Gaussian输出文件中直接获取GroupInfo
def readGroupFromOutFiles(listElementInfo, isNega):
    #选择out文件所在目录
    strEnergyPath = "D:\Data/posi-out-energy/"
    strPolarPath = "D:\Data/posi-out-polar/"
    if isNega:
        strEnergyPath = "D:\Data/nega-out-energy/"
        strPolarPath = "D:\Data/nega-out-polar/"

    listGroupInfo = readEnergyFromOut(strEnergyPath)
    listGroupInfo = readDiPoleAndPolarFromOut(strPolarPath, listGroupInfo)

    #计算每个Group的bondValencePara和Flexibility
    listGroupInfo = calFlexibility(listGroupInfo,listElementInfo, isNega)

    #计算每个Group的四极矩各向异性，极化率
    nIndex = -1
    for item in listGroupInfo:
        nIndex = nIndex + 1
        listGroupInfo[nIndex].check() #校正
        listGroupInfo[nIndex].calDerivatives()

    return listGroupInfo

#从完整的文本中直接获取GroupInfo
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


#从txt文件中读取所有Group相关信息
#可以先用readGroupFromOutFiles从out文件中读取，补充完整后现用此函数读取
def readAllGroupInfo():
    #选择out文件所在目录
    strPositiveFile = objects.BASISDIR + "positive-group-in.csv"
    listPosiGroups = readPostiveGroupInfo(strPositiveFile)
    strPositiveFile = objects.BASISDIR + "negative-group-in.csv"
    listNegaGroups = readPostiveGroupInfo(strPositiveFile)

    return listNegaGroups,listPosiGroups

#保存所有的信息
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


#单独调用时将结果直接导出csv文件
if __name__ == '__main__':
    print('getGroupInfo')




