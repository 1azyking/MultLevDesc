#根据所有的信息创建featureInfo
import math
import random
import objects
import method

#从完整的文本中直接获取
def readFeatureInfoFromTxt(strFile, listFeatureInfo):
    listFeatureInfo.clear()
    file = open(strFile,"r")
    while 1:
        strLine = file.readline()
        if not strLine:
            break

        strLine = strLine.strip();
        listInfos = strLine.split(",")
        if len(listInfos) < 13:
            continue

    file.close()


#从txt文件中读取所有Element相关信息
def readAllFeatureInfo(listFeatureInfo):
    #选择txt文件所在目录
    strFeatureFile = objects.BASISDIR + "feature-in.csv"
    readFeatureInfoFromTxt(strFeatureFile,listFeatureInfo)


#根据信息创建所有Feature
def createAllFeatureInfo(listPosiGroups, listNegaGroups, listElementInfo, listCrystalInfo):

    listFeatureInfo = []

    #根据过滤条件获得相应波长下的基团的极化率和超极化率，晶体的Gap和Dij
    dMinWaveLength = 999
    dMaxWaveLength = 1100
    # dMinWaveLength = -100
    # dMaxWaveLength = 100
    listGapLevelFilter = ["GGA","LDA","sXLDA","EXP","HSE06","PBE0","SNLO","MBJ","mBJ"]
    listBiRefLevelFilter = ["GGA","LDA","srLDA","EXP","HSE06","PBE0","SNLO","FIT"]
    listDijLevelFilter = ["GGA","LDA","srLDA","EXP","HSE06","PBE0","SNLO","FIT"]

    #阴离子中的所有中心原子
    listNegaAtom = ["Be","B",
                "Al","Si","P",
                "Ti","Ga","Ge","As",
                "Nb","Mo","Sn","Sb","Te",
                "Ta","W"]

    dictElements = {}
    for curCrystal in listCrystalInfo: #每个晶体
        curDescriptor = objects.CDescriptor(curCrystal.m_strID)

        #首先创建元素相关描述符
        dictElements.clear()
        method.getElementFromString(curCrystal.getFormula(),dictElements) #晶体中所有的元素,key为元素，value为个数
        for curElement in listElementInfo:#晶体中每个元素
            if curElement.m_strName not in dictElements.keys():
                continue
            curNum = dictElements[curElement.m_strName]

            if float(curDescriptor.m_strMinMass) > float(curElement.m_strMass):
                curDescriptor.m_strMinMass = curElement.m_strMass
            if float(curDescriptor.m_strMaxMass) < float(curElement.m_strMass):
                curDescriptor.m_strMaxMass = curElement.m_strMass
            curDescriptor.m_strAverMass = str(float(curDescriptor.m_strAverMass) + float(curElement.m_strMass) * curNum)

            if float(curDescriptor.m_strMinPeriodNum) > float(curElement.m_strPeriodNum):
                curDescriptor.m_strMinPeriodNum = curElement.m_strPeriodNum
            if float(curDescriptor.m_strMaxPeriodNum) < float(curElement.m_strPeriodNum):
                curDescriptor.m_strMaxPeriodNum = curElement.m_strPeriodNum
            curDescriptor.m_strAverPeriodNum = str(int(curDescriptor.m_strAverPeriodNum) + int(curElement.m_strPeriodNum) * curNum)

            if float(curDescriptor.m_strMinGroupNum) > float(curElement.m_strGroupNum):
                curDescriptor.m_strMinGroupNum = curElement.m_strGroupNum
            if float(curDescriptor.m_strMaxGroupNum) < float(curElement.m_strGroupNum):
                curDescriptor.m_strMaxGroupNum = curElement.m_strGroupNum
            curDescriptor.m_strAverGroupNum = str(int(curDescriptor.m_strAverGroupNum) + int(curElement.m_strGroupNum) * curNum)

            if float(curDescriptor.m_strMinRaduis) > float(curElement.m_strRadui):
                curDescriptor.m_strMinRaduis = curElement.m_strRadui
            if float(curDescriptor.m_strMaxRaduis) < float(curElement.m_strRadui):
                curDescriptor.m_strMaxRaduis = curElement.m_strRadui
            curDescriptor.m_strAverRaduis = str(float(curDescriptor.m_strAverRaduis) + float(curElement.m_strRadui) * curNum)

            if float(curDescriptor.m_strMinEleNegativity) > float(curElement.m_strEleNegativity):
                curDescriptor.m_strMinEleNegativity = curElement.m_strEleNegativity
            if float(curDescriptor.m_strMaxEleNegativity) < float(curElement.m_strEleNegativity):
                curDescriptor.m_strMaxEleNegativity = curElement.m_strEleNegativity
            curDescriptor.m_strAverEleNegativity = str(float(curDescriptor.m_strAverEleNegativity) + float(curElement.m_strEleNegativity) * curNum)

            curDescriptor.m_strAverSValenceEleNum = str(int(curDescriptor.m_strAverSValenceEleNum) + int(curElement.m_strSValenceEleNum) * curNum)
            curDescriptor.m_strAverPValenceEleNum = str(int(curDescriptor.m_strAverPValenceEleNum) + int(curElement.m_strPValenceEleNum) * curNum)
            curDescriptor.m_strAverDValenceEleNum = str(int(curDescriptor.m_strAverDValenceEleNum) + int(curElement.m_strDValenceEleNum) * curNum)
            curDescriptor.m_strAverFValenceEleNum = str(int(curDescriptor.m_strAverFValenceEleNum) + int(curElement.m_strFValenceEleNum) * curNum)
            curDescriptor.m_strAverSUnfilledState = str(int(curDescriptor.m_strAverSUnfilledState) + int(curElement.m_strSUnfilledState) * curNum)
            curDescriptor.m_strAverPUnfilledState = str(int(curDescriptor.m_strAverPUnfilledState) + int(curElement.m_strPUnfilledState) * curNum)
            curDescriptor.m_strAverDUnfilledState = str(int(curDescriptor.m_strAverDUnfilledState) + int(curElement.m_strDUnfilledState) * curNum)
            curDescriptor.m_strAverFUnfilledState = str(int(curDescriptor.m_strAverFUnfilledState) + int(curElement.m_strFUnfilledState) * curNum)

            if float(curDescriptor.m_strMinIonizationEner) > float(curElement.m_strionizationEner):
                curDescriptor.m_strMinIonizationEner = curElement.m_strionizationEner
            if float(curDescriptor.m_strMaxIonizationEner) < float(curElement.m_strionizationEner):
                curDescriptor.m_strMaxIonizationEner = curElement.m_strionizationEner
            curDescriptor.m_strAverIonizationEner = str(float(curDescriptor.m_strAverIonizationEner) + float(curElement.m_strionizationEner) * curNum)

            if float(curDescriptor.m_strMinEleAffinity) > float(curElement.m_strEleAffinity):
                curDescriptor.m_strMinEleAffinity = curElement.m_strEleAffinity
            if float(curDescriptor.m_strMaxEleAffinity) < float(curElement.m_strEleAffinity):
                curDescriptor.m_strMaxEleAffinity = curElement.m_strEleAffinity
            curDescriptor.m_strAverEleAffinity = str(float(curDescriptor.m_strAverEleAffinity) + float(curElement.m_strEleAffinity) * curNum)

            if float(curDescriptor.m_strMinMeltingPoint) > float(curElement.m_strMeltingPoint):
                curDescriptor.m_strMinMeltingPoint = curElement.m_strMeltingPoint
            if float(curDescriptor.m_strMaxMeltingPoint) < float(curElement.m_strMeltingPoint):
                curDescriptor.m_strMaxMeltingPoint = curElement.m_strMeltingPoint
            curDescriptor.m_strAverMeltingPoint = str(float(curDescriptor.m_strAverMeltingPoint) + float(curElement.m_strMeltingPoint) * curNum)

            if float(curDescriptor.m_strMinBoilingPoint) > float(curElement.m_strBoilingPoint):
                curDescriptor.m_strMinBoilingPoint = curElement.m_strBoilingPoint
            if float(curDescriptor.m_strMaxBoilingPoint) < float(curElement.m_strBoilingPoint):
                curDescriptor.m_strMaxBoilingPoint = curElement.m_strBoilingPoint
            curDescriptor.m_strAverBoilingPoint = str(float(curDescriptor.m_strAverBoilingPoint) + float(curElement.m_strBoilingPoint) * curNum)

            if float(curDescriptor.m_strMinDensity) > float(curElement.m_strDensity):
                curDescriptor.m_strMinDensity = curElement.m_strDensity
            if float(curDescriptor.m_strMaxDensity) < float(curElement.m_strDensity):
                curDescriptor.m_strMaxDensity = curElement.m_strDensity
            curDescriptor.m_strAverDensity = str(float(curDescriptor.m_strAverDensity) + float(curElement.m_strDensity) * curNum)

        #求一些平均值
        nTotalNum = 0
        for index in dictElements.values():
            nTotalNum = nTotalNum + index

        curDescriptor.m_strAverMass = str(float(curDescriptor.m_strAverMass) / nTotalNum)
        curDescriptor.m_strAverPeriodNum = str(float(curDescriptor.m_strAverPeriodNum) / nTotalNum)
        curDescriptor.m_strAverGroupNum = str(float(curDescriptor.m_strAverGroupNum) / nTotalNum)
        curDescriptor.m_strAverRaduis = str(float(curDescriptor.m_strAverRaduis) / nTotalNum)
        curDescriptor.m_strAverEleNegativity = str(float(curDescriptor.m_strAverEleNegativity) / nTotalNum)
        curDescriptor.m_strAverIonizationEner = str(float(curDescriptor.m_strAverIonizationEner) / nTotalNum)
        curDescriptor.m_strAverEleAffinity = str(float(curDescriptor.m_strAverEleAffinity) / nTotalNum)
        curDescriptor.m_strAverMeltingPoint = str(float(curDescriptor.m_strAverMeltingPoint) / nTotalNum)
        curDescriptor.m_strAverBoilingPoint = str(float(curDescriptor.m_strAverBoilingPoint) / nTotalNum)
        curDescriptor.m_strAverDensity = str(float(curDescriptor.m_strAverDensity) / nTotalNum)
        curDescriptor.m_strAverSValenceEleNum = str(int(curDescriptor.m_strAverSValenceEleNum)/nTotalNum)
        curDescriptor.m_strAverPValenceEleNum = str(int(curDescriptor.m_strAverPValenceEleNum)/nTotalNum)
        curDescriptor.m_strAverDValenceEleNum = str(int(curDescriptor.m_strAverDValenceEleNum)/nTotalNum)
        curDescriptor.m_strAverFValenceEleNum = str(int(curDescriptor.m_strAverFValenceEleNum)/nTotalNum)
        curDescriptor.m_strAverSUnfilledState = str(int(curDescriptor.m_strAverSUnfilledState)/nTotalNum)
        curDescriptor.m_strAverPUnfilledState = str(int(curDescriptor.m_strAverPUnfilledState)/nTotalNum)
        curDescriptor.m_strAverDUnfilledState = str(int(curDescriptor.m_strAverDUnfilledState)/nTotalNum)
        curDescriptor.m_strAverFUnfilledState = str(int(curDescriptor.m_strAverFUnfilledState)/nTotalNum)

        #创建基团相关描述符
        for curGroupName in curCrystal.m_listNegativeGroups:

            nGroupIndex = -1
            bIsOK = False
            for nGroupIndex in range(0,len(listNegaGroups)):
                strGroupName = listNegaGroups[nGroupIndex].m_strID.lstrip("0123456789") #需要去除数字
                if strGroupName == curGroupName:
                    bIsOK = True
                    break

            if not bIsOK:
                print(curGroupName,"阴离子基团不存在")
                continue

            curNegaGroup = listNegaGroups[nGroupIndex]

            if float(curDescriptor.m_strNegaMinGap) > float(curNegaGroup.m_strLumoHomo):
                curDescriptor.m_strNegaMinGap = curNegaGroup.m_strLumoHomo
            if float(curDescriptor.m_strNegaMaxGap) < float(curNegaGroup.m_strLumoHomo):
                curDescriptor.m_strNegaMaxGap = curNegaGroup.m_strLumoHomo
            curDescriptor.m_strNegaAverGap = str(float(curDescriptor.m_strNegaAverGap) + float(curNegaGroup.m_strLumoHomo))

            if float(curDescriptor.m_strMinNegaDipoleTotal) > float(curNegaGroup.m_strDipoleTotal):
                curDescriptor.m_strMinNegaDipoleTotal = curNegaGroup.m_strDipoleTotal
            if float(curDescriptor.m_strMaxNegaDipoleTotal) < float(curNegaGroup.m_strDipoleTotal):
                curDescriptor.m_strMaxNegaDipoleTotal = curNegaGroup.m_strDipoleTotal
            curDescriptor.m_strAverNegaDipoleTotal = str(float(curDescriptor.m_strAverNegaDipoleTotal) + float(curNegaGroup.m_strDipoleTotal))

            if float(curDescriptor.m_strMinNegaAnisoQuadrupole) > float(curNegaGroup.m_strAnisoQuadrupole):
                curDescriptor.m_strMinNegaAnisoQuadrupole = curNegaGroup.m_strAnisoQuadrupole
            if float(curDescriptor.m_strMaxNegaAnisoQuadrupole) < float(curNegaGroup.m_strAnisoQuadrupole):
                curDescriptor.m_strMaxNegaAnisoQuadrupole = curNegaGroup.m_strAnisoQuadrupole
            curDescriptor.m_strAverNegaAnisoQuadrupole = str(float(curDescriptor.m_strAverNegaDipoleTotal) + float(curNegaGroup.m_strAnisoQuadrupole))

            #根据过滤条件获得相应波长下的极化率和超极化率
            nIndex = -1
            for curFreq in curNegaGroup.m_listPolarFreq:
                nIndex = nIndex + 1
                if float(curFreq) >= dMinWaveLength and float(curFreq) <= dMaxWaveLength:
                    break

            if nIndex < len(curNegaGroup.m_listPolarFreq) and nIndex > -1:
                dCurAnisoPolar = float(curNegaGroup.m_listAnisoPolar[nIndex])
                dCurTotalHyperPolar = float(curNegaGroup.m_listTotalHyperPolar[nIndex])
                dCurVectorHyperPolar = float(curNegaGroup.m_listVectorHyperPolar[nIndex])

                if float(curDescriptor.m_strMinNegaAnisoPolar) > dCurAnisoPolar:
                    curDescriptor.m_strMinNegaAnisoPolar = str(dCurAnisoPolar)
                if float(curDescriptor.m_strMaxNegaAnisoPolar) < dCurAnisoPolar:
                    curDescriptor.m_strMaxNegaAnisoPolar = str(dCurAnisoPolar)
                curDescriptor.m_strAverNegaAnisoPolar = str(float(curDescriptor.m_strAverNegaAnisoPolar) + dCurAnisoPolar)

                if float(curDescriptor.m_strMinNegaTotalHyperPolar) > dCurTotalHyperPolar:
                    curDescriptor.m_strMinNegaTotalHyperPolar = str(dCurTotalHyperPolar)
                if float(curDescriptor.m_strMaxNegaTotalHyperPolar) < dCurTotalHyperPolar:
                    curDescriptor.m_strMaxNegaTotalHyperPolar = str(dCurTotalHyperPolar)
                curDescriptor.m_strAverNegaTotalHyperPolar = str(float(curDescriptor.m_strAverNegaTotalHyperPolar) + dCurTotalHyperPolar)

                if float(curDescriptor.m_strMinNegaVectorHyperPolar) > dCurVectorHyperPolar:
                    curDescriptor.m_strMinNegaVectorHyperPolar = str(dCurTotalHyperPolar)
                if float(curDescriptor.m_strMaxNegaVectorHyperPolar) < dCurVectorHyperPolar:
                    curDescriptor.m_strMaxNegaVectorHyperPolar = str(dCurTotalHyperPolar)
                curDescriptor.m_strAverNegaVectorHyperPolar = str(float(curDescriptor.m_strAverNegaVectorHyperPolar) + dCurVectorHyperPolar)

            if float(curDescriptor.m_strMinNegaFlexibility) > float(curNegaGroup.m_strAverFlexibility):
                curDescriptor.m_strMinNegaFlexibility = str(float(curNegaGroup.m_strAverFlexibility))
            if float(curDescriptor.m_strMaxNegaFlexibility) < float(curNegaGroup.m_strAverFlexibility):
                curDescriptor.m_strMaxNegaFlexibility = str(float(curNegaGroup.m_strAverFlexibility))
            curDescriptor.m_strTotalNegaFlexibility = str(float(curDescriptor.m_strTotalNegaFlexibility) + float(curNegaGroup.m_strAverFlexibility))

            curDescriptor.m_strTotalNegaCharge = str(float(curDescriptor.m_strTotalNegaCharge) + float(curNegaGroup.m_strCharge))
            curDescriptor.m_strTotalNegaMultiplicity = str(int(curDescriptor.m_strTotalNegaMultiplicity) + int(curNegaGroup.m_strMultiplicity))

            #最大体积
            if float(curDescriptor.m_strMaxNegaVolume) < float(curNegaGroup.m_strVolume):
                curDescriptor.m_strMaxNegaVolume = curNegaGroup.m_strVolume

        #设定值并求平均值
        curLen = len(curCrystal.m_listNegativeGroups)
        curDescriptor.m_strNegaAverGap = str(float(curDescriptor.m_strNegaAverGap) / curLen)
        curDescriptor.m_strAverNegaDipoleTotal = str(float(curDescriptor.m_strAverNegaDipoleTotal) / curLen)
        curDescriptor.m_strAverNegaAnisoQuadrupole = str(float(curDescriptor.m_strAverNegaAnisoQuadrupole) / curLen)
        curDescriptor.m_strAverNegaAnisoPolar = str(float(curDescriptor.m_strAverNegaAnisoPolar) / curLen)
        curDescriptor.m_strAverNegaTotalHyperPolar = str(float(curDescriptor.m_strAverNegaTotalHyperPolar) / curLen)
        curDescriptor.m_strAverNegaVectorHyperPolar = str(float(curDescriptor.m_strAverNegaVectorHyperPolar) / curLen)
        curDescriptor.m_strAverNegaFlexibility = str(float(curDescriptor.m_strTotalNegaFlexibility) / curLen)
        curDescriptor.m_strAverNegaCharge = str(float(curDescriptor.m_strTotalNegaCharge) / curLen)

       #创建基团相关描述符
        for curGroupName in curCrystal.m_listPostiveGroups:

            nGroupIndex = -1
            bIsOK = False
            for nGroupIndex in range(0,len(listPosiGroups)):
                strGroupName = listPosiGroups[nGroupIndex].m_strID.lstrip("0123456789") #需要去除数字
                if strGroupName == curGroupName:
                    bIsOK = True
                    break

            if not bIsOK:
                print(curGroupName,"金属氧化物不存在")
                continue

            curPosiGroup = listPosiGroups[nGroupIndex]

            if float(curDescriptor.m_strPosiMinGap) > float(curPosiGroup.m_strLumoHomo):
                curDescriptor.m_strPosiMinGap = curPosiGroup.m_strLumoHomo
            if float(curDescriptor.m_strPosiMaxGap) < float(curPosiGroup.m_strLumoHomo):
                curDescriptor.m_strPosiMaxGap = curPosiGroup.m_strLumoHomo
            curDescriptor.m_strPosiAverGap = str(float(curDescriptor.m_strPosiAverGap) + float(curPosiGroup.m_strLumoHomo))

            if float(curDescriptor.m_strMinPosiDipoleTotal) > float(curPosiGroup.m_strDipoleTotal):
                curDescriptor.m_strMinPosiDipoleTotal = curPosiGroup.m_strDipoleTotal
            if float(curDescriptor.m_strMaxPosiDipoleTotal) < float(curPosiGroup.m_strDipoleTotal):
                curDescriptor.m_strMaxPosiDipoleTotal = curPosiGroup.m_strDipoleTotal
            curDescriptor.m_strAverPosiDipoleTotal = str(float(curDescriptor.m_strAverPosiDipoleTotal) + float(curPosiGroup.m_strDipoleTotal))

            if float(curDescriptor.m_strMinPosiAnisoQuadrupole) > float(curPosiGroup.m_strAnisoQuadrupole):
                curDescriptor.m_strMinPosiAnisoQuadrupole = curPosiGroup.m_strAnisoQuadrupole
            if float(curDescriptor.m_strMaxPosiAnisoQuadrupole) < float(curPosiGroup.m_strAnisoQuadrupole):
                curDescriptor.m_strMaxPosiAnisoQuadrupole = curPosiGroup.m_strAnisoQuadrupole
            curDescriptor.m_strAverPosiAnisoQuadrupole = str(float(curDescriptor.m_strAverPosiDipoleTotal) + float(curPosiGroup.m_strAnisoQuadrupole))

            #根据过滤条件获得相应波长下的极化率和超极化率
            nIndex = -1
            for curFreq in curPosiGroup.m_listPolarFreq:
                nIndex = nIndex + 1
                if float(curFreq) >= dMinWaveLength and float(curFreq) <= dMaxWaveLength:
                    break

            if nIndex < len(curPosiGroup.m_listPolarFreq) and nIndex > -1:
                dCurAnisoPolar = float(curPosiGroup.m_listAnisoPolar[nIndex])
                dCurTotalHyperPolar = float(curPosiGroup.m_listTotalHyperPolar[nIndex])
                dCurVectorHyperPolar = float(curPosiGroup.m_listVectorHyperPolar[nIndex])

                if float(curDescriptor.m_strMinPosiAnisoPolar) > dCurAnisoPolar:
                    curDescriptor.m_strMinPosiAnisoPolar = str(dCurAnisoPolar)
                if float(curDescriptor.m_strMaxPosiAnisoPolar) < dCurAnisoPolar:
                    curDescriptor.m_strMaxPosiAnisoPolar = str(dCurAnisoPolar)
                curDescriptor.m_strAverPosiAnisoPolar = str(float(curDescriptor.m_strAverPosiAnisoPolar) + dCurAnisoPolar)

                if float(curDescriptor.m_strMinPosiTotalHyperPolar) > dCurTotalHyperPolar:
                    curDescriptor.m_strMinPosiTotalHyperPolar = str(dCurTotalHyperPolar)
                if float(curDescriptor.m_strMaxPosiTotalHyperPolar) < dCurTotalHyperPolar:
                    curDescriptor.m_strMaxPosiTotalHyperPolar = str(dCurTotalHyperPolar)
                curDescriptor.m_strAverPosiTotalHyperPolar = str(float(curDescriptor.m_strAverPosiTotalHyperPolar) + dCurTotalHyperPolar)

                if float(curDescriptor.m_strMinPosiVectorHyperPolar) > dCurVectorHyperPolar:
                    curDescriptor.m_strMinPosiVectorHyperPolar = str(dCurVectorHyperPolar)
                if float(curDescriptor.m_strMaxPosiVectorHyperPolar) < dCurVectorHyperPolar:
                    curDescriptor.m_strMaxPosiVectorHyperPolar = str(dCurVectorHyperPolar)
                curDescriptor.m_strAverPosiVectorHyperPolar = str(float(curDescriptor.m_strAverNegaVectorHyperPolar) + dCurVectorHyperPolar)


            if float(curDescriptor.m_strMinPosiFlexibility) > float(curPosiGroup.m_strAverFlexibility):
                curDescriptor.m_strMinPosiFlexibility = str(float(curPosiGroup.m_strAverFlexibility))
            if float(curDescriptor.m_strMaxPosiFlexibility) < float(curPosiGroup.m_strAverFlexibility):
                curDescriptor.m_strMaxPosiFlexibility = str(float(curPosiGroup.m_strAverFlexibility))
            curDescriptor.m_strTotalPosiFlexibility = str(float(curDescriptor.m_strTotalPosiFlexibility) + float(curPosiGroup.m_strAverFlexibility))

            curDescriptor.m_strTotalPosiCharge = str(float(curDescriptor.m_strTotalPosiCharge) + float(curPosiGroup.m_strCharge))
            curDescriptor.m_strTotalPosiMultiplicity = str(int(curDescriptor.m_strTotalPosiMultiplicity) + int(curPosiGroup.m_strMultiplicity))

            #最大体积
            if float(curDescriptor.m_strMaxPosiVolume) < float(curPosiGroup.m_strVolume):
                curDescriptor.m_strMaxPosiVolume = curPosiGroup.m_strVolume

        #设定值并求平均值
        curLen = len(curCrystal.m_listPostiveGroups)
        curDescriptor.m_strPosiAverGap = str(float(curDescriptor.m_strPosiAverGap) / curLen)
        curDescriptor.m_strAverPosiDipoleTotal = str(float(curDescriptor.m_strAverPosiDipoleTotal) / curLen)
        curDescriptor.m_strAverPosiAnisoQuadrupole = str(float(curDescriptor.m_strAverPosiAnisoQuadrupole) / curLen)
        curDescriptor.m_strAverPosiAnisoPolar = str(float(curDescriptor.m_strAverPosiAnisoPolar) / curLen)
        curDescriptor.m_strAverPosiTotalHyperPolar = str(float(curDescriptor.m_strAverPosiTotalHyperPolar) / curLen)
        curDescriptor.m_strAverPosiVectorHyperPolar = str(float(curDescriptor.m_strAverPosiVectorHyperPolar) / curLen)
        curDescriptor.m_strAverPosiFlexibility = str(float(curDescriptor.m_strTotalPosiFlexibility) / curLen)
        curDescriptor.m_strAverPosiCharge = str(float(curDescriptor.m_strTotalPosiCharge) / curLen)

        #添加阴阳基团的独自的扩展描述符
        curDescriptor.m_strNegaMaxMinDiffGap = str(float(curDescriptor.m_strNegaMaxGap) - float(curDescriptor.m_strNegaMinGap))
        curDescriptor.m_strPosiMaxMinDiffGap = str(float(curDescriptor.m_strPosiMaxGap) - float(curDescriptor.m_strPosiMinGap))
        curDescriptor.m_strNegaMaxMinDiffDipoleTotal = str(float(curDescriptor.m_strMaxNegaDipoleTotal) - float(curDescriptor.m_strMinNegaDipoleTotal))
        curDescriptor.m_strPosiMaxMinDiffDipoleTotal = str(float(curDescriptor.m_strMaxPosiDipoleTotal) - float(curDescriptor.m_strMinPosiDipoleTotal))
        curDescriptor.m_strNegaMaxMinDiffAnisoQuadrupole = str(float(curDescriptor.m_strMaxNegaAnisoQuadrupole) - float(curDescriptor.m_strMinNegaAnisoQuadrupole))
        curDescriptor.m_strPosiMaxMinDiffAnisoQuadrupole = str(float(curDescriptor.m_strMaxPosiAnisoQuadrupole) - float(curDescriptor.m_strMinPosiAnisoQuadrupole))
        curDescriptor.m_strNegaMaxMinDiffAnisoPolar = str(float(curDescriptor.m_strMaxNegaAnisoPolar) - float(curDescriptor.m_strMinNegaAnisoPolar))
        curDescriptor.m_strPosiMaxMinDiffAnisoPolar = str(float(curDescriptor.m_strMaxPosiAnisoPolar) - float(curDescriptor.m_strMinPosiAnisoPolar))
        curDescriptor.m_strNegaMaxMinDiffTotalHyperPolar = str(float(curDescriptor.m_strMaxNegaTotalHyperPolar) - float(curDescriptor.m_strMinNegaTotalHyperPolar))
        curDescriptor.m_strPosiMaxMinDiffTotalHyperPolar = str(float(curDescriptor.m_strMaxPosiTotalHyperPolar) - float(curDescriptor.m_strMinPosiTotalHyperPolar))
        curDescriptor.m_strNegaMaxMinDiffVectorHyperPolar = str(float(curDescriptor.m_strMaxNegaVectorHyperPolar) - float(curDescriptor.m_strMinNegaVectorHyperPolar))
        curDescriptor.m_strPosiMaxMinDiffVectorHyperPolar = str(float(curDescriptor.m_strMaxPosiVectorHyperPolar) - float(curDescriptor.m_strMinPosiVectorHyperPolar))
        curDescriptor.m_strNegaMaxMinDiffFlexibility = str(float(curDescriptor.m_strMaxNegaFlexibility) - float(curDescriptor.m_strMinNegaFlexibility))
        curDescriptor.m_strPosiMaxMinDiffFlexibility = str(float(curDescriptor.m_strMaxPosiFlexibility) - float(curDescriptor.m_strMinPosiFlexibility))
        curDescriptor.m_strNegaMaxMinSumGap = str(float(curDescriptor.m_strNegaMaxGap) + float(curDescriptor.m_strNegaMinGap))
        curDescriptor.m_strPosiMaxMinSumGap = str(float(curDescriptor.m_strPosiMaxGap) + float(curDescriptor.m_strPosiMinGap))
        curDescriptor.m_strNegaMaxMinSumDipoleTotal = str(float(curDescriptor.m_strMaxNegaDipoleTotal) + float(curDescriptor.m_strMinNegaDipoleTotal))
        curDescriptor.m_strPosiMaxMinSumDipoleTotal = str(float(curDescriptor.m_strMaxPosiDipoleTotal) + float(curDescriptor.m_strMinPosiDipoleTotal))
        curDescriptor.m_strNegaMaxMinSumAnisoQuadrupole = str(float(curDescriptor.m_strMaxNegaAnisoQuadrupole) + float(curDescriptor.m_strMinNegaAnisoQuadrupole))
        curDescriptor.m_strPosiMaxMinSumAnisoQuadrupole = str(float(curDescriptor.m_strMaxPosiAnisoQuadrupole) + float(curDescriptor.m_strMinPosiAnisoQuadrupole))
        curDescriptor.m_strNegaMaxMinSumAnisoPolar = str(float(curDescriptor.m_strMaxNegaAnisoPolar) + float(curDescriptor.m_strMinNegaAnisoPolar))
        curDescriptor.m_strPosiMaxMinSumAnisoPolar = str(float(curDescriptor.m_strMaxPosiAnisoPolar) + float(curDescriptor.m_strMinPosiAnisoPolar))
        curDescriptor.m_strNegaMaxMinSumTotalHyperPolar = str(float(curDescriptor.m_strMaxNegaTotalHyperPolar) + float(curDescriptor.m_strMinNegaTotalHyperPolar))
        curDescriptor.m_strPosiMaxMinSumTotalHyperPolar = str(float(curDescriptor.m_strMaxPosiTotalHyperPolar) + float(curDescriptor.m_strMinPosiTotalHyperPolar))
        curDescriptor.m_strNegaMaxMinSumVectorHyperPolar = str(float(curDescriptor.m_strMaxNegaVectorHyperPolar) + float(curDescriptor.m_strMinNegaVectorHyperPolar))
        curDescriptor.m_strPosiMaxMinSumVectorHyperPolar = str(float(curDescriptor.m_strMaxPosiVectorHyperPolar) + float(curDescriptor.m_strMinPosiVectorHyperPolar))
        curDescriptor.m_strNegaMaxMinSumFlexibility = str(float(curDescriptor.m_strMaxNegaFlexibility) + float(curDescriptor.m_strMinNegaFlexibility))
        curDescriptor.m_strPosiMaxMinSumFlexibility = str(float(curDescriptor.m_strMaxPosiFlexibility) + float(curDescriptor.m_strMinPosiFlexibility))

        #通过阴阳基团来构建虚拟晶体(VirtualCrystal,VC)的物理量,注意，此块处理的全部为阴阳基团之间的，不能为同类
        curDescriptor.m_strVCGap = float(curDescriptor.m_strNegaMinGap)
        dCurVCValue = float(curDescriptor.m_strPosiMinGap)
        if curDescriptor.m_strVCGap > dCurVCValue:
            curDescriptor.m_strVCGap = dCurVCValue
        curDescriptor.m_strVCGap = str(curDescriptor.m_strVCGap)
        curDescriptor.m_strMaxVCDiffDipoleTotal = math.fabs(float(curDescriptor.m_strMaxNegaDipoleTotal) - float(curDescriptor.m_strMinPosiDipoleTotal))
        dCurVCValue = math.fabs(float(curDescriptor.m_strMinNegaDipoleTotal) - float(curDescriptor.m_strMaxPosiDipoleTotal))
        if curDescriptor.m_strMaxVCDiffDipoleTotal < dCurVCValue:
            curDescriptor.m_strMaxVCDiffDipoleTotal = dCurVCValue
        curDescriptor.m_strMaxVCDiffDipoleTotal = str(curDescriptor.m_strMaxVCDiffDipoleTotal)
        curDescriptor.m_strMinVCDiffDipoleTotal = str(math.fabs(float(curDescriptor.m_strMinNegaDipoleTotal) - float(curDescriptor.m_strMinPosiDipoleTotal)))
        curDescriptor.m_strMaxVCDiffAnisoQuadrupole = math.fabs(float(curDescriptor.m_strMaxNegaAnisoQuadrupole) - float(curDescriptor.m_strMinPosiAnisoQuadrupole))
        dCurVCValue = math.fabs(float(curDescriptor.m_strMinNegaAnisoQuadrupole) - float(curDescriptor.m_strMaxPosiAnisoQuadrupole))
        if curDescriptor.m_strMaxVCDiffAnisoQuadrupole < dCurVCValue:
            curDescriptor.m_strMaxVCDiffAnisoQuadrupole = dCurVCValue
        curDescriptor.m_strMaxVCDiffAnisoQuadrupole = str(curDescriptor.m_strMaxVCDiffAnisoQuadrupole)
        curDescriptor.m_strMinVCDiffAnisoQuadrupole = str(math.fabs(float(curDescriptor.m_strMinNegaAnisoQuadrupole) - float(curDescriptor.m_strMinPosiAnisoQuadrupole)))
        curDescriptor.m_strMaxVCDiffTotalHyperPolar = math.fabs(float(curDescriptor.m_strMaxNegaTotalHyperPolar) - float(curDescriptor.m_strMinPosiTotalHyperPolar))
        dCurVCValue = math.fabs(float(curDescriptor.m_strMinNegaTotalHyperPolar) - float(curDescriptor.m_strMaxPosiTotalHyperPolar))
        if curDescriptor.m_strMaxVCDiffTotalHyperPolar < dCurVCValue:
            curDescriptor.m_strMaxVCDiffTotalHyperPolar = dCurVCValue
        curDescriptor.m_strMaxVCDiffTotalHyperPolar = str(curDescriptor.m_strMaxVCDiffTotalHyperPolar)
        curDescriptor.m_strMinVCDiffTotalHyperPolar = str(math.fabs(float(curDescriptor.m_strMinNegaTotalHyperPolar) - float(curDescriptor.m_strMinPosiTotalHyperPolar)))
        curDescriptor.m_strMaxVCDiffVectorHyperPolar = math.fabs(float(curDescriptor.m_strMaxNegaVectorHyperPolar) - float(curDescriptor.m_strMinPosiVectorHyperPolar))
        dCurVCValue = math.fabs(float(curDescriptor.m_strMinNegaVectorHyperPolar) - float(curDescriptor.m_strMaxPosiVectorHyperPolar))
        if curDescriptor.m_strMaxVCDiffVectorHyperPolar < dCurVCValue:
            curDescriptor.m_strMaxVCDiffVectorHyperPolar = dCurVCValue
        curDescriptor.m_strMaxVCDiffVectorHyperPolar = str(curDescriptor.m_strMaxVCDiffVectorHyperPolar)
        curDescriptor.m_strMinVCDiffVectorHyperPolar = str(math.fabs(float(curDescriptor.m_strMinNegaVectorHyperPolar) - float(curDescriptor.m_strMinPosiVectorHyperPolar)))
        curDescriptor.m_strMaxVCDiffFlexibility = math.fabs(float(curDescriptor.m_strMaxNegaFlexibility) - float(curDescriptor.m_strMinPosiFlexibility))
        dCurVCValue = math.fabs(float(curDescriptor.m_strMinNegaFlexibility) - float(curDescriptor.m_strMaxPosiFlexibility))
        if curDescriptor.m_strMaxVCDiffFlexibility < dCurVCValue:
            curDescriptor.m_strMaxVCDiffFlexibility = dCurVCValue
        curDescriptor.m_strMaxVCDiffFlexibility = str(curDescriptor.m_strMaxVCDiffFlexibility)
        curDescriptor.m_strMinVCDiffFlexibility = str(math.fabs(float(curDescriptor.m_strMinNegaFlexibility) - float(curDescriptor.m_strMinPosiFlexibility)))

        #基团相关基本描述符四则运算后
        curDescriptor.m_strDiffMinGap = str(float(curDescriptor.m_strNegaMinGap) - float(curDescriptor.m_strPosiMinGap))
        curDescriptor.m_strSumMinGap = str(float(curDescriptor.m_strNegaMinGap) + float(curDescriptor.m_strPosiMinGap))
        curDescriptor.m_strDiffAverGap = str(float(curDescriptor.m_strNegaAverGap) - float(curDescriptor.m_strPosiAverGap))
        curDescriptor.m_strSumAverGap = str(float(curDescriptor.m_strNegaAverGap) + float(curDescriptor.m_strPosiAverGap))
        curDescriptor.m_strDiffMaxGap = str(float(curDescriptor.m_strNegaMaxGap) - float(curDescriptor.m_strPosiMaxGap))
        curDescriptor.m_strSumMaxGap = str(float(curDescriptor.m_strNegaMaxGap) + float(curDescriptor.m_strPosiMaxGap))
        curDescriptor.m_strDiffMinDipoleTotal = str(float(curDescriptor.m_strMinNegaDipoleTotal) - float(curDescriptor.m_strMinPosiDipoleTotal))
        curDescriptor.m_strSumMinDipoleTotal = str(float(curDescriptor.m_strMinNegaDipoleTotal) + float(curDescriptor.m_strMinPosiDipoleTotal))
        curDescriptor.m_strDiffAverDipoleTotal = str(float(curDescriptor.m_strAverNegaDipoleTotal) - float(curDescriptor.m_strAverPosiDipoleTotal))
        curDescriptor.m_strSumAverDipoleTotal = str(float(curDescriptor.m_strAverNegaDipoleTotal) + float(curDescriptor.m_strAverPosiDipoleTotal))
        curDescriptor.m_strDiffMaxDipoleTotal = str(float(curDescriptor.m_strMaxNegaDipoleTotal) - float(curDescriptor.m_strMaxPosiDipoleTotal))
        curDescriptor.m_strSumMaxDipoleTotal = str(float(curDescriptor.m_strMaxNegaDipoleTotal) + float(curDescriptor.m_strMaxPosiDipoleTotal))
        curDescriptor.m_strDiffMinAnisoQuadrupole = str(float(curDescriptor.m_strMinNegaAnisoQuadrupole) - float(curDescriptor.m_strMinPosiAnisoQuadrupole))
        curDescriptor.m_strSumMinAnisoQuadrupole = str(float(curDescriptor.m_strMinNegaAnisoQuadrupole) + float(curDescriptor.m_strMinPosiAnisoQuadrupole))
        curDescriptor.m_strDiffAverAnisoQuadrupole = str(float(curDescriptor.m_strAverNegaAnisoQuadrupole) - float(curDescriptor.m_strAverPosiAnisoQuadrupole))
        curDescriptor.m_strSumAverAnisoQuadrupole = str(float(curDescriptor.m_strAverNegaAnisoQuadrupole) + float(curDescriptor.m_strAverPosiAnisoQuadrupole))
        curDescriptor.m_strDiffMaxAnisoQuadrupole = str(float(curDescriptor.m_strMaxNegaAnisoQuadrupole) - float(curDescriptor.m_strMaxPosiAnisoQuadrupole))
        curDescriptor.m_strSumMaxAnisoQuadrupole = str(float(curDescriptor.m_strMaxNegaAnisoQuadrupole) + float(curDescriptor.m_strMaxPosiAnisoQuadrupole))
        curDescriptor.m_strDiffMinAnisoPolar = str(float(curDescriptor.m_strMinNegaAnisoPolar) - float(curDescriptor.m_strMinPosiAnisoPolar))
        curDescriptor.m_strSumMinAnisoPolar = str(float(curDescriptor.m_strMinNegaAnisoPolar) + float(curDescriptor.m_strMinPosiAnisoPolar))
        curDescriptor.m_strDiffAverAnisoPolar = str(float(curDescriptor.m_strAverNegaAnisoPolar) - float(curDescriptor.m_strAverPosiAnisoPolar))
        curDescriptor.m_strSumAverAnisoPolar = str(float(curDescriptor.m_strAverNegaAnisoPolar) + float(curDescriptor.m_strAverPosiAnisoPolar))
        curDescriptor.m_strDiffMaxAnisoPolar = str(float(curDescriptor.m_strMaxNegaAnisoPolar) - float(curDescriptor.m_strMinNegaAnisoPolar))
        curDescriptor.m_strSumMaxAnisoPolar = str(float(curDescriptor.m_strMaxNegaAnisoPolar) + float(curDescriptor.m_strMinNegaAnisoPolar))
        curDescriptor.m_strDiffMinTotalHyperPolar = str(float(curDescriptor.m_strMinNegaTotalHyperPolar) - float(curDescriptor.m_strMinPosiTotalHyperPolar))
        curDescriptor.m_strSumMinTotalHyperPolar = str(float(curDescriptor.m_strMinNegaTotalHyperPolar) + float(curDescriptor.m_strMinPosiTotalHyperPolar))
        curDescriptor.m_strDiffAverTotalHyperPolar = str(float(curDescriptor.m_strAverNegaTotalHyperPolar) - float(curDescriptor.m_strAverPosiTotalHyperPolar))
        curDescriptor.m_strSumAverTotalHyperPolar = str(float(curDescriptor.m_strAverNegaTotalHyperPolar) + float(curDescriptor.m_strAverPosiTotalHyperPolar))
        curDescriptor.m_strDiffMaxTotalHyperPolar = str(float(curDescriptor.m_strMaxNegaTotalHyperPolar) - float(curDescriptor.m_strMaxPosiTotalHyperPolar))
        curDescriptor.m_strSumMaxTotalHyperPolar = str(float(curDescriptor.m_strMaxNegaTotalHyperPolar) + float(curDescriptor.m_strMaxPosiTotalHyperPolar))
        curDescriptor.m_strDiffMinVectorHyperPolar = str(float(curDescriptor.m_strMinNegaVectorHyperPolar) - float(curDescriptor.m_strMinPosiVectorHyperPolar))
        curDescriptor.m_strSumMinVectorHyperPolar = str(float(curDescriptor.m_strMinNegaVectorHyperPolar) + float(curDescriptor.m_strMinPosiVectorHyperPolar))
        curDescriptor.m_strDiffAverVectorHyperPolar = str(float(curDescriptor.m_strAverNegaVectorHyperPolar) + float(curDescriptor.m_strAverPosiVectorHyperPolar))
        curDescriptor.m_strSumAverVectorHyperPolar = str(float(curDescriptor.m_strAverNegaVectorHyperPolar) - float(curDescriptor.m_strAverPosiVectorHyperPolar))
        curDescriptor.m_strDiffMaxVectorHyperPolar = str(float(curDescriptor.m_strMaxNegaVectorHyperPolar) + float(curDescriptor.m_strMaxPosiVectorHyperPolar))
        curDescriptor.m_strSumMaxVectorHyperPolar = str(float(curDescriptor.m_strMaxNegaVectorHyperPolar) - float(curDescriptor.m_strMaxPosiVectorHyperPolar))
        curDescriptor.m_strDiffMinFlexibility = str(float(curDescriptor.m_strMinNegaFlexibility) + float(curDescriptor.m_strMinPosiFlexibility))
        curDescriptor.m_strSumMinFlexibility = str(float(curDescriptor.m_strMinNegaFlexibility) - float(curDescriptor.m_strMinPosiFlexibility))
        curDescriptor.m_strDiffTotalFlexibility = str(float(curDescriptor.m_strTotalNegaFlexibility) + float(curDescriptor.m_strTotalPosiFlexibility))
        curDescriptor.m_strSumTotalFlexibility = str(float(curDescriptor.m_strTotalNegaFlexibility) + float(curDescriptor.m_strTotalPosiFlexibility))
        curDescriptor.m_strDiffAverFlexibility = str(float(curDescriptor.m_strAverNegaFlexibility) - float(curDescriptor.m_strAverPosiFlexibility))
        curDescriptor.m_strSumAverFlexibility = str(float(curDescriptor.m_strAverNegaFlexibility) + float(curDescriptor.m_strAverPosiFlexibility))
        curDescriptor.m_strDiffMaxFlexibility = str(float(curDescriptor.m_strMaxNegaFlexibility) - float(curDescriptor.m_strMaxPosiFlexibility))
        curDescriptor.m_strSumMaxFlexibility = str(float(curDescriptor.m_strMaxNegaFlexibility) + float(curDescriptor.m_strMaxPosiFlexibility))

        #添加晶体相关描述符
        curDescriptor.m_strSpaceGroupName = curCrystal.m_strSpaceGroupName
        curDescriptor.m_strSpaceGroupID = curCrystal.m_strSpaceGroupID
        curDescriptor.m_strLengthMinDiff = curCrystal.m_strLengthMinDiff
        curDescriptor.m_strLengthMaxDiff = curCrystal.m_strLengthMaxDiff
        curDescriptor.m_strAngleMinDiff = curCrystal.m_strAngleMinDiff
        curDescriptor.m_strAngleMaxDiff = curCrystal.m_strAngleMaxDiff
        curDescriptor.m_strCellVolume = curCrystal.m_strVolume

        #处理SymmFunc
        dTotalSymmFunc = 0.0
        for curAtom in curCrystal.m_listNegaAtoms:
            dTotalSymmFunc = dTotalSymmFunc + float(curAtom.m_strSymFunc)
            if float(curDescriptor.m_strMinSymmFunc) > float(curAtom.m_strSymFunc):
                curDescriptor.m_strMinSymmFunc = str(float(curAtom.m_strSymFunc))
            if float(curDescriptor.m_strMaxSymmFunc) < float(curAtom.m_strSymFunc):
                curDescriptor.m_strMaxSymmFunc = str(float(curAtom.m_strSymFunc))
        if len(curCrystal.m_listNegaAtoms) < 1:
            curDescriptor.m_strAverSymmFunc = "0.0"
        else:
            curDescriptor.m_strAverSymmFunc = str(dTotalSymmFunc / len(curCrystal.m_listNegaAtoms))

        #处理Flexibility和键长分布
        nBondNum1T1_5 = 0
        nBondNum1_5T2 = 0
        nBondNum2T2_8 = 0
        dTotalFlexibility = 0.0
        for curBond in curCrystal.m_listNegaBonds:
            dTotalFlexibility = dTotalFlexibility + float(curBond.m_strFlexibility)
            dBondLen = float(curBond.m_strLength)
            if dBondLen >= 1.0 and dBondLen < 1.5:
                nBondNum1T1_5 = nBondNum1T1_5 + 1
            elif dBondLen >= 1.5 and dBondLen < 2.0:
                nBondNum1_5T2 = nBondNum1_5T2 + 1
            elif dBondLen >= 2.0 and dBondLen < 2.8:
                nBondNum2T2_8 = nBondNum2T2_8 + 1

        curDescriptor.m_strTotalFlexibility = str(dTotalFlexibility)
        if (len(curCrystal.m_listNegaBonds) == 0):
            curDescriptor.m_strAverFlexibility = "0"
        else:
            curDescriptor.m_strAverFlexibility = str(dTotalFlexibility / len(curCrystal.m_listNegaBonds))
        curDescriptor.m_strBondNum1T1_5 = str(nBondNum1T1_5)
        curDescriptor.m_strBondNum1_5T2 = str(nBondNum1_5T2)
        curDescriptor.m_strBondNum2T2_8 = str(nBondNum2T2_8)

        if curCrystal.m_strGapLevel in listGapLevelFilter:
            curDescriptor.m_strBandGap = curCrystal.m_strBandGap
            if curCrystal.m_strGapLevel == "GGA": #GGA计算的能隙要低50%左右，所以添加一个随机变量
                dCurBandGap = float(curCrystal.m_strBandGap)
                #dCurBandGap = dCurBandGap + dCurBandGap * ( random.random() * 0.5 ) #返回 [0.0, 1.0) 范围内的下一个随机浮点数
                #dCurBandGap = dCurBandGap + dCurBandGap * 0.5
                curDescriptor.m_strBandGap = str(dCurBandGap)
        if curCrystal.m_strBiRefLevel in listBiRefLevelFilter:
            nIndex = -1
            for curLength in curCrystal.m_listBeReflength:
                nIndex = nIndex + 1
                if float(curLength) >= dMinWaveLength and float(curLength) <= dMaxWaveLength:
                    curDescriptor.m_strBiRefracIndex = curCrystal.m_listBeRef[nIndex]
                    break
        if curCrystal.m_strMaxDijLevel in listBiRefLevelFilter:
            nIndex = -1
            for curLength in curCrystal.m_lisMaxDijlength:
                nIndex = nIndex + 1
                if float(curLength) >= dMinWaveLength and float(curLength) <= dMaxWaveLength:
                    curDescriptor.m_strMaxDij = curCrystal.m_listMaxDij[nIndex]
                    break

        listFeatureInfo.append(curDescriptor)

    return listFeatureInfo

#保存所有的信息
def saveAllFeatureInfo(listFeatureInfo):
    strFeatureFile = objects.BASISDIR + "feature-out.csv"
    file = open(strFeatureFile,"a+")
    file.write(objects.CDescriptor.joinTitle(","))
    for curFeature in listFeatureInfo:
        strCurLine = curFeature.joinToString(",")
        if "999.99" in strCurLine:
            print(curFeature.m_strID ,"存在999.99")
            continue

        file.write("\n")
        file.write(strCurLine)
    file.close()



