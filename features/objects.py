# 定义全局的常量、变量等
# 为了保持格式纯净，类中的所有成员变量均为string

import math
import method

BASISDIR = "D:/Data/" #数据存储的基本地址
#BASISDIR = "/public4/home/sc55809/mySoft/machineLearning/Data/" #数据存储的基本地址

#所有object的父类
class CBaseObject:
    def __init__(self,strID):
        self.m_strID=strID
        self.m_strName = ""

    @staticmethod
    def joinTitle(strBreak):
        strTitle = "strID" + strBreak
        strTitle = strTitle + "strName" + strBreak
        return  strTitle

    def joinToString(self, strBreak):
        strLine = self.m_strID + strBreak
        strLine = strLine + self.m_strName + strBreak
        #为了保险起见多添加一个分隔符，
        return strLine

# 计算基团Flexibility的bondValence值
class CBondValencePara(CBaseObject):
    def __init__(self,strID):
        CBaseObject.__init__(self,strID)
        self.m_strRo = ""
        self.m_strB = ""

    #检查当前Para是属于输入的两个原子
    def isBond(self,strElement1, strElement2):
        if strElement1 == self.m_strID and strElement2 == self.m_strName:
            return True
        elif strElement2 == self.m_strID and strElement1 == self.m_strName:
            return True

        return False



class CAtom(CBaseObject):
    def __init__(self,strID):
        CBaseObject.__init__(self,strID)
        self.m_strX = ""
        self.m_strY = ""
        self.m_strZ = ""
        self.m_strSymFunc = ""

    #计算当前键的自由度
    # listNearAtoms是在截断半径中的所有其他原子
    def calSymmetryFunction(self,listNearAtoms, listNearDist, dCutOffRadiu=6.0, dMiu = 0.0, dEta = 2.0):

        dTotalSymFunc = 0.0
        nAtomNum = len(listNearAtoms)
        for nIndex in range(0,nAtomNum):
            dCurSymFunc =  math.pow(listNearDist[nIndex] - dMiu, 2)
            dCurSymFunc = math.exp( -1 * dEta * dCurSymFunc)
            dCutOff = 0.5 * (( math.cos(math.pi * float(listNearDist[nIndex]) / dCutOffRadiu)) + 1)
            dCurSymFunc = dCurSymFunc * dCutOff
            dTotalSymFunc = dTotalSymFunc + dCurSymFunc

        #因为当前只是从一个单胞中获取的原子，所以只能求平均才能保证不随周围对称原子数改变而变化
        if nAtomNum < 1:
            self.m_strSymFunc = "0.0"
        else:
            self.m_strSymFunc = str(dTotalSymFunc/nAtomNum)


#键
class CBond(CBaseObject):
    def __init__(self,strID):
        CBaseObject.__init__(self,strID)
        self.m_CAtom1 = CAtom
        self.m_strLength = "999999999999.9999" #键长
        self.m_strFlexibility = "0.0"

    #计算当前键的自由度
    # dictAtomValEleNum原子名对应的价电子数
    def calFlexibility(self,listBondValParas,dictAtomValEleNum,dictAtomEleNegativity):
        #获得BondValence索引
        nParaIndex = -1
        for paraItem in listBondValParas:
            nParaIndex = nParaIndex + 1
            if paraItem.isBond(self.m_CAtom1.m_strName, self.m_CAtom2.m_strName):
                break

        if nParaIndex < 0 or nParaIndex >= len(listBondValParas):
            return 0.0

        #计算Flexibility
        dBondLengthDiff = float(listBondValParas[nParaIndex].m_strRo) - float(self.m_strLength)
        dBondValenece = math.exp(dBondLengthDiff/float(listBondValParas[nParaIndex].m_strB))
        dValenceEleNum = math.sqrt(dictAtomValEleNum[self.m_CAtom1.m_strName]) + math.sqrt(dictAtomValEleNum[self.m_CAtom2.m_strName])

        if float(self.m_strLength) - 0.0 < 1e-4:
            dValenceEleNum = 0
        else:
            dValenceEleNum = dValenceEleNum / float(self.m_strLength)

        dValenceEleNum = dValenceEleNum * dValenceEleNum
        dEleNegaDiff = math.fabs(dictAtomEleNegativity[self.m_CAtom1.m_strName] - dictAtomEleNegativity[self.m_CAtom2.m_strName])

        if dEleNegaDiff - 0.0 < 1e-4 or dValenceEleNum - 0.0 < 1e-4:
            dFlexibility = 0.0
        else:
            dFlexibility = dBondValenece / dValenceEleNum / dEleNegaDiff

        self.m_strFlexibility = str(dFlexibility)

#角
class CAngel(CBaseObject):
    def __init__(self,strID):
        CBaseObject.__init__(self,strID)
        self.m_strLength = "999999999999.9999" #键长


#字符串中括号类，m_strID表明左括号在字符串中的位置，m_strRightID表示右括号所在位置
# m_strName为括号符号，用于括号的匹配
class CBracketObject(CBaseObject):
    def __init__(self,strID):
        CBaseObject.__init__(self,strID)
        self.m_strRightID = ""
        self.m_strRightName = ""


#https://pubchem.ncbi.nlm.nih.gov/periodic-table/#view=list
class CElementObject(CBaseObject):
     def __init__(self,strID):
        CBaseObject.__init__(self,strID)
        self.m_strMass="0.0"
        self.m_strPeriodNum="0"
        self.m_strGroupNum="0"
        self.m_strRadui="0.0" #van der vaal, pm
        self.m_strEleNegativity="0.0" #pauling
        self.m_strSValenceEleNum="0"
        self.m_strPValenceEleNum="0"
        self.m_strDValenceEleNum="0"
        self.m_strFValenceEleNum="0"
        self.m_strSUnfilledState="0"
        self.m_strPUnfilledState="0"
        self.m_strDUnfilledState="0"
        self.m_strFUnfilledState="0"
        self.m_strionizationEner="0.0" #eV
        self.m_strEleAffinity="0.0" #eV
        self.m_strMeltingPoint="0.0" #K
        self.m_strBoilingPoint="0.0" #K
        self.m_strDensity="0.0" #g/cm-3

     def joinToString(self, strBreak):
        strLine = CBaseObject.joinToString(self,strBreak)
        strLine = strLine + self.m_strMass + strBreak
        strLine = strLine + self.m_strPeriodNum + strBreak
        strLine = strLine + self.m_strGroupNum + strBreak
        strLine = strLine + self.m_strRadui + strBreak
        strLine = strLine + self.m_strEleNegativity + strBreak
        strLine = strLine + self.m_strSValenceEleNum + strBreak
        strLine = strLine + self.m_strPValenceEleNum + strBreak
        strLine = strLine + self.m_strDValenceEleNum + strBreak
        strLine = strLine + self.m_strFValenceEleNum + strBreak
        strLine = strLine + self.m_strSUnfilledState + strBreak
        strLine = strLine + self.m_strPUnfilledState + strBreak
        strLine = strLine + self.m_strDUnfilledState + strBreak
        strLine = strLine + self.m_strFUnfilledState + strBreak
        strLine = strLine + self.m_strionizationEner + strBreak
        strLine = strLine + self.m_strEleAffinity + strBreak
        strLine = strLine + self.m_strMeltingPoint + strBreak
        strLine = strLine + self.m_strBoilingPoint + strBreak
        strLine = strLine + self.m_strDensity + strBreak
        return strLine


#基团信息,父类
class CGroupInfo(CBaseObject):
    #构造函数,只传入Group名称
    def __init__(self,strID):
        CBaseObject.__init__(self,strID)
        self.m_strHomo=""
        self.m_strLumo=""
        self.m_strLumoHomo=""
        self.m_strDipoleTotal="" #单位是Debye
        self.m_strQuadrupoleXX="" #单位是Debye/Ang
        self.m_strQuadrupoleYY=""
        self.m_strQuadrupoleZZ=""
        self.m_strQuadrupoleXY=""
        self.m_strQuadrupoleXZ=""
        self.m_strQuadrupoleYZ=""
        self.m_strAnisoQuadrupole="" #四极矩各向异性
        self.m_listPolarFreq = [] #外场频率,nm
        self.m_listIsoPolar= [] #极化率各向异性, 原子单位制au
        self.m_listAnisoPolar= [] #极化率各向异性
        self.m_listHyperPolarX = [] #第一超极化率X分量
        self.m_listHyperPolarY = [] #第一超极化率Y分量
        self.m_listHyperPolarZ = [] #第一超极化率
        self.m_listTotalHyperPolar= [] #总超极化率各向异性,
        self.m_listVectorHyperPolar= [] #总超极化率在偶极矩上的投影
        self.m_strAtomicNumber=""
        self.m_strCharge="" #Group电荷数
        self.m_strMultiplicity= "" #单电子个数
        self.m_strAverFlexibility = "" #平均共价键自由度
        self.m_strVolume = "" #体积，简单的通过(4/3)πr^3表示，r为平均键长。是否要添加范德华表面？？


    @staticmethod
    def joinTitle(strBreak):
        strLine = CBaseObject.joinTitle(strBreak)
        strLine = strLine + "strAtomicNumber" + strBreak
        strLine = strLine + "strCharge" + strBreak
        strLine = strLine + "strMultiplicity" + strBreak
        strLine = strLine + "strHomo" + strBreak
        strLine = strLine + "strLumo" + strBreak
        strLine = strLine + "strLumoHomo" + strBreak
        strLine = strLine + "strDipoleTotal" + strBreak
        strLine = strLine + "strQuadrupoleXX" + strBreak
        strLine = strLine + "strQuadrupoleYY" + strBreak
        strLine = strLine + "strQuadrupoleZZ" + strBreak
        strLine = strLine + "strQuadrupoleXY" + strBreak
        strLine = strLine + "strQuadrupoleXZ" + strBreak
        strLine = strLine + "strQuadrupoleYZ" + strBreak
        strLine = strLine + "strAnisoQuadrupole" + strBreak
        strLine = strLine + "listPolarFreq" + strBreak
        strLine = strLine + "listIsoPolar" + strBreak
        strLine = strLine + "listAnisoPolar" + strBreak
        strLine = strLine + "listHyperPolarX" + strBreak
        strLine = strLine + "listHyperPolarY" + strBreak
        strLine = strLine + "listHyperPolarZ" + strBreak
        strLine = strLine + "listTotalHyperPolar" + strBreak
        strLine = strLine + "listVectorHyperPolar" + strBreak
        strLine = strLine + "strAverFlexibility" + strBreak
        strLine = strLine + "strVolume" + strBreak
        return strLine

    #构建一个字符串
    def joinToString(self,strBreak):
        strLine = CBaseObject.joinToString(self,strBreak)
        strLine = strLine + self.m_strAtomicNumber + strBreak
        strLine = strLine + self.m_strCharge + strBreak
        strLine = strLine + self.m_strMultiplicity + strBreak
        strLine = strLine + self.m_strHomo + strBreak
        strLine = strLine + self.m_strLumo + strBreak
        strLine = strLine + self.m_strLumoHomo + strBreak
        strLine = strLine + self.m_strDipoleTotal + strBreak
        strLine = strLine + self.m_strQuadrupoleXX + strBreak
        strLine = strLine + self.m_strQuadrupoleYY + strBreak
        strLine = strLine + self.m_strQuadrupoleZZ + strBreak
        strLine = strLine + self.m_strQuadrupoleXY + strBreak
        strLine = strLine + self.m_strQuadrupoleXZ + strBreak
        strLine = strLine + self.m_strQuadrupoleYZ + strBreak
        strLine = strLine + self.m_strAnisoQuadrupole + strBreak
        strLine = strLine + ";".join(self.m_listPolarFreq) + strBreak
        strLine = strLine + ";".join(self.m_listIsoPolar) + strBreak
        strLine = strLine + ";".join(self.m_listAnisoPolar) + strBreak
        strLine = strLine + ";".join(self.m_listHyperPolarX) + strBreak
        strLine = strLine + ";".join(self.m_listHyperPolarY) + strBreak
        strLine = strLine + ";".join(self.m_listHyperPolarZ) + strBreak
        strLine = strLine + ";".join(self.m_listTotalHyperPolar) + strBreak
        strLine = strLine + ";".join(self.m_listVectorHyperPolar) + strBreak
        strLine = strLine + self.m_strAverFlexibility + strBreak
        strLine = strLine + self.m_strVolume + strBreak
        return strLine

    def calDerivatives(self):

        #计算四极矩各向异性
        dAnisoXY = float(self.m_strQuadrupoleXX) - float(self.m_strQuadrupoleYY)
        dAnisoXY = dAnisoXY*dAnisoXY
        dAnisoXZ = float(self.m_strQuadrupoleXX) - float(self.m_strQuadrupoleZZ)
        dAnisoXZ = dAnisoXZ*dAnisoXZ
        dAnisoYZ = float(self.m_strQuadrupoleYY) - float(self.m_strQuadrupoleZZ)
        dAnisoYZ = dAnisoYZ*dAnisoYZ
        dAnisoXYZ = float(self.m_strQuadrupoleXY) * float(self.m_strQuadrupoleXY)
        dAnisoXYZ = dAnisoXYZ + float(self.m_strQuadrupoleXZ) * float(self.m_strQuadrupoleXZ)
        dAnisoXYZ = dAnisoXYZ + float(self.m_strQuadrupoleYZ) * float(self.m_strQuadrupoleYZ)
        dAnisoXYZ = dAnisoXYZ * 6
        dAnisoXYZ = dAnisoXYZ + dAnisoXY + dAnisoXZ + dAnisoYZ
        self.m_strAnisoQuadrupole = str(math.sqrt(dAnisoXYZ/2))
        #计算总超极化率
        nIndex = -1
        for item in self.m_listPolarFreq:
            #计算总第一超极化率
            nIndex = nIndex + 1
            dBetaX = float(self.m_listHyperPolarX[nIndex])
            dBetaX = dBetaX * dBetaX
            dBetaY = float(self.m_listHyperPolarY[nIndex])
            dBetaY = dBetaY * dBetaY
            dBetaZ = float(self.m_listHyperPolarZ[nIndex])
            dBetaZ = dBetaZ * dBetaZ
            self.m_listTotalHyperPolar.append(str(math.sqrt(dBetaX + dBetaY + dBetaZ)))

    #检查格式是否符合要求
    def check(self):
        nLen = len(self.m_listPolarFreq)
        if len(self.m_listVectorHyperPolar) < 1:
            for nIndex in range(0,nLen):
                self.m_listVectorHyperPolar.append("0.0000")

#阳性基团信息
class CPostiveGroupInfo(CGroupInfo):
    #构造函数,只传入Group名称
    def __init__(self,strID):
        CGroupInfo.__init__(self,strID)

#阴性基团信息
class CNegativeGroupInfo(CGroupInfo):
    #构造函数,只传入Group名称
    def __init__(self,strID):
        CGroupInfo.__init__(self, strID)


#晶体信息
class CCrystalObject(CBaseObject):
    def __init__(self,strID):
        CBaseObject.__init__(self,strID)
        self.m_strSpaceGroupID=""
        self.m_strSpaceGroupName=""
        self.m_strLengthMinDiff = "" #晶胞参数最小的差值
        self.m_strLengthMaxDiff = "" #晶胞参数最大的差值
        self.m_strAngleMinDiff = "" #晶胞参数最小的差值
        self.m_strAngleMaxDiff = "" #晶胞参数最大的差值
        self.m_strLengthMin = "" #晶胞参数最小的值
        self.m_strLengthMiddle = "" #晶胞参数最小的值
        self.m_strLengthMax = "" #晶胞参数最大值
        self.m_strAngleMin = "" #晶胞参数最小的值
        self.m_strAngleMiddle = "" #晶胞参数最小的值
        self.m_strAngleMax = "" #晶胞参数最大的值
        self.m_strVolume="0.0"
        self.m_listPostiveGroups=[]
        self.m_listNegativeGroups=[]
        self.__m_strFormula="" #私有属性，必须通过成员函数设置，需要解析
        self.m_dictElements={} #键为元素符号，值为个数
        self.m_strGapLevel = "" #bandgap的水平
        self.m_strBandGap = ""
        self.m_strBiRefLevel = "" #BiRefracIndex的水平
        self.m_listBeReflength = [] #双折射率的所有波长
        self.m_listBeRef = [] #双折射率的值
        self.m_strMaxDijLevel = "" #二阶非线性系数的水平
        self.m_lisMaxDijlength = [] #二阶非线性系数的所有波长
        self.m_listMaxDij = [] #双折射率的值

        #不能直接从文本中得到
        self.m_listNegaAtoms = [] #晶体中所有原子
        self.m_listNegaBonds = [] #晶体中所有键
        self.m_listNegaAngels = []#晶体中所有角

    def joinToString(self, strBreak):
        strLine = CBaseObject.joinToString(self,strBreak)
        strLine = strLine + self.m_strSpaceGroupID + strBreak
        strLine = strLine + self.m_strSpaceGroupName + strBreak
        strLine = strLine + self.m_strLengthMinDiff + strBreak
        strLine = strLine + self.m_strLengthMaxDiff + strBreak
        strLine = strLine + self.m_strAngleMinDiff + strBreak
        strLine = strLine + self.m_strAngleMaxDiff + strBreak
        # strLine = strLine + self.m_strLengthMin + strBreak
        # strLine = strLine + self.m_strLengthMiddle + strBreak
        # strLine = strLine + self.m_strLengthMax + strBreak
        # strLine = strLine + self.m_strAngleMin + strBreak
        # strLine = strLine + self.m_strAngleMiddle + strBreak
        # strLine = strLine + self.m_strAngleMax + strBreak
        strLine = strLine + self.m_strVolume + strBreak
        #处理所有的阴阳基团
        strLine = strLine + ";".join(self.m_listPostiveGroups) + strBreak
        strLine = strLine + ";".join(self.m_listNegativeGroups) + strBreak

        strLine = strLine + self.__m_strFormula + strBreak
        strLine = strLine + self.m_strGapLevel + strBreak
        strLine = strLine + self.m_strBandGap + strBreak
        strLine = strLine + self.m_strBiRefLevel + strBreak
        strLine = strLine + ";".join(self.m_listBeReflength) + strBreak
        strLine = strLine + ";".join(self.m_listBeRef) + strBreak
        strLine = strLine + self.m_strMaxDijLevel + strBreak
        strLine = strLine + ";".join(self.m_lisMaxDijlength) + strBreak
        strLine = strLine + ";".join(self.m_listMaxDij) + strBreak
        return strLine

    def setFormula(self, strFormula):
        dictElements = {}
        self.__m_strFormula = method.getElementFromFormula(strFormula,dictElements)
        method.getElementFromString(self.__m_strFormula,self.m_dictElements) #晶体中所有的元素,key为元素，value为个数
    def getFormula(self):
        return self.__m_strFormula


#描述符
class CDescriptor(CBaseObject):
    def __init__(self,strID):
        CBaseObject.__init__(self, strID)
        #基本元素相关的，基本上都是取平均值
        self.m_strMinMass="9999999999999.99999"
        self.m_strAverMass="0.0" #所有元素的平均质量
        self.m_strMaxMass="-9999999999999.99999"
        self.m_strMinPeriodNum="999999"
        self.m_strAverPeriodNum="0"
        self.m_strMaxPeriodNum="-999999"
        self.m_strMinGroupNum="999999"
        self.m_strAverGroupNum="0"
        self.m_strMaxGroupNum="-999999"
        self.m_strMinRaduis="9999999999999.99999"
        self.m_strAverRaduis="0.0"
        self.m_strMaxRaduis="-9999999999999.99999"
        self.m_strMinEleNegativity="9999999999999.99999"
        self.m_strAverEleNegativity="0.0"
        self.m_strMaxEleNegativity="-9999999999999.99999"
        self.m_strAverSValenceEleNum="0"
        self.m_strAverPValenceEleNum="0"
        self.m_strAverDValenceEleNum="0"
        self.m_strAverFValenceEleNum="0"
        self.m_strAverSUnfilledState="0"
        self.m_strAverPUnfilledState="0"
        self.m_strAverDUnfilledState="0"
        self.m_strAverFUnfilledState="0"
        self.m_strMinIonizationEner="9999999999999.99999"
        self.m_strAverIonizationEner="0.0"
        self.m_strMaxIonizationEner="-9999999999999.99999"
        self.m_strMinEleAffinity="9999999999999.99999"
        self.m_strAverEleAffinity="0.0"
        self.m_strMaxEleAffinity="-9999999999999.99999"
        self.m_strMinMeltingPoint="9999999999999.99999"
        self.m_strAverMeltingPoint="0.0"
        self.m_strMaxMeltingPoint="-9999999999999.99999"
        self.m_strMinBoilingPoint="9999999999999.99999"
        self.m_strAverBoilingPoint="0.0"
        self.m_strMaxBoilingPoint="-9999999999999.99999"
        self.m_strMinDensity="9999999999999.99999"
        self.m_strAverDensity="0.0"
        self.m_strMaxDensity="-9999999999999.99999"

        #阴阳基团相关
        self.m_strNegaMinGap= "9999999999999.99999" #阴性基团中最小的Lumo-最大的Homo
        self.m_strNegaAverGap= "0.0" #阴性基团中平均的Lumo-最大的Homo
        self.m_strNegaMaxGap= "-9999999999999.99999" #阴性基团中最大的Lumo-最大的Homo
        self.m_strPosiMinGap="9999999999999.99999" #阴性基团中最小的Lumo-最大的Homo
        self.m_strPosiAverGap="0.0" #阴性基团中最小的Lumo-最大的Homo
        self.m_strPosiMaxGap="-9999999999999.99999" #阴性基团中最大的Lumo-最大的Homo

        self.m_strMinNegaDipoleTotal="9999999999999.99999" #阴性基团最小总偶极
        self.m_strAverNegaDipoleTotal="0.0" #阴性基团平均总偶极
        self.m_strMaxNegaDipoleTotal="-9999999999999.99999" #阴性基团最大总偶极
        self.m_strMinPosiDipoleTotal="9999999999999.99999"#阳性基团最小总偶极
        self.m_strAverPosiDipoleTotal="0.0"#阳性基团平均总偶极
        self.m_strMaxPosiDipoleTotal="-9999999999999.99999"#阳性基团最大总偶极

        self.m_strMinNegaAnisoQuadrupole="9999999999999.99999" #阴性基团最小四极矩各向异性
        self.m_strAverNegaAnisoQuadrupole="0.0" #阴性基团平均各向异性
        self.m_strMaxNegaAnisoQuadrupole="-9999999999999.99999" #阴性基团最大各向异性
        self.m_strMinPosiAnisoQuadrupole="9999999999999.99999" #阳性基团最小各向异性
        self.m_strAverPosiAnisoQuadrupole="0.0" #阳性基团平均各向异性
        self.m_strMaxPosiAnisoQuadrupole="-9999999999999.99999" #阳性基团最大各向异性

        self.m_strMinNegaAnisoPolar="9999999999999.99999" #阴性基团最小极化率各向异性
        self.m_strAverNegaAnisoPolar="0.0" #阴性基团平均极化率各向异性
        self.m_strMaxNegaAnisoPolar="-9999999999999.99999" #阴性基团最大极化率各向异性
        self.m_strMinPosiAnisoPolar="9999999999999.99999" #阳性基团最小极化率各向异性
        self.m_strAverPosiAnisoPolar="0.0" #阳性基团平均极化率各向异性
        self.m_strMaxPosiAnisoPolar="-9999999999999.99999" #阳性基团最大极化率各向异性

        self.m_strMinNegaTotalHyperPolar="9999999999999.99999" #阴性基团最小第一总超极化率
        self.m_strAverNegaTotalHyperPolar="0.0" #阴性基团平均第一总超极化率
        self.m_strMaxNegaTotalHyperPolar="-9999999999999.99999" #阴性基团最大第一总超极化率
        self.m_strMinPosiTotalHyperPolar="9999999999999.99999" #阳性基团最小第一总超极化率
        self.m_strAverPosiTotalHyperPolar="0.0" #阳性基团平均第一总超极化率
        self.m_strMaxPosiTotalHyperPolar="-9999999999999.99999" #阳性基团最大第一总超极化率

        self.m_strMinNegaVectorHyperPolar="9999999999999.99999" #阴性基团最小第一超极化率在偶极矩方向的分量
        self.m_strAverNegaVectorHyperPolar="0.0" #阴性基团平均第一超极化率在偶极矩方向的分量
        self.m_strMaxNegaVectorHyperPolar="-9999999999999.99999" #阴性基团最大第一超极化率在偶极矩方向的分量
        self.m_strMinPosiVectorHyperPolar="9999999999999.99999" #阳性基团最小第一超极化率在偶极矩方向的分量
        self.m_strAverPosiVectorHyperPolar="0.0" #阳性基团平均第一超极化率在偶极矩方向的分量
        self.m_strMaxPosiVectorHyperPolar="-9999999999999.99999" #阳性基团最大第一超极化率在偶极矩方向的分量

        self.m_strTotalNegaCharge="0"
        self.m_strTotalPosiCharge="0"
        self.m_strTotalNegaMultiplicity="0"
        self.m_strTotalPosiMultiplicity="0"

        self.m_strMinNegaFlexibility = "999999999999.99999" #基团的最小共价键自由度
        self.m_strTotalNegaFlexibility = "0.0" #总共价键自由度
        self.m_strAverNegaFlexibility = "0.0" #平均共价键自由度
        self.m_strMaxNegaFlexibility = "-999999999999.99999" #最大共价键自由度

        self.m_strMinPosiFlexibility = "999999999999.99999" #基团的最小共价键自由度
        self.m_strTotalPosiFlexibility = "0.0" #总共价键自由度
        self.m_strAverPosiFlexibility = "0.0" #平均共价键自由度
        self.m_strMaxPosiFlexibility = "-999999999999.99999" #最大共价键自由度

        #额外的基团结构信息,阴阳基团各自的
        self.m_strMaxNegaVolume = "-999999999.9999"
        self.m_strMaxPosiVolume = "-999999999.9999"
        self.m_strAverNegaCharge="0" #阴性基团的平均电荷
        self.m_strAverPosiCharge="0"
        self.m_strNegaMaxMinDiffGap= "" #阴性基团中最大gap与最小gap之差
        self.m_strPosiMaxMinDiffGap= "" #阴性基团中最大gap与最小gap之差
        self.m_strNegaMaxMinDiffDipoleTotal="" #阴性基团最小总偶极
        self.m_strPosiMaxMinDiffDipoleTotal="" #阴性基团最小总偶极
        self.m_strNegaMaxMinDiffAnisoQuadrupole=""
        self.m_strPosiMaxMinDiffAnisoQuadrupole=""
        self.m_strNegaMaxMinDiffAnisoPolar=""
        self.m_strPosiMaxMinDiffAnisoPolar=""
        self.m_strNegaMaxMinDiffTotalHyperPolar=""
        self.m_strPosiMaxMinDiffTotalHyperPolar=""
        self.m_strNegaMaxMinDiffVectorHyperPolar=""
        self.m_strPosiMaxMinDiffVectorHyperPolar=""
        self.m_strNegaMaxMinDiffFlexibility=""
        self.m_strPosiMaxMinDiffFlexibility=""
        self.m_strNegaMaxMinSumGap= "" #阴性基团中最大gap与最小gap之差
        self.m_strPosiMaxMinSumGap= "" #阴性基团中最大gap与最小gap之差
        self.m_strNegaMaxMinSumDipoleTotal="" #阴性基团最小总偶极
        self.m_strPosiMaxMinSumDipoleTotal="" #阴性基团最小总偶极
        self.m_strNegaMaxMinSumAnisoQuadrupole=""
        self.m_strPosiMaxMinSumAnisoQuadrupole=""
        self.m_strNegaMaxMinSumAnisoPolar=""
        self.m_strPosiMaxMinSumAnisoPolar=""
        self.m_strNegaMaxMinSumTotalHyperPolar=""
        self.m_strPosiMaxMinSumTotalHyperPolar=""
        self.m_strNegaMaxMinSumVectorHyperPolar=""
        self.m_strPosiMaxMinSumVectorHyperPolar=""
        self.m_strNegaMaxMinSumFlexibility=""
        self.m_strPosiMaxMinSumFlexibility=""

        #通过阴阳基团来构建虚拟晶体(VirtualCrystal,VC)的物理量,注意，此块处理的全部为阴阳基团之间的，不能为同类
        self.m_strVCGap="" #阴阳基团中gap的最小值
        self.m_strMaxVCDiffDipoleTotal="" #阴阳基团中DipoleTotal之差的最大值
        self.m_strMinVCDiffDipoleTotal="" #阴阳基团中DipoleTotal之差的最大值
        self.m_strMaxVCDiffAnisoQuadrupole="" #阴阳基团中Dipole之差的最大值
        self.m_strMinVCDiffAnisoQuadrupole="" #阴阳基团中Dipole之差的最大值
        self.m_strMaxVCDiffTotalHyperPolar="" #阴阳基团中TotalHyperPolar之差的最大值
        self.m_strMinVCDiffTotalHyperPolar="" #阴阳基团中TotalHyperPolar之差的最小值
        self.m_strMaxVCDiffVectorHyperPolar="" #阴阳基团中VectorHyperPolar之差的最大值
        self.m_strMinVCDiffVectorHyperPolar="" #阴阳基团中VectorHyperPolar之差的最小值
        self.m_strMaxVCDiffFlexibility="" #阴阳基团中VectorHyperPolar之差的最大值
        self.m_strMinVCDiffFlexibility="" #阴阳基团中VectorHyperPolar之差的最小值

        #四则运算后的特征，阴阳基团之间的
        self.m_strDiffMinGap = ""
        self.m_strSumMinGap = ""
        self.m_strDiffAverGap = ""
        self.m_strSumAverGap = ""
        self.m_strDiffMaxGap = ""
        self.m_strSumMaxGap = ""
        self.m_strDiffMinDipoleTotal = ""
        self.m_strSumMinDipoleTotal = ""
        self.m_strDiffAverDipoleTotal = ""
        self.m_strSumAverDipoleTotal = ""
        self.m_strDiffMaxDipoleTotal = ""
        self.m_strSumMaxDipoleTotal = ""
        self.m_strDiffMinAnisoQuadrupole = ""
        self.m_strSumMinAnisoQuadrupole = ""
        self.m_strDiffAverAnisoQuadrupole = ""
        self.m_strSumAverAnisoQuadrupole = ""
        self.m_strDiffMaxAnisoQuadrupole = ""
        self.m_strSumMaxAnisoQuadrupole = ""
        self.m_strDiffMinAnisoPolar = ""
        self.m_strSumMinAnisoPolar = ""
        self.m_strDiffAverAnisoPolar = ""
        self.m_strSumAverAnisoPolar = ""
        self.m_strDiffMaxAnisoPolar = ""
        self.m_strSumMaxAnisoPolar = ""
        self.m_strDiffMinTotalHyperPolar = ""
        self.m_strSumMinTotalHyperPolar = ""
        self.m_strDiffAverTotalHyperPolar = ""
        self.m_strSumAverTotalHyperPolar = ""
        self.m_strDiffMaxTotalHyperPolar = ""
        self.m_strSumMaxTotalHyperPolar = ""
        self.m_strDiffMinVectorHyperPolar = ""
        self.m_strSumMinVectorHyperPolar = ""
        self.m_strDiffAverVectorHyperPolar = ""
        self.m_strSumAverVectorHyperPolar = ""
        self.m_strDiffMaxVectorHyperPolar = ""
        self.m_strSumMaxVectorHyperPolar = ""
        self.m_strDiffMinFlexibility = ""
        self.m_strSumMinFlexibility = ""
        self.m_strDiffTotalFlexibility = ""
        self.m_strSumTotalFlexibility = ""
        self.m_strDiffAverFlexibility = ""
        self.m_strSumAverFlexibility = ""
        self.m_strDiffMaxFlexibility = ""
        self.m_strSumMaxFlexibility = "9999"

        #晶体相关，只计算但是不会输出，为了少改代码而己 20210122
        self.m_strSpaceGroupName=""
        self.m_strSpaceGroupID=""
        self.m_strLengthMinDiff = "" #晶胞参数最小的值
        self.m_strLengthMaxDiff = "" #晶胞参数最大的值
        self.m_strAngleMinDiff = "" #晶胞参数最小的值
        self.m_strAngleMaxDiff = "" #晶胞参数最大的值
        self.m_strCellVolume=""

        self.m_strMinSymmFunc="9999999999999.99999" #真实晶体中Group的最小原子的SymmFunc
        self.m_strAverSymmFunc="" #真实晶体中Group的总SymmFunc
        self.m_strMaxSymmFunc="-9999999999999.99999" #真实晶体中Group最大原小孩子的SymmFunc

        self.m_strAverFlexibility="" #真实晶体中Group的平均Flexibility
        self.m_strTotalFlexibility="" #真实晶体中Group的总Flexibility
        self.m_strBondNum1T1_5 = "0" #键长从1到1.5的个数
        self.m_strBondNum1_5T2 = "0" #键长从1.5到2的个数
        self.m_strBondNum2T2_8 = "0" #键长从2到2.8的个数

        self.m_strBandGap="" #与特定计算水平有关
        self.m_strBiRefracIndex=""#与特定计算水平和波长有关
        self.m_strMaxDij=""#与特定计算水平和波长有关

    @staticmethod
    def joinTitle(strBreak):
        strLine = CBaseObject.joinTitle(strBreak)
        strLine = strLine + "strMinMass" + strBreak
        strLine = strLine + "strAverMass" + strBreak
        strLine = strLine + "strMaxMass" + strBreak
        strLine = strLine + "strMinPeriodNum" + strBreak
        strLine = strLine + "strAverPeriodNum" + strBreak
        strLine = strLine + "strMaxPeriodNum" + strBreak
        strLine = strLine + "strMinGroupNum" + strBreak
        strLine = strLine + "strAverGroupNum" + strBreak
        strLine = strLine + "strMaxGroupNum" + strBreak
        strLine = strLine + "strMinRaduis" + strBreak
        strLine = strLine + "strAverRaduis" + strBreak
        strLine = strLine + "strMaxRaduis" + strBreak
        strLine = strLine + "strMinEleNegativity" + strBreak
        strLine = strLine + "strAverEleNegativity" + strBreak
        strLine = strLine + "strMaxEleNegativity" + strBreak
        strLine = strLine + "strAverSValenceEleNum" + strBreak
        strLine = strLine + "strAverPValenceEleNum" + strBreak
        strLine = strLine + "strAverDValenceEleNum" + strBreak
        strLine = strLine + "strAverFValenceEleNum" + strBreak
        strLine = strLine + "strAverSUnfilledState" + strBreak
        strLine = strLine + "strAverPUnfilledState" + strBreak
        strLine = strLine + "strAverDUnfilledState" + strBreak
        strLine = strLine + "strAverFUnfilledState" + strBreak
        strLine = strLine + "strMinIonizationEner" + strBreak
        strLine = strLine + "strAverIonizationEner" + strBreak
        strLine = strLine + "strMaxIonizationEner" + strBreak
        strLine = strLine + "strMinEleAffinity" + strBreak
        strLine = strLine + "strAverEleAffinity" + strBreak
        strLine = strLine + "strMaxEleAffinity" + strBreak
        strLine = strLine + "strMinMeltingPoint" + strBreak
        strLine = strLine + "strAverMeltingPoint" + strBreak
        strLine = strLine + "strMaxMeltingPoint" + strBreak
        strLine = strLine + "strMinBoilingPoint" + strBreak
        strLine = strLine + "strAverBoilingPoint" + strBreak
        strLine = strLine + "strMaxBoilingPoint" + strBreak
        strLine = strLine + "strMinDensity" + strBreak
        strLine = strLine + "strAverDensity" + strBreak
        strLine = strLine + "strMaxDensity" + strBreak

        #group基本的，属于同一类
        strLine = strLine + "strNegaMinGap" + strBreak
        strLine = strLine + "strNegaAverGap" + strBreak
        strLine = strLine + "strNegaMaxGap" + strBreak
        strLine = strLine + "strPosiMinGap" + strBreak
        strLine = strLine + "strPosiAverGap" + strBreak
        strLine = strLine + "strPosiMaxGap" + strBreak
        strLine = strLine + "strMinNegaDipoleTotal" + strBreak
        strLine = strLine + "strAverNegaDipoleTotal" + strBreak
        strLine = strLine + "strMaxNegaDipoleTotal" + strBreak
        strLine = strLine + "strMinPosiDipoleTotal" + strBreak
        strLine = strLine + "strAverPosiDipoleTotal" + strBreak
        strLine = strLine + "strMaxPosiDipoleTotal" + strBreak
        strLine = strLine + "strMinNegaAnisoQuadrupole" + strBreak
        strLine = strLine + "strAverNegaAnisoQuadrupole" + strBreak
        strLine = strLine + "strMaxNegaAnisoQuadrupole" + strBreak
        strLine = strLine + "strMinPosiAnisoQuadrupole" + strBreak
        strLine = strLine + "strAverPosiAnisoQuadrupole" + strBreak
        strLine = strLine + "strMaxPosiAnisoQuadrupole" + strBreak
        strLine = strLine + "strMinNegaAnisoPolar" + strBreak
        strLine = strLine + "strAverNegaAnisoPolar" + strBreak
        strLine = strLine + "strMaxNegaAnisoPolar" + strBreak
        strLine = strLine + "strMinPosiAnisoPolar" + strBreak
        strLine = strLine + "strAverPosiAnisoPolar" + strBreak
        strLine = strLine + "strMaxPosiAnisoPolar" + strBreak
        strLine = strLine + "strMinNegaTotalHyperPolar" + strBreak
        strLine = strLine + "strAverNegaTotalHyperPolar" + strBreak
        strLine = strLine + "strMaxNegaTotalHyperPolar" + strBreak
        strLine = strLine + "strMinPosiTotalHyperPolar" + strBreak
        strLine = strLine + "strAverPosiTotalHyperPolar" + strBreak
        strLine = strLine + "strMaxPosiTotalHyperPolar" + strBreak
        # strLine = strLine + "strMinNegaVectorHyperPolar" + strBreak
        # strLine = strLine + "strAverNegaVectorHyperPolar" + strBreak
        # strLine = strLine + "strMaxNegaVectorHyperPolar" + strBreak
        # strLine = strLine + "strMinPosiVectorHyperPolar" + strBreak
        # strLine = strLine + "strAverPosiVectorHyperPolar" + strBreak
        # strLine = strLine + "strMaxPosiVectorHyperPolar" + strBreak
        strLine = strLine + "strTotalNegaCharge" + strBreak
        strLine = strLine + "strTotalPosiCharge" + strBreak
        strLine = strLine + "strTotalNegaMultiplicity" + strBreak
        strLine = strLine + "strTotalPosiMultiplicity" + strBreak
        strLine = strLine + "strMinNegaFlexibility" + strBreak
        strLine = strLine + "strTotalNegaFlexibility" + strBreak
        strLine = strLine + "strAverNegaFlexibility" + strBreak
        strLine = strLine + "strMaxNegaFlexibility" + strBreak
        strLine = strLine + "strMinPosiFlexibility" + strBreak
        strLine = strLine + "strTotalPosiFlexibility" + strBreak
        strLine = strLine + "strAverPosiFlexibility" + strBreak
        strLine = strLine + "strMaxPosiFlexibility" + strBreak

        ##基团独自的扩展描述符
        # strLine = strLine + "strMaxNegaVolume" + strBreak
        # strLine = strLine + "strMaxPosiVolume" + strBreak
        strLine = strLine + "strAverNegaCharge" + strBreak
        # strLine = strLine + "strAverPosiCharge" + strBreak
        # strLine = strLine + "strNegaMaxMinDiffGap" + strBreak
        # strLine = strLine + "strPosiMaxMinDiffGap" + strBreak
        # strLine = strLine + "strNegaMaxMinDiffDipoleTotal" + strBreak
        # strLine = strLine + "strPosiMaxMinDiffDipoleTotal" + strBreak
        # strLine = strLine + "strNegaMaxMinDiffAnisoQuadrupole" + strBreak
        # strLine = strLine + "strPosiMaxMinDiffAnisoQuadrupole" + strBreak
        # strLine = strLine + "strNegaMaxMinDiffAnisoPolar" + strBreak
        # strLine = strLine + "strPosiMaxMinDiffAnisoPolar" + strBreak
        # strLine = strLine + "strNegaMaxMinDiffTotalHyperPolar" + strBreak
        # strLine = strLine + "strPosiMaxMinDiffTotalHyperPolar" + strBreak
        # strLine = strLine + "strNegaMaxMinDiffVectorHyperPolar" + strBreak
        # strLine = strLine + "strPosiMaxMinDiffVectorHyperPolar" + strBreak
        # strLine = strLine + "strNegaMaxMinDiffFlexibility" + strBreak
        # strLine = strLine + "strPosiMaxMinDiffFlexibility" + strBreak
        # strLine = strLine + "strNegaMaxMinSumGap" + strBreak
        # strLine = strLine + "strPosiMaxMinSumGap" + strBreak
        # strLine = strLine + "strNegaMaxMinSumDipoleTotal" + strBreak
        # strLine = strLine + "strPosiMaxMinSumDipoleTotal" + strBreak
        # strLine = strLine + "strNegaMaxMinSumAnisoQuadrupole" + strBreak
        # strLine = strLine + "strPosiMaxMinSumAnisoQuadrupole" + strBreak
        # strLine = strLine + "strNegaMaxMinSumAnisoPolar" + strBreak
        # strLine = strLine + "strPosiMaxMinSumAnisoPolar" + strBreak
        # strLine = strLine + "strNegaMaxMinSumTotalHyperPolar" + strBreak
        # strLine = strLine + "strPosiMaxMinSumTotalHyperPolar" + strBreak
        # strLine = strLine + "strNegaMaxMinSumVectorHyperPolar" + strBreak
        # strLine = strLine + "strPosiMaxMinSumVectorHyperPolar" + strBreak
        # strLine = strLine + "strNegaMaxMinSumFlexibility" + strBreak
        # strLine = strLine + "strPosiMaxMinSumFlexibility" + strBreak

        #虚拟晶体的信息
        # strLine = strLine + "strVCGap" + strBreak
        # strLine = strLine + "strMaxVCDiffDipoleTotal" + strBreak
        # strLine = strLine + "strMinVCDiffDipoleTotal" + strBreak
        # strLine = strLine + "strMaxVCDiffAnisoQuadrupole" + strBreak
        # strLine = strLine + "strMinVCDiffAnisoQuadrupole" + strBreak
        # strLine = strLine + "strMaxVCDiffTotalHyperPolar" + strBreak
        # strLine = strLine + "strMinVCDiffTotalHyperPolar" + strBreak
        # strLine = strLine + "strMaxVCDiffVectorHyperPolar" + strBreak
        # strLine = strLine + "strMinVCDiffVectorHyperPolar" + strBreak
        # strLine = strLine + "strMaxVCDiffFlexibility" + strBreak
        # strLine = strLine + "strMinVCDiffFlexibility" + strBreak


        #四则运算后的，所有尝试中均会拉低指标
        # strLine = strLine + "strDiffMinGap" + strBreak
        # strLine = strLine + "strSumMinGap" + strBreak
        # strLine = strLine + "strDiffAverGap" + strBreak
        # strLine = strLine + "strSumAverGap" + strBreak
        # strLine = strLine + "strDiffMaxGap" + strBreak
        # strLine = strLine + "strSumMaxGap" + strBreak
        # strLine = strLine + "strDiffMinDipoleTotal" + strBreak
        # strLine = strLine + "strSumMinDipoleTotal" + strBreak
        # strLine = strLine + "strDiffAverDipoleTotal" + strBreak
        # strLine = strLine + "strSumAverDipoleTotal" + strBreak
        # strLine = strLine + "strDiffMaxDipoleTotal" + strBreak
        # strLine = strLine + "strSumMaxDipoleTotal" + strBreak
        # strLine = strLine + "strDiffMinAnisoQuadrupole" + strBreak
        # strLine = strLine + "strSumMinAnisoQuadrupole" + strBreak
        # strLine = strLine + "strDiffAverAnisoQuadrupole" + strBreak
        # strLine = strLine + "strSumAverAnisoQuadrupole" + strBreak
        # strLine = strLine + "strDiffMaxAnisoQuadrupole" + strBreak
        # strLine = strLine + "strSumMaxAnisoQuadrupole" + strBreak
        # strLine = strLine + "strDiffMinAnisoPolar" + strBreak
        # strLine = strLine + "strSumMinAnisoPolar" + strBreak
        # strLine = strLine + "strDiffAverAnisoPolar" + strBreak
        # strLine = strLine + "strSumAverAnisoPolar" + strBreak
        # strLine = strLine + "strDiffMaxAnisoPolar" + strBreak
        # strLine = strLine + "strSumMaxAnisoPolar" + strBreak
        # strLine = strLine + "strDiffMinTotalHyperPolar" + strBreak
        # strLine = strLine + "strSumMinTotalHyperPolar" + strBreak
        # strLine = strLine + "strDiffAverTotalHyperPolar" + strBreak
        # strLine = strLine + "strSumAverTotalHyperPolar" + strBreak
        # strLine = strLine + "strDiffMaxTotalHyperPolar" + strBreak
        # strLine = strLine + "strSumMaxTotalHyperPolar" + strBreak
        # strLine = strLine + "strDiffMinVectorHyperPolar" + strBreak
        # strLine = strLine + "strSumMinVectorHyperPolar" + strBreak
        # strLine = strLine + "strDiffAverVectorHyperPolar" + strBreak
        # strLine = strLine + "strSumAverVectorHyperPolar" + strBreak
        # strLine = strLine + "strDiffMaxVectorHyperPolar" + strBreak
        # strLine = strLine + "strSumMaxVectorHyperPolar" + strBreak
        # strLine = strLine + "strDiffMinFlexibility" + strBreak
        # strLine = strLine + "strSumMinFlexibility" + strBreak
        # strLine = strLine + "strDiffTotalFlexibility" + strBreak
        # strLine = strLine + "strSumTotalFlexibility" + strBreak
        # strLine = strLine + "strDiffAverFlexibility" + strBreak
        # strLine = strLine + "strSumAverFlexibility" + strBreak
        # strLine = strLine + "strDiffMaxFlexibility" + strBreak
        # strLine = strLine + "strSumMaxFlexibility" + strBreak
        #
        #strLine = strLine + "strSpaceGroupName" + strBreak
        strLine = strLine + "strSpaceGroupID" + strBreak
        strLine = strLine + "strLengthMinDiff" + strBreak
        strLine = strLine + "strLengthMaxDiff" + strBreak
        strLine = strLine + "strAngleMinDiff" + strBreak
        strLine = strLine + "strAngleMaxDiff" + strBreak
        strLine = strLine + "strCellVolume" + strBreak
        #
        # strLine = strLine + "strMinSymmFunc" + strBreak
        # strLine = strLine + "strAverSymmFunc" + strBreak
        # strLine = strLine + "strMaxSymmFunc" + strBreak
        #
        # strLine = strLine + "strAverFlexibility" + strBreak
        # strLine = strLine + "strTotalFlexibility" + strBreak
        # strLine = strLine + "strBondNum1T1_5" + strBreak
        # strLine = strLine + "strBondNum1_5T2" + strBreak
        # strLine = strLine + "strBondNum2T2_8" + strBreak

        strLine = strLine + "strBandGap" + strBreak
        strLine = strLine + "strBiRefracIndex" + strBreak
        strLine = strLine + "strMaxDij" + strBreak
        return strLine

    def joinToString(self, strBreak):
        strLine = CBaseObject.joinToString(self,strBreak)
        strLine = strLine + self.m_strMinMass + strBreak
        strLine = strLine + self.m_strAverMass + strBreak
        strLine = strLine + self.m_strMaxMass + strBreak
        strLine = strLine + self.m_strMinPeriodNum + strBreak
        strLine = strLine + self.m_strAverPeriodNum + strBreak
        strLine = strLine + self.m_strMaxPeriodNum + strBreak
        strLine = strLine + self.m_strMinGroupNum + strBreak
        strLine = strLine + self.m_strAverGroupNum + strBreak
        strLine = strLine + self.m_strMaxGroupNum + strBreak
        strLine = strLine + self.m_strMinRaduis + strBreak
        strLine = strLine + self.m_strAverRaduis + strBreak
        strLine = strLine + self.m_strMaxRaduis + strBreak
        strLine = strLine + self.m_strMinEleNegativity + strBreak
        strLine = strLine + self.m_strAverEleNegativity + strBreak
        strLine = strLine + self.m_strMaxEleNegativity + strBreak
        strLine = strLine + self.m_strAverSValenceEleNum + strBreak
        strLine = strLine + self.m_strAverPValenceEleNum + strBreak
        strLine = strLine + self.m_strAverDValenceEleNum + strBreak
        strLine = strLine + self.m_strAverFValenceEleNum + strBreak
        strLine = strLine + self.m_strAverSUnfilledState + strBreak
        strLine = strLine + self.m_strAverPUnfilledState + strBreak
        strLine = strLine + self.m_strAverDUnfilledState + strBreak
        strLine = strLine + self.m_strAverFUnfilledState + strBreak
        strLine = strLine + self.m_strMinIonizationEner + strBreak
        strLine = strLine + self.m_strAverIonizationEner + strBreak
        strLine = strLine + self.m_strMaxIonizationEner + strBreak
        strLine = strLine + self.m_strMinEleAffinity + strBreak
        strLine = strLine + self.m_strAverEleAffinity + strBreak
        strLine = strLine + self.m_strMaxEleAffinity + strBreak
        strLine = strLine + self.m_strMinMeltingPoint + strBreak
        strLine = strLine + self.m_strAverMeltingPoint + strBreak
        strLine = strLine + self.m_strMaxMeltingPoint + strBreak
        strLine = strLine + self.m_strMinBoilingPoint + strBreak
        strLine = strLine + self.m_strAverBoilingPoint + strBreak
        strLine = strLine + self.m_strMaxBoilingPoint + strBreak
        strLine = strLine + self.m_strMinDensity + strBreak
        strLine = strLine + self.m_strAverDensity + strBreak
        strLine = strLine + self.m_strMaxDensity + strBreak

        #group同一类
        strLine = strLine + self.m_strNegaMinGap + strBreak
        strLine = strLine + self.m_strNegaAverGap + strBreak
        strLine = strLine + self.m_strNegaMaxGap + strBreak
        strLine = strLine + self.m_strPosiMinGap + strBreak
        strLine = strLine + self.m_strPosiAverGap + strBreak
        strLine = strLine + self.m_strPosiMaxGap + strBreak
        strLine = strLine + self.m_strMinNegaDipoleTotal + strBreak
        strLine = strLine + self.m_strAverNegaDipoleTotal + strBreak
        strLine = strLine + self.m_strMaxNegaDipoleTotal + strBreak
        strLine = strLine + self.m_strMinPosiDipoleTotal + strBreak
        strLine = strLine + self.m_strAverPosiDipoleTotal + strBreak
        strLine = strLine + self.m_strMaxPosiDipoleTotal + strBreak
        strLine = strLine + self.m_strMinNegaAnisoQuadrupole + strBreak
        strLine = strLine + self.m_strAverNegaAnisoQuadrupole + strBreak
        strLine = strLine + self.m_strMaxNegaAnisoQuadrupole + strBreak
        strLine = strLine + self.m_strMinPosiAnisoQuadrupole + strBreak
        strLine = strLine + self.m_strAverPosiAnisoQuadrupole + strBreak
        strLine = strLine + self.m_strMaxPosiAnisoQuadrupole + strBreak
        strLine = strLine + self.m_strMinNegaAnisoPolar + strBreak
        strLine = strLine + self.m_strAverNegaAnisoPolar + strBreak
        strLine = strLine + self.m_strMaxNegaAnisoPolar + strBreak
        strLine = strLine + self.m_strMinPosiAnisoPolar + strBreak
        strLine = strLine + self.m_strAverPosiAnisoPolar + strBreak
        strLine = strLine + self.m_strMaxPosiAnisoPolar + strBreak
        strLine = strLine + self.m_strMinNegaTotalHyperPolar + strBreak
        strLine = strLine + self.m_strAverNegaTotalHyperPolar + strBreak
        strLine = strLine + self.m_strMaxNegaTotalHyperPolar + strBreak
        strLine = strLine + self.m_strMinPosiTotalHyperPolar + strBreak
        strLine = strLine + self.m_strAverPosiTotalHyperPolar + strBreak
        strLine = strLine + self.m_strMaxPosiTotalHyperPolar + strBreak
        # strLine = strLine + self.m_strMinNegaVectorHyperPolar + strBreak
        # strLine = strLine + self.m_strAverNegaVectorHyperPolar + strBreak
        # strLine = strLine + self.m_strMaxNegaVectorHyperPolar + strBreak
        # strLine = strLine + self.m_strMinPosiVectorHyperPolar + strBreak
        # strLine = strLine + self.m_strAverPosiVectorHyperPolar + strBreak
        # strLine = strLine + self.m_strMaxPosiVectorHyperPolar + strBreak
        strLine = strLine + self.m_strTotalNegaCharge + strBreak
        strLine = strLine + self.m_strTotalPosiCharge + strBreak
        strLine = strLine + self.m_strTotalNegaMultiplicity + strBreak
        strLine = strLine + self.m_strTotalPosiMultiplicity + strBreak
        strLine = strLine + self.m_strMinNegaFlexibility + strBreak
        strLine = strLine + self.m_strTotalNegaFlexibility + strBreak
        strLine = strLine + self.m_strAverNegaFlexibility + strBreak
        strLine = strLine + self.m_strMaxNegaFlexibility + strBreak
        strLine = strLine + self.m_strMinPosiFlexibility + strBreak
        strLine = strLine + self.m_strTotalPosiFlexibility + strBreak
        strLine = strLine + self.m_strAverPosiFlexibility + strBreak
        strLine = strLine + self.m_strMaxPosiFlexibility + strBreak

        #阴阳基团各自的扩展描述符
        # strLine = strLine + self.m_strMaxNegaVolume + strBreak
        # strLine = strLine + self.m_strMaxPosiVolume + strBreak
        strLine = strLine + self.m_strAverNegaCharge + strBreak
        # strLine = strLine + self.m_strAverPosiCharge + strBreak
        # strLine = strLine + self.m_strNegaMaxMinDiffGap + strBreak
        # strLine = strLine + self.m_strPosiMaxMinDiffGap + strBreak
        # strLine = strLine + self.m_strNegaMaxMinDiffDipoleTotal + strBreak
        # strLine = strLine + self.m_strPosiMaxMinDiffDipoleTotal + strBreak
        # strLine = strLine + self.m_strNegaMaxMinDiffAnisoQuadrupole + strBreak
        # strLine = strLine + self.m_strPosiMaxMinDiffAnisoQuadrupole + strBreak
        # strLine = strLine + self.m_strNegaMaxMinDiffAnisoPolar + strBreak
        # strLine = strLine + self.m_strPosiMaxMinDiffAnisoPolar + strBreak
        # strLine = strLine + self.m_strNegaMaxMinDiffTotalHyperPolar + strBreak
        # strLine = strLine + self.m_strPosiMaxMinDiffTotalHyperPolar + strBreak
        # strLine = strLine + self.m_strNegaMaxMinDiffVectorHyperPolar + strBreak
        # strLine = strLine + self.m_strPosiMaxMinDiffVectorHyperPolar + strBreak
        # strLine = strLine + self.m_strNegaMaxMinDiffFlexibility + strBreak
        # strLine = strLine + self.m_strPosiMaxMinDiffFlexibility + strBreak
        # strLine = strLine + self.m_strNegaMaxMinSumGap + strBreak
        # strLine = strLine + self.m_strPosiMaxMinSumGap + strBreak
        # strLine = strLine + self.m_strNegaMaxMinSumDipoleTotal + strBreak
        # strLine = strLine + self.m_strPosiMaxMinSumDipoleTotal + strBreak
        # strLine = strLine + self.m_strNegaMaxMinSumAnisoQuadrupole + strBreak
        # strLine = strLine + self.m_strPosiMaxMinSumAnisoQuadrupole + strBreak
        # strLine = strLine + self.m_strNegaMaxMinSumAnisoPolar + strBreak
        # strLine = strLine + self.m_strPosiMaxMinSumAnisoPolar + strBreak
        # strLine = strLine + self.m_strNegaMaxMinSumTotalHyperPolar + strBreak
        # strLine = strLine + self.m_strPosiMaxMinSumTotalHyperPolar + strBreak
        # strLine = strLine + self.m_strNegaMaxMinSumVectorHyperPolar + strBreak
        # strLine = strLine + self.m_strPosiMaxMinSumVectorHyperPolar + strBreak
        # strLine = strLine + self.m_strNegaMaxMinSumFlexibility + strBreak
        # strLine = strLine + self.m_strPosiMaxMinSumFlexibility + strBreak

        #虚拟晶体，不同类
        # strLine = strLine + self.m_strVCGap + strBreak
        # strLine = strLine + self.m_strMaxVCDiffDipoleTotal + strBreak
        # strLine = strLine + self.m_strMinVCDiffDipoleTotal + strBreak
        # strLine = strLine + self.m_strMaxVCDiffAnisoQuadrupole + strBreak
        # strLine = strLine + self.m_strMinVCDiffAnisoQuadrupole + strBreak
        # strLine = strLine + self.m_strMaxVCDiffTotalHyperPolar + strBreak
        # strLine = strLine + self.m_strMinVCDiffTotalHyperPolar + strBreak
        # strLine = strLine + self.m_strMaxVCDiffVectorHyperPolar + strBreak
        # strLine = strLine + self.m_strMinVCDiffVectorHyperPolar + strBreak
        # strLine = strLine + self.m_strMaxVCDiffFlexibility + strBreak
        # strLine = strLine + self.m_strMinVCDiffFlexibility + strBreak

        #四则运算后
        # strLine = strLine + self.m_strDiffMinGap + strBreak
        # strLine = strLine + self.m_strSumMinGap + strBreak
        # strLine = strLine + self.m_strDiffAverGap + strBreak
        # strLine = strLine + self.m_strSumAverGap + strBreak
        # strLine = strLine + self.m_strDiffMaxGap + strBreak
        # strLine = strLine + self.m_strSumMaxGap + strBreak
        # strLine = strLine + self.m_strDiffMinDipoleTotal + strBreak
        # strLine = strLine + self.m_strSumMinDipoleTotal + strBreak
        # strLine = strLine + self.m_strDiffAverDipoleTotal + strBreak
        # strLine = strLine + self.m_strSumAverDipoleTotal + strBreak
        # strLine = strLine + self.m_strDiffMaxDipoleTotal + strBreak
        # strLine = strLine + self.m_strSumMaxDipoleTotal + strBreak
        # strLine = strLine + self.m_strDiffMinAnisoQuadrupole + strBreak
        # strLine = strLine + self.m_strSumMinAnisoQuadrupole + strBreak
        # strLine = strLine + self.m_strDiffAverAnisoQuadrupole + strBreak
        # strLine = strLine + self.m_strSumAverAnisoQuadrupole + strBreak
        # strLine = strLine + self.m_strDiffMaxAnisoQuadrupole + strBreak
        # strLine = strLine + self.m_strSumMaxAnisoQuadrupole + strBreak
        # strLine = strLine + self.m_strDiffMinAnisoPolar + strBreak
        # strLine = strLine + self.m_strSumMinAnisoPolar + strBreak
        # strLine = strLine + self.m_strDiffAverAnisoPolar + strBreak
        # strLine = strLine + self.m_strSumAverAnisoPolar + strBreak
        # strLine = strLine + self.m_strDiffMaxAnisoPolar + strBreak
        # strLine = strLine + self.m_strSumMaxAnisoPolar + strBreak
        # strLine = strLine + self.m_strDiffMinTotalHyperPolar + strBreak
        # strLine = strLine + self.m_strSumMinTotalHyperPolar + strBreak
        # strLine = strLine + self.m_strDiffAverTotalHyperPolar + strBreak
        # strLine = strLine + self.m_strSumAverTotalHyperPolar + strBreak
        # strLine = strLine + self.m_strDiffMaxTotalHyperPolar + strBreak
        # strLine = strLine + self.m_strSumMaxTotalHyperPolar + strBreak
        # strLine = strLine + self.m_strDiffMinVectorHyperPolar + strBreak
        # strLine = strLine + self.m_strSumMinVectorHyperPolar + strBreak
        # strLine = strLine + self.m_strDiffAverVectorHyperPolar + strBreak
        # strLine = strLine + self.m_strSumAverVectorHyperPolar + strBreak
        # strLine = strLine + self.m_strDiffMaxVectorHyperPolar + strBreak
        # strLine = strLine + self.m_strSumMaxVectorHyperPolar + strBreak
        # strLine = strLine + self.m_strDiffMinFlexibility + strBreak
        # strLine = strLine + self.m_strSumMinFlexibility + strBreak
        # strLine = strLine + self.m_strDiffTotalFlexibility + strBreak
        # strLine = strLine + self.m_strSumTotalFlexibility + strBreak
        # strLine = strLine + self.m_strDiffAverFlexibility + strBreak
        # strLine = strLine + self.m_strSumAverFlexibility + strBreak
        # strLine = strLine + self.m_strDiffMaxFlexibility + strBreak
        # strLine = strLine + self.m_strSumMaxFlexibility + strBreak
        #
        # strLine = strLine + self.m_strSpaceGroupName + strBreak
        strLine = strLine + self.m_strSpaceGroupID + strBreak
        strLine = strLine + self.m_strLengthMinDiff + strBreak
        strLine = strLine + self.m_strLengthMaxDiff + strBreak
        strLine = strLine + self.m_strAngleMinDiff + strBreak
        strLine = strLine + self.m_strAngleMaxDiff + strBreak
        strLine = strLine + self.m_strCellVolume + strBreak
        #
        # strLine = strLine + self.m_strMinSymmFunc + strBreak
        # strLine = strLine + self.m_strAverSymmFunc + strBreak
        # strLine = strLine + self.m_strMaxSymmFunc + strBreak
        #
        # strLine = strLine + self.m_strAverFlexibility + strBreak
        # strLine = strLine + self.m_strTotalFlexibility + strBreak
        # strLine = strLine + self.m_strBondNum1T1_5 + strBreak
        # strLine = strLine + self.m_strBondNum1_5T2 + strBreak
        # strLine = strLine + self.m_strBondNum2T2_8 + strBreak

        strLine = strLine + self.m_strBandGap + strBreak
        strLine = strLine + self.m_strBiRefracIndex + strBreak
        strLine = strLine + self.m_strMaxDij + strBreak
        return strLine

#vdW表面ESP的描述符
class CvdWESPDescriptor(CBaseObject):
    def __init__(self,strID):
        CBaseObject.__init__(self, strID)

        self.m_strTotalVolume = "999999.99999"
        self.m_strTotalArea = "999999.99999"
        self.m_strPositiveArea = "999999.99999"
        self.m_strNegativeArea = "999999.99999"
        self.m_strPolarArea = "999999.99999"
        self.m_strNonPolarArea = "999999.99999"

        self.m_strGlobalMin = "-999999.99999"
        self.m_strGlobalMax = "999999.99999"
        self.m_strTotalVar = "999999.99999"
        self.m_strMinimaNum = "0" #小于某个值的极小点个数
        self.m_strMaximaNum = "0" #大于某个值的极大点个数
        self.m_strMinimaSum = "0.0" #小于某个值的极小点的和
        self.m_strMaximaSum = "0.0" #大于某个值的极大点的和

    @staticmethod
    def joinTitle(strBreak):
        strLine = CBaseObject.joinTitle(strBreak)
        # strLine = strLine + "strvdWTotalVolume" + strBreak
        # strLine = strLine + "strvdWTotalArea" + strBreak
        # strLine = strLine + "strvdWPositiveArea" + strBreak
        # strLine = strLine + "strvdWNegativeArea" + strBreak
        # strLine = strLine + "strvdWPolarArea" + strBreak
        # strLine = strLine + "strvdWNonPolarArea" + strBreak

        strLine = strLine + "strvdWGlobalMin" + strBreak
        strLine = strLine + "strvdWGlobalMax" + strBreak
        strLine = strLine + "strvdWTotalVar" + strBreak
        strLine = strLine + "strvdWMinimaNum" + strBreak
        strLine = strLine + "strvdWMaximaNum" + strBreak
        strLine = strLine + "strvdWMinimaSum" + strBreak
        strLine = strLine + "strvdWMaximaSum" + strBreak
        return strLine

    def joinToString(self, strBreak):
        strLine = CBaseObject.joinToString(self,strBreak)

        # strLine = strLine + self.m_strTotalVolume + strBreak
        # strLine = strLine + self.m_strTotalArea + strBreak
        # strLine = strLine + self.m_strPositiveArea + strBreak
        # strLine = strLine + self.m_strNegativeArea + strBreak
        # strLine = strLine + self.m_strPolarArea + strBreak
        # strLine = strLine + self.m_strNonPolarArea + strBreak

        strLine = strLine + self.m_strGlobalMin + strBreak
        strLine = strLine + self.m_strGlobalMax + strBreak
        strLine = strLine + self.m_strTotalVar + strBreak
        strLine = strLine + self.m_strMinimaNum + strBreak
        strLine = strLine + self.m_strMaximaNum + strBreak
        strLine = strLine + self.m_strMinimaSum + strBreak
        strLine = strLine + self.m_strMaximaSum + strBreak

        return strLine

#vdW表面ESP的QSAR相关描述符, J. Phys. Chem. A 1999, 103, 1853-1856
class CvdWESPQSARDescriptor(CBaseObject):
    def __init__(self,strID):
        CBaseObject.__init__(self, strID)

        self.m_strNegativeArea = "999999.99999" #AS-
        self.m_strNegativeAver = "999999.99999" #VS-
        self.m_strGlobalMin = "-999999.99999" #VS,min
        self.m_strGlobalMax = "999999.99999" #VS,max

    @staticmethod
    def joinTitle(strBreak):
        strLine = CBaseObject.joinTitle(strBreak)
        strLine = strLine + "strGlobalMax" + strBreak
        strLine = strLine + "strGlobalMin" + strBreak
        strLine = strLine + "strNegativeArea" + strBreak
        strLine = strLine + "strNegativeAver" + strBreak
        return strLine

    def joinToString(self, strBreak):
        strLine = CBaseObject.joinToString(self,strBreak)
        strLine = strLine + self.m_strGlobalMax + strBreak
        strLine = strLine + self.m_strGlobalMin + strBreak
        strLine = strLine + self.m_strNegativeArea + strBreak
        strLine = strLine + self.m_strNegativeAver + strBreak

        return strLine

#小分子溶剂化能描述符
class CSolvDescriptor(CBaseObject):
    def __init__(self,strID):
        CBaseObject.__init__(self, strID)

        self.m_listEleRatio=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0] #不同原子所占比例，H,C,N,O,F,P,S,Cl,Br,I

        self.m_strSCFEner = "999999.99999" #SCF Done:  E(RB3LYP)
        self.m_strSCFKE = "999999.99999" #KE= 2.660105974559D+02 PE=-1.728616267669D+03
        self.m_strSCFPE = "999999.99999"
        self.m_strSCFEE = "999999.99999"

        self.m_strDipole="999999.99999"
        self.m_strAnisoPolar="999999.99999"

        self.m_listHarmonicFreq=[] #Harmonic frequencies，最高的3个
        self.m_strZeroPointEner="999999.999" # Zero-point vibrational energy

        self.m_nShellNum = 20 #总共20个Shell
        self.m_listMaxESPShells=[]
        self.m_listMinESPShells=[]
        self.m_listAverESPShells=[]
        self.m_listVarESPShells=[] #方差
        # self.m_listTMomESPShells=[] #三阶中心距
        # self.m_listFMomESPShells=[] #四阶中心距
        self.m_listSkewESPShells=[] #偏度
        self.m_listPeakESPShells=[] #峰度
        # self.m_listTotDenESPShells=[] #在球面上的整体密度
        self.m_listSmallESPRatioShells=[] #在球面上的ESP小于0.001的点的比例
        self.m_listPowerSpectrumShells=[] #在当前球面上的平方谱，二维数组
        self.m_strVdWSurfaceArea="999999.99999"

    @staticmethod
    def joinTitle(strBreak):
        strLine = CBaseObject.joinTitle(strBreak)
        strLine = strLine + "strHRatio" + strBreak
        strLine = strLine + "strCRatio" + strBreak
        strLine = strLine + "strNRatio" + strBreak
        strLine = strLine + "strORatio" + strBreak
        strLine = strLine + "strFRatio" + strBreak
        strLine = strLine + "strPRatio" + strBreak
        strLine = strLine + "strSRatio" + strBreak
        strLine = strLine + "strClRatio" + strBreak
        strLine = strLine + "strBrRatio" + strBreak
        strLine = strLine + "strIRatio" + strBreak

        strLine = strLine + "strSCFEner" + strBreak
        strLine = strLine + "strSCFKE" + strBreak
        strLine = strLine + "strSCFPE" + strBreak
        strLine = strLine + "strSCFEE" + strBreak

        strLine = strLine + "strDipole" + strBreak
        strLine = strLine + "strAnisoPolar" + strBreak

        for nIndex in range(1, 4):
            strLine = strLine + "strHarmFreq%d"%(nIndex) + strBreak

        strLine = strLine + "strZeroPointEner" + strBreak

        for nIndex in range(1, 21):
            strLine = strLine + "strMaxESPShell%d"%(nIndex) + strBreak
            strLine = strLine + "strMinESPShell%d"%(nIndex) + strBreak
            strLine = strLine + "strAverESPShell%d"%(nIndex) + strBreak
            strLine = strLine + "strVarESPShell%d"%(nIndex) + strBreak
            # strLine = strLine + "strTMomESPShell%d"%(nIndex) + strBreak
            # strLine = strLine + "strFMomESPShell%d"%(nIndex) + strBreak
            strLine = strLine + "strSkewESPShell%d"%(nIndex) + strBreak
            strLine = strLine + "strPeakESPShell%d"%(nIndex) + strBreak
            # strLine = strLine + "strTotDenESPShell%d"%(nIndex) + strBreak
            strLine = strLine + "strSmallESPRatioShell%d"%(nIndex) + strBreak

            #20为球谐函数展开时的lmax
            lmax = 7
            for nSubIndex in range(0, lmax+1):
                strLine = strLine + "strPowerSpectrum%dShells%d"%(nSubIndex,nIndex) + strBreak

        strLine = strLine + "strSASA" + strBreak
        strLine = strLine + "strUSRScore" + strBreak
        return strLine

    def joinToString(self, strBreak):
        strLine = CBaseObject.joinToString(self,strBreak)

        for nIndex in range(0, len(self.m_listEleRatio)):
            strLine = strLine + self.m_listEleRatio[nIndex] + strBreak

        strLine = strLine + self.m_strSCFEner + strBreak
        strLine = strLine + self.m_strSCFKE + strBreak
        strLine = strLine + self.m_strSCFPE + strBreak
        strLine = strLine + self.m_strSCFEE + strBreak

        strLine = strLine + self.m_strDipole + strBreak
        strLine = strLine + self.m_strAnisoPolar + strBreak

        nFreqNum = len(self.m_listHarmonicFreq)
        if nFreqNum < 3:
            print("频率个数小于3")

        for nIndex in range(nFreqNum-1, nFreqNum-4, -1):
            if nFreqNum < 3:
                strLine = strLine + "999.999" + strBreak
            else:
                strLine = strLine + self.m_listHarmonicFreq[nIndex] + strBreak

        strLine = strLine + self.m_strZeroPointEner + strBreak

        for nIndex in range(0, self.m_nShellNum):
            strLine = strLine + self.m_listMaxESPShells[nIndex] + strBreak
            strLine = strLine + self.m_listMinESPShells[nIndex] + strBreak
            strLine = strLine + self.m_listAverESPShells[nIndex] + strBreak
            strLine = strLine + self.m_listVarESPShells[nIndex] + strBreak
            # strLine = strLine + self.m_listTMomESPShells[nIndex] + strBreak
            # strLine = strLine + self.m_listFMomESPShells[nIndex] + strBreak
            strLine = strLine + self.m_listSkewESPShells[nIndex] + strBreak
            strLine = strLine + self.m_listPeakESPShells[nIndex] + strBreak
            # strLine = strLine + self.m_listTotDenESPShells[nIndex] + strBreak
            strLine = strLine + self.m_listSmallESPRatioShells[nIndex] + strBreak

            #20为球谐函数展开时的lmax
            for nSubIndex in range(0, len(self.m_listPowerSpectrumShells[nIndex])):
                strLine = strLine + self.m_listPowerSpectrumShells[nIndex][nSubIndex] + strBreak

        strLine = strLine + self.m_strVdWSurfaceArea + strBreak
        strLine = strLine + self.m_strVdWSurfaceArea + strBreak
        return strLine


#训练的测试集结果
#strID表示二值化拆分点
#strName表示不平衡分类时处理方程
#m_strThreshold表示进行特征筛选时的准确率门槛
class CPredictResultObject(CBaseObject):
    def __init__(self,strID):
        CBaseObject.__init__(self, strID)
        self.m_strThreshold = ""
        self.m_arrRawTestLabel = ""
        self.m_arrPredictLabel = ""
        self.m_arrPredictProba = ""

#--------------------------------------------------------------------------------------------------#
#qm4d中所有的原子信息，可以从pdb，psf，sxyz中获取信息
class CQm4dAtom(CBaseObject):
    def __init__(self,strID):
        CBaseObject.__init__(self,strID)
        self.m_strStart = "ATOM"
        self.m_strChain = ""
        self.m_strResidueID = ""
        self.m_strResidueName = ""
        self.m_strPrmType = "" #参数文件中的类型
        self.m_strX = ""
        self.m_strY = ""
        self.m_strZ = ""
        self.m_strCharge = ""
        self.m_strOccupancy = "1.00"
        self.m_strTempFactor = "0.00"
        self.m_strElement = ""
        self.m_strAtomMass = ""
        self.m_strTinkerClass = "" #Tinker中的classID

    def addPDBInfo(self, strInputLine):
        strLine = strInputLine.strip()
        strLine = strLine.strip("\n")
        strLine = ' '.join(strLine.split()) #处理多个空格间隔的情况
        listItems = strLine.split()
        if (len(listItems) < 11): #change
            print("Error：%s"%(strLine))
            return

        self.m_strStart = listItems[0]
        self.m_strID = listItems[1]
        self.m_strName = listItems[2]
        self.m_strResidueName = listItems[3][0:3]
        self.m_strChain = listItems[4][0]
        self.m_strResidueID = listItems[5]
        if self.m_strResidueID.find(".") > -1:
            print("注意，列数可能不正确")
        self.m_strX = listItems[6]
        self.m_strY = listItems[7]
        self.m_strZ = listItems[8]
        self.m_strOccupancy = "1.00" #listItems[9]
        self.m_strTempFactor = "0.00" #listItems[10]
        #self.m_strElement = listItems[11]

    #注意，忽略了链的结尾
    def joinToPDBString(self):
        strLine = "%-6s"%(self.m_strStart) #1-6
        strLine += "%5s"%(self.m_strID) #7-11
        strLine += " " #12
        strLine += "%-4s"%(self.m_strName) #13-16
        strLine += " " #17
        strLine += "%3s"%(self.m_strResidueName[0:3]) #18-20
        strLine += " "#21
        strLine += "%1s"%(self.m_strChain) #22
        strLine += "%4s"%(self.m_strResidueID[0:3]) #23-26
        strLine += "%4s"%("    ")#27-30
        strLine += "%7.3f "%(float(self.m_strX)) #31-38
        strLine += "%7.3f "%(float(self.m_strY)) #39-46
        strLine += "%7.3f "%(float(self.m_strZ)) #47-54
        strLine += "%6s"%(self.m_strOccupancy) #55-60
        strLine += "%6s"%(self.m_strTempFactor) #61-66
        strLine += "%10s"%("          ")#67-76
        strLine += "%2s"%(self.m_strElement) #77-78
        strLine += "%2s"%("  \n") #78-79
        return strLine

    def addPsfInfo(self, strInputLine):
        strLine = strInputLine.strip()
        strLine = strLine.strip("\n")
        strLine = ' '.join(strLine.split()) #处理多个空格间隔的情况
        listItems = strLine.split()
        if (len(listItems) < 8):
            print("Error input in addPsfInfo!")
            return

        self.m_strID = listItems[0]
        self.m_strChain = listItems[1][0]
        self.m_strResidueID = listItems[2]
        self.m_strResidueName = listItems[3]
        self.m_strName = listItems[4]
        self.m_strPrmType = listItems[5]
        self.m_strCharge = listItems[6]
        self.m_strAtomMass = listItems[7]

    def joinToPsfString(self):
        strLine = "%8s%2s    %-5s%-5s%-5s%-5s %-15s %-18s0\n"\
                      %(self.m_strID,self.m_strChain,self.m_strResidueID,self.m_strResidueName,
                        self.m_strName,self.m_strPrmType,self.m_strCharge,self.m_strAtomMass)
        return strLine

    def addSxyzInfo(self, strInputLine):
        strLine = strInputLine.strip()
        strLine = strLine.strip("\n")
        strLine = ' '.join(strLine.split()) #处理多个空格间隔的情况
        listItems = strLine.split()
        if (len(listItems) < 8):
            print("Error input in addSxyzInfo!")
            return

        self.m_strStart = listItems[0]
        self.m_strID = listItems[1]
        self.m_strName = listItems[2]
        self.m_strResidueName = listItems[3]
        self.m_strResidueID = listItems[4]
        self.m_strX = listItems[5]
        self.m_strY = listItems[6]
        self.m_strZ = listItems[7]
        # self.m_strX = "%.4f"%(float(listItems[5]))
        # self.m_strY = "%.4f"%(float(listItems[6]))
        # self.m_strZ = "%.4f"%(float(listItems[7]))

    def joinToSxyzString(self):
        strLine = "%-6s"%(self.m_strStart) #
        strLine += "%5s"%(self.m_strID) #
        strLine += " " #
        strLine += "%-4s"%(self.m_strName) #
        strLine += " " #
        strLine += "%4s"%(self.m_strResidueName) #
        strLine += "    "#
        strLine += "%5s"%(self.m_strResidueID) #
        strLine += "  " #
        strLine += "%12s"%(self.m_strX) #
        strLine += "  "
        strLine += "%12s"%(self.m_strY) #
        strLine += "  "
        strLine += "%12s"%(self.m_strZ) #
        strLine += "\n"
        return strLine

    def addXyzInfo(self, strInputLine):
        strLine = strInputLine.strip()
        strLine = strLine.strip("\n")
        strLine = ' '.join(strLine.split()) #处理多个空格间隔的情况
        listItems = strLine.split()
        if (len(listItems) < 4):
            print("Error input in addSxyzInfo!")
            return

        self.m_strName = listItems[0]
        self.m_strX = listItems[1]
        self.m_strY = listItems[2]
        self.m_strZ = listItems[3]

    def addGroInfo(self, strInputLine):
        strLine = strInputLine.strip()
        strLine = strLine.strip("\n")
        strLine = ' '.join(strLine.split()) #处理多个空格间隔的情况
        listItems = strLine.split()
        if (len(listItems) < 6):
            print("Error input in addGroInfo!")
            return

        #self.m_strStart = listItems[0]
        self.m_strName = listItems[1]
        self.m_strID = listItems[2]
        self.m_strX = "%.4f"%(float(listItems[3])*10)
        self.m_strY = "%.4f"%(float(listItems[4])*10)
        self.m_strZ = "%.4f"%(float(listItems[5])*10)

    def joinToGroString(self):
        strLine = "%8s"%(self.m_strResidueID + self.m_strResidueName) #
        strLine += "%7s"%(self.m_strName) #
        strLine += "%5s"%(self.m_strID) #

        strX = "%.3f"%(float(self.m_strX)/10)
        strY = "%.3f"%(float(self.m_strY)/10)
        strZ = "%.3f"%(float(self.m_strZ)/10)
        strLine += "%8s"%(strX) #
        strLine += "%8s"%(strY) #
        strLine += "%8s"%(strZ) #
        strLine += "\n"
        return strLine

    def getTinkerClass(self):

        #目前功能还有问题

        # charmm22cmap.prm中原子类型
        # 部分类型如H,HB,HR3,CT1等有多个值，在此简化为直接取第一个值。目标只是计算GB能
        listCharmm = ["HA","HP","H","HB","HC","H","HR1","HR2","HR3","HS",
                      "C","CA","CC","CT1","CT2","CT3","CP1","CP2","CP3","CH1",
                      "CH2","CY","CPT","CT","NH1","NH2","NH3","N","NP","NR1",
                      "NR2","NR3","NY","NC2","O","OH1","OC","S","SM"]

        listTinker = [1,2,3,4,6,8,12,13,14,16,
                      20,21,22,23,26,27,30,31,32,42,
                      45,50,52,62,63,64,65,66,67,68,
                      69,70,71,72,74,76,78,80,82]

        #小分子部分
        dictOther = {}
        dictOther["OR"] = "OH1"
        dictOther["O2CM"] = "OC"
        dictOther["CR"] = "CA"
        dictOther["C=C"] = "CA"
        dictOther["CO2M"] = "CC"
        dictOther["HCMM"] = "HP"
        dictOther["HOR"] = "H"
        dictOther["CPH1"] = "CP1"
        dictOther["CPH2"] = "CP2"

        if self.m_strPrmType in dictOther.keys():
            self.m_strPrmType = dictOther[self.m_strPrmType]

        nIndex = listCharmm.index(self.m_strPrmType)
        if (nIndex >= 0):
            self.m_strTinkerClass = str(listTinker[nIndex])
