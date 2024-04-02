#常见的公共方法
import math
import numpy as np
import objects

#判断是否为数字
def is_number(inputstring):
    s = str(inputstring)
    if s.count('.') == 1:#小数
        new_s = s.split('.')
        left_num = new_s[0]
        right_num = new_s[1]
        if right_num.isdigit():
            if left_num.isdigit():
               return True
            elif left_num.count('-') == 1 and left_num.startswith('-'):#负小数
                tmp_num = left_num.split('-')[-1]
                if tmp_num.isdigit():
                    return True
    elif s.count(".") == 0:#整数
        if s.isdigit():
            return True
        elif s.count('-') == 1 and s.startswith('-'):#负整数
            ss = s.split('-')[-1]
            if ss.isdigit():
                return True
    return False


#括号匹配，只处理()，返回括号匹配的范围，()()类型的括号依次存储即可，使用时直接展开
#aa((xxxx)yyy(zzz))bb类型的括号本函数处理完之后，最里层的括号最后入栈listBackets
def checkBrackets(stringInput, listBackets):
    listBackets.clear()
    listLeftIndex = []
    for index in range(len(stringInput)):
        if stringInput[index] == "(":
            listLeftIndex.append(index)
        if stringInput[index] == ")" and len(listLeftIndex) > 0: #左括号出栈并构建括号对象
            leftIndex = listLeftIndex.pop()
            curObject = objects.CBracketObject(str(leftIndex))
            curObject.m_strName = "("
            curObject.m_strRightID = index
            curObject.m_strRightName = ")"
            listBackets.append(curObject)

    return len(listBackets) > 0


#从字符串中得到元素的组成和个数，返回dict
#只处理形如Al2O3,CaO,Na2O，有可能出现Na3B5O9H2O等
def getElementFromString(stringInput, dictElements):
    dictElements.clear()
    strElement = ""
    strNum = ""
    for curChar in stringInput:
        if curChar.isupper(): #旧元素的结束，新元素开始
            if len(strNum) < 1:
                strNum = "1"
            if strElement in dictElements.keys():
                dictElements[strElement] = dictElements[strElement] + int(strNum)
            elif len(strElement) > 0:
                dictElements[strElement] = int(strNum)
            strElement = curChar
            strNum = ""
        elif curChar.islower():#元素结束，没有3个字符的元素符号
            strElement = strElement + curChar
        elif curChar.isdigit():
            if len(strElement) > 0: #有可能存在11Na2O等形式，要排除前面的数字
                strNum = strNum + curChar

    #最后有可能是数字所以需要处理最后一个元素
    if len(strElement) > 0:
        if len(strNum) < 1:
            strNum = "1"
        if strElement in dictElements.keys():
            dictElements[strElement] = dictElements[strElement] + int(strNum)
        else:
            dictElements[strElement] = int(strNum)


#从化学式中获取元素的组成和个数，修改dict并返回是否处理成功
#可处理带括号的形式，如(Al2O3)10(B2O3)2、Ca(B3O4(OH)3)(H2O)
# 注意，可能会修改strFormula进行括号展开
def getElementFromFormula(strFormula, dictElements):
    dictElements.clear()
    listBackets = []
    dictSub = {}#用来临时存储括号内的元素
    strNum = ""

    #循环直到所有括号都被展开
    while (checkBrackets(strFormula, listBackets)):
        leftIndex = int(listBackets[0].m_strID)
        rightIndex = int(listBackets[0].m_strRightID)

        #得到子串
        strSubFormula = strFormula[leftIndex+1:rightIndex]
        dictSub.clear()
        getElementFromString(strSubFormula,dictSub)

        #获得括号后面的个数，如果后面是字母或其他符号则break
        #还需要将个数移除掉
        nDigitNum = 0
        strNum = ""
        for index in range(rightIndex+1,len(strFormula)):
            if not strFormula[index].isdigit():
                break
            strNum = strNum + strFormula[index]
            nDigitNum += 1

        #更新括号中元素个数
        #去除原字符串中括号，并将dictSub中元素展开
        strLeft = ""
        if leftIndex > 0:
            strLeft = strFormula[0:int(leftIndex)]
        strRight = ""
        if rightIndex < len(strFormula):
            strRight = strFormula[int(rightIndex)+1+nDigitNum:len(strFormula)]
        strDictSub = ""
        for key in dictSub.keys():
            if len(strNum) > 0:
                dictSub[key] = int(dictSub[key]) * int(strNum)
            strDictSub = strDictSub + key + str(dictSub[key])
        strFormula = strLeft + strDictSub + strRight #输入字符串已经被改变了

    #最后获得所有的元素和数量
    getElementFromString(strFormula,dictElements)
    strOutFormula = ""
    for key in dictElements.keys():
        strOutFormula = strOutFormula + key + str(dictElements[key])

    return strOutFormula


#根据输入将Label二分化,arrLabel1D(adarray)
def biDvidedDataVectorByPoint(arrLabel1D, dDivPoint):
    nIndex = 0
    nNum = arrLabel1D.size
    while(nIndex < nNum):
        if arrLabel1D[nIndex] < dDivPoint :
            arrLabel1D[nIndex] = 0
        else:
            arrLabel1D[nIndex] = 1
        nIndex += 1

    return dDivPoint


#将输入的一维数组，根据中间值二分化，常用于分类任务中标签的二分化
#一维数组应该是数值类型的，float，int，而不能是string
#返回分值点
def biDvidedDataVectorByMid(arrLabel1D):
    dMid = 0.1
    for item in arrLabel1D:
        dMid += item

    nNum = arrLabel1D.size
    dMid = dMid / nNum

    biDvidedDataVectorByPoint(arrLabel1D, dMid)

    return dMid


#计算数组arrData1D中的的信息熵，一维数组，且内容只有0和1
def calInfoEntropy(arrData1D):
    dEntropy = 0.0
    dictValue = {}
    nNum = arrData1D.size
    nIndex = 0
    while nIndex < nNum:
        nKey = int(arrData1D[nIndex])
        if nKey in dictValue.keys():
            dictValue[nKey] += 1
        else:
            dictValue[nKey] = 1

        nIndex += 1

    for value in dictValue.values():
        dEntropy = (value * math.log(value/nNum)) / nNum

    return dEntropy


#计算两点距离，注意有可能传入的是cif中的fractional坐标，所以
# CAtom是objects.CAtom类，listLengths是晶体三个方向的边长，listAngels是晶体常数中的角度
# http://pd.chem.ucl.ac.uk/pdnn/refine2/bonds.htm
def calDistance(CAtom1, CAtom2, listLengths = [1.0,1.0,1.0], listAngels = [90.0,90.0,90.0]):

        dDistX = float(CAtom1.m_strX) - float(CAtom2.m_strX)
        dDistX2 = dDistX * dDistX * listLengths[0] * listLengths[0]
        dDistY = float(CAtom1.m_strY) - float(CAtom2.m_strY)
        dDistY2 = dDistY * dDistY * listLengths[1] * listLengths[1]
        if CAtom1.m_strZ == "-." or CAtom2.m_strZ == "-.":
            print("ddd")
        dDistZ = float(CAtom1.m_strZ) - float(CAtom2.m_strZ)
        dDistZ2 = dDistZ * dDistZ * listLengths[2] * listLengths[2]
        dDistance2 = dDistX2 + dDistY2 + dDistZ2
        dDistance2 = dDistance2 + 2 * listLengths[1] * listLengths[2] * dDistY * dDistZ * math.cos(listAngels[0]*math.pi/180)
        dDistance2 = dDistance2 + 2 * listLengths[0] * listLengths[2] * dDistX * dDistZ * math.cos(listAngels[1]*math.pi/180)
        dDistance2 = dDistance2 + 2 * listLengths[0] * listLengths[1] * dDistX * dDistY * math.cos(listAngels[2]*math.pi/180)

        return math.sqrt(dDistance2)

#计算两点距离
def calDistanceSimple(listCoords1, listCoords2):

    dDistX = listCoords1[0] - listCoords2[0]
    dDistY = listCoords1[1] - listCoords2[1]
    dDistZ = listCoords1[2] - listCoords2[2]
    dR2 = dDistX*dDistX + dDistY*dDistY + dDistZ*dDistZ
    return math.sqrt(dR2)


#读取计算基团Flexibility的bondValence值
# https://www.iucr.org/resources/data/data-sets/bond-valence-parameters
#  https://doi.org/10.1038/s41598-020-60410-x
def readBondValenceParas():
    strFile = objects.BASISDIR + "public/bondValence2020.cif"
    listBondValParas = []
    file = open(strFile,"r")
    while 1:
        strLine = file.readline()
        if not strLine:
            break

        #处理多个空格分隔的情况
        strLine = ' '.join(strLine.split())
        strLine = strLine.strip()
        listInfos = strLine.split(" ")
        if len(listInfos) < 6:
            continue

        curPara = objects.CBondValencePara(listInfos[0])
        curPara.m_strName = listInfos[2]
        curPara.m_strRo = listInfos[4]
        curPara.m_strB = listInfos[5]
        listBondValParas.append(curPara)

    file.close()

    return listBondValParas


# 根据一系列原子坐标和晶胞参数来过滤所有的键，并计算Flexibility
# 注意，要保证listAtoms对象中ID是元素序号
# 如果Name为空则此函数会修改listAtoms中的Name
def getBondsFromAtoms(listAtoms, listElementInfo, listLengths = [1.0,1.0,1.0], listAngels = [90.0,90.0,90.0]):

    listBonds = []

    #获取不同原子对的bondValence值
    listBondValParas = readBondValenceParas()

    #获取不同原子的价电子数和电负性
    dictAtomIndex = {}
    dictAtomName = {}
    dictAtomValEleNum = {}
    dictAtomEleNegativity = {}
    for elem in listElementInfo:
        dictAtomIndex[elem.m_strID] = elem.m_strName
        dictAtomName[elem.m_strName] = elem.m_strID
        dValEleNum = int(elem.m_strSValenceEleNum)
        dValEleNum = dValEleNum + int(elem.m_strPValenceEleNum)
        dValEleNum = dValEleNum + int(elem.m_strDValenceEleNum)
        dValEleNum = dValEleNum + int(elem.m_strFValenceEleNum)
        dictAtomValEleNum[elem.m_strName] = dValEleNum
        dictAtomEleNegativity[elem.m_strName] = float(elem.m_strEleNegativity)


    #根据ID处理Name
    for item in listAtoms:
        if item.m_strID in dictAtomIndex.keys():
            item.m_strName = dictAtomIndex[item.m_strID]


    #两两之间不能成键(或当前不存在)的元素
    listExcludeAtoms1 = ['N', 'F', 'Cl', 'Br', 'I']
    listExcludeAtoms2 = ['O', 'F', 'Cl', 'Br']
    listExcludeAtoms3 = ['Se', 'N', 'F', 'Cl', 'Br', 'I']
    listExcludeAtoms4 = ['S', 'N', 'I']

    #计算所有键,通过原子类型和距离进行过滤
    nAtomNum = len(listAtoms)
    for nIndex in range(0,nAtomNum):
        for nSubIndex in range(nIndex+1,nAtomNum):

            #不考虑两个相同的原子，#不考虑H原子
            strAtom1 = listAtoms[nIndex].m_strName
            strAtom2 = listAtoms[nSubIndex].m_strName
            bIsNotValid = (strAtom1 == strAtom2 or strAtom1 == "H" or strAtom2 == "H")
            bIsNotValid = bIsNotValid or (strAtom1 in listExcludeAtoms1 and strAtom2 in listExcludeAtoms1)
            bIsNotValid = bIsNotValid or (strAtom1 in listExcludeAtoms2 and strAtom2 in listExcludeAtoms2)
            bIsNotValid = bIsNotValid or (strAtom1 in listExcludeAtoms3 and strAtom2 in listExcludeAtoms3)
            bIsNotValid = bIsNotValid or (strAtom1 in listExcludeAtoms4 and strAtom2 in listExcludeAtoms4)
            if bIsNotValid:
                continue

            dDistance = calDistance(listAtoms[nIndex], listAtoms[nSubIndex], listLengths, listAngels)

            #普通键长不会小于1；
            if dDistance < 1.0:
                continue

            nAtomNum1 = int(dictAtomName[strAtom1])
            nAtomNum2 = int(dictAtomName[strAtom2])
            if nAtomNum1 < 18 and nAtomNum2 < 18 and dDistance > 2.18:
                continue
            elif dDistance > 2.8:
                continue

            curBond = objects.CBond("")
            curBond.m_CAtom1 = listAtoms[nIndex]
            curBond.m_CAtom2 = listAtoms[nSubIndex]
            curBond.m_strLength = str(dDistance)
            curBond.calFlexibility(listBondValParas,dictAtomValEleNum,dictAtomEleNegativity)
            listBonds.append(curBond)

    return listBonds

# 根据一系列原子坐标和晶胞参数来过滤所有的键，并计算SymmetryFunction
# 注意，要保证listAtoms对象中ID是元素序号
# 如果Name为空则此函数会修改listAtoms中的Name
def getSymmeFuncFromAtoms(listAtoms, listElementInfo, listLengths = [1.0,1.0,1.0], listAngels = [90.0,90.0,90.0]):

    listBonds = []

    #获取不同原子对的bondValence值
    listBondValParas = readBondValenceParas()

    #获取不同原子的价电子数和电负性
    dictAtomIndex = {}
    dictAtomName = {}
    dictAtomValEleNum = {}
    dictAtomEleNegativity = {}
    for elem in listElementInfo:
        dictAtomIndex[elem.m_strID] = elem.m_strName
        dictAtomName[elem.m_strName] = elem.m_strID
        dValEleNum = int(elem.m_strSValenceEleNum)
        dValEleNum = dValEleNum + int(elem.m_strPValenceEleNum)
        dValEleNum = dValEleNum + int(elem.m_strDValenceEleNum)
        dValEleNum = dValEleNum + int(elem.m_strFValenceEleNum)
        dictAtomValEleNum[elem.m_strName] = dValEleNum
        dictAtomEleNegativity[elem.m_strName] = float(elem.m_strEleNegativity)


    #根据ID处理Name
    for item in listAtoms:
        if item.m_strID in dictAtomIndex.keys():
            item.m_strName = dictAtomIndex[item.m_strID]

    #计算所有键,通过原子类型和距离进行过滤
    dCutOffRadui = 6.0
    nAtomNum = len(listAtoms)
    for nIndex in range(0,nAtomNum):
        for nSubIndex in range(nIndex+1,nAtomNum):

            #不考虑是否真正成键
            # International Journal of Quantum Chemistry 2015, 115, 1032–1050
            dDistance = calDistance(listAtoms[nIndex], listAtoms[nSubIndex], listLengths, listAngels)
            if dDistance > dCutOffRadui:
                continue

            curBond = objects.CBond("")
            curBond.m_CAtom1 = listAtoms[nIndex]
            curBond.m_CAtom2 = listAtoms[nSubIndex]
            curBond.m_strLength = str(dDistance)
            curBond.calSymmetryFunction(listBondValParas,dCutOffRadui)
            listBonds.append(curBond)

    return listBonds


#获得两个list中的不同项
def compareTwoListDiff(listData1, listData2):
    listNotIn2 = []
    for item in listData1:
        if item not in listData2:
            listNotIn2.append(item)

    listNotIn1 = []
    for item in listData2:
        if item not in listData1:
            listNotIn1.append(item)

    return listNotIn1,listNotIn2

#对字典按其keys/values进行排序，并返回排序后的所有keys
#只处理values是float的类型
def sortDictions(dictTarget, strType):

    if strType == "keys":
        listKeys = list(dictTarget.keys())
        listKeys.sort()
        return listKeys

    listValues = []
    for value in dictTarget.values():
        listValues.append(value)

    listValues.sort()
    listKeys = []
    for value in listValues:
        for key in dictTarget.keys():
            if math.isclose(value, dictTarget[key], rel_tol=1e-6) and key not in listKeys:
                listKeys.append(key)
                break

    return listKeys


#获得list中的最小项对应的索引
def getMinimalIndex(listData1):

    nMiniIndex = -1
    dMinValue = 9999999.99999
    for nSubIndex in range(0,len(listData1)):
        if (listData1[nSubIndex] < dMinValue):
            dMinValue = listData1[nSubIndex]
            nMiniIndex = nSubIndex

    return nMiniIndex

#获得list中的最小项若干项
def getMinimalSubArr(listData1, nNum):

    dictData = {}
    listValues = []
    for nIndex in range(0,len(listData1)):
        dictData[nIndex] = listData1[nIndex]
        listValues.append(listData1[nIndex])

    listValues.sort()
    listRetunIndexs = []
    for value in listValues:
        for key in dictData.keys():
            if math.isclose(value, dictData[key], rel_tol=1e-6) and key not in listRetunIndexs:
                listRetunIndexs.append(key)
                break
        if len(listRetunIndexs) > nNum-1:
            break

    return listRetunIndexs

#获得list中的小于Threshold若干项
def getMinimalSubArrByThreshold(listData1, dThreshold):

    listRetunIndexs = []
    for nIndex in range(0,len(listData1)):
        if (listData1[nIndex] < dThreshold):
            listRetunIndexs.append(nIndex)

    return listRetunIndexs

#获得list中的最大项对应的索引
def getMaximalIndex(listData1):

    nMaxIndex = -1
    dMaxValue = -9999999.99999
    for nSubIndex in range(0,len(listData1)):
        if (listData1[nSubIndex] > dMaxValue):
            dMaxValue = listData1[nSubIndex]
            nMaxIndex = nSubIndex

    return nMaxIndex

#计算质心
def CalMassCenter(listAtoms, listCoords):

    listMassCenter = [0.0,0.0,0.0]

    nAtoms = len(listAtoms)
    if nAtoms != len(listCoords):
        print("原子数与坐标数不匹配！")
        return listMassCenter

    dictAtomMass = {}
    dictAtomMass["H"] = 1.008
    dictAtomMass["C"] = 12.011
    dictAtomMass["N"] = 14.007
    dictAtomMass["O"] = 15.999
    dictAtomMass["F"] = 18.998
    dictAtomMass["P"] = 30.974
    dictAtomMass["S"] = 32.065
    dictAtomMass["Cl"] = 35.453
    dictAtomMass["Br"] = 79.904
    dictAtomMass["I"] = 126.904

    dTotMass = 0.0
    for nIndex in range(0,nAtoms):

        if listAtoms[nIndex] not in dictAtomMass.keys():
            print("%s的质量不存在\n"%(listAtoms[nIndex]))
        else:
            dTotMass += dictAtomMass[listAtoms[nIndex]]

            listMassCenter[0] += dictAtomMass[listAtoms[nIndex]] * listCoords[nIndex][0]
            listMassCenter[1] += dictAtomMass[listAtoms[nIndex]] * listCoords[nIndex][1]
            listMassCenter[2] += dictAtomMass[listAtoms[nIndex]] * listCoords[nIndex][2]

    listMassCenter[0] /= dTotMass
    listMassCenter[1] /= dTotMass
    listMassCenter[2] /= dTotMass
    return dTotMass, listMassCenter

# 求二面角
def calDihedral():

    C = [-1.123, 1.383, 0.000]
    C1 = [-0.609, -0.069,0.000]
    C2 = [0.931, -0.069,-0.000]
    C3 = [-2.663, 1.383, 0.000]

    A = np.array(C2)
    B = np.array(C)
    C = np.array(C1)
    D = np.array(C3)

    # 需要注意减法的顺序，因为这决定了向量的方向，从而决定了平面法向量的方向
    a, b, c = B-A, C-B, D-C
    # 使用数学上向量的方法叉乘，n1=a×b，n2=b×c；（注意顺序，决定了法向量的方向，方向错了二面角会错）
    n1, n2 = np.cross(a, b), np.cross(b, c)
    # 得出的两个法向量用余弦定理求角度 θ = arccos{（n1·n2）/（|n1|*|n2|)}
    normal1 = np.linalg.norm(n1)
    normal2 = np.linalg.norm(n2)
    dotProduct = np.dot(n1, n2)
    dihedral = math.acos(dotProduct/(normal1*normal2)) * 180 / math.pi
    print(dihedral)

#单独调用时进行测试
if __name__ == '__main__':

    calDihedral()


