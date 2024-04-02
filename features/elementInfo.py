#从外部文件中获得所有元素信息

import os
import sys
sys.path.append("..") #引入上级目录
import objects
import method

#从完整的文本中直接获取ElementInfo
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

        curElement = objects.CElementObject(listInfos[0]); #以原子序数
        curElement.m_strName = listInfos[1]#以元素符号
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


#从txt文件中读取所有Element相关信息
def readAllElementInfo():
    #选择txt文件所在目录
    strElementFile = objects.BASISDIR + "public/element-in.csv"
    listElementInfo = readElementInfoFromTxt(strElementFile)
    return  listElementInfo

#保存所有的信息
def saveAllElementInfo(listElementInfo):
    strElementFile = objects.BASISDIR + "element-out.csv"
    file = open(strElementFile,"a+")
    file.write("\n--------Elements---------------\n")
    for curElement in listElementInfo:
        file.write(curElement.joinToString(","))
        file.write("\n")
    file.close()


