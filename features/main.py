# 数据处理的主函数
import elementInfo, featureInfo, crystalInfo, groupInfo

from objects import BASISDIR
from collections import Counter
import numpy as np


if __name__ == '__main__':
    # 获取所有元素信息
    listElementInfo = elementInfo.readAllElementInfo()

    # 获得Group
    listNegaGroups = []
    listPosiGroups = []
    # listNegaGroups = groupInfo.readGroupFromOutFiles(listElementInfo, True)  # 从out文件中获取所有Group的信息
    # listPosiGroups = groupInfo.readGroupFromOutFiles(listElementInfo, False)  # 从out文件中获取所有Group的信息
    listNegaGroups, listPosiGroups = groupInfo.readAllGroupInfo()
    # groupInfo.saveAllGroupInfo(listNegaGroups, True)
    # groupInfo.saveAllGroupInfo(listPosiGroups, False)

    # 获取所有NLO晶体信息
    listCrystalInfo = crystalInfo.readAllCrystalInfo(listElementInfo)
    #crystalInfo.saveAllCrystalInfo(listCrystalInfo)

    # 构建晶体对应的描述符
    listFeatureInfo = featureInfo.createAllFeatureInfo(listPosiGroups, listNegaGroups, listElementInfo, listCrystalInfo)
    featureInfo.saveAllFeatureInfo(listFeatureInfo)

    # 预处理，二阶非线性系数与dij关系确定

    # 训练，参见featureInfo等
