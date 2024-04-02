# Common public methods
import math
import numpy as np
import objects

# Determine if it is a number
def is_number(inputstring):
    s = str(inputstring)
    if s.count('.') == 1:  # decimal
        new_s = s.split('.')
        left_num = new_s[0]
        right_num = new_s[1]
        if right_num.isdigit():
            if left_num.isdigit():
               return True
            elif left_num.count('-') == 1 and left_num.startswith('-'):  # nagative decimal
                tmp_num = left_num.split('-')[-1]
                if tmp_num.isdigit():
                    return True
    elif s.count(".") == 0:  # integer
        if s.isdigit():
            return True
        elif s.count('-') == 1 and s.startswith('-'):  # negative integer
            ss = s.split('-')[-1]
            if ss.isdigit():
                return True
    return False


# Bracket matching. Only handle '()', return the range of matching brackets. For '()()' type of brackets, store them sequentially. When used, they can be expanded directly.
# After processing brackets like 'aa((xxxx)yyy(zzz))bb', the innermost brackets are pushed to the stack listBrackets
def checkBrackets(stringInput, listBackets):
    listBackets.clear()
    listLeftIndex = []
    for index in range(len(stringInput)):
        if stringInput[index] == "(":
            listLeftIndex.append(index)
        if stringInput[index] == ")" and len(listLeftIndex) > 0:  # Popping the left bracket off the stack and constructing a bracket object
            leftIndex = listLeftIndex.pop()
            curObject = objects.CBracketObject(str(leftIndex))
            curObject.m_strName = "("
            curObject.m_strRightID = index
            curObject.m_strRightName = ")"
            listBackets.append(curObject)

    return len(listBackets) > 0


# Extract the composition and count of elements from the string, and return a dictionary
# Only handle strings like Al2O3, CaO, Na2O. It's possible to encounter strings like Na3B5O9H2O
def getElementFromString(stringInput, dictElements):
    dictElements.clear()
    strElement = ""
    strNum = ""
    for curChar in stringInput:
        if curChar.isupper():  # The end of the old element, the start of the new element
            if len(strNum) < 1:
                strNum = "1"
            if strElement in dictElements.keys():
                dictElements[strElement] = dictElements[strElement] + int(strNum)
            elif len(strElement) > 0:
                dictElements[strElement] = int(strNum)
            strElement = curChar
            strNum = ""
        elif curChar.islower():  # Element ended, no element symbol with 3 characters
            strElement = strElement + curChar
        elif curChar.isdigit():
            if len(strElement) > 0: # It's possible to have forms like 11Na2O, we need to exclude the preceding numbers
                strNum = strNum + curChar

    # It's possible that the last element is a number, so we need to handle the last element
    if len(strElement) > 0:
        if len(strNum) < 1:
            strNum = "1"
        if strElement in dictElements.keys():
            dictElements[strElement] = dictElements[strElement] + int(strNum)
        else:
            dictElements[strElement] = int(strNum)


# Extract the composition and count of elements from the chemical formula, modify the dictionary, and return whether the operation was successful
# Formulas with brackets, such as (Al2O3)10(B2O3)2, Ca(B3O4(OH)3)(H2O) can also be handled
# Note that the strFormula may be modified for bracket expansion
def getElementFromFormula(strFormula, dictElements):
    dictElements.clear()
    listBackets = []
    dictSub = {}  # Used to temporarily store elements within brackets
    strNum = ""

    # Loop until all brackets are expanded
    while (checkBrackets(strFormula, listBackets)):
        leftIndex = int(listBackets[0].m_strID)
        rightIndex = int(listBackets[0].m_strRightID)

        # Obtain the substring
        strSubFormula = strFormula[leftIndex+1:rightIndex]
        dictSub.clear()
        getElementFromString(strSubFormula,dictSub)

        # Retrieve the count after the bracket. 
        # If there are letters or other symbols following, break
        # Also need to remove the count
        nDigitNum = 0
        strNum = ""
        for index in range(rightIndex+1,len(strFormula)):
            if not strFormula[index].isdigit():
                break
            strNum = strNum + strFormula[index]
            nDigitNum += 1

        # Update the count of elements within brackets
        # Remove brackets from the original string and expand the elements in dictSub
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
        strFormula = strLeft + strDictSub + strRight  # The input string has been modified

    # Finally, obtain all elements and their quantities
    getElementFromString(strFormula,dictElements)
    strOutFormula = ""
    for key in dictElements.keys():
        strOutFormula = strOutFormula + key + str(dictElements[key])

    return strOutFormula


# Divide Label into two based on the input, arrLabel1D (array)
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


# Divide the input one-dimensional array into two based on the median value, commonly used for binary classification tasks to binaryize labels
# The one-dimensional array should be of numeric types, such as float or int, and not string
# Return the divided site
def biDvidedDataVectorByMid(arrLabel1D):
    dMid = 0.1
    for item in arrLabel1D:
        dMid += item

    nNum = arrLabel1D.size
    dMid = dMid / nNum

    biDvidedDataVectorByPoint(arrLabel1D, dMid)

    return dMid


# Calculate the information entropy of the array arrData1D. The array is one-dimensional and contains only 0 or 1
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


# Calculate the distance between two points. Note that fractional coordinates from CIF files might be passed as input
# So CAtom is an objects.CAtom class, listLengths are the lengths of the three directions in the crystal, and listAngles are the angles in the crystal constants
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

# Calculate the distance between two points
def calDistanceSimple(listCoords1, listCoords2):

    dDistX = listCoords1[0] - listCoords2[0]
    dDistY = listCoords1[1] - listCoords2[1]
    dDistZ = listCoords1[2] - listCoords2[2]
    dR2 = dDistX*dDistX + dDistY*dDistY + dDistZ*dDistZ
    return math.sqrt(dR2)


# Read the bondValence values used for calculating the flexibility of functional groups
# https://www.iucr.org/resources/data/data-sets/bond-valence-parameters
# https://doi.org/10.1038/s41598-020-60410-x
def readBondValenceParas():
    strFile = objects.BASISDIR + "public/bondValence2020.cif"
    listBondValParas = []
    file = open(strFile,"r")
    while 1:
        strLine = file.readline()
        if not strLine:
            break

        # Handle cases where multiple spaces are used as separators
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


# Filter all bonds based on a series of atomic coordinates and unit cell parameters, and calculate Flexibility
# Note: Ensure that the IDs in the listAtoms object correspond to the atomic indices
# If the Name is empty, this function will modify the Name in the listAtoms
def getBondsFromAtoms(listAtoms, listElementInfo, listLengths = [1.0,1.0,1.0], listAngels = [90.0,90.0,90.0]):

    listBonds = []

    # Retrieve the bondValence values for different atomic pairs
    listBondValParas = readBondValenceParas()

    # Retrieve the valence electron count and electronegativity for different atoms
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


    # Process Name based on ID
    for item in listAtoms:
        if item.m_strID in dictAtomIndex.keys():
            item.m_strName = dictAtomIndex[item.m_strID]


    # Elements that cannot form bonds with each other (or currently do not exist)
    listExcludeAtoms1 = ['N', 'F', 'Cl', 'Br', 'I']
    listExcludeAtoms2 = ['O', 'F', 'Cl', 'Br']
    listExcludeAtoms3 = ['Se', 'N', 'F', 'Cl', 'Br', 'I']
    listExcludeAtoms4 = ['S', 'N', 'I']

    # Calculate all bonds, filtered by atomic types and distance
    nAtomNum = len(listAtoms)
    for nIndex in range(0,nAtomNum):
        for nSubIndex in range(nIndex+1,nAtomNum):

            # Exclude pairs of identical atoms
            # Exclude hydrogen atoms
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

            # The regular bond length cannot be smaller than 1
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

# Filter all bonds and calculate SymmetryFunction based on a series of atomic coordinates and unit cell parameters
# Note: Ensure that the IDs in the listAtoms object correspond to the atomic indices
# If the Name is empty, this function will modify the Name in the listAtoms
def getSymmeFuncFromAtoms(listAtoms, listElementInfo, listLengths = [1.0,1.0,1.0], listAngels = [90.0,90.0,90.0]):

    listBonds = []

    # Retrieve the bondValence values for different atomic pairs
    listBondValParas = readBondValenceParas()

    # Retrieve the valence electron count and electronegativity for different atoms
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


    # Process Name based on ID
    for item in listAtoms:
        if item.m_strID in dictAtomIndex.keys():
            item.m_strName = dictAtomIndex[item.m_strID]

    # Calculate all bonds, filtered by atomic types and distance
    dCutOffRadui = 6.0
    nAtomNum = len(listAtoms)
    for nIndex in range(0,nAtomNum):
        for nSubIndex in range(nIndex+1,nAtomNum):

            # Ignore whether they are actually bonded
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


# Retrieve the different items between two lists
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

# Sort the dictionary by its keys/values and return all keys after sorting
# Only handle values of type float
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


# Retrieve the index corresponding to the minimum item in the list
def getMinimalIndex(listData1):

    nMiniIndex = -1
    dMinValue = 9999999.99999
    for nSubIndex in range(0,len(listData1)):
        if (listData1[nSubIndex] < dMinValue):
            dMinValue = listData1[nSubIndex]
            nMiniIndex = nSubIndex

    return nMiniIndex

# Retrieve several minimum items from the list
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

# Retrieve several items less than the Threshold from the list
def getMinimalSubArrByThreshold(listData1, dThreshold):

    listRetunIndexs = []
    for nIndex in range(0,len(listData1)):
        if (listData1[nIndex] < dThreshold):
            listRetunIndexs.append(nIndex)

    return listRetunIndexs

# Retrieve the index corresponding to the maximum item in the lis
def getMaximalIndex(listData1):

    nMaxIndex = -1
    dMaxValue = -9999999.99999
    for nSubIndex in range(0,len(listData1)):
        if (listData1[nSubIndex] > dMaxValue):
            dMaxValue = listData1[nSubIndex]
            nMaxIndex = nSubIndex

    return nMaxIndex

# Calculate the centroid
def CalMassCenter(listAtoms, listCoords):

    listMassCenter = [0.0,0.0,0.0]

    nAtoms = len(listAtoms)
    if nAtoms != len(listCoords):
        print("The number of atoms does not match the number of coordinates!")
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
            print("The mass of %s does not exist\n"%(listAtoms[nIndex]))
        else:
            dTotMass += dictAtomMass[listAtoms[nIndex]]

            listMassCenter[0] += dictAtomMass[listAtoms[nIndex]] * listCoords[nIndex][0]
            listMassCenter[1] += dictAtomMass[listAtoms[nIndex]] * listCoords[nIndex][1]
            listMassCenter[2] += dictAtomMass[listAtoms[nIndex]] * listCoords[nIndex][2]

    listMassCenter[0] /= dTotMass
    listMassCenter[1] /= dTotMass
    listMassCenter[2] /= dTotMass
    return dTotMass, listMassCenter

# Calculate the dihedral angle
def calDihedral():

    C = [-1.123, 1.383, 0.000]
    C1 = [-0.609, -0.069,0.000]
    C2 = [0.931, -0.069,-0.000]
    C3 = [-2.663, 1.383, 0.000]

    A = np.array(C2)
    B = np.array(C)
    C = np.array(C1)
    D = np.array(C3)

    # Pay attention to the order of subtraction, as it determines the direction of the vectors, 
    # and thus the direction of the normal vector of the plane
    a, b, c = B-A, C-B, D-C
    # Use the mathematical cross product of vectors: n1 = a × b, n2 = b × c 
    # Note the order, which determines the direction of the normal vectors 
    # If the direction is wrong, the dihedral angle will be incorrect
    n1, n2 = np.cross(a, b), np.cross(b, c)
    # Calculate the angle between the two normal vectors using the law of cosines
    # θ = arccos{（n1·n2）/（|n1|*|n2|)}
    normal1 = np.linalg.norm(n1)
    normal2 = np.linalg.norm(n2)
    dotProduct = np.dot(n1, n2)
    dihedral = math.acos(dotProduct/(normal1*normal2)) * 180 / math.pi
    print(dihedral)

# Test individually when calling
if __name__ == '__main__':

    calDihedral()


