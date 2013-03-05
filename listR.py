#                                                                version 19 --- Zhi Qiu
# level 1
"""
	Provide useful functions dealing with lists and strings.
	If a function ends with "I", it returns an iterator, "D" for dictionary,
	and "L" for list.
"""

import re

def toList(qvar):
    """ Make qvar a list if not."""
    if type(qvar) != type([]): qvar = [qvar]
    return qvar


def flatten(nested):
    """
        Flatten a list using iterator and exceptional control.
    """
    try:
        for sublist in nested:
          if type(sublist)==type(""):
            yield sublist
          else:
            for element in flatten(sublist):
                yield element
    except TypeError:
        yield nested


def FL(nested):
    """
        Flatten a list using iterator.
    """
    for sublist in nested:
        if type(sublist) == type([]):
            for element in FL(sublist): yield element
        else:
            yield sublist

FLI = FL # use FLI in future functions

def FLL(nested):
	""" Return a flattend list. """
	return list(FL(nested))


def totalLen(nested):
	""" Return the total number of elements in a list. """
	return len(FLL(nested))


def intersect(lists):
	"""
		Return the intersection of all lists in "lists".
	"""
	if len(lists) == 0: return lists
	if len(lists) == 1: return lists[0]
	finalList = set(lists[0])
	for aList in lists[1:]:
		finalList = finalList & set(aList)
	return list(finalList)


def union(lists):
	"""
		Return the union of all lists in "lists".
	"""
	if len(lists) == 0: return lists
	if len(lists) == 1: return lists[0]
	finalList = set(lists[0])
	for aList in lists[1:]:
		finalList = finalList | set(aList)
	return list(finalList)

def difference(lists):
	"""
		Return the first set minus the rest.
	"""
	if len(lists) == 0: return lists
	if len(lists) == 1: return lists[0]
	finalList = set(lists[0])
	for aList in lists[1:]:
		finalList = finalList - set(aList)
	return list(finalList)


def outer(lists, cuList=[]):
	"""
		Return all combination of elements in the sublist
		of lists.
	"""
	if lists == []:
		yield cuList
		return
	for element in lists[0]:
		newList = list(cuList)
		newList.append(element)
		for elem in outer(lists[1:], newList):
			yield elem



def strZip(list1, list2, string):
	""" Return a list of strings of the form x1stringx2 where x1
	and x2 are elements of list1 and list2 respectively.
	"""
	result = []
	for x1, x2 in zip(list1, list2):
		result.append(str(x1)+string+str(x2))
	return result


def listToStr(list, seperator=" "):
	""" Return a string consists of elements in list seperated by
	seperator.
	"""
	return seperator.join(map(str,list))


def applyOrderList(order, aList):
	""" Apply order to aList.
	An order of the form [2,0,1] means that the first element
	should be the 3rd element of aList, and so on.
	"""
	if order==[]:
		return aList
	else:
		return map(lambda v:aList[v], order)


def applyOrderDic(order, aDic):
	""" Apply order to aList.
	An order of the form ["a","b","c"] means that the first element
	should be aDic["a"], and so on.
	"""
	if order==[]:
		return aDic.value()
	else:
		return map(lambda v:aDic[v], order)


def createOrderList(wantedOrder, currentOrder):
	""" Create an order list that can transform currentOrder to
	wantedOrder by applying applyOrderList function.
	An order list is a list that specifies the position of the
	desired element in a list in the correct order, e.g:
	order of [3,1,2,4,6]->[1,2,3] which is got by using
	createOrderList([1,2,3],[3,1,2,4,6]) is [1,2,0].
	"""
	return map(lambda x:currentOrder.index(x), wantedOrder)


def firstOccurenceInStr(aList, aString):
	""" Return the first element in aList that is contained in the
	string aString.
	"""
	for elem in aList:
		if elem in aString:
			return elem
	else:
		return None


def getTailNumber(string):
    """
        Extracts the last occurence of a number from a string.
    """
    num = re.findall(r"[0-9][0-9\.]*", string)
    if num != []:
        return num[len(num)-1]
    else:
        return ""


def areDefined(names, dic):
	"""
		Return True if all keys that have "names" are defined and are not "None".
	"""
	key_list = dic.keys()
	names = toList(names)
	for aName in names:
		if aName not in key_list: return False
		if dic[aName] == None: return False
	else:
		return True


def itemsList(dic):
	""" Return dic.items() using list instead of tuple. """
	return map(list, dic.items())


def split(a_string, seperator):
	""" Split a string using seperator. """
	return a_string.split(seperator)




def _mimicCore(patternList, flatList):
	""" Core program for function mimic. """
	for elem in patternList:
		if type(elem) != type([]):
			yield flatList[0]
			flatList = flatList[1:]
			continue
		else:
			yield list(_mimicCore(elem, flatList[:len(list(FL(elem)))]))
			flatList = flatList[len(list(FL(elem))):]

def mimic(patternList, flatList):
	""" Make the flattened list (flatList) to have the same
	structure as patternList. (List only, no tuples)
	"""
	if len(list(FL(patternList))) != len(flatList):
		print("patternList must have the same number of total elements as flatList!")
		return None
	return list(_mimicCore(patternList, flatList))


def containedIn(smaller, larger):
	""" Return True if smaller (list) is contained in larger (list). """
	for elem in toList(smaller):
		if elem not in larger: return False
	else:
		return True

def biDifference(larger, smaller):
	""" Remove smaller from larger (not necessary to contain smaller). """
	largerCopy = list(larger)
	for elem in toList(smaller):
		if elem in larger:
			largerCopy.remove(elem)
	return largerCopy

def biSetDifference(larger, smaller):
  """ Remove smaller from larger (not necessary to contain smaller) using set operation. """
  return list(set(flatten(larger))-set(smaller))

def removeListFromDict(aDict, aList):
	""" Remove those items whose keys in aList from aDict. """
	aList = toList(aList)
	if aList == []: return aDict
	for elem in aList:
		aDict.pop(elem)
	return aDict

def biIntersectI(list1, list2):
	""" Return the intersection of the two lists. """
	for elem in toList(list1):
		if elem in list2: yield elem

def subDict(keys, aDict):
	""" Return a sub dictionary according to the list keys. """
	preDict = []
	allKeys = aDict.keys()
	for aKey in toList(keys):
		if aKey in allKeys: preDict.append([aKey, aDict[aKey]])
	return dict(preDict)

def removeDuplicatesSimple(aList):
	""" Remove duplicates in a simple list. """
	return list(set(aList))

def _removeDuplicatesOneLevel(aList):
	""" Remove first level duplicates. """
	result = []
	if aList == []: return
	if type(aList) != type([]): return aList
	for elem in aList:
		if elem not in result: result.append(elem)
	return result

def removeDuplicates(aList):
	""" Remove duplicates in a list. """
	if type(aList) != type([]): return aList
	result = []
	for elem in aList:
		result.append(removeDuplicates(elem))
	return _removeDuplicatesOneLevel(result)

def strEqual(str1, str2, ignoreCase=False):
	""" Return true if two strings are the same. """
	if ignoreCase==False:
		return str1.strip()==str2.strip()
	else:
		return str1.strip().upper()==str2.strip().upper()

def getValueListFromDict(keys, aDict):
    """ Return a list of values corresponding to the specified keys. """
    if type(keys)==type([]):
        return map(lambda x:aDict[x], keys)
    else:
        return []

def addItemsToDict(aList, aDict):
	""" Add a list of items to a dictionary. """
	return dict(removeDuplicates(itemsList(aDict)+aList))

def floatizeL(aList):
	""" Convert elements in aList to float number. """
	return map(float, aList)

def floatizeItemInDict(aDict, keyList):
    """ Convert items in aDict with keys in keyList to float number."""
    keyList = toList(keyList)
    result = aDict
    for elem in keyList:
        result[elem] = float(aDict[elem])
    return result


def stringizeL(aList):
    """ Convert elements in aList to strings. """
    if type(aList)==type([]):
        return map(str, flatten(aList))
    else:
        return str(aList)

def transpose(lists):
	""" Transpose a list of lists. """
	if not lists: return []
	return map(lambda *row: list(row), *lists)

def transpose2(lists, defval=0):
	""" Transpose a list of list. defval is used instead of None for uneven lengthed lists. """
	if not lists: return []
	return map(lambda *row: [elem or defval for elem in row], *lists)

def getColumn(data, colN):
	""" Return the column colN (counted starting from 0) in the data. """
	return transpose(data)[colN]

def getColumns(data, col_list):
	"""
		Take sublist given from all rows and columns specified by col_list from
		a doublely iterated list.
		The convention for the index of the rows and columns are the same as in slicing.
	"""
	result = []
	for aRow in data:
		result.append(map(lambda column: aRow[column], col_list))
	return result

def seperateStr(strV, seperationSymbols=[" ", "-", "\n"]):
	""" Split string according to seperationSymbols. """
	if not strV:
		return []
	else:
		strings = [strV]
		for aSeperator in seperationSymbols:
			strings = FLL(map(lambda x:split(x, aSeperator), strings))
	return strings

def readCSESD(strV, connectionSymbol="=", seperationSymbols=[",", " ", "-", "\n", "/", "\\"]): #CSE: Comma Seperated Equations
	""" Return a dic of the form {arg1:value1, ...} if with
	connectionSymbol="=" and seperationSymbol=",", strV is
	like arg1=value1,arg2=value2,...
	Values are in string form.
	"""
	result = []
	strings = seperateStr(strV, seperationSymbols) # total seperation
	for aStr in strings:
		# print(aStr)
		if connectionSymbol not in aStr: continue
		result.append(split(aStr, connectionSymbol))
	result = removeDuplicates(result) # remove duplicates
	# print(result)
	if not result: result=[] # make it "empty" instead of "None"
	if ([""] in result):
		result.remove([""]) # remove empty lists
	return dict(result)

def connectCSES(aList, connectionSymbol="=", seperationSymbol=","):
	""" Return a string represents a list of the form [name, value] using the form name=value. """
	if not aList: return ""
	return seperationSymbol.join(map(lambda x: connectionSymbol.join(x), aList))

def removeTailReturn(aStr):
	""" Return a string without the tail "\n" if it has one. """
	if aStr[-1:] == "\n":
		return aStr[:-1]
	else:
	   return aStr

def takeBlock(aList, row_l,row_r,col_l,col_r):
	""" Take sublist given from row row_l to row_r and column col_l to col_r from a double list.
	The convention for the index of the rows and columns are the same as in slicing.
	"""
	result = []
	for aRow in aList[row_l:row_r]:
		result.append(aRow[col_l:col_r])
	return result

def takeBlock2(aList, row_list, col_list):
	"""
		Take sublist given from rows specified by row_list and column specified by col_list from
		a doublely iterated list.
		The convention for the index of the rows and columns are the same as in slicing.
	"""
	result = []
	for row in row_list:
		result.append(map(lambda column: aList[row][column], col_list))
	return result

def intStr(i, total=3):
	""" Return a sting of the integer i begin with 0. """
	return '0'*(total-len(str(i)))+str(i)

def isNested(data):
	""" Return true if the data is a nested list. """
	if not data: return False # if empty then return False
	if type(data[0]) != type([]): return False # if not a nested list then return False
	return True

def sliceMatrixData(data, columnStep=0, centralLargeness=0):
	""" Slice data into smaller nested lists of specified vertical size (columnStep).
	If centralLargeness is not 0, only a smaller central block of specified size is used. """
	Ny, Nx = len(data), len(data[0])
	if columnStep==0: columnStep = Nx
	if Ny % columnStep != 0: # check if Ny is dividable by Nx
		print("Total length of data is not dividable by columnStep (or did you use float number for colStep?)!")
		return []

	if centralLargeness == 0: centralLargeness = min(columnStep, Nx) # set visible area
	y_left = int((columnStep-centralLargeness)/2) # set block size in y direction (row direction)
	y_right = int((columnStep+centralLargeness)/2)
	x_left = int((Nx-centralLargeness)/2) # set block in size x direction (column direction)
	x_right = int((Nx+centralLargeness)/2)

	result = []
	for i in range(Ny/columnStep):
	    result.append(take(data[i*columnStep:(i+1)*columnStep][:], y_left, y_right, x_left, x_right))

	return result

def next(aList, index):
	""" Return the index to the next element (compared to the element
	with index "index") in aList, or 0 if it already is the last one.
	Useful to make a list of loop.
	"""
	return (index+1) % len(aList)

def isFloat(string):
	""" Return a true is string is convertable to a float, otherwise false.
	"""
	try:
		float(string)
		return True
	except ValueError:
		return False

def zeros(m,n):
	""" Return a zero "matrix" (iterated list) with m rows (1st index) and n columns. """
	return map(lambda var: map(lambda var: 0, range(n)), range(m));

