###################    Last edited on Oct. 26, 2010   #################################
#                                                              version 22 --- Zhi Qiu
# level 2
"""
	Provide functions related to file operations or data file manipulations,
	and functions that are closely related to them.
"""

import dirR
import listR

import os
import sys
import shutil
import subprocess

__q_debug__ = False


# --- Used in all functions that call external programs ---
_xterm_path_directory = {"darwin":"/usr/X11/bin/xterm", "linux2":"/usr/bin/xterm"}
if sys.platform in _xterm_path_directory.keys():
    _default_xterm = _xterm_path_directory[sys.platform]
else:
	_default_xterm = _xterm_path_directory["linux2"]

_konsole_path_directory = {"darwin":"/usr/X11/bin/konsole", "linux2":"/usr/bin/konsole"}
if sys.platform in _konsole_path_directory.keys():
    _default_konsole = _konsole_path_directory[sys.platform]
else:
	_default_konsole = _konsole_path_directory["linux2"]

_gnome_terminal_path_directory = {"darwin":"/usr/bin/gnome-terminal", "linux2":"/usr/bin/gnome-terminal"}
if sys.platform in _gnome_terminal_path_directory.keys():
    _default_gnome_terminal = _gnome_terminal_path_directory[sys.platform]
else:
	_default_gnome_terminal = _gnome_terminal_path_directory["linux2"]

_default_terminal = _default_xterm



def runCom(command, cwd=os.getcwd(), terminal=_default_terminal):
	""" Invoke a command and wait for it to stop. """
	proc = subprocess.Popen(terminal+" -e "+command, shell=True, cwd=cwd)
	while proc.wait() != 0:
		pass

def runComUnlimited(command, cwd=os.getcwd(), sleepTime=0, terminal=_default_terminal):
    """ Invode a command and wait for it to stop.
    ulimit is set to unlimited before the execution of
    the program.
    """
    tmpsh = open(os.path.join(cwd,"QZTEMP.sh"), "w")
    tmpsh.write("ulimit -s unlimited\n")
    tmpsh.write(command+"\n")
    tmpsh.write("sleep "+str(sleepTime))
    if __q_debug__ == True:
        tmpsh.write("sleep 10")
    tmpsh.close()
    proc = subprocess.Popen(terminal+" -e bash QZTEMP.sh", shell=True, cwd=cwd)
    while proc.wait() != 0:
        pass
    #os.remove(os.path.join(cwd,"QZTEMP.sh"))

runCommand = runComUnlimited
execute = runCom



def copyIn_calculate_copyOut(source_dir, target_dir, exe_file, para_str, copyIn, copyOut, silent=True):
	"""
		Produce a copy-into, calculate, then copy-out procedure. Files
		given in "copyIn" (list of file names) are copied into
		"target_dir", then "exe_file" (path+name) with "para_str"
		(string of ALL parameters) is invoked, then files given in
		"copyOut" are copied into "source_dir". This is a low-level function
		and you should not run this directly.
	"""

	# copyIn:
	copyFiles(copyIn, source_dir, target_dir, silent)

	# invoke executable file and wait for it to stop:
	runCommand(exe_file+" "+para_str, os.path.dirname(exe_file))

	# copyOut:
	copyFiles(copyOut, target_dir, source_dir, silent)




def preExePost(argValueList, preFncs=[], exeFile="", postFncs=[], order=[], sleepTime=0): # exeFile="" to do nothing
	"""
		Execute exeFile with a choice of arguments specified in
		argNameList. argValueList is a nested list, elements from
		each sublist consist a "choice" of arguments. order is applied
		before this choice of arguments is used. For each possible
		combination of values of arguments, preFncs are executed first,
		then exeFile is executed, then postFncs are executed.

		preFnc is used to do some possible clean work between each
		run, it has one parameter: a list of choice of parameters.
		postFun is similar, only executed after each run.

		Note that exeFile is a string that contains the executable
		file, and the string can contain other arguments, as long as
		these arguments are positioned in front of those automatically
		generated arguments list of the form argNameList. The list
		argValueList is very free, for example, the list
		[["1 2","3 4"],...] will use "1 2" or "3 4" as the first
		generated parameter. This is especially useful if a program
		reads argument for a loop to generate one file.
	"""
	for valueCmb in listR.outer(map(listR.toList, listR.toList(argValueList))): # take a perticular list of possible combination of values
		valueCmb = list(listR.FL(valueCmb))
		# execute preFnc
		if preFncs!=[]:
			for pre_func in listR.toList(preFncs):
				pre_func(valueCmb)

		# execute exeFile with given choice of arguments
		if exeFile!="":
			runCommand(exeFile + " " + listR.listToStr(listR.applyOrderList(order,valueCmb)), os.path.dirname(exeFile), sleepTime)

		# execute postFnc
		if postFncs!=[]:
			for post_func in listR.toList(postFncs):
				post_func(valueCmb)


def _argListToDir(argNameList, valueCmb, baseDir, connectionSymbol="=", seperationSymbol=","):
	""" Return a associated directory (path+name) to
	argNameList and the corresponding valueCmb.
		argNameList can be one level nested, e.g:
		[[a,b,c],[d,e],f,...]
		The corresponding valueCmb is flattened:
		[1,2,3,4,5,6,...]
		And the generated dir name is (w/ default "=" & ",")
		a=1,b=2,c=3/d=4,e=5/f=6/...
	"""
	# make a list of type [[arg1=xx, arg2=xx],arg3=xx,...]:
	tmp = listR.mimic(argNameList, listR.strZip(listR.FL(argNameList), listR.FL(valueCmb), connectionSymbol))
	# join sublist using seperationSymbol:
	tmp = map(lambda x:seperationSymbol.join(listR.toList(x)), tmp)
	# join them with os.path.join
	return os.path.join(baseDir, os.path.sep.join(tmp))


def preExeMkdirPost(argNameList, argValueList, baseDir, preFncs=[], exeFile="", callOrderList=None, postFncs=[], sleepTime=0, connectionSymbol="=", seperationSymbol=","):
	"""
		Execute exeFile with a choice of arguments specified in
		argNameList. argValueList is a nested list, elements from
		each sublist consist a "choice" of arguments. order is applied
		before this choice of arguments is used. For each possible
		combination of values of arguments, preFncs are executed first,
		then exeFile is executed, then postFncs are executed.

		Note that exeFile is a string that contains the executable
		file, and the string can contain other arguments, as long as
		these arguments are positioned in front of those automatically
		generated arguments list of the form argNameList. The list
		argValueList is very free, for example, the list
		[["1 2","3 4"],...] will use "1 2" or "3 4" as the first
		generated parameter. This is especially useful if a program
		reads argument for a loop to generate one file.

		A directory tree is constructed, with type (with connectionSymbol="=")
		arg1=xx/arg2=xx/...
		If argNameList and argValueList are nested (once), such as
		argNameList=[[arg11,arg12],arg2,...] with argValueList as
		[[["a","e"],["b","f"]],["c","d"]...], then directories are
		like arg11=a, arg12=e/arg2=c/...
		argNameList and argValueList must be consistent, i.e. one
		name corresponds to one list of values; a group of names
		correspond to one list of group of values.

		preFnc is used to do some possible clean work between each
		run, it has two parameters: a list of choice of parameters, and
		a path of the associated directory. postFun is similar, only
		executed after each run.

		callOrderList specifies the order of parameters when the external
		program is called.
	"""
	if len(argNameList) != len(argValueList):
		print("The lists argNameList and argValueList have different length!\n")
		return False

	def mkdirHook(valueCmb):
		# make a directory corresponding to a value combination:
		makeDir(_argListToDir(argNameList, valueCmb, baseDir, connectionSymbol, seperationSymbol))

	# add mkdirHook to preFncs:
	preFncs=listR.toList(preFncs)
	preFncs.insert(0,mkdirHook)

	if callOrderList == None: callOrderList = argNameList
	order = listR.createOrderList(listR.FLL(callOrderList), listR.FLL(argNameList))

	preExePost(argValueList, preFncs, exeFile, postFncs, order, sleepTime)




def makeDirTree(baseDir, argNameList, argValueList, connectionSymbol="=", seperationSymbol=","):
	""" Construct a tree directories struction of the form
		baseDir/arg1=xx/arg2=xx/arg3=xx/...
	"""
	preExeMkdirPost(argNameList, argValueList, baseDir)



# need change:
def copyinPreExeMkdirCopyoutPost(exeFile, toBeCopiedDir, copyToBaseDir, toBeCopiedFileNames,
	argNameList, argValueList, callOrderList=None, sleepTime=0, preFncs=[], postFncs=[],
	connectionSymbol="=",seperationSymbol=","):
	"""
		Generate data files using exeFile. argNameList is a list of
		variable names. argValueList is a nested list, whose sublist
		at position i is a list of possible values of variable given at
		position i in argNameList. For each possible combination of
		values of arguments, exeFile is executed, then toBeCopiedFileNames
		in toBeCopiedDir are copied to generated directory under
		copyToBaseDir.The structure of the directories under
		copyToBaseDir are like (when connectionSymbol is "="):
		copyToBaseDir/arg1=xx/arg2=xx/arg3=xx/...
		the data file corresponding to a certain combination of
		values of argument is copied to the bottom directory with
		which the name of the path indicates the choice of arguments.

		callOrderList specifies the order of argument when being
		called by exeFile. For example, if argNameList is like
		["a1","a2","a3",...], callOrderList can be like
		["a3","a1","a2",...]

		preFnc is used to do some possible clean work between each
		run, it has two parameters: a list of choice of parameters
		and a path to the generated directory corresponding to this
		choice of paramters.

		Note that exeFile is a string that contains the executable
		file, and the string can contain other arguments, as long as
		these arguments are positioned in front of those automatically
		generated arguments list of the form argNameList. The list
		argValueList is very free, for example, the list
		[["1 2","3 4"],...] will use "1 2" or "3 4" as the first
		generated parameter. This is especially useful if a program
		reads argument for a loop to generate one file.

		See help for preExeMkdirPost for additional information.
	"""
	def copyinHook(valueCmb):
		# make a list of type [arg1=xx, arg2=xx,...], join them with os.path.join
		copyToFullDir = _argListToDir(argNameList, valueCmb, copyToBaseDir)
		copyFiles(toBeCopiedFileNames, toBeCopiedDir, copyToFullDir)
	# add copyinHook to preFncs:
	preFncs=listR.toList(preFncs)
	preFncs.insert(0,copyinHook)

	preExeMkdirPost(argNameList, argValueList, copyToBaseDir, preFncs,
		exeFile, callOrderList, postFncs, sleepTime, connectionSymbol, seperationSymbol)


# listR.readCSESD is more powerful than readCSED and readCSEFullpathD
def readCSED(dir_name, connectionSymbol="=", seperationSymbol=","): #CSE: Comma Seperated Equations
	""" Return a dic of the form {arg1:value1, ...} if with
	connectionSymbol="=" and seperationSymbol=",", dir_name is
	like arg1=value1,arg2=value2,...
	Values are in string form.
	"""
	if connectionSymbol not in dir_name: return {}
	return dict(map(lambda x:listR.split(x, connectionSymbol), listR.toList(dir_name.split(seperationSymbol))))


# listR.readCSESD is more powerful than readCSED and readCSEFullpathD
def readCSEFullpathD(fullpath, baseDir="", connectionSymbol="=", seperationSymbol=","):
	""" Return a dic of the form {arg1:value1, ...} for the full path. """
	dictList = []
	while fullpath!="":
		aDic = os.path.basename(fullpath)
		if connectionSymbol not in aDic: break
		dictList.extend(listR.itemsList(readCSED(aDic, connectionSymbol, seperationSymbol)))
		fullpath = os.path.dirname(fullpath)
	return dict(dictList)




def descendDirTree(baseDir, mustBeDefined, connectionSymbol="=", seperationSymbol=","):
	"""
		Return a list of directories and a list of the corresponding var_name:var_value
		dictionaries under baseDir. The valuse are given using string, like "1" instead of 1.
		The subdirectaries are so chosen that all variables in mustBeDefined
		must be defined.
	"""
	mustBeDefined = listR.toList(mustBeDefined)
	tmp_result = dirR.listDir(baseDir)
	def qf(var): # used in the filter, any result that does not fulled define vars in mustBeDefined will be filterd
		return listR.containedIn(mustBeDefined, readCSEFullpathD(var, connectionSymbol="=", seperationSymbol=",").keys())
	return filter(qf, tmp_result)


def descendDirTreeSharp(baseDir, mustBeDefined, connectionSymbol="=", seperationSymbol=","):
	"""
		Return a list of directories and a list of the corresponding var_name:var_value
		dictionaries under baseDir. The valuse are given using string, like "1" instead of 1.
		The subdirectaries are so chosen that all variables in mustBeDefined
		must be defined.

		Only those directories that barely defined mustBeDefined (sharp) are returned.
	"""
	tmp_result = descendDirTree(baseDir, mustBeDefined, connectionSymbol, seperationSymbol)
	if tmp_result==[]: return tmp_result
	tmp_list = []
	for aPath in tmp_result: # for each dir in result
		while list(listR.biIntersectI(mustBeDefined,readCSED(os.path.basename(aPath),connectionSymbol,seperationSymbol).keys()))==[]:
			# if the lowerest dir definds a variable in the mustBeDefined list
			aPath = os.path.dirname(aPath) # cut the lowerest dir
			if aPath == baseDir: break # if it's the shortest dir, return it
		tmp_list.append(aPath)
	tmp_list = listR.removeDuplicates(tmp_list)
	return tmp_list



def sortByColumn(data_file, column_to_sort=1):
    """
        Sort the data file by the specified column. The specified
        column must be numerical.
    """
    cmp = lambda qvar: float(qvar.split()[column_to_sort-1]) # used in "sorted" method
    in_file = open(data_file, "r")
    out_file = open(data_file+".TEMP", "w")
    for a_line in sorted(in_file.readlines(),key=cmp):
        out_file.write(a_line)
    in_file.close()
    out_file.close()
    shutil.copy(data_file+".TEMP", data_file)
    os.remove(data_file+".TEMP")


def switchColumn(data_file, column1, column2):
	"""
		Switch two columns specified by column1 and column2 in data_file.
	"""
	data = []
	for dataLine in readData(data_file):
		tmp = dataLine[column1-1]
		dataLine[column1-1] = dataLine[column2-1]
		dataLine[column2-1] = tmp
		data.append(dataLine)
	writeData(data_file, data)

def copyFile(filename, sourceDir, targetDir, renameTo=None, silent=True):
	""" Copy file from sourceDir to targetDir.
	"""
	if renameTo == None: renameTo = filename
	fullname_source = os.path.join(sourceDir, filename)
	fullname_target = os.path.join(targetDir, renameTo)
	shutil.copy(fullname_source, fullname_target)
	if silent==False:
		print("File "+fullname_source+" copied to "+source_dir)

def copy(source, target):
	""" Copy source file to target. """
	shutil.copy(source, target)


def copyFiles(fileNames, sourceDir, targetDir, silent=True):
	""" Copy files whose names are given in the list filenames
	from sourceDir to targetDir.
	"""
	fileNames = listR.toList(fileNames)
	for filename in fileNames:
		copyFile(filename, sourceDir, targetDir, None, silent)


def xcopy(namePatterns, sourceDir, targetDir, renameTo=None, flags=None):
	"""
		Copy files that match namePatterns in sourceDir to targetDir.
		All files match namePatterns are copied if renameTo is set to None;
		otherwise the first one matches is copied and renamed to renameTo.
		flags are used in listFilesMatch function (indirectly in re module).
	"""
	nameL = dirR.listFilesMatch(sourceDir, namePatterns, flags)
	if len(nameL) == 0: return
	if not os.path.exists(targetDir): makeDir(targetDir)
	if renameTo == None:
		for name in nameL:
			full_source_path = os.path.join(sourceDir, name)
			full_target_path = os.path.join(targetDir, name)
			shutil.copy(full_source_path, full_target_path)
	else:
		full_source_path = os.path.join(sourceDir, nameL[0])
		full_target_path = os.path.join(targetDir, renameTo)
		shutil.copy(full_source_path, full_target_path)



def nestedXcopy(namePatterns, sourceDir, targetDir, renameTo=None, flags=None):
	"""
		Copy files that match namePatterns under sourceDir to the corresponding
		directories under targetDir (i.e. they have the same relative path).
		All files match namePatterns are copied if renameTo is set to None;
		otherwise the first one matches is copied and renamed to renameTo.
		flags are used in listFilesMatch function (indirectly in re module).
	"""
	for aDir in dirR.listNestedDirContainsOneOfFilesM(sourceDir, namePatterns, flags):
		xcopy(namePatterns, aDir, os.path.join(targetDir, dirR._relativePathString(sourceDir, aDir)), renameTo, flags)



def readData(filename, commentSymbol="#"):
	"""
		Read a data file and return a nested list (data block).
		Each line contains data from each row (sub-list).
		All lines contain the commentSymbol are ignored.
	"""
	inFile = open(filename, "r")
	data = []
	for aLine in inFile.readlines():
		if aLine.find(commentSymbol)!=-1: continue
		lineData = []
		for piece in aLine.split():
			if listR.isFloat(piece): # is numerical
				lineData.append(float(piece))
			else: # is a string
				lineData.append(piece)
		data.append(lineData)
	inFile.close()
	return data



def writeData(filename, data, seperator="   "):
	"""
		Write a nested list (data block) into a file.
		Each line contains data from each row (sub-list).
	"""
	outFile = open(filename, "w")
	for dataLine in data:
		outFile.write(seperator.join(map(repr, dataLine))+"\n")
	outFile.close()



def nestedRenameFiles(dir_path, old_filenames, new_filenames, silent_level=0):
    """
        Rename all files in "old_filenames" under "dir_path" to
        "new_filenames". The first file in "old_filenames" will be
        renamed to the first file in "new_filenames", and similar for
        the rest. If "leaf_only" is specified as "True" then only
        files in leaf subdirectories are modified. "silent_level"=0:
        no output on screen; 1: short output; 2: full output
    """
    old_filenames = listR.toList(old_filenames)
    new_filenames = listR.toList(new_filenames)
    for old_name, new_name in zip(old_filenames, new_filenames):
        dirL = dirR.listNestedDirContainsFiles(dir_path, old_name)
        for aDir in dirL:
            if os.path.exists(os.path.join(aDir, new_name)):
                print("File "+os.path.join(aDir, new_name)+" already exists! skipped.")
            else:
                os.rename(os.path.join(aDir, old_name), os.path.join(aDir, new_name))
                if silent_level==0:
                    pass
                elif silent_level==1:
                    print("File "+_relativePath(dir_path, full_path)
                          +" renamed to "+_relativePath(dir_path,new_full_path))
                else:
                    print("File "+full_path+" renamed to "+new_full_path)





def nestedRenameFilesAdd(dir_path, filenames, str_add, add_to_front=True):
	"""
		Rename all files in "filenames" under "dir_path" by adding
		the string "str_add" to the front. If "add_to_front" is
		specified as "False", changes will be made to the end of the
		file name, otherwise (by default) it will be add to the front
		of the file name.
	"""
	filenames = listR.toList(filenames)
	for name in filenames:
		if add_to_front == True:
			str_mode = str_add+"%s%s"
		else:
			str_mode = "%s" + str_add + "%s"
		new_filename = str_mode % (os.path.splitext(name))
		nestedRenameFiles(dir_path, name, new_filename)


def nestedDeleteFiles(dir_path, filenames, silence_level=0, leaf_only=False):
    """
        Delete all files in the list "filenames" under "dir_path".
    """
    filenames = listR.toList(filenames)
    for name in filenames:
        dirL = dirR.listNestedDirContainsFiles(dir_path, name)
        for aDir in dirL:
            os.remove(os.path.join(aDir, name))
            if silence_level>0: print("File "+full_path+" deleted.")


def groupingDataOneL(dir_path, data_filename, frontAdd=None, isValid=None, change_name_to=""): # low level function
	"""
		Group data files by copying them into one large file with the
		same name placed in the subdirectory one level up. If
		"frontAdd" is given, a string suggested by "frontAdd" function
		using directory name as argument will be added to the file. If
		"isValid" is given, only for those directories that it
		returns True (the input of isValid is the full path of the dir), the
		data file will be combined, and the string "change_name_to"
		will be added to the tail of the combined data
		file name.

	"""
	if isValid == None: # build a trivial "isValid" function
		isValid = lambda qvar: True
		cmb_data_filename = data_filename
	else: # modify file name for the combined data file
		if change_name_to == "": change_name_to = data_filename
		cmb_data_filename = change_name_to

	untreated = []
	for aHLDir in dirR.nested_oneL_oneSubDir_hasAll(dir_path, data_filename):
		if __q_debug__: print(aHLDir)
		toWrite = open(os.path.join(aHLDir, cmb_data_filename), "w") # the combined data file is created here
		is_empty = True # combined data file is empty <==> no treatment
		for aDir in os.listdir(aHLDir):
			full_path = os.path.normpath(os.path.join(aHLDir, aDir))
			if os.path.isdir(full_path) == False: continue # not a directory
			# if hasNoSubDir(full_path) == False: continue # not a bottom directory
			if os.path.exists(os.path.join(full_path, data_filename)) == False: continue # does not contain the data file
			if not isValid(full_path): continue # skip this directory according to "isValid"
			is_empty = False
			if frontAdd != None: # add a string (usually a number) suggested by the name of dir in the combined data file
				toWrite.write(frontAdd(aDir)+" ")
			toRead = open(os.path.join(full_path, data_filename), "r")
			textBuffer = toRead.read()
			if textBuffer[-1:] != "\n": textBuffer = textBuffer + "\n" # force one and only one return after each data file
			toWrite.write(textBuffer)
			toRead.close()
		toWrite.close()
		if is_empty == True:
			os.remove(os.path.join(aHLDir, cmb_data_filename))
			untreated.append(aHLDir)

	return untreated



def makeDir(dir_path):
    """ Make directory at dir_path. If parent directory does not exist, it is created too. """
    if os.path.exists(dir_path): return
    dir_path = os.path.realpath(dir_path)
    dir_path = os.path.normpath(dir_path)
    if os.path.exists(os.path.dirname(dir_path)):
        os.mkdir(dir_path)
    else:
        makeDir(os.path.dirname(dir_path))
        os.mkdir(dir_path)


def removeDir(dir_path):
	""" Remove a directory. """
	if not os.path.exists(dir_path): return
	for name in os.listdir(dir_path):
		full_path = os.path.join(dir_path, name)
		if os.path.isfile(full_path):
			os.remove(full_path)
			continue
		if os.path.isdir(full_path):
			removeDir(full_path)
			continue
	os.rmdir(dir_path)


def delete(file):
	""" Delete a file. """
	if not os.path.exists(file): return
	os.remove(file)


def extractToken(filename, token, numOfLines=2):
	""" Return a list of strings consist of numOfLines lines
	in filename that follow immediately after the first line
	that contains string token.
	"""
	file = open(filename, "r")
	lines = file.readlines()
	file.close()
	for aLine in lines:
		if aLine.find(token) != -1:
			return lines[lines.index(aLine)+1:lines.index(aLine)+1+numOfLines]
	return ""

def collectFile(pathDir, filenames, targetDir=None):
	""" Collect files of name filenames under pathDir to the targetDir,
	then rename them according to the parameters.
	"""
	filenames = listR.toList(filenames)
	if targetDir == None: targetDir=os.getcwd()
	for aDir in dirR.listNestedDirContainsFiles(pathDir, filenames):
		for aFile in filenames:
			paraList = readCSEFullpathD(aDir).items() # read directory structure
			toAdd="-" + ",".join(map(lambda x:"=".join(x), paraList)) # transform the directory strunction into a string
			changedName = os.path.splitext(aFile)[0] + toAdd + os.path.splitext(aFile)[1] # add the string to the name of the file
			copyFile(aFile, aDir, targetDir, changedName) # copy file


def addColumnsToFile(filename, columns, add_before_original=True):
	""" Add columns of data into a file with filename, before or after the original data
	in each line. The variable "columns" should be a list of strings. Each string will be
	inserted accordingly into the file. This list will be re-used if it is shorter than the
	length of the file.
	"""
	tempFile = "TEMP.tmp"
	outFile = open(tempFile, "w")
	inFile = open(filename, "r")
	columns = listR.toList(columns)
	index = 0
	aLine = inFile.readline()
	while aLine:
		if add_before_original:
			outFile.write(columns[index] + aLine)
		else:
			outFile.write(aLine + columns[index])
		index = listR.next(columns, index)
		aLine = inFile.readline()
	inFile.close
	outFile.close()
	copy(tempFile, filename)


def combineFilesWithParas(dir_path, filename, targetDir=None, connector="="):
    """ Combine all files with name filename. Also, parameters will be added as seperated
    columns in the combined file. The order of parameters of corresponding to the inserted
    columns is indicated by the file name. """
    filenameL = listR.toList(filename)

    for filename in filenameL: # do the following for each file:
        # get list of common parameters:
        allParas = []
        for aDir in dirR.listNestedDirContainsFiles(dir_path, filename):
            allParas.append(listR.readCSESD(os.path.join(aDir, filename)).keys())
        allParas = listR.removeDuplicates(listR.intersect(allParas))

        # generate modified file name:
        if targetDir == None: targetDir=os.getcwd()
        new_filename = os.path.join(targetDir,os.path.splitext(filename)[0]+"("+",".join(listR.stringizeL(allParas))+")"+os.path.splitext(filename)[1])

        # generate the combined data file:
        outFile = open(new_filename, "w")
        for aDir in dirR.listNestedDirContainsFiles(dir_path, filename):
            to_read = os.path.join(aDir, filename) # target file: path+name
            parasD = listR.readCSESD(to_read) # get a dictionary of parameter values
            para_str = " ".join(listR.stringizeL(listR.getValueListFromDict(allParas, parasD))) # get the associated parameter list in string format
            inFile = open(to_read, "r")
            buffer = inFile.readlines() # read target file
            for aLine in buffer: # concarnate target file with parameter list:
                outFile.write(para_str+" "+aLine)
            inFile.close()
        outFile.close()

def takeRatioWithFirstLine(filename, columns=[2], newfilename=""):
		""" For all specified columns, convert all elements to the ratio of it over
		the 1st one in this column. Note that columns are indexed from 1. """
		# initialize
		if newfilename == "": newfilename = filename;
		columns = listR.toList(columns);

		buffer = readData(filename); # read data from disk
		first = list(buffer[0]); # a list of first elements in each column

		for i in range(len(buffer)): # loop through rows
				for j in columns: # loop through columns
						buffer[i][j-1] = buffer[i][j-1] / float(first[j-1]); # take ratio
		writeData(newfilename, buffer); # write data back to disk

def takeRatioWithFirstLineForAll(dir_path, filename, columns=[2], newfilename=""):
		""" For all files under dir_path, apply takeRatioWithFirstLine.
		Note that columns are indexed from 1. """
		if newfilename=="": newfilename = filename;
		for aDir in dirR.listNestedDirContainsFiles(dir_path, filename):
				fullpath_in = os.path.join(aDir, filename);
				fullpath_out = os.path.join(aDir, newfilename);
				takeRatioWithFirstLine(fullpath_in, columns, fullpath_out);


def getRowWithGivenValue(filename, value, column=1):
		""" Return the row in which the element in the specified column has a
		value closest to the specified one. Column is specified starting from 1."""
		buffer = readData(filename);
		min_diff = 99999999; # hopefully this is enough for real use
		target_row=[]; # the result
		for aLine in buffer:
				if abs(value-aLine[column-1]) < min_diff:
						min_diff = abs(value-aLine[column-1]);
						target_row = aLine;
		return target_row;


if __name__ == "__main__":
	print("Morning!")
