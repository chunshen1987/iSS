###################    Last edited on Sep 24, 2009 (8:58 pm)   ####################
#                                            version 22 --- Zhi Qiu
# level 1
"""
	A collection of functions helping finding or testing directories under variece condition.
	Any function that modifies the file structure is in fileR module.
"""

import os
import sys
import re


_q_debug_ = False

def _toList(qvar):
    """ Make qvar a list if not."""
    if type(qvar) != type([]): qvar = [qvar]
    return qvar

def __obsoleted_listNestedDir(dir_path): #--- low level function --- obsoleted ---
# **slightly** slower than _listNestedDir function
    """
        Make a list of all the nested subdirectories of "dir_path".
    """
    dirL=[dir_path]
    for name in os.listdir(dir_path):
        full_path = os.path.join(dir_path, name)
        if os.path.isdir(full_path):
            dirL.extend(__obsoleted_listNestedDir(full_path))
    return dirL



def _listNestedDir(dir_path): #--- low level function
    """
        Make a list of all the nested subdirectories of "dir_path".
    """
    return list(map(lambda qvar: qvar[0], os.walk(dir_path)))


def _listNestedDir2(dir_path): #--- low level function
    """
        Make an iterator of all the nested subdirectories of
        "dir_path", faster than _listNestedDir.
    """
    return map(lambda qvar: qvar[0], os.walk(dir_path))


def _listSpecificNestedDir(dir_path, judgeFunction): #--- low level function
    """
        Returns a list of subdirectories under "dir_path" specified by
        "judgeFunction". The judge function should have one argument,
        directory path, and return True if the directory is wanted.
    """
    return filter(judgeFunction, _listNestedDir(dir_path))




def hasNoSubDir(dir_path):
    """
        Returns True if "dir_path" has no subdirectories.
    """
    for name in os.listdir(dir_path):
        full_path = os.path.join(dir_path, name)
        if os.path.isdir(full_path): break
    else:
        return True
    return False


def _listNestedLeafDir(dir_path):
    """
        Make a list of all the nested leaf subdirectories of
        "dir_path".
    """
    return _listSpecificNestedDir(dir_path, hasNoSubDir)


def _q_match(pattern, target, flags=None):
	""" Return True if match, False if not."""
	if flags == None:
		return re.findall(pattern, target) != []
	else:
		return re.findall(pattern, target, flags) != []


def hasFiles(dir_path, filenames):
    """
        Returns True if "dir_path" contains all the files listed in
        "filenames".
    """
    filenames = _toList(filenames)
    for name in filenames:
        full_path = os.path.join(dir_path, name)
        if os.path.exists(full_path) == False: break
    else:
        return True
    return False


def listNestedDirContainsFiles(dir_path, filenames):
    """
        Make a list of all the nested subdirectories which contain all
        the file in "filenames".
    """
    def hasFilesHook(dir_path):
        return hasFiles(dir_path, filenames)
    return _listSpecificNestedDir(dir_path, hasFilesHook)




def hasFilesM(dir_path, filePatterns, flags=None): # M stands for "match"
	"""
		Returns True if every pattern defined in "filePatterns" matches some file
		in "dir_path". Pattens are matched using regular expression module.
		"flags" are standard flags argument used in re module.
	"""
	filenames = _toList(filePatterns)
	listNames = os.listdir(dir_path)
	for namePattern in filenames:
		for existedName in listNames:
			if _q_debug_ == True: print(existedName)
			if _q_match(namePattern, existedName, flags): # successfully find a match
				break # check next pattern
		else: # no match for one pattern
			return False
	else: # all patterns match
		return True


def listNestedDirContainsFilesM(dir_path, filePatterns, flags=None):
    """
        Make a list of all the nested subdirectories which contain all
        the file specified in "filePatterns". (see "hasFilesM")
    """
    def hasFilesHook(dir_path):
        return hasFilesM(dir_path, filePatterns, flags)
    return _listSpecificNestedDir(dir_path, hasFilesHook)



def hasOneOfFiles(dir_path, filenames):
    """
        Returns True if "dir_path" contains one of the files listed in
        "filenames".
    """
    filenames = _toList(filenames)
    for name in filenames:
        full_path = os.path.join(dir_path, name)
        if os.path.exists(full_path) == True: break
    else:
        return False
    return True


def listNestedDirContainsOneOfFiles(dir_path, filenames):
	"""
		Make a list of all the nested subdirectories which contain at least one of the
		file in "filenames".
	"""
	def hasOneOfFilesHook(dir_path):
		return hasOneOfFiles(dir_path, filenames)
	return _listSpecificNestedDir(dir_path, hasOneOfFilesHook)




def hasOneOfFilesM(dir_path, filePatterns, flags=None):
	"""
		Returns True if some pattern defined in "filePatterns" matches some
		file in "dir_path". Pattens are matched using regular expression module.
		"flags" are standard flags used in re module.
	"""
	filenames = _toList(filePatterns)
	listNames = os.listdir(dir_path)
	for namePattern in filenames:
		for existedName in listNames:
			if _q_debug_ == True: print(existedName)
			if _q_match(namePattern, existedName, flags): # successfully find a match
				return True # already satisfy the criteria
	else: # no pattern matches
		return False


def listNestedDirContainsOneOfFilesM(dir_path, filePatterns, flags=None):
    """
        Make a list of all the nested subdirectories which contain at least one ofl
        the files specified in "filePatterns". (see "hasOneOfFilesM")
    """
    def hasFilesHook(dir_path):
        return hasOneOfFilesM(dir_path, filePatterns, flags)
    return _listSpecificNestedDir(dir_path, hasFilesHook)


def _oneLevelAbove_allSubDirSatisfies(dir_path, judgeFunction): #--- low level function
    """
        Returns True if "judgeFunction" return True on all the
        subdirectories of "dir_path".
    """
    if hasNoSubDir(dir_path): return False
    for name in os.listdir(dir_path):
        full_path = os.path.join(dir_path, name)
        if os.path.isdir(full_path) == False: continue
        if judgeFunction(full_path) == False: break
    else:
        return True
    return False


def nested_oneL_allSubDir_hasAll(dir_path, filenames):
    """
        Returns a list of subdirectories under "dir_path" which is one
        level above its subdirectories that all contain all the files
        in "filenames".
    """
    def hook(dpath):
        return _oneLevelAbove_allSubDirSatisfies(dpath,
                                                lambda qvar: hasFiles(qvar, filenames))
    return _listSpecificNestedDir(dir_path, hook)


def nested_oneL_allSubDir_hasOneOf(dir_path, filenames):
    """
        Returns a list of subdirectories under "dir_path" which is one
        level above its subdirectories that all contain at least one
        of the files in "filenames".
    """
    def hook(dpath):
        return _oneLevelAbove_allSubDirSatisfies(dpath,
                                                lambda qvar: hasOneOfFiles(qvar, filenames))
    return _listSpecificNestedDir(dir_path, hook)



def _oneLevelAbove_oneSubDirSatisfies(dir_path, judgeFunction): #--- low level function
    """
        Returns True if "judgeFunction" return True on one of the
        subdirectories of "dir_path".
    """
    for name in os.listdir(dir_path):
        full_path = os.path.join(dir_path, name)
        if os.path.isdir(full_path) == False: continue
        if judgeFunction(full_path) == True: break
    else:
        return False
    return True


def nested_oneL_oneSubDir_hasAll(dir_path, filenames):
    """
        Returns a list of subdirectories under "dir_path" which is one
        level above its subdirectories at least one of which contains
        all the files in "filenames".
    """
    def hook(dpath):
        return _oneLevelAbove_oneSubDirSatisfies(dpath,
                                                lambda qvar: hasFiles(qvar, filenames))
    return _listSpecificNestedDir(dir_path, hook)




def nested_oneL_oneSubDir_hasOneOf(dir_path, filenames):
    """
        Returns a list of subdirectories under "dir_path" which is one
        level above its subdirectories at least one of which contains
        one of the files in "filenames".
    """
    def hook(dpath):
        return _oneLevelAbove_oneSubDirSatisfies(dpath,
                                                lambda qvar: hasOneOfFiles(qvar, filenames))
    return _listSpecificNestedDir(dir_path, hook)



def listDir(dir_path, leaf_only=False):
    """
        Make a list of all the nested leaf subdirectories of
        "dir_path". If "leaf_only" is True only leaf directories are
        given.
    """
    if leaf_only == False:
        functionToUse = _listNestedDir
    else:
        functionToUse = _listNestedLeafDir
    dirL = functionToUse(dir_path)
    return dirL




def listFilesMatch(dir_path, patterns, flags=None):
	"""
		Return a list of file names (directly) in dir_path that match one of patterns.
		flags are used in re module.
	"""
	patterns = _toList(patterns)
	result = []
	for name in os.listdir(dir_path):
		# file names loop is the ourside loop to ensure that only one file name
		# is in the list if more patterns match the same file name
		if os.path.isdir(os.path.join(dir_path, name)): continue
		for pattern in patterns:
			if _q_match(pattern, name, flags): # match one of the patterns
				result.append(name)
				break # to next name
	return result



def _relativePath(dir_path, full_path): #--- low level function
    """
        Give a string represent the relative path of "full_path"
        relative to "dir_path", beginning with "./".
    """
    return os.path.join(".",full_path[len(os.path.normpath(dir_path))+1:])


def _relativePathString(dir_path, full_path): #--- low level function
	"""
		Give a string represent the relative path of "full_path"
		relative to "dir_path", beginning with the name of the
		next level directory.
	"""
	return full_path[len(os.path.normpath(dir_path))+1:]



def lookUpForFiles(dir_path, filenames, level=2, notExact=False):
	"""
		Give a full path to the file filename in the parent directories (including
		self, level=2 means with one direct parent directory) of dir_path within
		"level". Return None if not found with levels. If notExact is True, pattern
		match searching method will be used.
	"""
	if notExact == False:
		matchFnc = hasFiles
	else:
		matchFnc = hasFilesM

	currentPath = os.path.normpath(dir_path)
	for ii in range(level):
		if matchFnc(currentPath, filenames):
			return currentPath
		else:
			currentPath = os.path.dirname(currentPath)
	return None



def lookUpForOneOfFiles(dir_path, filenames, level=2, notExact=False):
	"""
		Give a full path to one of the file filename in the parent directories (including
		self, level=2 means with one direct parent directory) of dir_path within
		"level". Return None if not found with levels. If notExact is True, pattern
		match searching method will be used.
	"""
	if notExact == False:
		matchFnc = hasOneOfFiles
	else:
		matchFnc = hasOneOfFilesM

	currentPath = os.path.normpath(dir_path)
	for ii in range(level):
		if matchFnc(currentPath, filenames):
			return currentPath
		else:
			currentPath = os.path.dirname(currentPath)
	return None



def listSubDirectories(dir_path):
	"""
		Return a list of subdirectories (full path) in the folder specified by dir_path.
	"""
	if not os.path.exists(dir_path):
		print("Path "+dir_path+" does not exist.");
		return [];
	return filter(os.path.isdir, map(lambda x: os.path.join(dir_path,x), os.listdir(dir_path)));



if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("Welcome to dirR module. Use help(dirR) for information. (Zhi Qiu, 2009)")
    else:
        print("Executing: "+sys.argv[1]+"('"+"','".join(sys.argv[2:])+"')")
        exec(sys.argv[1]+"('"+"','".join(sys.argv[2:])+"')")
