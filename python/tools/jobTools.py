import os

def getInputFileNames(inputFilePath, inputFile_matcher = None):
    inputFileNames = []
    files = os.listdir(inputFilePath)
    for file in files:
        if os.path.isdir(os.path.join(inputFilePath, file)):
            inputFileNames.extend(getInputFileNames(os.path.join(inputFilePath, file)))
        else:
            # add inputFile to list if name of inputFile matches regular expression
            # or if no regular expression has been given as argument to getInputFileNames function
            if inputFile_matcher is None or inputFile_matcher.match(file):
                inputFileNames.append("file:%s" % os.path.join(inputFilePath, file))
    return inputFileNames
