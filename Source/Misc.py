import os
import sys

rv_dict = {'A':'T',
           'G':'C',
           'C':'G',
           'T':'A',
           'N':'N'}

class Misc():

    @staticmethod
    def checkDir(dir, dirName):
    # Check if dir exists, and is writeable, else create new dir.

        if os.path.exists(dir) and os.path.isdir(dir):
            # Check if write permission is granted for dir
            if os.access(dir, os.W_OK):
                sys.stdout.write("Set {} at: {}\n".format(dirName, dir))
            else:
                raise FileNotFoundError("{} exists but cannot be written.\n".format(dirName))
        else:
            os.makedirs(dir)
            sys.stdout.write("Set {} at: {}\n".format(dirName, dir))

    @staticmethod
    def rcSeq(seq):
    #Return Reverse Complement of Sequence
        rcSeq = "".join([rv_dict[nuc] for nuc in seq[::-1]])
        return rcSeq

    @staticmethod
    def checkProg(progPath, progName):
    # Check if Software Path exists and can be executed

        if os.path.exists(progPath) and os.access(progPath, os.X_OK):
            sys.stdout.write("{} Accessible at {}\n".format(progName, progPath))
        else:
            raise FileNotFoundError("{} cannot be accessed at {}\n".format(progName, progPath))

    @staticmethod
    def checkFile(filePath, fileName):
    # Check if file exists at path and can be read

        if os.path.exists(filePath) and os.access(filePath, os.R_OK):
            sys.stdout.write("{} Accessible at {}\n".format(fileName, filePath))
        else:
            raise FileNotFoundError("{} cannot be accessed at {}\n".format(fileName, filePath))
